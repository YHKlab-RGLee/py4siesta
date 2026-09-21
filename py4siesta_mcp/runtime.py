"""Call cataloged APIs and translate Python objects at the MCP boundary."""

import contextlib
import importlib
import inspect
import io
import itertools
import math
import multiprocessing
import os
import sys
import tempfile
import threading
from collections.abc import ItemsView, Iterator
from pathlib import Path

import numpy as np


class Runtime:
    def __init__(self, catalog, workdir):
        self.catalog = catalog
        self.workdir = Path(workdir).resolve()
        self.objects = {}
        self.next_id = 1

    def get_object(self, object_id):
        if object_id not in self.objects:
            raise ValueError('Unknown or expired object_id: %s' % object_id)
        return self.objects[object_id]

    def store(self, value):
        for object_id, existing in self.objects.items():
            if existing is value: break
        else:
            object_id = 'obj_%d' % self.next_id
            self.next_id += 1
            self.objects[object_id] = value
        result = {'$object': object_id, 'type': '%s.%s' % (type(value).__module__, type(value).__name__)}
        if isinstance(value, np.ndarray): result.update(shape=list(value.shape), dtype=str(value.dtype))
        return result

    def decode(self, value):
        if isinstance(value, dict):
            if '$object' in value: return self.get_object(value['$object'])
            if '$array' in value:
                return np.array(self.decode(value['$array']), dtype=value.get('dtype'))
            if '$tuple' in value: return tuple(self.decode(v) for v in value['$tuple'])
            if '$float' in value: return float(value['$float'])
            if '$complex' in value: return complex(*value['$complex'])
            if '$bytes' in value: return bytes.fromhex(value['$bytes'])
            return {k: self.decode(v) for k, v in value.items()}
        if isinstance(value, list): return [self.decode(v) for v in value]
        return value

    def encode(self, value):
        if value is None or isinstance(value, (str, bool, int)): return value
        if isinstance(value, float):
            return value if math.isfinite(value) else {'$float': str(value)}
        if isinstance(value, np.generic): return self.encode(value.item())
        if isinstance(value, Path): return str(value)
        if isinstance(value, complex): return {'$complex': [value.real, value.imag]}
        if isinstance(value, bytes): return {'$bytes': value.hex()}
        if isinstance(value, np.ndarray):
            if value.size > 4096: return self.store(value)
            return {'$array': self.encode(value.tolist()), 'shape': list(value.shape), 'dtype': str(value.dtype)}
        if isinstance(value, dict): return {str(k): self.encode(v) for k, v in value.items()}
        if isinstance(value, (list, tuple, ItemsView)):
            if len(value) > 4096: return self.store(value)
            return [self.encode(v) for v in value]
        return self.store(value)

    def _arguments(self, function, parameters):
        signature = inspect.signature(function)
        positional = []; keyword = {}
        parameters = dict(parameters)
        varargs = any(p.kind == p.VAR_POSITIONAL for p in signature.parameters.values())
        for name, p in signature.parameters.items():
            if name not in parameters:
                if varargs and p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD) and p.default is not p.empty:
                    positional.append(p.default)
                continue
            value = self.decode(parameters.pop(name))
            if p.annotation is Path and isinstance(value, str): value = Path(value)
            if p.kind == p.VAR_POSITIONAL: positional.extend(value)
            elif p.kind == p.VAR_KEYWORD:
                if any(k in signature.parameters for k in value):
                    raise ValueError('kwargs duplicates a named parameter')
                keyword.update(value)
            elif p.kind == p.POSITIONAL_ONLY or varargs and p.kind == p.POSITIONAL_OR_KEYWORD:
                positional.append(value)
            else: keyword[name] = value
        if parameters: raise ValueError('Unexpected parameters: %s' % ', '.join(parameters))
        signature.bind(*positional, **keyword)
        return positional, keyword

    def invoke(self, tool, arguments):
        if tool == 'py4siesta_api_catalog': return self.catalog.search(**arguments)
        if tool == 'py4siesta_object_read':
            value = self.get_object(arguments['object_id'])
            if 'attribute' in arguments:
                attribute = arguments['attribute']
                if attribute.startswith('_'): raise ValueError('Only public data attributes can be read')
                value = getattr(value, attribute)
                if callable(value): raise ValueError('Use the cataloged method tool')
                return self.encode(value)
            start = arguments.get('offset', 0); count = arguments.get('limit', 100)
            if isinstance(value, Iterator):
                return {'items': self.encode(list(itertools.islice(value, count))), 'consumed': True}
            if isinstance(value, (np.ndarray, list, tuple)):
                part = value.reshape(-1)[start:start+count] if isinstance(value, np.ndarray) else value[start:start+count]
                return {'items': self.encode(part), 'offset': start,
                        'total': int(value.size) if isinstance(value, np.ndarray) else len(value)}
            # Inspect data objects such as BandStructureData without calling methods.
            return {'type': self.store(value)['type'], 'attributes': {
                k: self.encode(v) for k, v in vars(value).items() if not k.startswith('_')
            }} if hasattr(value, '__dict__') else self.store(value)
        if tool == 'py4siesta_object_release':
            value = self.get_object(arguments['object_id'])
            if isinstance(value, io.IOBase): value.close()
            del self.objects[arguments['object_id']]
            return {'released': arguments['object_id']}
        entry = self.catalog.tools.get(tool)
        if entry is None: raise ValueError('Unknown or unavailable API tool: %s' % tool)
        parameters = arguments.get('parameters', {})
        if entry['kind'] == 'cli':
            from py4siesta.tool_cli import execute
            argv = [entry['member']]
            for action in entry['actions']:
                if action['dest'] not in parameters: continue
                value = parameters[action['dest']]
                if action['boolean']:
                    if value != action['default']: argv.append(action['flag'])
                    continue
                values = value if action['append'] else [value]
                for group in values:
                    if action['flag']: argv.append(action['flag'])
                    argv.extend(str(v) for v in (group if isinstance(group, list) else [group]))
            return execute(argv)
        module = importlib.import_module(entry['module'])
        if entry['kind'] in ('method', 'property'):
            obj = self.get_object(arguments['object_id'])
            cls = getattr(module, entry['owner'])
            if not isinstance(obj, cls): raise TypeError('object_id must refer to %s' % cls.__name__)
            function = getattr(obj, entry['member'])
            if entry['kind'] == 'property': return self.encode(function)
        elif entry['kind'] == 'constructor': function = getattr(module, entry['owner'])
        elif entry['kind'] == 'class_method': function = getattr(getattr(module, entry['owner']), entry['member'])
        else: function = getattr(module, entry['member'])
        positional, keyword = self._arguments(function, parameters)
        return self.encode(function(*positional, **keyword))


@contextlib.contextmanager
def captured_output():
    """Capture Python and subprocess output inside the dedicated worker only."""
    with tempfile.TemporaryFile(mode='w+b') as stream:
        sys.stdout.flush(); sys.stderr.flush()
        saved = [os.dup(1), os.dup(2)]
        try:
            os.dup2(stream.fileno(), 1); os.dup2(stream.fileno(), 2)
            yield stream
        finally:
            sys.stdout.flush(); sys.stderr.flush()
            os.dup2(saved[0], 1); os.dup2(saved[1], 2)
            for fd in saved: os.close(fd)


def worker_main(connection, catalog, workdir):
    runtime = Runtime(catalog, workdir)
    while True:
        try: request = connection.recv()
        except EOFError: break
        if request is None: break
        tool, arguments = request
        previous = Path.cwd()
        with captured_output() as output:
            try:
                target = Path(arguments.get('workdir', runtime.workdir)).expanduser()
                if not target.is_absolute(): target = runtime.workdir / target
                if not target.is_dir(): raise ValueError('workdir must be an existing directory')
                os.chdir(target)
                value = runtime.invoke(tool, arguments)
                failed = isinstance(value, dict) and value.get('ok') is False
                result = dict(ok=not failed, tool=tool, result=value, workdir=str(target.resolve()))
            except (Exception, SystemExit) as exc:
                result = dict(ok=False, tool=tool, error=dict(type=type(exc).__name__, message=str(exc)))
            finally:
                os.chdir(previous)
            sys.stdout.flush(); sys.stderr.flush()
            output.seek(0)
            log = output.read(65537)
            if log: result['log'] = log[:65536].decode('utf-8', errors='replace')
            if len(log) > 65536: result['log_truncated'] = True
        connection.send(result)
    for value in runtime.objects.values():
        if isinstance(value, io.IOBase): value.close()
    connection.close()


class Worker:
    """One serialized worker keeps cwd changes and legacy stdout off the wire."""
    def __init__(self, catalog, workdir):
        context = multiprocessing.get_context('spawn')
        self.connection, child = context.Pipe()
        self.process = context.Process(target=worker_main, args=(child, catalog, str(workdir)), daemon=True)
        self.process.start()
        child.close()
        self.lock = threading.Lock()

    def call(self, tool, arguments):
        with self.lock:
            if not self.process.is_alive(): raise RuntimeError('API worker stopped; restart the MCP session. Do not retry submissions blindly.')
            self.connection.send((tool, arguments))
            try: return self.connection.recv()
            except EOFError: raise RuntimeError('API worker exited during the call; side effects may already have occurred.')

    def close(self):
        if self.process.is_alive():
            self.connection.send(None)
            self.process.join(timeout=2)
            if self.process.is_alive(): self.process.terminate(); self.process.join(timeout=2)
        self.connection.close()
