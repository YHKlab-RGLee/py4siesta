"""Discover local public APIs and describe them without running calculations."""

import argparse
import ast
import hashlib
import importlib
import inspect
import json
import re
import io
import textwrap
import tokenize
from pathlib import Path

from .metadata import EXCLUDED, PARAMETERS, describe


OPERATORS = {
    '__getitem__': 'operator_getitem', '__len__': 'operator_length', '__add__': 'operator_add',
    '__sub__': 'operator_subtract', '__mul__': 'operator_multiply', '__truediv__': 'operator_divide',
    '__neg__': 'operator_negate', '__abs__': 'operator_absolute', '__eq__': 'operator_equals',
    '__repr__': 'operator_describe',
}
REFERENCE = {'type': 'object', 'properties': {'$object': {'type': 'string'}},
             'required': ['$object'], 'additionalProperties': True}


def tool_name(api):
    name = re.sub(r'[^a-zA-Z0-9_]', '_', api).strip('_')
    if len(name) > 64:
        name = name[:53] + '_' + hashlib.sha256(api.encode()).hexdigest()[:10]
    return name


def parameter_schema(parameter):
    """Describe permissive legacy arguments without guessing hidden constraints."""
    annotation = parameter.annotation
    schema = {}
    if annotation in (int, float, bool, str):
        schema['type'] = {int: 'integer', float: 'number', bool: 'boolean', str: 'string'}[annotation]
    elif annotation is Path:
        schema['type'] = 'string'
    elif inspect.isclass(annotation) and annotation.__module__.startswith(('nanocore', 'py4siesta')):
        schema = dict(REFERENCE)
    default = parameter.default
    if default is not inspect.Parameter.empty:
        try:
            json.dumps(default, allow_nan=False)
            schema['default'] = default
        except (TypeError, ValueError):
            pass
    # Defaults are not type restrictions: nanocore commonly accepts multiple types.
    schema['description'] = PARAMETERS.get(
        parameter.name, 'Argument %s; accepted values follow the API documentation below.' % parameter.name)
    if annotation is not inspect.Parameter.empty:
        schema['description'] += ' Python annotation: %s.' % str(annotation)
    if parameter.kind == inspect.Parameter.VAR_POSITIONAL:
        schema.update(type='array', items={})
    elif parameter.kind == inspect.Parameter.VAR_KEYWORD:
        schema.update(type='object', additionalProperties=True)
    return schema


def input_schema(signature, instance=False):
    properties = {}; required = []
    for p in signature.parameters.values():
        properties[p.name] = parameter_schema(p)
        if p.default is inspect.Parameter.empty and p.kind not in (
                inspect.Parameter.VAR_KEYWORD, inspect.Parameter.VAR_POSITIONAL):
            required.append(p.name)
    arguments = {'type': 'object', 'properties': properties,
                 'required': required, 'additionalProperties': False}
    result = {'type': 'object', 'properties': {
        'parameters': arguments,
        'workdir': {'type': 'string', 'description': 'Existing calculation directory; defaults to the server workdir.'}},
        'required': ['parameters'] if required else [], 'additionalProperties': False}
    if instance:
        result['properties']['object_id'] = {
            'type': 'string', 'description': 'Session handle returned when this object was created/read.'}
        result['required'].append('object_id')
    return result


def public_nodes(tree):
    """Source declarations only: imported libraries are never exposed as APIs."""
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if not node.name.startswith('_'):
                yield node


class Catalog:
    def __init__(self):
        self.entries = []
        self.tools = {}
        self.modules = []
        for package in ('nanocore', 'py4siesta'):
            root = Path(importlib.import_module(package).__file__).parent
            for path in sorted(root.rglob('*.py')):
                if path.name.startswith('__'): continue
                module = package + '.' + '.'.join(path.relative_to(root).with_suffix('').parts)
                self._module(module, path)
        self._cli_commands()
        self.entries.sort(key=lambda item: item['api'])
        for entry in self.entries:
            if entry['available']:
                if entry['tool'] in self.tools:
                    raise ValueError('Duplicate MCP tool name: %s' % entry['tool'])
                self.tools[entry['tool']] = entry

    def _unavailable(self, api, reason, module, line=None):
        self.entries.append(dict(api=api, tool=tool_name(api), available=False,
                                 reason=reason, module=module, line=line,
                                 category='unavailable'))

    def _module(self, module, path):
        source = path.read_text()
        try:
            tree = ast.parse(source)
            loaded = importlib.import_module(module)
        except Exception as exc:
            reason = '%s: %s' % (type(exc).__name__, str(exc))
            self.modules.append(dict(module=module, available=False, reason=reason))
            # Keep old Python-2 declarations visible in the inventory, too.
            owner = None
            ignored = set()
            try:
                for token in tokenize.generate_tokens(io.StringIO(source).readline):
                    if token.type == tokenize.STRING and token.end[0] > token.start[0]:
                        ignored.update(range(token.start[0], token.end[0]+1))
            except (tokenize.TokenError, IndentationError):
                pass
            for lineno, line in enumerate(source.splitlines(), 1):
                if lineno in ignored: continue
                match = re.match(r'^(class|def)\s+([A-Za-z]\w*)', line)
                if match:
                    owner = match[2] if match[1] == 'class' else None
                    self._unavailable(module + '.' + match[2], reason, module, lineno)
                else:
                    match = re.match(r'^    def\s+([A-Za-z]\w*)', line)
                    if owner and match:
                        self._unavailable(module + '.' + owner + '.' + match[1], reason, module, lineno)
            return
        self.modules.append(dict(module=module, available=True))
        for node in public_nodes(tree):
            if isinstance(node, ast.ClassDef):
                cls = getattr(loaded, node.name)
                self._add(module, node.name, '', cls, 'constructor', node.lineno)
                # Include inherited project methods (operation.run, etc.), not object/FileIO internals.
                members = {}
                for base in reversed(cls.__mro__):
                    if not base.__module__.startswith(('nanocore', 'py4siesta')): continue
                    members.update(vars(base))
                for name, raw in sorted(members.items()):
                    if name.startswith('_') and name not in OPERATORS: continue
                    if isinstance(raw, property):
                        self._add(module, node.name, name, raw.fget, 'property', node.lineno)
                    elif callable(raw) or isinstance(raw, (staticmethod, classmethod)):
                        kind = 'class_method' if isinstance(raw, (staticmethod, classmethod)) else 'method'
                        self._add(module, node.name, name, getattr(cls, name), kind, node.lineno)
            else:
                self._add(module, '', node.name, getattr(loaded, node.name), 'function', node.lineno)

    def _add(self, module, owner, name, function, kind, line):
        api = '.'.join(part for part in (module, owner, name) if part)
        reason = EXCLUDED.get(api)
        doc = inspect.getdoc(function) or ''
        try:
            signature = inspect.signature(function)
            if kind in ('method', 'property'):
                signature = signature.replace(parameters=list(signature.parameters.values())[1:])
            source = inspect.getsource(function) if kind != 'constructor' else ''
            if re.search(r'^\s*raise NotImplementedError\b', source, re.M):
                reason = 'Abstract operation hook; use an implemented subclass.'
            if '@contextmanager' in source:
                reason = 'Context-manager lifecycle is not a standalone MCP call; use the operation.run API.'
            if 'return working_dir(' in source:
                reason = 'Returns a directory context manager; use the operation.run API.'
        except (TypeError, ValueError, OSError) as exc:
            self._unavailable(api, 'Cannot inspect public signature: %s' % exc, module, line)
            return
        metadata = describe(module, owner, name, kind, doc)
        # Keep output names available even when legacy docstrings lack Returns.
        outputs = []
        try: source_tree = ast.parse(textwrap.dedent(source.expandtabs()))
        except SyntaxError: source_tree = ast.Module(body=[], type_ignores=[])
        if source and hasattr(ast, 'unparse'):
            for node in ast.walk(source_tree):
                if isinstance(node, ast.Return) and node.value is not None:
                    expression = ast.unparse(node.value)
                    if len(expression) <= 200 and expression not in outputs:
                        outputs.append(expression)
        entry = dict(api=api, tool=tool_name(api.replace(name, OPERATORS.get(name, name)) if name in OPERATORS else api),
                     module=module, owner=owner, member=name, kind=kind,
                     signature=str(signature), documentation=doc, line=line,
                     return_expressions=outputs,
                     available=reason is None, reason=reason,
                     input_schema=input_schema(signature, kind in ('method', 'property')),
                     **metadata)
        # Operation.run forwards **kwargs to the subclass's case_parameters.
        if owner and name == 'run' and 'kwargs' in signature.parameters:
            cls = getattr(importlib.import_module(module), owner)
            case_parameters = getattr(cls, 'case_parameters', None)
            if case_parameters:
                case_signature = inspect.signature(case_parameters)
                case_signature = case_signature.replace(parameters=list(case_signature.parameters.values())[1:])
                kwargs_schema = input_schema(case_signature)['properties']['parameters']
                entry['input_schema']['properties']['parameters']['properties']['kwargs'] = kwargs_schema
                entry['purpose'] += ' kwargs follow case_parameters%s.' % case_signature
        self.entries.append(entry)

    def _cli_commands(self):
        from py4siesta.tool_cli import build_parser
        parser = build_parser()
        subparsers = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))
        for name, command in sorted(subparsers.choices.items()):
            properties = {}; required = []; actions = []
            for action in command._actions:
                if action.dest == 'help': continue
                schema = {'description': action.help or action.dest}
                schema['type'] = 'integer' if action.type is int or getattr(action.type, '__name__', '') == '_positive_int' else 'number' if action.type is float else 'string'
                if action.choices is not None: schema['enum'] = list(action.choices)
                if action.nargs in ('+', '*') or isinstance(action.nargs, int):
                    schema = {'type': 'array', 'items': schema}
                    if action.nargs == '+': schema['minItems'] = 1
                    if isinstance(action.nargs, int): schema.update(minItems=action.nargs, maxItems=action.nargs)
                if isinstance(action, argparse._AppendAction):
                    schema = {'type': 'array', 'items': schema}
                if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
                    schema = {'type': 'boolean', 'description': action.help or action.dest}
                if action.default is not None and action.default != argparse.SUPPRESS:
                    schema['default'] = action.default
                properties[action.dest] = schema
                if action.required: required.append(action.dest)
                actions.append(dict(dest=action.dest, flag=action.option_strings[-1] if action.option_strings else None,
                                    append=isinstance(action, argparse._AppendAction),
                                    boolean=isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)),
                                    default=action.default))
            api = 'py4siesta-tool.' + name
            self.entries.append(dict(
                api=api, tool=tool_name(api), module='py4siesta.tool_cli', owner='', member=name,
                kind='cli', available=True, reason=None, category='submission' if name == 'submit' else 'cli-tools',
                purpose=next((a.help for a in subparsers._choices_actions if a.dest == name), name),
                documentation=command.format_help(), signature=command.format_usage().strip(),
                units='Geometry/displacements: angstrom; energy windows: eV; ratios and counts: dimensionless.',
                returns='Existing tool_cli envelope: ok, command, result or error. Generated paths follow workdir.',
                side_effects=['May write/overwrite calculation files; submit runs sbatch and must not be blindly retried.'],
                read_only=False, actions=actions,
                input_schema={'type': 'object', 'properties': {
                    'parameters': {'type': 'object', 'properties': properties, 'required': required, 'additionalProperties': False},
                    'workdir': {'type': 'string', 'description': 'Existing calculation root, containing origin/ where required.'}},
                    'required': ['parameters'] if required else [], 'additionalProperties': False}))

    def search(self, query='', category=None, include_unavailable=True, offset=0, limit=30, details=False):
        entries = [e for e in self.entries if
                   (include_unavailable or e['available']) and
                   (not category or e['category'] == category) and
                   query.lower() in ('%s %s' % (e['api'], e.get('purpose', ''))).lower()]
        fields = ('api', 'tool', 'available', 'reason', 'category', 'purpose', 'signature')
        page = entries[offset:offset+limit]
        if not details: page = [{k: e[k] for k in fields if k in e} for e in page]
        return dict(total=len(entries), offset=offset, entries=page,
                    categories=sorted({e['category'] for e in self.entries}))

    def document(self):
        return dict(modules=self.modules, entries=self.entries,
                    note='Available means importable and callable through MCP, not scientifically validated for every input. Legacy failures are returned unchanged.')
