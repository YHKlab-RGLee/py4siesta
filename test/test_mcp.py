import ast
import asyncio
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from py4siesta_mcp.catalog import Catalog, tool_name
from py4siesta_mcp.runtime import Runtime


class McpCatalogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = Catalog()

    def test_every_public_source_function_and_class_is_accounted_for(self):
        names = {entry['api'] for entry in self.catalog.entries}
        root = Path(__file__).resolve().parent.parent
        for package in ('nanocore', 'py4siesta'):
            for path in (root / package).rglob('*.py'):
                if path.name.startswith('__'): continue
                try: tree = ast.parse(path.read_text())
                except SyntaxError: continue
                module = '.'.join(path.relative_to(root).with_suffix('').parts)
                for node in tree.body:
                    if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and not node.name.startswith('_'):
                        self.assertIn(module + '.' + node.name, names)
                        if isinstance(node, ast.ClassDef):
                            for member in node.body:
                                if isinstance(member, ast.FunctionDef) and not member.name.startswith('_'):
                                    self.assertIn(module + '.' + node.name + '.' + member.name, names)

    def test_catalog_has_contracts_cli_commands_and_no_imported_library_tools(self):
        self.assertGreater(len(self.catalog.tools), 300)
        for entry in self.catalog.tools.values():
            for field in ('purpose', 'units', 'returns', 'side_effects', 'input_schema', 'signature'):
                self.assertTrue(entry[field], (entry['api'], field))
            self.assertLessEqual(len(entry['tool']), 64)
        self.assertIn('py4siesta_tool_band', self.catalog.tools)
        self.assertIn('nanocore_atoms_AtomsSystem_calc_connectivity', self.catalog.tools)
        self.assertNotIn('nanocore_atoms_np_load', self.catalog.tools)
        self.assertFalse(self.catalog.search(query='nanocore.quest.read_seqquest')['entries'][0]['available'])

    def test_runtime_can_compose_original_atom_and_system_apis(self):
        runtime = Runtime(self.catalog, Path.cwd())
        def call(api, parameters=None, object_id=None):
            args = {'parameters': parameters or {}}
            if object_id: args['object_id'] = object_id
            return runtime.invoke(tool_name(api), args)
        first = call('nanocore.atoms.Atom', {'symbol': 'C', 'position': [0, 0, 0]})
        second = call('nanocore.atoms.Atom', {'symbol': 'H', 'position': [1, 0, 0]})
        atoms = call('nanocore.atoms.AtomsSystem', {'atoms': [first, second]})
        call('nanocore.atoms.AtomsSystem.calc_connectivity', object_id=atoms['$object'])
        atom = call('nanocore.atoms.AtomsSystem.operator_getitem', {'i': 0}, atoms['$object'])
        result = call('nanocore.atoms.Atom.get_connectivity', {'details': True}, atom['$object'])
        self.assertEqual(result, [[2, ['C', 'H'], 1., [0, 0, 0]]])
        self.assertEqual(call('nanocore.atoms.AtomsSystem.distance', {'selected': [1, 2]}, atoms['$object']), 1.)
        with self.assertRaises(TypeError):
            call('nanocore.atoms.AtomsSystem.select_all', object_id=first['$object'])
        with self.assertRaises(ValueError):
            runtime.invoke('unregistered_function', {})

    def test_arrays_attributes_and_released_handles(self):
        runtime = Runtime(self.catalog, Path.cwd())
        values = np.arange(5000.)
        handle = runtime.encode(values)
        page = runtime.invoke('py4siesta_object_read', dict(object_id=handle['$object'], offset=10, limit=3))
        self.assertEqual(page['items']['$array'], [10., 11., 12.])
        self.assertIs(runtime.decode(handle), values)
        runtime.invoke('py4siesta_object_release', {'object_id': handle['$object']})
        with self.assertRaises(ValueError): runtime.decode(handle)

    def test_cli_adapter_preserves_repeated_arguments_and_failure_envelope(self):
        runtime = Runtime(self.catalog, Path.cwd())
        with mock.patch('py4siesta.tool_cli.execute', return_value={'ok': True}) as execute:
            result = runtime.invoke('py4siesta_tool_pdos', {'parameters': {
                'pdos_path': 'example.PDOS', 'orbital': [['C_0', 'H_0'], ['C_2_1']], 'emin': -4.}})
        self.assertTrue(result['ok'])
        argv = execute.call_args[0][0]
        self.assertEqual(argv.count('--orbital'), 2)
        self.assertIn('C_2_1', argv)
        self.assertIn('-4.0', argv)


@unittest.skipUnless(importlib.util.find_spec('mcp'), 'Optional MCP SDK is not installed')
class McpProtocolTests(unittest.TestCase):
    def test_stdio_discovery_structure_files_errors_and_logging(self):
        asyncio.run(asyncio.wait_for(self._exercise_stdio(), timeout=45))

    async def _exercise_stdio(self):
        from mcp import ClientSession, StdioServerParameters
        from mcp.client.stdio import stdio_client
        repository = Path(__file__).resolve().parent.parent
        environment = os.environ.copy()
        environment['PYTHONPATH'] = str(repository)
        environment['PYTHONDONTWRITEBYTECODE'] = '1'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'molecule.xyz').write_text('2\nCH\nC 0 0 0\nH 1 0 0\n')
            (root / 'stdout.txt').write_text('siesta:         Total = -12.5\n')
            parameters = StdioServerParameters(command=sys.executable,
                args=['-m', 'py4siesta_mcp', '--workdir', directory],
                cwd=str(repository), env=environment)
            async with stdio_client(parameters) as (read, write):
                async with ClientSession(read, write) as client:
                    await client.initialize()
                    listed = await client.list_tools()
                    names = {tool.name for tool in listed.tools}
                    self.assertGreater(len(names), 300)
                    self.assertIn('py4siesta_tool_submit', names)
                    resources = await client.list_resources()
                    self.assertEqual(len(resources.resources), 2)
                    catalog = await client.call_tool('py4siesta_api_catalog', {
                        'query': 'calc_connectivity', 'details': True})
                    self.assertEqual(catalog.structuredContent['result']['total'], 1)
                    env = await client.read_resource('py4siesta://runtime')
                    self.assertEqual(json.loads(env.contents[0].text)['workdir'], directory)

                    async def call(tool, arguments):
                        result = await client.call_tool(tool, arguments)
                        self.assertFalse(result.isError, result.content)
                        return result.structuredContent

                    loaded = await call('nanocore_io_read_xyz', {'parameters': {'file_name': 'molecule.xyz'}})
                    atoms = loaded['result']['$object']
                    await call('nanocore_atoms_AtomsSystem_calc_connectivity', {'object_id': atoms})
                    first = await call('nanocore_atoms_AtomsSystem_operator_getitem', {
                        'object_id': atoms, 'parameters': {'i': 0}})
                    details = await call('nanocore_atoms_Atom_get_connectivity', {
                        'object_id': first['result']['$object'], 'parameters': {'details': True}})
                    self.assertEqual(details['result'], [[2, ['C', 'H'], 1., [0, 0, 0]]])
                    await call('nanocore_atoms_AtomsSystem_set_cell', {
                        'object_id': atoms, 'parameters': {'cell_vector': [[8,0,0],[0,8,0],[0,0,8]]}})
                    await call('nanocore_siestaio_write_struct', {'parameters': {
                        'atoms': {'$object': atoms}, 'file_path': 'copy.fdf'}})
                    self.assertTrue((root / 'copy.fdf').is_file())
                    energy = await call('nanocore_siesta_get_total_energy', {'parameters': {'output_file': 'stdout.txt'}})
                    self.assertEqual(energy['result'], -12.5)
                    self.assertFalse((root / 'OUT').exists())
                    # A native function printing to stdout must not corrupt MCP framing.
                    description = await call('nanocore_atoms_AtomsSystem_operator_describe', {'object_id': atoms})
                    self.assertIn('Number of atoms = 2', description['log'])
                    error = await client.call_tool('py4siesta_tool_band', {
                        'parameters': {'bands_path': 'missing.bands'}})
                    self.assertTrue(error.isError)
                    error = await client.call_tool('nanocore_atoms_Atom_get_symbol', {'object_id': 'missing'})
                    self.assertTrue(error.isError)
                    invalid = await client.call_tool('nanocore_atoms_Atom', {'parameters': {'symbol': 'C'}})
                    self.assertTrue(invalid.isError)
                    alive = await client.call_tool('py4siesta_api_catalog', {'query': 'read_xyz'})
                    self.assertFalse(alive.isError)


if __name__ == '__main__':
    unittest.main()
