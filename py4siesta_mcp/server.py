"""MCP stdio transport; scientific work stays in the existing APIs."""

import asyncio
import json
import os
import sys
from pathlib import Path

from .catalog import Catalog
from .runtime import Worker


def support_tools():
    from mcp import types
    return [
        types.Tool(name='py4siesta_api_catalog', description=(
            'Search the complete public API inventory by name/category. details=True '
            'includes inputs, return conventions, units, effects and original documentation. '
            'Unavailable legacy APIs include their exclusion reason.'),
            inputSchema={'type': 'object', 'properties': {
                'query': {'type': 'string', 'default': ''},
                'category': {'type': 'string'},
                'include_unavailable': {'type': 'boolean', 'default': True},
                'details': {'type': 'boolean', 'default': False},
                'offset': {'type': 'integer', 'minimum': 0, 'default': 0},
                'limit': {'type': 'integer', 'minimum': 1, 'maximum': 100, 'default': 30}},
                'additionalProperties': False},
            annotations=types.ToolAnnotations(readOnlyHint=True, openWorldHint=False)),
        types.Tool(name='py4siesta_object_read', description=(
            'Read public data attributes of a session object (for example BandStructureData.energies), '
            'or page a large array/list. Arrays are flattened in row-major order. '
            'Iterator reads consume up to limit items. Use API method tools for calculations.'),
            inputSchema={'type': 'object', 'properties': {
                'object_id': {'type': 'string'}, 'attribute': {'type': 'string', 'pattern': '^[^_]'},
                'offset': {'type': 'integer', 'minimum': 0, 'default': 0},
                'limit': {'type': 'integer', 'minimum': 1, 'maximum': 4096, 'default': 100}},
                'required': ['object_id'], 'additionalProperties': False},
            annotations=types.ToolAnnotations(readOnlyHint=False, destructiveHint=False, openWorldHint=False)),
        types.Tool(name='py4siesta_object_release', description='Release a session object handle and close it if it is a file.',
            inputSchema={'type': 'object', 'properties': {'object_id': {'type': 'string'}},
                         'required': ['object_id'], 'additionalProperties': False},
            annotations=types.ToolAnnotations(readOnlyHint=False, destructiveHint=False, openWorldHint=False)),
    ]


def api_tool(entry):
    from mcp import types
    description = '%s\nAPI: %s%s\nCategory: %s\nUnits: %s\nReturns: %s\nEffects: %s' % (
        entry['purpose'], entry['api'], entry['signature'], entry['category'],
        entry['units'], entry['returns'], ' '.join(entry['side_effects']))
    if entry.get('return_expressions'):
        description += '\nNative return expressions: ' + '; '.join(entry['return_expressions'])
    description += ('\nPass objects as {"$object":"obj_N"}; use object_id for instance methods. '
                    'Use {"$array":[...],"dtype":"float64"} where a numpy array is required. '
                    'Full documentation: py4siesta_api_catalog(details=true, query="%s").' % entry['api'])
    return types.Tool(name=entry['tool'], description=description,
        inputSchema=entry['input_schema'],
        outputSchema={'type': 'object', 'properties': {'ok': {'type': 'boolean'}, 'result': {}}, 'required': ['ok']},
        annotations=types.ToolAnnotations(readOnlyHint=entry['read_only'],
            destructiveHint=not entry['read_only'], idempotentHint=entry['read_only'],
            openWorldHint=entry['category'] in ('execution', 'submission', 'visualization')))


async def serve(workdir, categories=None):
    from mcp import types
    from mcp.server.lowlevel import Server
    from mcp.server.lowlevel.helper_types import ReadResourceContents
    from mcp.server.stdio import stdio_server

    catalog = Catalog()
    server = Server('py4siesta', version='0.1.0')
    selected = [entry for entry in catalog.tools.values()
                if not categories or entry['category'] in categories]
    tools = support_tools() + [api_tool(entry) for entry in selected]
    allowed = {tool.name for tool in tools}
    worker = Worker(catalog, workdir)

    @server.list_tools()
    async def list_tools():
        return tools

    @server.call_tool()
    async def call_tool(name, arguments):
        if name not in allowed:
            payload = {'ok': False, 'error': {'type': 'ValueError', 'message': 'Tool is not enabled in this server.'}}
        else:
            try:
                if name == 'py4siesta_api_catalog':
                    payload = {'ok': True, 'result': catalog.search(**arguments)}
                else:
                    payload = await asyncio.to_thread(worker.call, name, arguments)
            except Exception as exc:
                payload = {'ok': False, 'error': {'type': type(exc).__name__, 'message': str(exc)}}
        return types.CallToolResult(content=[types.TextContent(type='text', text=json.dumps(payload, allow_nan=False))],
                                    structuredContent=payload, isError=not payload['ok'])

    @server.list_resources()
    async def list_resources():
        return [types.Resource(uri='py4siesta://api/catalog', name='Public API catalog', mimeType='application/json'),
                types.Resource(uri='py4siesta://runtime', name='Runtime environment', mimeType='application/json')]

    @server.read_resource()
    async def read_resource(uri):
        if str(uri) == 'py4siesta://api/catalog': payload = catalog.document()
        elif str(uri) == 'py4siesta://runtime':
            import nanocore
            import py4siesta
            payload = dict(python=sys.executable, python_version=sys.version,
                           nanocore=nanocore.__file__, py4siesta=py4siesta.__file__,
                           workdir=str(workdir), available_apis=len(catalog.tools),
                           enabled_apis=len(selected), unavailable_apis=len(catalog.entries)-len(catalog.tools))
        else: raise ValueError('Unknown resource: %s' % uri)
        return [ReadResourceContents(json.dumps(payload, allow_nan=False), mime_type='application/json')]

    try:
        async with stdio_server() as (read_stream, write_stream):
            await server.run(read_stream, write_stream, server.create_initialization_options())
    finally:
        worker.close()


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description='Expose existing nanocore and py4siesta APIs through MCP stdio.')
    parser.add_argument('--workdir', default='.', help='Existing default calculation directory.')
    parser.add_argument('--categories', nargs='+', help='Optionally expose only selected catalog categories; default: all.')
    parser.add_argument('--catalog', action='store_true', help='Print the full API inventory as JSON and exit.')
    args = parser.parse_args(argv)
    # Headless plotting is confined to this optional server process.
    os.environ.setdefault('MPLBACKEND', 'Agg')
    workdir = Path(args.workdir).expanduser().resolve()
    if not workdir.is_dir(): parser.error('--workdir must be an existing directory')
    if args.catalog:
        print(json.dumps(Catalog().document(), indent=2, allow_nan=False))
        return 0
    if sys.version_info < (3, 10): parser.error('MCP serving requires Python 3.10 or newer.')
    try:
        import mcp
    except ImportError:
        parser.error('Install the optional MCP dependencies: python -m pip install ".[mcp]"')
    asyncio.run(serve(workdir, args.categories))
    return 0
