# py4siesta MCP

This optional stdio server exposes the public functions, constructors, properties
and methods defined in NanoCore and py4siesta, plus every existing
`py4siesta-tool` subcommand. It calls their original implementations.

## Start and discover

Use Python 3.10 or newer for the MCP server:

```bash
python -m pip install '.[mcp]'
python -m py4siesta_mcp --workdir /absolute/calculation/directory
```

In an MCP client, configure the command as the absolute path to that Python
executable and arguments as `-m`, `py4siesta_mcp`, `--workdir`, and the calculation
directory. For an uninstalled checkout, include the repository in `PYTHONPATH`.
The installed `py4siesta-mcp` command is equivalent.

The server uses standard input/output for MCP messages. It does not accept an LLM
prompt, run a workflow automatically, or change existing menus and CLI commands.

Discover APIs without reading source code:

- `tools/list`: one tool per available API, with its signature, parameter schema,
  role, units, return conventions and effects.
- `py4siesta_api_catalog`: search/paginate by `query` and `category`; use
  `details=true` for the full original docstring and schema. Unavailable APIs
  remain visible with a reason.
- Resource `py4siesta://api/catalog`: full inventory.
- Resource `py4siesta://runtime`: actual Python/module paths and enabled counts.
- `python -m py4siesta_mcp --catalog`: export the same inventory as JSON.

All API categories are enabled by default. To limit the tools offered to a
client, use e.g. `--categories structure file-io post-processing`. This is a
discovery filter, not a security sandbox.

## Classification

| Category | Existing APIs |
| --- | --- |
| `structure` | Vector, Atom, AtomsSystem, Trajectory: geometry, selection, connectivity, transforms |
| `structure-builders` | carbonlab and surflab constructors |
| `file-io` | Structure/FDF/binary readers and writers, Fortran records, VASP files |
| `siesta` | Siesta configuration, input generation, result readers and XML objects |
| `case-workflows` | Calculation contexts, k-point/EOS/distance/interpolation preparation and analysis |
| `post-processing` | Band, PDOS, PLDOS and planar averages, including numerical and figure outputs |
| `cli-tools` | Existing deterministic `py4siesta-tool` commands and `execute()` |
| `submission` | Existing Slurm submission APIs |
| `execution` | Existing `Siesta.run()` |
| `visualization` | Existing XCrySDen helpers |
| `utilities` | Element lookup, file-copy and other public helpers |

Imported third-party functions, private helpers and Python implementation
details are not separate tools. Public inherited project methods are included.
Useful operators are named `operator_getitem`, `operator_multiply`, etc.
Interactive entry points, abstract hooks, context-manager-only operations and
unimportable legacy modules are recorded but not callable tools.

## Inputs and outputs

Each API tool takes:

```json
{
  "parameters": {"the_original_argument_name": "value"},
  "object_id": "obj_1",
  "workdir": "/absolute/calculation/directory"
}
```

`object_id` is required only for instance methods and properties. `parameters`
may be omitted when all original parameters are optional. `workdir` defaults
to the server's calculation directory. No directory is implicitly created.
Original `*args` and `**kwargs` use array/object fields with those parameter names.
Legacy unannotated parameters intentionally allow multiple JSON types; defaults
are shown without incorrectly making the default's type a restriction.

Constructors and object-returning readers return a handle such as:

```json
{"$object": "obj_1", "type": "nanocore.atoms.AtomsSystem"}
```

Pass that object to another API using `{"$object":"obj_1"}`. Handles are local
to the current server process and expire when the session/server ends. They do
not add fields or ownership references to NanoCore objects. Instance mutations
apply to that existing handle. A native `.copy()` still creates a separate object.

Small NumPy arrays are encoded as `{"$array":[...],"shape":[...],"dtype":"float64"}`;
the same marker can supply an ndarray input. Plain JSON lists remain lists.
Large arrays/lists and iterators receive handles: `py4siesta_object_read` pages
their values. Iterator reads consume items. Use the same tool's `attribute`
argument for public fields such as `BandStructureData.energies`.
`py4siesta_object_release` releases handles and closes file objects.

Successful calls return `{ "ok": true, "result": ... }`. Python errors or
failed CLI envelopes set `ok=false` and MCP `isError=true`. Native Python and
subprocess output is captured in `log`, limited to 64 KiB. Returned file paths
are interpreted relative to the returned `workdir` unless absolute.

## Example: connectivity

1. Call `nanocore_siestaio_read_fdf` with
   `{"parameters":{"file_name":"STRUCT.fdf"}}`; retain the returned object ID.
2. Call `nanocore_atoms_AtomsSystem_set_pbc` with that `object_id` and
   `{"parameters":{"pbc":[true,true,false]}}`.
3. Call `nanocore_atoms_AtomsSystem_calc_connectivity` with that `object_id`
   and `{"parameters":{"scale":1.1}}`.
4. Call `nanocore_atoms_AtomsSystem_operator_getitem` with the system ID
   and `{"parameters":{"i":0}}`; this returns the existing API's Atom copy.
5. Call `nanocore_atoms_Atom_get_connectivity` and
   `nanocore_atoms_Atom_get_iconnectivity` on the Atom ID, each with
   `{"parameters":{"details":true}}`.

Atom serials are ordinarily one-based; `operator_getitem` is zero-based.
Connectivity distances are in angstrom, and image shifts are integer lattice
translations. PBC is supplied by the user, never inferred by the MCP layer.

## Existing behavior and limitations

The adapter preserves file formats, scientific algorithms, selection semantics,
overwriting behavior and external executable requirements. MCP annotations are
descriptions, not permission enforcement. This is a local trusted-code interface
running with its launching user's filesystem and execution permissions.

Some original methods have incomplete docstrings, mixed accepted argument types,
or legacy runtime bugs. `available` means the API can be resolved and called,
not that every input or scientific result has been validated. The catalog includes
native documentation and explicitly marks raw-file units that are not converted.
Python-2 syntax errors in legacy modules are reported rather than repaired here.

Calls are serialized in a separate worker so process-wide cwd changes and
subprocess stdout cannot interfere with other calls or the MCP protocol. This
is not a background job manager. Submission may partially succeed before an
error; never automatically repeat it. A successful API return does not prove
SCF convergence, physical validity or scheduler completion.

The existing skill can discover and compose these tools. Workflow decisions,
user preferences and personal memory remain outside this server.
