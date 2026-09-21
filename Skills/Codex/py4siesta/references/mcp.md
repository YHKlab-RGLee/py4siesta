# MCP prerequisites and calling contract

## User-prepared environment

Before scientific execution, confirm the connected py4siesta MCP server exposes the
required tools. Read `py4siesta://runtime` for its actual Python, module locations,
default workdir, and enabled counts when resource access is available. If unavailable,
use equivalent runtime evidence supplied by the user; do not claim verification of
facts that could not be checked. A catalog entry being `available` does not prove
it is enabled in this connection or valid for the requested inputs.

Search `py4siesta_api_catalog` by task/operation/category with `details=true` only
for relevant candidates. Match returned tool names to the connected tool list
(client-added namespace prefixes may differ). Inspect parameters, units, return
conventions, side effects, and missing external programs. Check only prerequisites
needed by this task, not every API or dependency.

The user prepares a compatible Python environment, py4siesta with MCP dependencies,
MCP connection, and required SIESTA utilities/scheduler. A dedicated environment is
recommended, not mandatory. Do not create environments, install packages, change
PATH/PYTHONPATH, edit MCP settings, or restart the server as part of this skill.
An already-running server does not change Python when a shell activates a venv.
If a prerequisite fails, report the specific missing item and refer to the
repository's Skills/Codex/README.md setup guidance. Stop dependent execution;
continue independent document work if useful. Never fall back to CLI/direct Python.

## Calls and result bindings

A calculation API tool receives `parameters`, optional `object_id` for an instance,
and `workdir`. Use an explicit absolute calculation workdir (already existing), never
the installed skill directory. Relative returned paths resolve against the returned
workdir. Input/output inspection and state-file edits can use host filesystem tools;
calculation, submission, and scheduler-query operations must use available MCP tools.
If scheduler querying is unavailable, retain IDs and request external status evidence;
do not invent a scheduler tool or switch to shell submission/monitoring.

Example API mapping: `py4siesta-tool.kpoint-bulk` is exposed by this checkout as
`py4siesta_tool_kpoint_bulk`, accepting `parameters: {"kpoints": [2, 4, 6]}`.
These samples illustrate schema, not recommended scientific settings. Confirm the
connected schema before calling. PDOS's repeatable orbital argument is nested:
`parameters: {"orbital": [["Mg_0", "O_0"]]}` for this implementation.

Constructors/readers may return `{"$object":"obj_N"}`. Pass object references in
parameters, or use the live ID in `object_id` for instance methods/properties.
Workflow documents bind logical results (e.g. `read_structure.result`) rather than
persisting concrete IDs. IDs expire on server restart; resume from durable files and
reconstruct objects, without repeating completed mutations or submissions.
Use `py4siesta_object_read` for public fields or paging; iterator reads consume items.
Use `py4siesta_object_release` when finished. `$array` carries ndarray values;
plain lists remain lists. Check native units and indexing conventions per API.

Check MCP `isError`, structured `ok`, nested CLI envelopes, and captured `log`.
Logs are limited to 64 KiB; absence from the returned log is not proof an action did
not happen. Retain necessary evidence from actual calculation outputs and job IDs.
The worker and workdir do not enforce filesystem isolation, and categories are a
selection filter, not a sandbox.

## File and calculation constraints

The migrated bundled workflows wrap existing CLI operations through MCP. Their
origin-dependent steps initialize from `workdir/origin/input/STRUCT.fdf`; band,
PDOS, and PLDOS do not require origin. Apply each API's own contract to other tools.
Case generation copies origin, including scientific inputs and job scripts. A
readable STRUCT is sufficient only for structure generation, not a runnable case.

BaseOperation may delete/recreate its output tree before validating all inputs.
Validate paths and protect existing results before calling. Keep source files outside
recreated output trees. Preparation tools do not themselves complete a SIESTA run.
Use each workflow's file checks and scientific completion criteria, not merely `ok`.
Package version alone may not distinguish source revisions: record relevant verified
module/version/source evidence when needed for compatibility; do not invent a commit
when the runtime does not expose one. Shell inspection of known source metadata is
permitted, but importing/executing calculation APIs outside MCP is not.
