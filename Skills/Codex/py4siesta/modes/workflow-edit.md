# Workflow registration and editing mode

Use this mode only for a request to register, revise, or validate a workflow.
A combined calculation-and-registration request covers both activities.

1. Determine the goal and scope, existing workflow ID if any, and intended storage.
   Default to personal `<state>/workflows/<id>/<revision>.md`; write public
   `workflows/<id>/<revision>.md` only when public/repository registration is requested.
   Verify the writable checkout; never infer it from a copied installed skill.
   Follow [state.md](../references/state.md), preserving prior revisions and evidence.
2. Use the [workflow template](../assets/workflow.md). Start from user instructions,
   a selected prior procedure, or actual free-mode calls and results. Generalize
   variable inputs only within supported scientific and environment conditions.
3. Resolve each required MCP tool, schema, units, outputs, and effects through the
   connected catalog/tool list. Define explicit result bindings, branches, stopping
   conditions, per-step validation, output locations, and recovery behavior.
   A missing operation is a blocked step, not permission to invent a tool or algorithm.
4. Without MCP, author a `draft` and label tool contracts and execution as unverified.
   Separate schema/document review from execution validation. Registration alone
   does not authorize running calculations, submitting jobs, or overwriting results.
   When execution validation is requested, apply the common MCP prerequisites and
   trial the declared procedure within that scope.
5. Save a first draft or meaningful new revision. Promote only when evidence meets
   its declared completion criteria; preparation evidence validates preparation only.
   Link existing outputs/logs, note remaining untested branches, and limit validated
   applicability accordingly. Read back the saved document and report its path/status.

Local Python paths and current server object IDs belong outside portable procedures.
Use symbolic step-result bindings and reconstruct objects from persistent inputs on
resume. Modifying core/tool code, installing dependencies, or changing the MCP server
is outside this mode. Do not update unrelated workflows or the skill's mode rules.
