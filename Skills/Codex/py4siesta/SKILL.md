---
name: py4siesta
description: Use when explicitly requested to perform py4siesta calculation work or register/edit py4siesta workflows. Solve tasks through connected py4siesta MCP tools in free or registered-workflow mode. Exclude general DFT advice, environment setup, and core/tool source development.
---

# py4siesta MCP workflows

Use the user's prepared py4siesta MCP environment to solve calculation tasks and
manage reusable workflows. Existing menu recipes are examples, not the boundary of
available operations. This skill does not require `py4siesta-agent`.

## Select a mode

Identify the goal, calculation directory, inputs, and scope (preparation, submission,
analysis, or complete calculation). Honor an explicitly requested mode. Otherwise:

| Request | Mode and instructions |
| --- | --- |
| Register, revise, or validate a reusable procedure | [Workflow registration/editing](modes/workflow-edit.md) |
| Calculation with a clearly matching compatible validated workflow | [Workflow execution](modes/workflow-run.md) |
| Other calculation work | [Free mode](modes/free.md) |

Announce the selected mode and read only its instructions and relevant references.
Use [state.md](references/state.md) to search relevant workflow metadata and personal
facts; do not load all history. Public candidates are in `workflows/`; the optional
[menu map](references/menu-map.md) indexes their original menu mappings.

When a registered procedure cannot handle a case, explain the gap before switching
to free mode. Proceed within existing authorization; ask only when a changed goal,
scientific choice, or action needs information or authorization not already supplied.
An explicit request to follow only the workflow forbids an automatic free-mode switch.
Free and execution modes do not register, revise, or promote workflows automatically;
enter editing mode only when requested. A combined execute-and-register request
already authorizes that transition.

## Common execution contract

- Before calculation calls, follow [mcp.md](references/mcp.md). A connected,
  compatible py4siesta MCP server and the required tools are prerequisites.
  The user owns installation, environment preparation, and MCP configuration.
  Do not create environments, install packages, reconfigure/restart servers, or
  bypass unavailable MCP tools with shell CLI calls or direct Python API calls.
- Non-calculation file inspection, artifact viewing, and workflow/state document
  edits may use host tools. They must not implement missing domain operations or
  substitute for unavailable calculation, submission, or monitoring MCP functions.
- Use only existing MCP operations for scientific/deterministic work. Glue may
  connect their inputs, outputs, branches, and checks; it must not implement or
  modify core/tool functionality. Report missing operations for separate development.
- Inspect inputs and existing outputs before mutations. Preparation does not
  authorize submission. Honor actual overwrite behavior and the user's scope.
  Never silently change scientific settings or repeat an ambiguous submission.
- Validate outputs against task criteria: `ok=true` does not establish calculation
  completion, SCF convergence, or physical validity. Retry only after a testable
  changed condition; stop the step if the same cause survives a targeted retry.
- Preserve necessary evidence and resume state under [state.md](references/state.md).
  Inspect actual outputs/status on resume; waiting jobs remain `waiting`.

Before finishing, apply the independent save conditions in state.md and read back
any changed records. Report mode, workflow/revision if used, actual results and
paths, validation, pending work, and records only when saved. Disclose necessary
persistence failures. There is no background learning or monitoring after the session.
