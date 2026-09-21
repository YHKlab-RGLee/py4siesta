# Free mode

Use this mode to solve the requested problem by composing connected MCP tools.

1. Establish the inputs, scientific choices, allowed side effects, and completion
   criteria. Read the common [MCP contract](../references/mcp.md).
2. Search `py4siesta_api_catalog` for relevant operations and inspect their schemas,
   units, outputs, and effects. Confirm they are enabled in the connected tool list.
   A compatible workflow may inform the approach but does not constrain the sequence.
3. Form a task-local sequence connecting verified inputs and outputs. Execute and
   inspect results before deciding the next call; adapt within the authorized goal.
4. Stop a blocked step when a required tool, environment capability, input, or
   validation is missing. State the attempted operation and observed limitation.
   Do not install dependencies, bypass MCP, or implement the missing operation.
5. Report verified results and unfinished work. Save only qualifying preferences,
   environment facts, observations, or necessary resume evidence under
   [state.md](../references/state.md).

Do not save or improve a reusable workflow just because the task succeeded. If the
user also requested registration, transition to [editing](workflow-edit.md) using
the actual calls and evidence; otherwise leave the procedure task-local.
