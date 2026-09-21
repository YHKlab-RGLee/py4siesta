# Workflow execution mode

1. Select the requested workflow/revision, or a clearly matching latest compatible
   validated revision. Read the entire selected procedure, its linked operational
   constraints, inputs, branches, and completion criteria. Do not routinely read
   older revisions. See [state.md](../references/state.md).
2. Check the connected runtime and required tools under [mcp.md](../references/mcp.md).
   Bind task inputs to declared parameters; resolve scientific choices from the
   current request rather than illustrative recipe values or another material's defaults.
3. Execute the defined MCP sequence and permitted branches. Carry outputs forward
   using named step results, with live session handles resolved at execution time.
   Validate each declared checkpoint before continuing.
4. If tools or environment requirements are unavailable, stop the affected step and
   report the missing prerequisite. If the procedure itself does not fit, explain
   the gap before a permitted free-mode transition. Never silently substitute tools,
   add branches, alter scientific parameters, or rewrite the workflow.
5. Check the final criteria and report the exact revision, results, and pending work.
   Save resume information only when needed; routine repeats create no revision.

Drafts are not automatic execution candidates. If the user explicitly requests a
trial of a draft, disclose its unvalidated scope and execute only the authorized
trial. Registering validation evidence or promoting it is editing-mode work and
requires a registration/edit/validation request; execution alone does not promote it.
Legacy CLI/API-only personal recipes need an explicitly requested MCP adaptation
before registered execution. They may be relevant evidence in free mode, never a
reason to bypass MCP.
