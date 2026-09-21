> Historical CLI recipe, preserved as migration evidence. Do not execute its
> commands under the current skill. Use the [MCP draft](../../workflows/submission/r2.md)
> and the current mode and MCP rules instead.

# Job submission and resume — revision 1

Menus: `9` (k-points), `10` (optimization). Read runtime.md before execution.

Use only when submission is included in the user's request. Inspect origin, scripts,
resources, and actual target list first. `submit --mode kpt` chooses the FIRST sorted
`01.*` directory; `submit --mode opt` chooses the FIRST sorted `02.*` directory. It
then calls sbatch for EVERY `slm_*` in every immediate child directory. It does not
filter completed jobs, choose a requested scan, or deduplicate submissions. An
`optimized_structure` child containing a copied script can therefore also be submitted.

If that exact target set does not match the request, do not run the batch command.
Explain the mismatch and use an explicitly scoped submission method only if authorized.
No matching directory can return successfully without submitting anything.

```bash
py4siesta-tool submit --mode kpt
# Or, after verifying the optimization target set:
py4siesta-tool submit --mode opt
```

Capture each sbatch response, job ID, script, and case path immediately in history.
On partial failure, reconcile successful IDs and scheduler status before any retry;
never rerun the full command blindly. If IDs were lost, inspect scheduler/accounting
records and script logs. Unknown submission state must not trigger automatic resubmission.

Submission success means submitted, not calculated. Use the environment's available
scheduler monitoring mechanism when requested. Track pending/running/failed/completed
states and inspect SIESTA convergence separately. Save the next step and IDs when the
session stops. Resume from these records without regenerating cases or repeating jobs.
