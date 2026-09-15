# Register the py4siesta Codex skill

Run from the repository root in your terminal (Linux/macOS):

```bash
bash Skills/Codex/register.sh --dry-run
bash Skills/Codex/register.sh
```

Requires Python 3.11+, or an older Python 3 with `tomli` or `toml` installed.
Select an interpreter with `PYTHON=/path/to/python bash Skills/Codex/register.sh`.
Registration does not install py4siesta or its calculation dependencies.

The script links `~/.agents/skills/py4siesta` to this checkout, creates
`~/.local/state/py4siesta-skill`, and adds that absolute state path to
`[sandbox_workspace_write].writable_roots` in `${CODEX_HOME:-~/.codex}/config.toml`.
Existing settings and writable paths are retained; changed configuration is backed
up beside the original. Repeated registration is safe. An existing skill at another
location is not replaced. Keep this checkout in place while using the link.

For another personal state location, export `PY4SIESTA_SKILL_STATE` before running
the script and keep it set when starting Codex. The location must be outside the
repository and installed skill tree. No workflow or history records are created
during registration.

Start a new Codex session, check `/skills`, then invoke `$py4siesta`. Additional
writable roots apply to `workspace-write`; the script does not change your sandbox
mode. If the host manages session permissions separately, it must also allow the
state path. Registration alone cannot override that policy. When a task needs persistent state, verify that the skill reports a saved record
in this directory.

## Selective memory

Routine unchanged success creates no new record. The skill saves meaningful
workflow improvements, new preferences/environment facts, and short observations
of actual friction only when useful. It retains minimal history for resumption,
uncertain submissions, necessary validation, or explicitly requested auditing.
Existing calculation logs are referenced rather than copied.

Reusable procedures belong in workflows; observations describe what happened,
including inefficient successful work, without proposing implementation changes.
Only relevant state is read. Prior revisions and legacy lessons are preserved,
but are not loaded routinely. Registration and sandbox setup are unchanged.
