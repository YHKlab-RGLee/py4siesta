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

Start a new Codex session, check `/skills`, then invoke `$py4siesta` explicitly.
Automatic invocation is disabled; SIESTA-related requests alone do not activate
the skill. Additional
writable roots apply to `workspace-write`; the script does not change your sandbox
mode. If the host manages session permissions separately, it must also allow the
state path. Registration alone cannot override that policy. When a task needs persistent state, verify that the skill reports a saved record
in this directory.

## MCP prerequisites

Follow the [Conda setup, Codex MCP registration, and connection check](../../README.MD#register-once-for-all-projects)
for the complete setup sequence. Skill registration and MCP registration are separate steps.

Calculation execution requires a connected, compatible py4siesta MCP server. Prepare
Python 3.10+ with py4siesta and its MCP dependencies, configure the server in your
client, and provide the external programs required by your task (for example SIESTA
utilities or Slurm). A dedicated Python environment is recommended, but a working
existing environment is acceptable. See the [MCP setup guide](../../py4siesta_mcp/README.md)
for installation and connection details. Start the server with the intended Python's
absolute path; activating another shell environment does not change a running server.

The skill checks the connection and task-specific prerequisites. It does not create
environments, install packages, change MCP settings, or fall back to CLI/direct Python
execution. If a prerequisite is missing, it reports the missing item and stops the
dependent work. Skill registration above does not install or connect the MCP server.

## Three modes

| Mode | Usage |
| --- | --- |
| Free | `$py4siesta 자유 모드로 이 구조의 연결성을 분석해줘.` — discover and compose available MCP tools for the goal. |
| Workflow execution | `$py4siesta 등록된 workflow로 k-point 수렴을 분석해줘.` — follow a compatible procedure and its permitted branches. |
| Workflow registration/editing | `$py4siesta 방금 수행한 절차를 개인 workflow로 등록해줘.` — write or revise a reusable MCP procedure. |

An explicit mode takes precedence. Otherwise the skill selects execution mode when
a validated workflow clearly matches, and free mode for other calculation requests;
registration/editing requests select editing mode. The skill announces its mode.
Execution does not automatically register or change a workflow. A draft trial can
be requested explicitly, but it is not automatically treated as validated.

Public candidates live in [py4siesta/workflows/](py4siesta/workflows/). The five migrated
MCP recipes are drafts: tool names/argument schemas were reviewed, but their MCP
calculation paths have not been execution-validated. Historical CLI recipes remain
as migration evidence and are not current execution instructions. Editing can produce
an unverified draft without MCP; actual calculation validation requires it. A request
to register a procedure alone does not submit jobs or run calculations.

## Selective memory

Personal workflows, preferences, environment facts, observations, and necessary resume
records stay outside the repository, by default in `~/.local/state/py4siesta-skill`.
Workflow registration defaults to personal storage; request public/repository
registration explicitly to add a distributable procedure to `workflows/`.

Only requested registration/editing changes workflow revisions. Other modes save
useful new preferences/environment facts, short nonduplicate observations of actual
friction, and minimal history when needed for safe resumption, uncertain submissions,
validation, or requested auditing. Existing calculation logs are referenced, not copied.
Unchanged success creates no new record. Prior revisions and legacy evidence are
preserved and read only when relevant. Changed state files are read back to verify
persistence. No background learning or monitoring continues after the session ends.
