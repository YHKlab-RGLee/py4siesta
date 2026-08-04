# DFT setup agent implementation

## Scope

The first `py4siesta-agent` workflow routes a natural-language DFT setup
request to `DFTSetupAgent`. It prepares a structure, creates `origin/`, runs the
existing k-point and lattice workflows, prepares a conjugate-gradient geometry
calculation, and generates `final_input/`. Workflow state is checkpointed
throughout the run.

The agent does not recover failed calculations or change scientific inputs
automatically.

## Package structure and boundaries

- `py4siesta_agent/shared/model_client.py` is the only model initialization
  location. It returns a provider-neutral `ChatModelClient`.
- `py4siesta_agent/shared/types.py` contains validated routing, request,
  structure-selection, job, and workflow-state values.
- `py4siesta_agent/router/agent_router.py` selects from declared agent
  capabilities. It contains no DFT operations.
- `py4siesta_agent/agents/dft_setup/dft_setup_agent.py` defines and executes
  the LangGraph state graph.
- `py4siesta_agent/agents/dft_setup/workflow_operations.py` translates agent
  configuration into calls to existing deterministic scientific operations.
- `py4siesta_agent/scheduler.py` owns the scheduler backend, managed-job
  accounting, pending queues, and the per-workflow node limit.
- `py4siesta_agent/structure_sources.py` owns the database integration
  introduced for structure selection by this agent.
- `py4siesta/operations.py` exposes the existing deterministic scientific
  operations reused by the agent without depending on agent packages.

The model is used only for routing, request parsing, and selection among a
small recorded candidate set. It cannot submit commands, parse outputs,
calculate convergence, change resource limits, or invent scientific data.

## Model configuration

Choose one supported LangChain chat-model integration:

```bash
export PY4SIESTA_LLM_PROVIDER="openai"
export PY4SIESTA_LLM_MODEL="<model-name>"
export OPENAI_API_KEY="<credential>"
```

Supported provider values are `openai`, `anthropic`, and `ollama`. Install the
matching optional extra:

```bash
pip install '.[openai]'
pip install '.[anthropic]'
pip install '.[ollama]'
```

Anthropic uses `ANTHROPIC_API_KEY`. Ollama uses the endpoint conventions of
`langchain-ollama` and normally needs no remote API credential. Credentials
are never copied into calculation inputs, scheduler scripts, summaries, or
LangGraph state.

`create_model_client()` validates the provider, model, integration package,
and required credential. `ChatModelClient.structured()` first uses the common
model `with_structured_output()` interface and validates a Pydantic schema.
Models without native structured output may use a validated plain-JSON
fallback. A model that supports neither path fails before deterministic
workflow execution.

## Structure selection

A user-supplied local FDF path is used directly. A 2D or monolayer request
queries the public C2DB OPTIMADE endpoint directly:

```text
https://c2db.fysik.dtu.dk/optimade/v1/structures
```

No ASE installation, local C2DB file, or C2DB credential is required.

A 3D bulk request uses the official Materials Project client:

```bash
pip install '.[materials-project]'
export MP_API_KEY="<credential>"
```

Database queries and exact formula, dimensionality, structure-type, and
available-metadata filtering are deterministic. Candidate identifiers and
metadata are recorded. If more than one filtered candidate remains, the model
may choose only one of those identifiers.

## Scientific configuration and reused operations

The workflow does not invent k-point defaults. Configure the initial mesh and
convergence cases:

```bash
export PY4SIESTA_INITIAL_KPOINTS="2 2 1"
export PY4SIESTA_KPOINT_SERIES="2 4 6 8"
```

`DFTWorkflowOperations.initialize_origin()` calls the existing
`py4siesta.operations.initialize_origin()` Python entry point used by
`py4siesta-tool init`; it never invokes a model-generated shell command.
K-point generation and deterministic convergence analysis reuse
`SiestaWorkflow.kpoint_sampling()` and `kpoint_analysis()`. Lattice cases and
curve analysis reuse `eos_slab()`/`eos_bulk()` and
`find_optimized_lattice()`.

The optimized lattice case becomes `03.geometry_optimization/`. Existing
NanoCore SIESTA inputs already configure conjugate-gradient optimization.
Completion requires both `OUT/0_NORMAL_EXIT` and a parseable
`OUT/*.STRUCT_OUT`. The latter replaces `final_input/input/STRUCT.fdf`.

## Scheduler configuration and resource control

The repository uses Slurm directives and `sbatch`, so the initial backend is
Slurm:

```bash
export PY4SIESTA_SLURM_TEMPLATE="/path/to/slm_siesta_run"
export PY4SIESTA_SCHEDULER_BACKEND="slurm"
export PY4SIESTA_MAX_TOTAL_NODES="10"
```

An explicit scheduler path parsed from the request takes precedence over
`PY4SIESTA_SLURM_TEMPLATE`. The script must be readable, contain Slurm
directives, and include `#SBATCH --nodes=<count>` or `#SBATCH -N <count>`.

`SchedulerManager` constructs only fixed `sbatch`, `squeue`, `sacct`, and
`scancel` argument lists. It never accepts arbitrary shell strings. Before
submission it sums nodes for all managed queued and running jobs. Jobs that
would exceed the default aggregate limit of ten remain `pending`. Released
capacity is filled in recorded workflow order.

## LangGraph workflow and persistence

The explicit graph stages are:

```text
parse_request → select_structure → initialize_origin
→ prepare/submit/monitor/analyze k-points
→ prepare/submit/monitor/analyze lattice cases
→ prepare/submit/monitor/validate CG geometry
→ generate_final_input → write_summary
```

Monitor nodes have explicit transitions for waiting, completion, and failure.
They do not call the model. A failed or cancelled job advances only to the
summary and never to numerical analysis.

`DFTSetupAgent.run()` remains active when a monitor node returns `waiting`.
After `PY4SIESTA_POLL_INTERVAL_SECONDS` (60 seconds by default), it reloads the
checkpoint and invokes the graph again. This repeats until the graph records
`completed` or `failed`.

SQLite checkpoints are stored under:

```text
<workdir>/.py4siesta-agent/checkpoints.sqlite
```

Start a workflow:

```bash
python py4siesta-agent \
  "can you setup the dft calculation with lda for ReS2 monolayer"
```

Use `--workdir <path>` or `PY4SIESTA_WORKDIR` when the calculation root differs
from the current directory. The initial command polls recorded jobs, releases
node capacity, submits pending jobs, runs deterministic analysis, and advances
the graph without another user command.

The final Markdown report is
`<workdir>/.py4siesta-agent/<workflow-id>.md`. It is generated solely from
recorded state, including incomplete stages and errors.

## Testing

Run the focused mocked suite:

```bash
python -m unittest discover -s test -p 'test_agent_workflow.py' -v
```

The tests fake the common model boundary, structure source, scheduler backend,
and scientific outputs. They cover unsupported routing, invalid model output,
missing model configuration, scheduler validation, aggregate node accounting,
pending jobs, autonomous one-command continuation, scheduler failure,
incomplete geometry output, and rejection of model-provided command strings.
No live model, database, scheduler, SIESTA executable, or DFT calculation is
required.
