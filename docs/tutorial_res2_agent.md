# ReS2 monolayer agent tutorial

This tutorial runs the complete ReS2 workflow with one command. It performs
real database access, submits real Slurm jobs, and runs the configured SIESTA
executable.

## 1. Configure an LLM

Choose one provider. For OpenAI:

```bash
python -m pip install -e '.[openai]'

export PY4SIESTA_LLM_PROVIDER="openai"
export PY4SIESTA_LLM_MODEL="<model-name>"
export OPENAI_API_KEY="<your-api-key>"
```

For Anthropic, install `.[anthropic]` and set `ANTHROPIC_API_KEY`. For a local
Ollama model, install `.[ollama]`; no remote API key is normally required.
Never put an API key in the repository or a scheduler script.

## 2. Configure the calculation

The following paths match the current machine:

```bash
export PY4SIESTA_WORKDIR="$HOME/py4siesta-runs/ReS2-LDA"
mkdir -p "$PY4SIESTA_WORKDIR"

export PY4SIESTA_INITIAL_KPOINTS="2 2 1"
export PY4SIESTA_KPOINT_SERIES="2 4 6 8"

export PY4SIESTA_SLURM_TEMPLATE="$PWD/test/agent_res2/slm_siesta_run"
export PY4SIESTA_MAX_TOTAL_NODES="10"
export PY4SIESTA_POLL_INTERVAL_SECONDS="60"
export PY4SIESTA_SIESTA_BIN="/home2/rong/bin/siesta_4.1b4"
export PY4SIESTA_MPIRUN="/opt/intel/oneapi/mpi/2021.12/bin/mpirun"
```

`NanoCore/env.py` already points to `/home2/rong/bin/psf`. The required LDA
pseudopotentials are therefore found automatically:

```text
/home2/rong/bin/psf/LDA/Re.psf
/home2/rong/bin/psf/LDA/S.psf
```

The supplied scheduler file, `test/agent_res2/slm_siesta_run`, requests one X2
node and 12 MPI tasks. Change its Slurm directives if those resources are not
appropriate for the cluster.

Use a new calculation directory. Existing `origin/` and calculation-stage
directories are not overwritten.

## 3. Run

From the repository root:

```bash
python py4siesta-agent \
  "can you setup the dft calculation with lda for ReS2 monolayer"
```

Keep this command running. It automatically:

```text
gets the ReS2 monolayer structure from C2DB
→ creates origin/
→ submits and monitors k-point jobs
→ analyzes k-point convergence
→ submits and monitors lattice-optimization jobs
→ submits and monitors geometry optimization
→ creates final_input/
```

No separate preflight, status, or continue command is required.

The C2DB request is necessary because the command names a material but does
not provide atomic coordinates. It occurs once when the real workflow selects
the starting structure; the tutorial does not make an additional validation
request. C2DB requires no API key or ASE installation.

To avoid database access, provide an existing local FDF structure path in the
request instead.

## 4. Results

The calculation directory contains:

```text
origin/
01.kpoint_sampling/
02.slab_eos/
03.geometry_optimization/
final_input/
.py4siesta-agent/
```

Raw SIESTA output remains in each case’s `OUT/` directory. The final Markdown
summary is stored in `.py4siesta-agent/`.

This workflow automates execution, but the resulting structure,
pseudopotentials, basis, convergence range, and optimization criteria still
require scientific review.
