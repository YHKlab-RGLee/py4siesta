# py4siesta

`py4siesta` provides interactive, deterministic, and AI-orchestrated command-line interfaces for preparing, organizing, submitting, and post-processing SIESTA calculation cases.

## Table of Contents

- [Important Notice](#important-notice)
- [Quick Start](#quick-start)
- [Requirements](#requirements)
- [Configure SIESTA Paths](#configure-siesta-paths)
- [Working Directory Setup](#working-directory-setup)
- [Menu Reference](#menu-reference)
  - [K-point Sampling](#k-point-sampling)
  - [Structure Optimization](#structure-optimization)
  - [Job Submission](#job-submission)
  - [Post-Process](#post-process)
  - [Utility](#utility)
- [Running Modes](#running-modes)
- [Direct JSON Tools](#direct-json-tools)
- [AI-Agent Interface](#ai-agent-interface)
  - [Install Agent Support](#install-agent-support)
  - [Configure the Agent](#configure-the-agent)
  - [Run the Agent](#run-the-agent)
  - [Agent Outputs](#agent-outputs)
- [Bundled NanoCore](#bundled-nanocore)

The interactive and deterministic interfaces are designed around this basic
workflow:

1. Prepare a reference SIESTA input in an `origin/` directory.
2. Start the CLI.
3. Choose a numbered menu option.
4. Let `py4siesta` generate or analyze the calculation cases.

## Important Notice

This project includes the bundled `NanoCore` code developed by Dr. Hu Sung Kim, KISTI, and related original work from **OORINano / YHKLab**.

- Keep the original attribution and copyright context for `NanoCore` and OORINano.
- The local bundled `NanoCore/` package is used directly from this repository.
- You do not need to download `NanoCore` from an external source.

## Quick Start

Clone and enter the repository:

```bash
git clone https://github.com/YHKlab-RGLee/py4siesta.git
cd py4siesta
```

Install the package:

```bash
pip install .
```

For development, use editable install instead:

```bash
pip install -e .
```

Run the interactive numbered-menu interface:

```bash
py4siesta
```

Or run it directly from the repository root:

```bash
python -m py4siesta
```

The project has three separate user interfaces:

| Interface | Purpose | Status |
| --- | --- | --- |
| `python py4siesta` | Interactive GUI/menu-based workflow. | Available |
| `py4siesta-tool <command> [options]` | Direct, deterministic CLI tool invocation with JSON output. | Available |
| `py4siesta-agent "<request>"` | Natural-language AI-agent orchestration over existing tools. | Initial DFT setup workflow available |

After installation, the equivalent entry points are `py4siesta`, `py4siesta-tool`, and `py4siesta-agent`.

## Requirements

Core runtime dependencies are declared in `pyproject.toml`:

- Python 3.7 or newer for the non-agent interfaces
- `numpy`
- `scipy`
- `matplotlib`

The repository also installs the bundled `NanoCore` package.

Agent frameworks, model-provider integrations, and the Materials Project
client are optional extras. The agent workflow requires Python 3.10 or newer.
C2DB monolayer structures use its public OPTIMADE endpoint and require no
additional Python package.

## Configure SIESTA Paths

The bundled `NanoCore` module reads local executable paths from `NanoCore/env.py`. Update this file before using features that call SIESTA utilities.

Common values to edit:

```python
siesta_calculator = "/path/to/siesta"
siesta_psf_location = "/path/to/psf"
siesta_util_location = "/path/to/siesta/utils"
siesta_util_band = "gnubands"
siesta_util_dos = "Eig2DOS"
siesta_util_pdos = "fmpdos"
siesta_util_rho = "rho2xsf"
siesta_util_vh = "macroave"
siesta_util_pldos = "/path/to/PyProjection.py"
```

For PDOS extraction, make sure `siesta_util_location` points to the directory containing `fmpdos`, or set `siesta_util_pdos` to the full executable path.

## Working Directory Setup

Create a complete `origin/` directory non-interactively from a SIESTA FDF
structure, exchange-correlation functional, initial k-point mesh, and SLURM
script:

```bash
py4siesta-tool init \
    --structure C_diamond.fdf \
    --xc lda \
    --kpt 2 2 2 \
    --slurm ./slm_siesta_run
```

The command obtains the required pseudopotentials from the database configured
by `siesta_psf_location` in `NanoCore/env.py`. It refuses to overwrite an
existing `origin/` directory.

Expected layout:

```text
origin/
├── slm_siesta_run
└── input/
    ├── RUN.fdf
    ├── STRUCT.fdf
    ├── BASIS.fdf
    ├── KPT.fdf
    └── *.psf
```

You may also prepare this layout manually. Required inputs:

| Path | Purpose |
| --- | --- |
| `origin/input/STRUCT.fdf` | Reference SIESTA structure/input file copied into generated cases. |
| `origin/slm_*` | One or more job scripts copied into generated cases. |

## Menu Reference

Start the program, then select a menu number:

```bash
py4siesta
```

### K-point Sampling

| Menu | Feature | Purpose |
| --- | --- | --- |
| `1` | Bulk | Generate bulk cases using isotropic `k x k x k` Monkhorst-Pack meshes. |
| `2` | Slab | Generate slab cases using in-plane `k x k x 1` meshes. |
| `3` | Analysis | Run k-point convergence analysis from completed k-point calculations. |

### Structure Optimization

| Menu | Feature | Purpose |
| --- | --- | --- |
| `4` | Bulk | Create bulk equation-of-state structures and choose scaled directions (`x y z`). |
| `5` | Slab | Create slab equation-of-state structures with in-plane-focused scaling. |
| `6` | Sliding | Translate selected atoms laterally to scan stacking or interface registry. |
| `7` | Distance | Change the vertical distance between selected atoms and the remaining structure. |
| `8` | Analysis | Analyze completed EOS or distance calculations and generate an optimized structure. |

### Job Submission

| Menu | Feature | Purpose |
| --- | --- | --- |
| `9` | K-point Sampling | Submit generated k-point sampling jobs. |
| `10` | Structure Optimization | Submit generated structure-optimization jobs. |

### Post-Process

| Menu | Feature | Inputs | Outputs |
| --- | --- | --- | --- |
| `11` | Plot Band Structure | A `.bands` file and energy window. | `band.png`, plus band-path helper files. |
| `12` | Generate PDOS CSV | SIESTA `.PDOS`, matching `.EIG` and `.DOS`, energy window, and `fmpdos` atom/orbital selections. | `PDOS.csv` and `pdos.png`. |
| `13` | Plot PLDOS | A `.PDOS` file, matching `.xyz` and `.EIG`, and energy window. | `pldos.png`, plus PLDOS helper files. |

Notes:

- Menu `11` automatically selects a local `.bands` file when one is present.
- Menu `12` uses a matching `<label>.bands` file for `E_v` when available; otherwise it estimates `E_v` from `<label>.EIG`.
- Menu `12` removes existing same-named `fmpdos` output files before regenerating PDOS data.
- Menu `13` automatically selects a local `.PDOS` file when one is present.

### Utility

| Menu | Feature | Purpose |
| --- | --- | --- |
| `01` | Generate Geometries | Open the geometry-generation submenu. |

Geometry-generation submenu:

| Submenu | Feature | Purpose |
| --- | --- | --- |
| `1` | Move structure | Translate the full structure by a user-defined displacement vector and write a new `STRUCT.fdf`. |
| `2` | Interpolate structure | Generate interpolated structures between two `STRUCT.fdf` files, with optional extrapolated points. |

## Running Modes

Use one of these commands depending on how you installed or checked out the project:

```bash
# From an installed package
py4siesta

# From the repository root
python -m py4siesta

# Direct script usage
python {py4siesta_path}/run.py
```

## Direct JSON Tools

For scripts and other callers that need direct non-interactive access to deterministic operations, use:

```bash
py4siesta-tool <command> [options]
```

Without installation, the same interface is available from the repository root as:

```bash
python -m py4siesta.tool_cli <command> [options]
```

The command prints one JSON object. Successful commands include `"ok": true`; failures return a nonzero exit status and a JSON error object.

Examples:

```bash
py4siesta-tool kpoint-bulk --kpoints 2 4 6
py4siesta-tool eos-sliding --selection 20-30 --mode fractional --vector 0.25 0.50
py4siesta-tool band --bands-path Solid.bands --emin -2 --emax 4
py4siesta-tool pdos --pdos-path MgO.PDOS --orbital Mg_0 O_0 --emin -4 --emax 12
py4siesta-tool pldos --pdos-path Device.PDOS --emin -4 --emax 2
```

For `py4siesta-tool pdos`, pass one or more `--orbital` selections in `atom_or_species[_n[_l[_m]]]` format, for example `Mg_0`, `O_0`, or `C_2_1_0`. Multiple selections may be space-separated, comma-separated, or provided with repeated `--orbital` options.

Run `py4siesta-tool --help` or `py4siesta-tool <command> --help` for the full command list and required inputs.

`py4siesta-tool` is a deterministic CLI, not an AI agent. It performs only the command and options explicitly supplied by the caller.

## AI-Agent Interface

The agent interface currently supports preparing a complete SIESTA DFT setup,
including k-point convergence, lattice optimization, geometry optimization,
and final-input generation.

### Install Agent Support

Install one model provider and any optional bulk-structure client you need:

```bash
# OpenAI
pip install '.[openai]'

# Anthropic with Materials Project
pip install '.[anthropic,materials-project]'

# Local Ollama model
pip install '.[ollama]'
```

### Configure the Agent

Set the model provider, model name, k-points, scheduler script, and calculation
directory.

Example for OpenAI and a 2D material from online C2DB:

```bash
export PY4SIESTA_LLM_PROVIDER="openai"
export PY4SIESTA_LLM_MODEL="<model-name>"
export OPENAI_API_KEY="<your-api-key>"

export PY4SIESTA_INITIAL_KPOINTS="2 2 1"
export PY4SIESTA_KPOINT_SERIES="2 4 6 8"

export PY4SIESTA_SLURM_TEMPLATE="/absolute/path/to/slm_siesta_run"
export PY4SIESTA_MAX_TOTAL_NODES="10"
export PY4SIESTA_POLL_INTERVAL_SECONDS="60"
export PY4SIESTA_WORKDIR="/absolute/path/to/calculation"
```

For Anthropic, use:

```bash
export PY4SIESTA_LLM_PROVIDER="anthropic"
export PY4SIESTA_LLM_MODEL="<model-name>"
export ANTHROPIC_API_KEY="<your-api-key>"
```

For a local Ollama model, use:

```bash
export PY4SIESTA_LLM_PROVIDER="ollama"
export PY4SIESTA_LLM_MODEL="<local-model-name>"
```

For a bulk structure from Materials Project, set:

```bash
export MP_API_KEY="<your-materials-project-api-key>"
```

Also configure `siesta_psf_location` in `NanoCore/env.py` as described in
[Configure SIESTA Paths](#configure-siesta-paths).

The scheduler script must be readable and contain either:

```bash
#SBATCH --nodes=1
```

or:

```bash
#SBATCH -N 1
```

Keep API keys private. Do not add them to source files, scheduler scripts, or
calculation directories.

### Run the Agent

From the repository root:

```bash
python py4siesta-agent \
    "can you setup the dft calculation with lda for ReS2 monolayer"
```

After installation:

```bash
py4siesta-agent \
    "can you setup the dft calculation with lda for ReS2 monolayer"
```

To specify the calculation directory directly:

```bash
python py4siesta-agent \
    --workdir /absolute/path/to/calculation \
    "can you setup the dft calculation with lda for ReS2 monolayer"
```

The command stays active while scheduler jobs are queued or running. It
automatically submits each calculation stage, monitors the jobs, analyzes
completed outputs, and advances to the next stage. No additional agent command
is required.

For a complete ReS2 example covering model configuration, local
pseudopotentials, the supplied Slurm template, and the one-command workflow,
see the
[ReS2 agent tutorial](docs/tutorial_res2_agent.md).

To use an existing local structure, include its absolute FDF path in the
request.

### Agent Outputs

The calculation directory contains:

```text
origin/
01.kpoint_sampling/
02.slab_eos/ or 02.volume_eos/
03.geometry_optimization/
final_input/
.py4siesta-agent/checkpoints.sqlite
.py4siesta-agent/<workflow-id>.md
```

If a scheduler job fails or required SIESTA output is missing, the workflow
stops and records the error in its JSON state and Markdown summary. Failed
calculations are not automatically modified or resubmitted.

## Bundled NanoCore

This project is based on `NanoCore`, an older version of the ongoing [`OORINano`](https://github.com/yhkimlab/OORINano) project. `NanoCore` is imported locally from this repository.

That means:

- `python -m py4siesta` works from this repository because `NanoCore/` is present locally.
- `pip install .` or `pip install -e .` installs `py4siesta`, the
  `py4siesta_agent` package, and the bundled `NanoCore` package defined in
  `pyproject.toml`.
- There is no need to fetch `NanoCore` from an external site or package index.
