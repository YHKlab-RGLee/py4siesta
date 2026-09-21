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
- [MCP Tools](#mcp-tools)
- [Codex Skill: Setup and Usage](#codex-skill-setup-and-usage)
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
- The local bundled `nanocore/` package is used directly from this repository.
- You do not need to download `NanoCore` from an external source.

The bundled Python package is named `nanocore`: use `import nanocore`
(previously `import NanoCore`). After updating an existing checkout, rerun
`python -m pip install -e .` to refresh the installed package names.

## Quick Start

Clone and enter the repository:

```bash
git clone https://github.com/YHKlab-RGLee/py4siesta.git
cd py4siesta
```

Prepare a Python environment for all interfaces. Reuse a compatible environment,
or create a dedicated Conda environment:

```bash
conda create -n py4siesta python=3.11 pip
conda activate py4siesta
```

Install from the repository root:

```bash
python -m pip install . 'numpy==1.26.4'
```

Use `-e` for an editable installation that follows changes in this checkout.
For MCP tools and the Codex skill, replace `.` with `'.[mcp]'` in the installation
command, for example:

```bash
python -m pip install -e '.[mcp]' 'numpy==1.26.4'
python -m pip check
```

NumPy 1.26.4 is the tested compatibility version; NumPy 2.4.6 caused a structure
rotation failure. SIESTA executables, utilities, and scheduler access are configured
separately under [Configure SIESTA Paths](#configure-siesta-paths).

Run the interactive numbered-menu interface:

```bash
py4siesta
```

Or run it directly from the repository root:

```bash
python -m py4siesta
```

The project has three program interfaces and an additional LLM assistant/skill mode:

| Interface / Mode | Purpose | Workflow scope | Status |
| --- | --- | --- | --- |
| `python py4siesta` | Interactive GUI/menu-based usage. | The user selects an existing numbered-menu workflow and supplies its inputs. | Available |
| `py4siesta-tool <command> [options]` | Direct, deterministic CLI tools with JSON output. | The user or a script invokes specific operations with explicit arguments. | Available |
| `py4siesta-agent "<request>"` | AI execution of predefined tasks requested in natural language. | Runs tools, checks results, and chooses the next step within the task's defined branches and stopping conditions. | Under development; initial DFT setup workflow available |
| LLM assistant with the `py4siesta` skill | Solves calculation tasks through py4siesta MCP tools. | Free execution, registered workflow execution, or requested workflow registration/editing. | Codex skill available; [separate registration required](#codex-skill-setup-and-usage) |

After installation, the equivalent entry points are `py4siesta`, `py4siesta-tool`, and `py4siesta-agent`.

## Requirements

Core runtime dependencies are declared in `pyproject.toml`:

- Python 3.7 or newer for the non-agent interfaces
- `numpy`
- `scipy`
- `matplotlib`

The repository also installs the bundled `nanocore` package.

Agent frameworks, model-provider integrations, and the Materials Project
client and MCP support are optional extras. MCP and the agent workflow require
Python 3.10 or newer.
C2DB monolayer structures use its public OPTIMADE endpoint and require no
additional Python package.

## Configure SIESTA Paths

The bundled `NanoCore` module reads local executable paths from `nanocore/env.py`. Update this file before using features that call SIESTA utilities.

Common values to edit:

```python
siesta_calculator = "/path/to/siesta"
siesta_psf_location = "/path/to/psf"
siesta_util_location = "/path/to/siesta/utils"
siesta_util_band = "gnubands"
siesta_util_dos = "Eig2DOS"
siesta_util_pdos = "fmpdos"
siesta_util_rho = "rho2xsf"
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
by `siesta_psf_location` in `nanocore/env.py`. It refuses to overwrite an
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
| `14` | Planeaverage grid | Select VH, VT, RHO, or DRHO and axis x/y/z (0/1/2); a matching grid file. | `planeaverage_<TARGET>_<axis>.png` and `.txt` beside the input file. |

Notes:

- Menu `11` automatically selects a local `.bands` file when one is present.
- Menu `12` uses a matching `<label>.bands` file for `E_v` when available; otherwise it estimates `E_v` from `<label>.EIG`.
- Menu `12` removes existing same-named `fmpdos` output files before regenerating PDOS data.
- Menu `13` automatically selects a local `.PDOS` file when one is present.
- Menu `14` accepts case-insensitive targets and axes, selects the first matching local file or asks for its path, and overwrites same-named outputs. TXT columns are distance (Å) and average value: eV for VH/VT, e/Å³ for RHO/DRHO. Spin components are averaged.

### Utility

| Menu | Feature | Purpose |
| --- | --- | --- |
| `01` | Generate Geometries | Open the geometry-generation submenu. |

Geometry-generation submenu:

| Submenu | Feature | Purpose |
| --- | --- | --- |
| `1` | Move structure | Translate the full structure by a user-defined displacement vector and write a new `STRUCT.fdf`. |
| `2` | Interpolate structure | Generate interpolated structures between two `STRUCT.fdf` files, with optional extrapolated points. |
| `3` | Job submission | Submit `slm_*` scripts in generated cases under `11.interpolate_structure` using `sbatch`. |

In `python py4siesta` → `01` → `2`, case folders use polar mode labels:
initial = `-1`, final = `+1`, and `0` is their midpoint (the para phase when
the input pair has the corresponding symmetry). The division count includes
both endpoints; choose an odd count to include mode `0`. Extrapolate `N` adds
**N points on each side**, below `-1` and above `+1`, for `2N` extra cases.
For example, division `5` and extrapolate `2` produce modes
`-2, -1.5, -1, -0.5, 0, +0.5, +1, +1.5, +2`, with folders such as
`01-mode_-2.0000` and `05-mode_+0.0000`.

Both this menu and `python py4siesta-tool interpolate-structure` also write
`11.interpolate_structure/interpolate.ANI`: concatenated XYZ frames in case
order, including both endpoints and any extrapolated structures.

Use `01` → `3` to submit these generated jobs. The equivalent direct CLI
command is `python py4siesta-tool submit --mode geometry`.

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
py4siesta-tool planeaverage-grid --file-path Device.VH --target vh --axis z
```

For `py4siesta-tool pdos`, pass one or more `--orbital` selections in `atom_or_species[_n[_l[_m]]]` format, for example `Mg_0`, `O_0`, or `C_2_1_0`. Multiple selections may be space-separated, comma-separated, or provided with repeated `--orbital` options.

For explicit Slurm jobs, use any script filename; `origin/` is not required:

```bash
py4siesta-tool job-submit --case /work/case01 --script run.slurm
py4siesta-tool job-status --job-id 12345
py4siesta-tool job-cancel --job-id 12345
```

Relative script paths resolve inside `--case`. Submission returns the job ID and
optional cluster name; pass `--cluster NAME` to status/cancel if returned. Status
reports the Slurm state and accounting exit code; cancellation reports a request,
not confirmed termination. Array task IDs such as `12345_2` are accepted, but array
parent status is not aggregated. Missing records are `unknown`; query errors and
uncertain submissions return errors with details. Never blindly retry submission.
These commands reuse your configured Slurm environment and do not edit scripts.
Existing `submit --mode ...` behavior is unchanged.

Run `py4siesta-tool --help` or `py4siesta-tool <command> --help` for the full command list and required inputs.

`py4siesta-tool` is a deterministic CLI, not an AI agent. It performs only the command and options explicitly supplied by the caller.

## MCP Tools

An optional MCP stdio server exposes existing NanoCore and py4siesta APIs with
searchable descriptions, input schemas, units and side effects. Install the `mcp`
extra as described in [Quick Start](#quick-start), then launch the server with:

```bash
python -m py4siesta_mcp --workdir /absolute/calculation/directory
```

Configure this command in your MCP client. `py4siesta_api_catalog` searches the
API inventory; `python -m py4siesta_mcp --catalog` exports it as JSON. The server
can read structures, manipulate session objects, generate files, analyze outputs
and invoke existing submission APIs when called. Legacy APIs that cannot be
exposed are listed with reasons. Existing menu and direct CLI usage remain available.
See [MCP usage and API categories](py4siesta_mcp/README.md) for arguments,
object handles, output conventions and examples.

## Codex Skill: Setup and Usage

The [py4siesta skill](Skills/Codex/py4siesta/SKILL.md) lets Codex solve calculation
tasks through MCP and manage reusable workflows. It is separate from
`py4siesta-agent` and does not require agent-framework dependencies.

### Register once for all projects

Complete [Quick Start](#quick-start) with the `mcp` extra first. The following steps
connect that existing environment to Codex; the skill does not install or repair it.

**1. Register the MCP server.** In the prepared environment, find the Python path:

```bash
python -c "import sys; print(sys.executable)"
```

Use that absolute path and an existing calculation directory:

```bash
codex mcp add py4siesta -- \
  /absolute/path/to/envs/py4siesta/bin/python \
  -m py4siesta_mcp --workdir /absolute/calculation/directory
codex mcp list
```

If already registered, update the existing entry to the intended Python and workdir.
The absolute Python path avoids needing to activate Conda whenever Codex starts.
See the [official MCP guide](https://learn.chatgpt.com/docs/extend/mcp?surface=cli).

**2. Register the skill.** From the repository root:

```bash
bash Skills/Codex/register.sh --dry-run
bash Skills/Codex/register.sh
```

This links the skill to the checkout and configures its external personal state
location. Keep the checkout in place. See [registration details](Skills/Codex/README.md)
for custom paths and permissions.

**3. Restart Codex and verify.** Check `/mcp` and `/skills`, then ask:

```text
Use py4siesta MCP to read py4siesta://runtime and search
py4siesta_api_catalog for read_xyz. Report the Python path and tools.
```

A successful resource read and tool call confirm the connection. Registration or
local `--catalog` output alone does not. If connection fails, check the configured
Python, workdir, server errors, and client permissions.

### Give Codex a calculation task

Invoke `$py4siesta` explicitly with the calculation directory, inputs, and scope.
Generation usually needs `origin/input/STRUCT.fdf` and matching calculation inputs;
post-processing uses the relevant output files.

| Mode | Example request |
| --- | --- |
| Free | `$py4siesta In free mode, prepare bulk k-point grids 2, 4, 6 in /work/diamond. Do not submit.` |
| Workflow execution | `$py4siesta Run my validated k-point analysis workflow in /work/diamond with a 0.01 eV tolerance.` |
| Workflow registration/editing | `$py4siesta Register the procedure we just executed as a personal workflow with validation evidence.` |

Without an explicit mode, the skill uses a clearly matching validated workflow or
free mode. The five bundled [workflow examples](Skills/Codex/py4siesta/workflows/)
are drafts pending execution validation. Registration alone does not run calculations.
Submission and overwriting existing results must be within the requested scope;
prepared inputs or submitted jobs are not completed calculations.

### Personal history and reusable workflows

Personal workflows and necessary memory/history are stored outside the repository,
by default in `~/.local/state/py4siesta-skill`. Set `PY4SIESTA_SKILL_STATE` before
registration and when starting Codex to choose another external location.

Workflows are registered or revised only when requested; public registration must
be explicit. Unchanged success creates no new record. The skill retains useful
facts and necessary resume evidence, reports persistence failures, and performs
no background learning or monitoring. See [state details](Skills/Codex/README.md#selective-memory).

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

Also configure `siesta_psf_location` in `nanocore/env.py` as described in
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

This project is based on `NanoCore`, an older version of the ongoing [`OORINano`](https://github.com/yhkimlab/OORINano) project. `nanocore` is imported locally from this repository.

That means:

- `python -m py4siesta` works from this repository because `nanocore/` is present locally.
- `pip install .` or `pip install -e .` installs `py4siesta`, the
  `py4siesta_agent` package, and the bundled `nanocore` package defined in
  `pyproject.toml`.
- There is no need to fetch `NanoCore` from an external site or package index.
