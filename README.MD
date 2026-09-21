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

The project has three program interfaces and an additional LLM assistant/skill mode:

| Interface / Mode | Purpose | Workflow scope | Status |
| --- | --- | --- | --- |
| `python py4siesta` | Interactive GUI/menu-based usage. | The user selects an existing numbered-menu workflow and supplies its inputs. | Available |
| `py4siesta-tool <command> [options]` | Direct, deterministic CLI tools with JSON output. | The user or a script invokes specific operations with explicit arguments. | Available |
| `py4siesta-agent "<request>"` | AI execution of predefined tasks requested in natural language. | Runs tools, checks results, and chooses the next step within the task's defined branches and stopping conditions. | Under development; initial DFT setup workflow available |
| LLM assistant with the `py4siesta` skill | Composes and improves workflows for the user's goal using existing NanoCore and py4siesta-tool operations. | Uses menu workflows as initial examples, can create new combinations, and accumulates validated revisions in personal memory and workflow registries. AI composes existing operations; it does not implement core functionality. | Codex skill available; [separate registration required](#codex-skill-setup-and-usage) |

After installation, the equivalent entry points are `py4siesta`, `py4siesta-tool`, and `py4siesta-agent`.

## Requirements

Core runtime dependencies are declared in `pyproject.toml`:

- Python 3.7 or newer for the non-agent interfaces
- `numpy`
- `scipy`
- `matplotlib`

The repository also installs the bundled `nanocore` package.

Agent frameworks, model-provider integrations, and the Materials Project
client are optional extras. The agent workflow requires Python 3.10 or newer.
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

Run `py4siesta-tool --help` or `py4siesta-tool <command> --help` for the full command list and required inputs.

`py4siesta-tool` is a deterministic CLI, not an AI agent. It performs only the command and options explicitly supplied by the caller.

## MCP Tools

An optional MCP stdio server exposes existing NanoCore and py4siesta APIs with
searchable descriptions, input schemas, units and side effects. Use Python 3.10+:

```bash
python -m pip install '.[mcp]'
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

The bundled [py4siesta skill](Skills/Codex/py4siesta/SKILL.md) uses connected
py4siesta MCP tools to prepare, run, and analyze calculation cases, or register/edit
reusable workflows. It supports free and registered-workflow execution; numbered
menus are optional references. The skill is separate from `py4siesta-agent` and does
not require its optional agent-framework dependencies.

### Register once for all projects

Prepare the Python environment, register the MCP server, and register the skill.
These are separate steps. Run the following commands in your terminal (Linux/macOS).
The skill itself does not install packages, change environments or MCP settings,
or fall back to CLI/direct Python when MCP is unavailable.

**1. Prepare a Python environment.** A dedicated Conda environment is recommended;
reuse a compatible environment if you already have one. For a new environment:

```bash
conda create -n py4siesta python=3.11 pip
conda activate py4siesta
cd /absolute/path/to/py4siesta
python -m pip install -e '.[mcp]' 'numpy==1.26.4'
python -m pip check
python -c "import sys; print(sys.executable)"
```

NumPy 1.26.4 is used here because a structure-rotation compatibility failure was
observed with NumPy 2.4.6. This installs the current checkout in editable mode;
keep it at that path. SIESTA, external post-processing utilities, and scheduler
access must be provided separately for tasks that need them.

**2. Register the MCP server with Codex.** Use the absolute Python path printed
above and an existing calculation directory; replace both example paths below:

```bash
codex mcp add py4siesta -- \
  /absolute/path/to/envs/py4siesta/bin/python \
  -m py4siesta_mcp --workdir /absolute/calculation/directory
codex mcp list
```

If `py4siesta` is already registered, check its command and update that entry to use
the intended Python and workdir. The absolute interpreter path means you do not
need to activate Conda every time Codex starts. `--workdir` is the server's default;
provide the actual calculation directory in each task. See the
[official Codex MCP documentation](https://learn.chatgpt.com/docs/extend/mcp?surface=cli).

**3. Register the skill.** From the repository root:

```bash
bash Skills/Codex/register.sh --dry-run
bash Skills/Codex/register.sh
```

The registration script links `~/.agents/skills/py4siesta` to this checkout, creates
the external personal state directory, and adds that path to Codex's writable roots.
It retains existing settings, backs up changed configuration, and does not replace
an existing skill at another location. It does not install or register the MCP server.
See [skill registration details](Skills/Codex/README.md) for custom state paths and
managed filesystem permissions.

**4. Restart Codex and verify the connection.** Check `/mcp` for the server and
`/skills` for the skill. Then ask:

```text
Use the py4siesta MCP to read py4siesta://runtime and search
py4siesta_api_catalog for read_xyz. Report the Python path and available tools.
Do not generate files or submit calculations.
```

`codex mcp list` confirms registration; a successful MCP resource read and tool call
confirm the connection. `python -m py4siesta_mcp --catalog` only checks local catalog
generation. If connection fails, check the configured Python, existing workdir,
server errors, and client permissions before starting a calculation. Passing a
standalone MCP test does not establish that the Codex connection works.

### Give Codex a calculation task

Invoke `$py4siesta` explicitly and supply the absolute calculation directory, inputs,
and desired scope. Generation workflows usually require `origin/input/STRUCT.fdf`
and matching calculation inputs; post-processing uses the relevant output files.

| Mode | Behavior |
| --- | --- |
| Free | Discover and compose MCP tools for the user's goal. |
| Workflow execution | Follow a registered procedure and its allowed branches. |
| Workflow registration/editing | Create or revise a reusable procedure when requested. |

An explicit mode takes precedence. Otherwise a clearly matching compatible validated
workflow is selected; other calculation requests use free mode. The skill announces
its mode. The five bundled [MCP workflow examples](Skills/Codex/py4siesta/workflows/)
cover k-points, geometry scans/fitting, geometry utilities, submission, and
post-processing. They remain drafts pending workflow execution validation; they
are not automatically selected as validated procedures.

Example prompts (replace paths and scientific parameters with your own):

```text
$py4siesta In free mode, prepare bulk k-point cases in /work/diamond
for grids 2, 4, and 6. Check generated inputs; do not submit jobs.

$py4siesta Run my validated k-point analysis workflow in /work/diamond,
using a 0.01 eV tolerance. Report missing or unconverged calculations.

$py4siesta Register the procedure we just executed as a personal workflow,
including its MCP tools, inputs, steps, and validation evidence.

$py4siesta Resume my previous task in /work/diamond. Inspect existing
outputs and submission status before continuing; do not resubmit blindly.
```

Case generation can recreate existing output directories; job submission must be
part of your request. A generated structure or submitted job is not yet a completed,
converged SIESTA calculation. Missing MCP tools or environment prerequisites stop
the dependent work. Workflow drafting can proceed offline as unverified.

### Personal history and reusable workflows

Personal workflows, memory, and necessary history stay outside this repository,
by default in `~/.local/state/py4siesta-skill`. Set `PY4SIESTA_SKILL_STATE` before
registration and keep it set when starting Codex to use another external location.

Workflow registration/editing occurs only when requested and defaults to personal
storage. Request public/repository registration explicitly to add a distributable
workflow. New or unverified procedures remain drafts; meaningful revisions preserve
prior versions and evidence. Registration alone does not authorize calculations.

The skill saves only useful preferences, verified environment facts, nonduplicate
observations, and history needed for resumption, uncertain submissions, investigation,
validation, or requested auditing. Existing logs are referenced rather than copied;
unchanged repeats create no new record. Saved records are read back, and unavailable
persistence is reported. Current instructions override memory. No background learning
or monitoring continues after the session ends.

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
