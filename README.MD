# py4siesta

`py4siesta` provides interactive and deterministic command-line interfaces for preparing, organizing, submitting, and post-processing SIESTA calculation cases. A separate AI-agent orchestration interface is planned.

It is designed around one workflow:

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
| `py4siesta-agent "<request>"` | Natural-language AI-agent orchestration over existing tools. | Skeleton only; under development |

After installation, the equivalent entry points are `py4siesta`, `py4siesta-tool`, and `py4siesta-agent`.

## Requirements

Runtime Python dependencies are declared in `pyproject.toml`:

- `numpy`
- `scipy`
- `matplotlib`

The repository also installs the bundled `NanoCore` package.

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

The planned natural-language orchestration entry point is:

```bash
py4siesta-agent "Can you set up the input parameters for ReS2?"
```

After installation, it can be invoked as `py4siesta-agent`. This interface is currently an under-development skeleton: it accepts the future prompt shape but does not make LLM decisions or run autonomous scientific workflows.

The intended layer hierarchy is:

```text
NanoCore
    ↓
py4siesta / py4siesta-tool
    ↓
py4siesta-agent
```

The dependency direction is one-way: `py4siesta-agent` may depend on `py4siesta-tool` and `py4siesta`, which may depend on `NanoCore`, never the reverse. `NanoCore` owns reusable scientific functionality. `py4siesta` and `py4siesta-tool` provide deterministic application workflows over it. The agent layer may orchestrate those existing tools through thin adapters, but must not duplicate their scientific or workflow implementations. NanoCore and the deterministic interfaces do not depend on the agent package or any future agent-framework dependencies.

## Bundled NanoCore

This project is based on `NanoCore`, an older version of the ongoing [`OORINano`](https://github.com/yhkimlab/OORINano) project. `NanoCore` is imported locally from this repository.

That means:

- `python -m py4siesta` works from this repository because `NanoCore/` is present locally.
- `pip install .` or `pip install -e .` installs `py4siesta`, the `py4siesta_agent` skeleton, and the bundled `NanoCore` package defined in `pyproject.toml`.
- There is no need to fetch `NanoCore` from an external site or package index.
