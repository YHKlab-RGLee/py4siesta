# Runtime and execution contract

The command map is based on the bundled `py4siesta/tool_cli.py`, `cli.py`,
`operations.py`, and `post_process.py`. Package version 1.0.0 alone does not identify
all behavior: record a source commit and dirty status when available, plus the
imported module location. Recheck relevant help/source if an installed version differs.

## Locate the program

Use the user-selected environment first, then the external state's saved environment,
then an executable on PATH. Verify `py4siesta-tool --help` and the selected subcommand's
`--help`. With a known Python interpreter, `python -m py4siesta.tool_cli --help`
is an equivalent entrypoint. Verify `python -c "import py4siesta; print(py4siesta.__file__)"`
in that same environment when diagnosing conflicting installs.

For an uninstalled checkout, use its absolute path in PYTHONPATH for the selected
Python process while retaining the calculation directory as cwd. Never infer the
source path from a globally copied skill's ancestors. If imports fail, identify the
missing environment/dependency; skill registration does not install dependencies.
Do not change environments or install software unless included in the user's task.

## Inputs and outputs

- All workflow-backed commands, including analysis, submission, move, and
  interpolation, initialize from `cwd/origin/input/STRUCT.fdf`. Post-processing
  commands `band`, `pdos`, and `pldos` do not require `origin`.
- Case generation copies the entire `origin/` tree, then replaces the generated
  input file. Actual calculations need appropriate RUN, BASIS, KPT, pseudopotentials,
  and job scripts; a readable STRUCT alone only suffices for generating structures.
- BaseOperation deletes and recreates its named output directory BEFORE validating
  all inputs. Resolve and validate inputs first. Keep source files outside that output
  tree. If results already exist, use a separate agreed calculation directory or
  obtain overwrite authorization before running. Do not silently move/delete results.
- Capture stdout, stderr, exit status, cwd, and exact arguments. Some operations print
  progress before the terminal JSON object; do not parse all stdout as a single JSON
  document. Preserve logs and inspect the terminal payload and exit status together.
  Argument-parser failures may be plain text rather than JSON.
- Verify file existence and contents using the specific recipe, and report missing
  external utilities. Local SIESTA utility configuration is in `NanoCore/env.py`.
- None of the preparation commands itself executes SIESTA. Job execution depends on
  the user's scheduler and scripts. Never present prepared inputs as optimized results.
