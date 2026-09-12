---
name: py4siesta
description: Prepare, run, and analyze SIESTA cases using existing py4siesta numbered-menu recipes, including geometry optimization, k-point convergence, structure interpolation, band structure, PDOS, and PLDOS. Use for py4siesta calculation work and workflow reuse, not general DFT advice or unrelated code changes.
---

# py4siesta calculation workflows

Use the installed py4siesta program as the execution engine. This skill supplies
menu-based procedures and persistent experience, not a replacement calculator.
Follow the user's requested scope and the working project's instructions.

## Start a task

1. Identify the calculation directory, requested menu/function, inputs, and whether
   the user wants preparation, submission, analysis, or a complete calculation.
2. Read [runtime.md](references/runtime.md). Verify the executable/Python environment
   and actual command help before writing calculation files.
3. Read [state.md](references/state.md). For execution tasks, select the external
   state directory, load only relevant memories and workflow revisions, and create
   a run record. Explanation-only requests do not create persistent records.
4. Select the recipe from [menu-map.md](references/menu-map.md) and read it fully.
   Record the recipe revision, inputs, planned steps, side effects, and completion
   criteria before execution. Current user instructions override remembered defaults.

## Execute and learn

- Prefer `py4siesta-tool`; preserve the numbered-menu mapping in plans and reports.
  Run from the calculation directory, never from the installed skill directory.
- Check existing outputs before generation, fitting, plotting, or resubmission.
  Follow the recipe's actual overwrite behavior. Preparation does not authorize
  submitting jobs. Do not repeat a batch submission after an ambiguous failure.
- Record each attempt and verify its outputs. A successful process or JSON result
  is not proof of SCF convergence, job completion, or physical validity.
- On failure, preserve evidence, identify a testable cause, and retry only when a
  changed condition supports it and the action remains in scope. If the same cause
  persists after a targeted retry, stop that step and report the missing requirement.
- For interrupted work, read history and inspect files/job status before continuing;
  do not regenerate completed cases. Record waiting jobs as `waiting`, not complete.
- Follow the promotion and version rules in [state.md](references/state.md) to
  revise personal workflows using verified lessons. Keep unknown causes in history.
- If no recipe fits, compose existing menu steps using the
  [workflow template](assets/workflow.md), save a personal draft, and validate it.
  Do not invent commands, silently alter scientific settings, or modify py4siesta
  source to make the recipe work. Report unsupported steps.

Finish with the menu/function used, generated paths, validation, outstanding work,
run-record location, and any memory/workflow revision made. State clearly when
persistence was unavailable. The skill performs no background learning or monitoring
after the agent session ends.
