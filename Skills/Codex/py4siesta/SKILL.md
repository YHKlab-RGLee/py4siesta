---
name: py4siesta
description: Use only when the user explicitly requests the py4siesta skill or calculation work using py4siesta. Prepare, run, and analyze cases by composing existing NanoCore and py4siesta-tool operations and reusing validated personal workflows. SIESTA-related requests alone do not trigger this skill; exclude general DFT advice and source-code development.
---

# py4siesta calculation workflows

Compose existing public NanoCore and py4siesta-tool operations to fulfill the user's
goal. Menu workflows are initial examples, not required paths or a closed catalog.
Use personal memory and workflow registries to reuse and improve validated compositions.
Do not implement or modify core scientific or deterministic tool functionality.
Workflow definitions and glue code may connect existing operations, inputs, outputs,
branches, and validation steps without reproducing their domain logic.

This LLM assistant/skill has broader workflow freedom than `py4siesta-agent`, which
executes predefined tasks in loops within allowed branches and stopping conditions.
It may invoke that agent for a matching task through a verified public interface;
using the skill does not require the agent or its optional dependencies.
Follow the user's requested scope and the working project's instructions.

## Start a task

1. Identify the calculation directory, goal/function, inputs, and whether
   the user wants preparation, submission, analysis, or a complete calculation.
2. Read [state.md](references/state.md). For execution tasks, select the external
   state directory and load only relevant memories and workflow revisions. Persist
   only useful new information or required resume/validation evidence. Repetition
   or explanation alone creates no record; independently useful new facts still qualify.
3. Read [runtime.md](references/runtime.md). Verify the executable/Python environment
   and actual command help before writing calculation files.
4. Reuse a compatible personal workflow, adapt an initial example from
   [menu-map.md](references/menu-map.md), or compose a new workflow using verified
   public operations. Read any selected example fully to retain its operational
   constraints. A menu mapping or base recipe is optional. Establish the applicable
   revision, inputs, steps, side effects, and completion criteria before execution;
   do not copy an unchanged workflow into a new persistent plan. Current user instructions override remembered defaults.

## Execute and learn

- Prefer `py4siesta-tool` for non-interactive execution; direct public NanoCore
  calls are also available after verifying their interface and environment.
  Include numbered-menu mappings only where applicable.
  Run from the calculation directory, never from the installed skill directory.
- Check existing outputs before generation, fitting, plotting, or resubmission.
  Follow the recipe's actual overwrite behavior. Preparation does not authorize
  submitting jobs. Do not repeat a batch submission after an ambiguous failure.
- Verify each attempt; retain evidence only under the state retention rules. A successful process or JSON result
  is not proof of SCF convergence, job completion, or physical validity.
- On failure, preserve evidence, identify a testable cause, and retry only when a
  changed condition supports it and the action remains in scope. If the same cause
  persists after a targeted retry, stop that step and report the missing requirement.
- For interrupted work, read history and inspect files/job status before continuing;
  do not regenerate completed cases. Record waiting jobs as `waiting`, not complete.
- Follow the promotion and version rules in [state.md](references/state.md) to
  improve reusable procedures directly in workflows when a meaningful change is
  validated. Do not create separate lesson files or revisions for repeated success.
- Use the [workflow template](assets/workflow.md) for new or adapted compositions,
  save a draft only when it merits reuse or further validation, and validate the
  executed steps before promotion.
  If a step cannot be composed from existing public operations, report the missing
  observed limitation, attempted operation, and actual workaround or failure.
  Record new nonduplicate friction in memory/observations.md, including inefficient
  successful work. Do not propose implementation, API design, layer ownership, or
  priority as part of these observations. Core/tool implementation is a
  separate, explicitly requested development task, not workflow improvement.
  Do not invent commands or APIs, silently alter scientific settings, or modify
  source to make a workflow work.

Before finishing, including after a failed or blocked task, evaluate each save condition
in [state.md](references/state.md), complete qualifying writes, and read back changed
files to verify their location and contents. Do not create records merely to prove
that this check occurred.

Finish with the workflow and operations used, generated paths, validation, outstanding work,
record locations only when records were needed, and any memory/workflow change made. State clearly when
persistence was unavailable. The skill performs no background learning or monitoring
after the agent session ends.
