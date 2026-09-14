---
name: py4siesta
description: Prepare, run, and analyze SIESTA cases by composing existing NanoCore and py4siesta-tool operations and progressively improving personal workflows from validated experience. Use for py4siesta calculation work and workflow reuse, not general DFT advice or unrelated code changes.
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
2. Read [runtime.md](references/runtime.md). Verify the executable/Python environment
   and actual command help before writing calculation files.
3. Read [state.md](references/state.md). For execution tasks, select the external
   state directory, load only relevant memories and workflow revisions, and create
   a run record. Explanation-only requests do not create persistent records.
4. Reuse a compatible personal workflow, adapt an initial example from
   [menu-map.md](references/menu-map.md), or compose a new workflow using verified
   public operations. Read any selected example fully to retain its operational
   constraints. A menu mapping or base recipe is optional. Record workflow revision,
   tool interfaces, inputs, planned steps, side effects, and completion criteria
   before execution. Current user instructions override remembered defaults.

## Execute and learn

- Prefer `py4siesta-tool` for non-interactive execution; direct public NanoCore
  calls are also available after verifying their interface and environment.
  Include numbered-menu mappings only where applicable.
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
- Use the [workflow template](assets/workflow.md) for new or adapted compositions,
  save a personal draft, and validate the executed steps before promotion.
  If a step cannot be composed from existing public operations, report the missing
  operation and required inputs/outputs/interface. Core/tool implementation is a
  separate, explicitly requested development task, not workflow improvement.
  Do not invent commands or APIs, silently alter scientific settings, or modify
  source to make a workflow work.

Finish with the workflow and operations used, generated paths, validation, outstanding work,
run-record location, and any memory/workflow revision made. State clearly when
persistence was unavailable. The skill performs no background learning or monitoring
after the agent session ends.
