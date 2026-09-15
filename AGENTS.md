## Project Overview

`py4siesta` is a utility for preparing, organizing, and managing SIESTA calculation cases.

The project provides the following user-facing interfaces:

* `python py4siesta` provides the GUI/menu-based user interface.
* `python py4siesta-tool` provides direct command-line access to individual tools and deterministic functionality.
* `python py4siesta-agent` provides the AI-agent-based interface and is currently under development.

The project follows a layered architecture:

* `NanoCore` contains reusable core scientific functionality.
* `py4siesta-tool` exposes reusable deterministic operations using NanoCore functionality.
* `py4siesta` provides user-facing numbered-menu workflows that compose these operations.
* `py4siesta-agent` composes existing operations in loops for predefined tasks, with decisions restricted to the task's allowed steps and branches.
* An LLM assistant/skill has broader responsibility for composing and progressively improving workflows using existing operations and accumulated personal memory and workflow registries. It may use `py4siesta-agent` to execute a predefined task.

This is a hierarchy of responsibilities, not a requirement to call through the GUI/menu interface. Both AI modes may use explicit public NanoCore and py4siesta-tool interfaces directly. Existing menu workflows are initial examples, not limits on possible compositions. Neither AI mode implements or modifies core scientific or deterministic tool functionality as part of workflow execution or improvement.

All GUI/menu-based functionality in `py4siesta` is accessed through numbered menu entries. Each feature should be implemented as a clearly separated functionality and exposed through the menu system using a menu number.

Direct CLI tools should be exposed through `py4siesta-tool`, while agent-specific orchestration should remain isolated under `py4siesta-agent`.

## Critical File Protection Rule

This `AGENTS.md` file must never be modified by agents.

Do not edit, rewrite, reformat, rename, move, delete, or automatically update this file under any circumstance. Any changes to this file will be made manually by the project owner.

## Design Rules

1. Preserve the existing user-interface designs.

   * All `py4siesta` GUI/menu functionality must remain accessible through numbered menu options.
   * New GUI/menu functionality should be added as a new menu entry.
   * Direct CLI functionality should be exposed through `py4siesta-tool`.
   * AI-agent functionality should be exposed through `py4siesta-agent`.
   * Existing user workflows should remain stable unless the project owner explicitly requests otherwise.

2. Preserve the layered architecture.

   * `NanoCore` must contain reusable core scientific functionality and must remain independent of py4siesta and py4siesta-agent.
   * `py4siesta` and `py4siesta-tool` may use NanoCore functionality but must remain independent of agent-specific dependencies.
   * `py4siesta-agent` may use NanoCore, py4siesta, and py4siesta-tool functionality through explicit interfaces.
   * NanoCore, py4siesta, and py4siesta-tool must not depend on py4siesta-agent.
   * Agent-specific prompts, model clients, workflow state, and orchestration logic must not be introduced into core modules.

3. Keep functionality isolated.

   * When adding or updating a specific functionality, do not modify unrelated features.
   * Avoid broad refactoring unless it is explicitly required for the requested change.
   * Do not change the behavior of existing functions while implementing a new feature.
   * Core scientific or deterministic workflow functionality must not be duplicated inside py4siesta-agent.
   * First distinguish a new composition of existing operations from a genuinely missing core/tool operation. AI modes may compose existing operations within their task scope. If an operation is missing, report the gap and required interface; lower-layer implementation is a separate, explicitly requested development task, not autonomous workflow improvement.

4. Maintain backward compatibility.

   * Existing menu options, CLI commands, arguments, and workflows should continue to work.
   * If menu numbering must be updated, only adjust the menu-number assignment and related display text.
   * Do not remove, rename, or alter existing functionality unless explicitly instructed.
   * Refactoring for the layered architecture must not change existing py4siesta or py4siesta-tool behavior.

5. Minimize the scope of changes.

   * Implement only the requested functionality.
   * Avoid unrelated cleanup, formatting changes, style changes, or dependency changes.
   * Do not rewrite large sections of code when a small, localized update is sufficient.
   * Keep agent-related dependencies isolated from NanoCore, py4siesta, and py4siesta-tool.

6. Keep the project structure organized.

   * New code should follow the existing project layout and naming conventions.
   * Feature-specific logic should be placed near related functionality.
   * Avoid mixing unrelated logic in the same function or module.
   * Core scientific functionality should be placed in NanoCore.
   * SIESTA-specific application and deterministic workflow functionality should be placed in py4siesta or py4siesta-tool.
   * Agent-specific orchestration, prompts, state management, and tool adapters should be placed in py4siesta-agent.
   * Agent tool adapters should remain thin and should call existing deterministic functionality rather than implement scientific operations directly.

## Update Rules

1. Do not update `AGENTS.md`.

   * This file is owned and maintained manually by the project owner.
   * Agents must treat this file as read-only.

2. When adding a new GUI/menu feature:

   * Add the feature as a new numbered menu option in `py4siesta`.
   * Update menu numbering only as needed.
   * Ensure the new option does not interfere with existing menu behavior.
   * Do not modify unrelated functionality.

3. When adding a new direct CLI tool:

   * Add the functionality through `py4siesta-tool`.
   * Preserve existing CLI commands and arguments.
   * Reuse NanoCore or py4siesta functionality rather than duplicating it.
   * Do not describe direct CLI tools as AI-agent functionality.

4. When adding or updating agent functionality:

   * Keep agent-specific code isolated within `py4siesta-agent`.
   * Use existing NanoCore, py4siesta, and py4siesta-tool interfaces.
   * Keep `py4siesta-agent` loops within predefined tasks, allowed branches, validation rules, and stopping conditions. Open-ended workflow creation and registry improvement currently belong to the LLM assistant/skill, not this agent interface.
   * Do not implement core scientific functionality inside the agent layer.
   * Do not introduce agent-framework or LLM dependencies into NanoCore, py4siesta, or py4siesta-tool.
   * Preserve the independent operation of all non-agent interfaces.
   * Until the agent is implemented, retain only the required skeleton and entry-point structure.

5. When updating an existing feature:

   * Restrict changes to the requested feature only.
   * Do not change other menu options, CLI tools, helper functions, file formats, or workflows unless strictly necessary.
   * If a shared utility must be changed, confirm that existing behavior remains compatible.
   * Place reusable functionality in the appropriate lower layer rather than implementing it only for a higher-level interface.

6. When updating documentation:

   * Do not rewrite the entire `README.md`.
   * Add only the user-facing documentation necessary to describe the new or updated functionality.
   * Preserve the existing README structure, wording, and sections whenever possible.
   * Clearly distinguish:

     * `python py4siesta` as GUI/menu-based usage,
     * `python py4siesta-tool` as direct CLI tool usage,
     * `python py4siesta-agent` as AI-agent-based usage.
   * Do not describe `py4siesta-tool` as an AI agent.
   * If `py4siesta-agent` is not yet implemented, describe it as planned or under development.
   * Keep internal architecture, module boundaries, dependency rules, refactoring details, and implementation-specific design decisions out of the README unless they are directly relevant to users.

## README Update Policy

The README is user-facing documentation and should be updated incrementally.

When a new feature is added, include only a concise user-facing description, such as:

* the relevant command or menu entry,
* the purpose of the feature,
* basic usage,
* required user inputs,
* generated outputs or visible side effects,
* whether the feature is available, experimental, or under development.

The README must consistently distinguish between:

* `python py4siesta` for GUI/menu-based usage,
* `python py4siesta-tool` for direct CLI usage,
* `python py4siesta-agent` for AI-agent-based usage.

Do not include detailed internal architecture, package dependency rules, module ownership, adapter design, or maintenance strategy in the README unless explicitly requested by the project owner.

Do not rewrite the full README unless explicitly requested by the project owner.


## Expected Agent Behavior

Before making changes, agents should identify:

* which specific functionality is being added or updated,
* which architectural layer should own the functionality,
* which files are directly relevant,
* whether menu numbering or CLI entry points need to be adjusted,
* whether the README requires a small incremental update,
* whether the requested change can reuse existing NanoCore, py4siesta, or py4siesta-tool functionality.

After making changes, agents should verify:

* the `py4siesta` GUI/menu interface still works,
* existing menu options still behave as before,
* existing `py4siesta-tool` commands and arguments still behave as before,
* the new or updated functionality is reachable through the correct interface,
* no core functionality was duplicated inside py4siesta-agent,
* no agent-specific dependency was introduced into NanoCore, py4siesta, or py4siesta-tool,
* unrelated files were not modified,
* `AGENTS.md` was not changed.

## Agent Skills Development

* Distribute reusable agent skills within `Skills/Codex/` in this repository, suitable for global registration by users.
* Treat existing numbered-menu workflows as verified initial examples, not mandatory execution paths or a closed set of workflows. Document menu mappings where applicable, required inputs, tool interfaces, execution steps, and output validation.
* Use existing py4siesta interfaces, preferring `py4siesta-tool` for non-interactive execution. Verify the executable or source path and Python environment before use; skill registration does not install py4siesta.
* Let the LLM assistant/skill retrieve, reuse, compose, and progressively improve workflows for the user's goal using existing public NanoCore and py4siesta-tool operations. Workflow definitions and glue code may connect operations, inputs, outputs, branches, and validation steps; they must not implement or modify core scientific or deterministic tool functionality.
* Keep skill development isolated from existing code, menu behavior, and project layout outside `Skills/Codex/`.
* Keep personal state outside the public repository and installed skill tree.
  Persist only information that changes future decisions, validates a meaningful
  workflow improvement, or enables safe resumption. Repeated unchanged success
  requires no new persistent record, unless the user requests an audit trail.
* Let workflows own reusable execution procedures and validation conditions.
  Create revisions only for meaningful changes; preserve previous revisions and
  evidence, but retrieve only the compatible current procedure for routine work.
  Do not maintain a separate duplicate collection of procedural lessons.
* Keep memory limited to explicit user preferences, verified environment facts,
  and concise nonduplicate observations of friction. Observations describe the
  attempted task, limitation, actual workaround or failure, and minimal evidence;
  successful but inefficient work also qualifies. Do not propose implementations,
  assign architectural ownership, design APIs, or prioritize development here.
* Retain history only for requested auditing, safe resumption, uncertain submission
  status, continuing investigation, or necessary validation evidence. Reference
  existing calculation outputs/logs rather than copying them. Preserve evidence
  already retained; do not automatically delete prior records.
* Retrieve only task-relevant state. Do not read all historical revisions, logs,
  memories, or observations at startup. Current user instructions take precedence
  over remembered defaults; unverified attempts must not become validated workflows.

## Non-Negotiable Constraints

* Never modify `AGENTS.md`.
* Never rewrite the entire README for a small feature update.
* Never change unrelated functionality when adding or updating one feature.
* Preserve the numbered-menu structure of `py4siesta`.
* Preserve the existing commands and behavior of `py4siesta-tool`.
* Keep py4siesta-agent functionality separate from core scientific and deterministic workflow functionality.
* Never duplicate NanoCore, py4siesta, or py4siesta-tool functionality inside py4siesta-agent.
* Never introduce agent-framework or LLM dependencies into NanoCore, py4siesta, or py4siesta-tool.
* Menu-number changes must be limited to menu assignment and documentation consistency.
