## Project Overview

`py4siesta` is a utility for preparing, organizing, and managing SIESTA calculation cases.

The project provides the following user-facing interfaces:

* `python py4siesta` provides the GUI/menu-based user interface.
* `python py4siesta-tool` provides direct command-line access to individual tools and deterministic functionality.
* `python py4siesta-agent` provides the AI-agent-based interface and is currently under development.

The project follows a layered architecture:

* `NanoCore` contains reusable core scientific functionality.
* `py4siesta` and `py4siesta-tool` apply NanoCore functionality through GUI/menu-based and CLI-based interfaces.
* `py4siesta-agent` orchestrates existing NanoCore and py4siesta tools and must not reimplement their underlying functionality.

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
   * If an agent requires a missing reusable operation, implement that operation first in NanoCore, py4siesta, or py4siesta-tool, as appropriate, and expose it to the agent through a thin interface.

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
