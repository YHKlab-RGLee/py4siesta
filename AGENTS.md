## Project Overview

`py4siesta` is a command-line interface utility for preparing, organizing, and managing SIESTA calculation cases.

All user-facing functionality is accessed through numbered menu entries in the CLI. Each feature should be implemented as a clearly separated functionality and exposed through the menu system using a menu number.

## Critical File Protection Rule

This `AGENTS.md` file must never be modified by agents.

Do not edit, rewrite, reformat, rename, move, delete, or automatically update this file under any circumstance. Any changes to this file will be made manually by the project owner.

## Design Rules

1. Preserve the menu-driven CLI design.

   * All functionality must remain accessible through numbered menu options.
   * New functionality should be added as a new menu entry.
   * Existing user workflows should remain stable unless the project owner explicitly requests otherwise.

2. Keep functionality isolated.

   * When adding or updating a specific functionality, do not modify unrelated features.
   * Avoid broad refactoring unless it is explicitly required for the requested change.
   * Do not change the behavior of existing functions while implementing a new feature.

3. Maintain backward compatibility.

   * Existing menu options and workflows should continue to work.
   * If menu numbering must be updated, only adjust the menu-number assignment and related display text.
   * Do not remove, rename, or alter existing functionality unless explicitly instructed.

4. Minimize the scope of changes.

   * Implement only the requested functionality.
   * Avoid unrelated cleanup, formatting changes, style changes, or dependency changes.
   * Do not rewrite large sections of code when a small, localized update is sufficient.

5. Keep the project structure organized.

   * New code should follow the existing project layout and naming conventions.
   * Feature-specific logic should be placed near related functionality.
   * Avoid mixing unrelated logic in the same function or module.

## Update Rules

1. Do not update `AGENTS.md`.

   * This file is owned and maintained manually by the project owner.
   * Agents must treat this file as read-only.

2. When adding a new CLI feature:

   * Add the feature as a new numbered menu option.
   * Update menu numbering only as needed.
   * Ensure the new option does not interfere with existing menu behavior.
   * Do not modify unrelated functionality.

3. When updating an existing feature:

   * Restrict changes to the requested feature only.
   * Do not change other CLI options, helper functions, file formats, or workflows unless strictly necessary.
   * If a shared utility must be changed, confirm that existing behavior remains compatible.

4. When updating documentation:

   * Do not rewrite the entire `README.md`.
   * Add only the documentation necessary to describe the new or updated functionality.
   * Preserve existing README structure, wording, and sections whenever possible.
   * If a new feature is added, document the corresponding menu option and a short usage description.

5. When changing menu numbers:

   * Only update the menu display and the corresponding selection mapping.
   * Do not change the underlying functionality.
   * Ensure all menu references in the README remain consistent.

## README Update Policy

The README should be updated incrementally.

When a new feature is added, only add a concise description of that feature, including:

* the menu number,
* the feature name,
* the purpose of the feature,
* any required inputs,
* any generated outputs or side effects.

Do not rewrite the full README unless explicitly requested by the project owner.

## Expected Agent Behavior

Before making changes, agents should identify:

* which specific functionality is being added or updated,
* which files are directly relevant,
* whether menu numbering needs to be adjusted,
* whether the README requires a small incremental update.

After making changes, agents should verify:

* the CLI menu still works,
* existing menu options still behave as before,
* the new or updated functionality is reachable by menu number,
* unrelated files were not modified,
* `AGENTS.md` was not changed.

## Agent Skills Development

* Distribute reusable agent skills within `Skills/Codex/` in this repository, suitable for global registration by users.
* Base skill workflows on existing numbered CLI menu functionality. Document menu mappings, required inputs, execution steps, and output validation.
* Use existing py4siesta interfaces, preferring `py4siesta-tool` for non-interactive execution. Verify the executable or source path and Python environment before use; skill registration does not install py4siesta.
* Keep skill development isolated from existing code, menu behavior, and project layout outside `Skills/Codex/`.
* Keep personal history, memory, and adapted workflows outside the public repository. Incorporate only verified lessons into workflows, recording their scope, applicable version, and supporting evidence.

## Non-Negotiable Constraints

* Never modify `AGENTS.md`.
* Never rewrite the entire README for a small feature update.
* Never change unrelated functionality when adding or updating one feature.
* Preserve the numbered-menu CLI structure.
* Menu-number changes must be limited to menu assignment and documentation consistency.
