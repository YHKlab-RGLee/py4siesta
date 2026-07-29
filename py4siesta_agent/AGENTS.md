# py4siesta-agent Development Instructions

## Scope

These instructions apply only to files under `py4siesta_agent/`.

Implement all AI-agent-specific functionality inside this package. Do not place
agent prompts, model clients, routing logic, workflow state, tool adapters, or
agent-framework dependencies in `NanoCore`, `py4siesta`, or `py4siesta-tool`.

Do not modify the project-level `AGENTS.md`.

## Package Responsibility

`py4siesta_agent` provides the AI-agent-based interface for the project. It
orchestrates existing scientific and deterministic functionality exposed by
`NanoCore`, `py4siesta`, and `py4siesta-tool`.

The agent package must not reimplement calculation, structure-processing,
file-management, or other deterministic domain functionality that already
exists in a lower layer. Reuse such functionality through stable, clearly
defined public interfaces.

Keep the non-agent interfaces independently usable. Agent code must not become
a dependency of `NanoCore`, `py4siesta`, or `py4siesta-tool`.

## Missing Lower-Layer Functionality

Before implementing a new agent or workflow, verify that every required
scientific and deterministic operation is already available through a suitable
public interface in `py4siesta` or `py4siesta-tool`.

If the available code or public interfaces are insufficient, stop the agent
implementation. Do not work around the limitation, reproduce the missing
behavior in `py4siesta_agent`, or modify a lower layer as an unrequested part
of the agent task.

Report the gap to the user and recommend the specific functionality that must
first be added to `py4siesta` or exposed through `py4siesta-tool`. The
recommendation should identify:

* the missing deterministic operation;
* the expected inputs, outputs, and failure behavior;
* the appropriate lower-layer module or public interface, when identifiable;
* why the agent requires the operation;
* how the future agent should call it once it is available.

Resume agent implementation only after the required lower-layer functionality
has been implemented or the user explicitly expands the task to include it.

## Initial Architecture

Use a router-based execution flow:

```text
User request
    ↓
AgentRouter
    ↓
Selected specialized agent
    ↓
Tool or workflow execution
    ↓
Result
```

`AgentRouter` is responsible only for interpreting a request, classifying it,
and delegating it to the appropriate specialized agent. It must not perform
domain-specific work.

Each specialized agent owns one well-defined responsibility. For example,
`DFTSetupAgent` handles requests for preparing an initial DFT calculation by
calling existing `py4siesta` or `py4siesta-tool` functionality through explicit
interfaces.

Keep routing logic explicit and easy to inspect. Prefer a small, direct workflow
for the first implementation. Do not introduce a complex multi-agent framework
when ordinary Python control flow is sufficient.

## Modularity and Package Structure

Add specialized agents incrementally as separate subpackages under
`py4siesta_agent/agents/`. Keep each agent independently understandable,
testable, and maintainable.

The following is a guideline for the initial organization:

```text
py4siesta_agent/
├── AGENTS.md
├── __init__.py
├── router/
│   ├── __init__.py
│   └── agent_router.py
├── agents/
│   ├── __init__.py
│   └── dft_setup/
│       ├── __init__.py
│       └── dft_setup_agent.py
└── shared/
    ├── __init__.py
    └── types.py
```

Treat this layout as an initial direction, not a rigid final architecture.
Create only the modules and directories required by implemented functionality.
Do not add speculative abstraction layers, base classes, registries,
dependency-injection mechanisms, or placeholders for agents that do not yet
exist.

Keep feature-specific logic close to the specialized agent that owns it. Put
shared code in `shared/` only after there is a demonstrated need for reuse.

## Naming Conventions

Use conventional names that clearly describe each component's role:

```text
Concept or class name: AgentRouter
Python module name:    agent_router.py
Instance name:         agent_router

Concept or class name: DFTSetupAgent
Python module name:    dft_setup_agent.py
Instance name:         dft_setup_agent
```

Use the same class/module/instance pattern for future agents and major
components. Avoid vague names such as `manager`, `handler`, or `processor` when
a domain-specific name is available.

## Implementation Principles

1. Keep each specialized agent focused on one well-defined responsibility.
2. Separate language-model decisions from deterministic Python operations.
3. Keep prompts concerned with interpretation and decisions, not duplicated
   scientific or workflow logic.
4. Reuse existing calculation and file-management functions through stable
   public interfaces.
5. Keep agent tool adapters thin; adapters translate inputs and outputs but do
   not implement scientific operations.
6. Make agent inputs, outputs, validation rules, and failure states explicit.
7. Prefer simple data structures and direct control flow until greater
   abstraction is justified by working features.
8. Preserve existing `py4siesta` and `py4siesta-tool` behavior and interfaces.
9. Restrict changes to the agent feature currently being implemented.
10. Do not add speculative support for agents or workflows that have not been
    requested.

LangChain or LangGraph may be used when a concrete implementation benefits
from their capabilities. Do not add either framework merely to characterize
the package as an agent system. Keep agent-framework and language-model
dependencies confined to `py4siesta_agent/`.

## Testing

Add tests as each router, specialized agent, tool adapter, or workflow is
implemented.

Test language-model-dependent decisions separately from deterministic
execution. Routing tests should verify request classification and delegation
without running domain workflows. Deterministic workflow tests should verify
the underlying execution independently of a language model.

Tests should cover explicit inputs, outputs, unsupported requests, validation
errors, tool failures, and other defined failure states. Prefer deterministic
fakes or stubs at model boundaries so the test suite does not require a live
model service.

Verify that agent changes do not alter the existing behavior of the
`py4siesta` menu interface or the `py4siesta-tool` command-line interface.

## Developer Documentation

Whenever a new specialized agent or major agent-system component is
implemented, add corresponding developer documentation under `docs/`.

Documentation must describe the actual implementation in enough detail for a
developer learning agent-system development. Include:

* the component's purpose and responsibility;
* its directory and module structure;
* the execution flow from user request to result;
* its main classes, functions, and data structures;
* how the language model interacts with deterministic Python code;
* which existing `py4siesta` or `py4siesta-tool` interfaces are reused;
* the reasons for important architectural decisions;
* how to run and test the implementation;
* known limitations and appropriate next steps.

Do not provide only a high-level feature summary. Keep documentation aligned
with implemented behavior and do not document speculative components as if
they exist.

## Change Checklist

Before implementing an agent feature:

* identify the specialized agent or router responsibility being added;
* identify existing public functionality that can be reused;
* confirm that all required deterministic operations are available through
  suitable public interfaces, and stop with a recommendation if they are not;
* determine the smallest package structure required;
* define inputs, outputs, and failure states;
* plan separate routing and deterministic-execution tests;
* determine what developer documentation is required.

After implementing an agent feature:

* verify that routing contains no domain-specific execution logic;
* verify that deterministic domain functionality was reused rather than
  duplicated;
* verify that agent dependencies remain inside `py4siesta_agent/`;
* verify that existing non-agent interfaces still operate independently;
* run relevant routing, workflow, and interface regression tests;
* confirm that documentation describes the implementation accurately;
* confirm that unrelated files and the project-level `AGENTS.md` were not
  modified.
