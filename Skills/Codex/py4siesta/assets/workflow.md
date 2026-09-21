# Workflow title

- Workflow ID / unique revision ID / date:
- Status: draft | validated | superseded
- Goal, applicability, exclusions, and menu mappings (if any):
- Parent revision / source recipe and fingerprint (if any):
- Required MCP tools and original APIs / compatible runtime conditions:
- Review status: tool schemas checked | unverified; execution checked | unverified
- Validation evidence: existing output/log paths and applicable source/version:
- Change and reason:

## Inputs and prerequisites

Define required inputs, units, default policy, scientific assumptions, required
external utilities, and calculation workdir selection. Keep personal installation
paths in external environment state. Identify unresolved requirements explicitly.

## Steps and bindings

For each step specify the MCP tool, parameters, symbolic prior-result references,
workdir, expected output, and substantive validation. Specify live object bindings
without fixed session IDs. Define allowed branches, retry and stopping conditions.
Distinguish preparation, submission, monitoring, and analysis. Mark missing tools
and untested steps. Glue only connects existing operations.

## Outputs and completion criteria

Declare artifacts and actual success criteria for the stated scope. Schema review
or successful preparation does not validate an end-to-end calculation. Limit
validated applicability to the branches/material/method/environment actually checked.

## Side effects and resumption

List overwritten paths, authorized submission scope, job-ID capture, ambiguous
submission handling, and how to reconstruct transient objects from durable inputs.
Never regenerate completed cases on resume. Record unavailable monitoring explicitly.
