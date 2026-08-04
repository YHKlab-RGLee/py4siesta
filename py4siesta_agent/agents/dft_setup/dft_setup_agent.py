"""LangGraph orchestration for the first persistent DFT setup workflow."""

import json
import os
import sqlite3
import time
import uuid
from pathlib import Path
from typing import Any, Dict, TypedDict

from langgraph.checkpoint.sqlite import SqliteSaver
from langgraph.graph import END, START, StateGraph

from py4siesta_agent.scheduler import SchedulerManager
from py4siesta_agent.structure_sources import StructureSource
from py4siesta_agent.agents.dft_setup.workflow_operations import (
    DFTWorkflowOperations,
)
from py4siesta_agent.shared.types import (
    DFTRequest,
    JobRecord,
    StructureSelectionDecision,
    WorkflowState,
)


PARSE_INSTRUCTIONS = """Extract only explicitly stated DFT setup requirements.
Use 2D and monolayer for a monolayer request, and 3D and bulk for a bulk request.
Normalize the exchange-correlation functional to LDA or GGA. Do not invent
scientific defaults. Paths and initial k-points must be null unless explicit."""


class GraphState(TypedDict, total=False):
    original_request: str
    workflow_id: str
    current_stage: str
    routing_decision: Dict[str, Any]
    parsed_request: Dict[str, Any]
    structure_query: Dict[str, Any]
    structure_candidates: list
    selected_structure: Dict[str, Any]
    calculation_root: str
    scheduler_script: str
    jobs: list
    kpoint_result: Dict[str, Any]
    lattice_result: Dict[str, Any]
    geometry_result: Dict[str, Any]
    final_input_location: str
    summary_location: str
    warnings: list
    errors: list
    status: str


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return value


class DFTSetupAgent:
    name = "DFTSetupAgent"
    capability = (
        "Prepare a SIESTA DFT setup from a material or local structure, including "
        "k-point convergence, lattice optimization, CG geometry optimization, "
        "and final input generation."
    )

    def __init__(
        self,
        model_client,
        calculation_root=".",
        scheduler_manager=None,
        scientific_workflow=None,
        structure_source=None,
        checkpoint_path=None,
        poll_interval=None,
        sleeper=None,
        environment=None,
    ):
        self.model_client = model_client
        self.environment = os.environ if environment is None else environment
        self.root = Path(calculation_root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.scheduler = scheduler_manager or SchedulerManager(
            environment=self.environment
        )
        self.science = scientific_workflow or DFTWorkflowOperations(
            self.root, environment=self.environment
        )
        self.structures = structure_source or StructureSource(
            self.root / ".py4siesta-agent" / "structures",
            c2db_endpoint=self.environment.get(
                "PY4SIESTA_C2DB_OPTIMADE_URL"
            ),
            materials_project_api_key=self.environment.get("MP_API_KEY"),
        )
        self.poll_interval = float(
            poll_interval
            if poll_interval is not None
            else self.environment.get("PY4SIESTA_POLL_INTERVAL_SECONDS", "60")
        )
        if self.poll_interval < 0:
            raise ValueError("PY4SIESTA_POLL_INTERVAL_SECONDS cannot be negative.")
        self._sleep = sleeper or time.sleep
        checkpoint_path = checkpoint_path or (
            self.root / ".py4siesta-agent" / "checkpoints.sqlite"
        )
        checkpoint_path = Path(checkpoint_path)
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
        self._connection = sqlite3.connect(str(checkpoint_path), check_same_thread=False)
        self._checkpointer = SqliteSaver(self._connection)
        self.graph = self._build_graph()

    def close(self):
        self._connection.close()

    def _config(self, workflow_id):
        return {"configurable": {"thread_id": str(workflow_id)}}

    def run(self, request):
        workflow_id = uuid.uuid4().hex
        state = WorkflowState(
            original_request=request,
            workflow_id=workflow_id,
            calculation_root=str(self.root),
            routing_decision={
                "agent": "dft_setup",
                "reason": "Selected by AgentRouter.",
            },
        ).to_dict()
        result = self.graph.invoke(state, config=self._config(workflow_id))
        while result.get("status") == "waiting":
            self._sleep(self.poll_interval)
            result = self._resume(workflow_id)
        return self._result(result)

    def _resume(self, workflow_id):
        snapshot = self.graph.get_state(self._config(workflow_id))
        if not snapshot.values:
            raise KeyError("Workflow was not found: %s" % workflow_id)
        return self.graph.invoke(
            dict(snapshot.values), config=self._config(workflow_id)
        )

    @staticmethod
    def _result(state):
        return {
            "ok": state.get("status") not in {"failed"},
            "status": state.get("status"),
            "workflow_id": state.get("workflow_id"),
            "current_stage": state.get("current_stage"),
            "selected_agent": "DFTSetupAgent",
            "state": _jsonable(state),
        }

    def _build_graph(self):
        builder = StateGraph(GraphState)
        nodes = {
            "parse_request": self._parse_request,
            "select_structure": self._select_structure,
            "initialize_origin": self._initialize_origin,
            "prepare_kpoint_cases": self._prepare_kpoint_cases,
            "submit_kpoint_jobs": self._submit_kpoint_jobs,
            "monitor_kpoint_jobs": self._monitor_kpoint_jobs,
            "analyze_kpoint_convergence": self._analyze_kpoints,
            "prepare_lattice_cases": self._prepare_lattice_cases,
            "submit_lattice_jobs": self._submit_lattice_jobs,
            "monitor_lattice_jobs": self._monitor_lattice_jobs,
            "analyze_lattice_results": self._analyze_lattice,
            "prepare_geometry_optimization": self._prepare_geometry,
            "submit_geometry_job": self._submit_geometry,
            "monitor_geometry_job": self._monitor_geometry,
            "validate_geometry_optimization": self._validate_geometry,
            "generate_final_input": self._generate_final_input,
            "write_summary": self._write_summary,
        }
        for name, function in nodes.items():
            builder.add_node(name, function)

        builder.add_conditional_edges(
            START,
            lambda state: state.get("current_stage", "parse_request"),
            {name: name for name in nodes},
        )
        sequence = [
            "parse_request",
            "select_structure",
            "initialize_origin",
            "prepare_kpoint_cases",
            "submit_kpoint_jobs",
        ]
        for current, following in zip(sequence, sequence[1:]):
            builder.add_edge(current, following)
        self._monitor_edges(
            builder,
            "monitor_kpoint_jobs",
            "analyze_kpoint_convergence",
        )
        builder.add_edge("submit_kpoint_jobs", "monitor_kpoint_jobs")
        builder.add_conditional_edges(
            "analyze_kpoint_convergence",
            self._next_or_summary,
            {
                "prepare_lattice_cases": "prepare_lattice_cases",
                "write_summary": "write_summary",
            },
        )
        builder.add_edge("prepare_lattice_cases", "submit_lattice_jobs")
        builder.add_edge("submit_lattice_jobs", "monitor_lattice_jobs")
        self._monitor_edges(
            builder, "monitor_lattice_jobs", "analyze_lattice_results"
        )
        builder.add_conditional_edges(
            "analyze_lattice_results",
            self._next_or_summary,
            {
                "prepare_geometry_optimization": "prepare_geometry_optimization",
                "write_summary": "write_summary",
            },
        )
        builder.add_edge("prepare_geometry_optimization", "submit_geometry_job")
        builder.add_edge("submit_geometry_job", "monitor_geometry_job")
        self._monitor_edges(
            builder,
            "monitor_geometry_job",
            "validate_geometry_optimization",
        )
        builder.add_conditional_edges(
            "validate_geometry_optimization",
            self._next_or_summary,
            {
                "generate_final_input": "generate_final_input",
                "write_summary": "write_summary",
            },
        )
        builder.add_edge("generate_final_input", "write_summary")
        builder.add_edge("write_summary", END)
        return builder.compile(checkpointer=self._checkpointer)

    @staticmethod
    def _monitor_edges(builder, monitor, completed_node):
        builder.add_conditional_edges(
            monitor,
            lambda state: (
                completed_node
                if state.get("current_stage") == completed_node
                else "write_summary"
                if state.get("current_stage") == "write_summary"
                else END
            ),
            {completed_node: completed_node, "write_summary": "write_summary", END: END},
        )

    @staticmethod
    def _next_or_summary(state):
        return state["current_stage"]

    def _parse_request(self, state):
        parsed = self.model_client.structured(
            PARSE_INSTRUCTIONS, state["original_request"], DFTRequest
        )
        return {
            "parsed_request": parsed.model_dump(),
            "current_stage": "select_structure",
        }

    def _select_structure(self, state):
        request = state["parsed_request"]
        if request.get("local_structure"):
            candidates = self.structures.local(
                request["local_structure"],
                request["material"],
                request["dimensionality"],
                request["structure_type"],
            )
        else:
            candidates = self.structures.query(
                request["material"],
                request["dimensionality"],
                request["structure_type"],
            )
        if len(candidates) == 1:
            selected = dict(candidates[0])
            selected["selection_reason"] = "Only compatible candidate after deterministic filtering."
        else:
            decision = self.model_client.structured(
                "Select only one identifier from the recorded compatible candidates. "
                "Do not invent metadata or identifiers.",
                json.dumps(
                    {"request": request, "candidates": candidates}, sort_keys=True
                ),
                StructureSelectionDecision,
            )
            matches = [
                candidate
                for candidate in candidates
                if candidate["identifier"] == decision.candidate_id
            ]
            if not matches:
                raise ValueError(
                    "Model selected an identifier outside the candidate set."
                )
            selected = dict(matches[0])
            selected["selection_reason"] = decision.reason
        return {
            "structure_query": getattr(
                self.structures,
                "last_query_record",
                {
                    "source": selected.get("source"),
                    "retrieved_candidate_ids": [
                        candidate["identifier"] for candidate in candidates
                    ],
                    "deterministic_filtering": [],
                },
            ),
            "structure_candidates": candidates,
            "selected_structure": selected,
            "current_stage": "initialize_origin",
        }

    def _initialize_origin(self, state):
        request = state["parsed_request"]
        scheduler_script = self.science.resolve_scheduler_script(
            request.get("scheduler_script")
        )
        self.scheduler.validate_script(scheduler_script)
        kpoints = self.science.initial_kpoints(request)
        result = self.science.initialize_origin(
            state["selected_structure"]["structure_path"],
            request["exchange_correlation_functional"],
            kpoints,
            scheduler_script,
        )
        return {
            "scheduler_script": scheduler_script,
            "warnings": state.get("warnings", [])
            + ["origin initialized at %s" % _jsonable(result["origin_dir"])],
            "current_stage": "prepare_kpoint_cases",
        }

    def _jobs_for(self, stage, case_directories):
        jobs = []
        for case in case_directories:
            script = self.science.scheduler_script_for_case(case)
            validation = self.scheduler.validate_script(script)
            jobs.append(
                JobRecord(
                    workflow_stage=stage,
                    calculation_directory=str(case),
                    scheduler_script=validation["path"],
                    requested_nodes=validation["requested_nodes"],
                    relevant_input=str(Path(case) / "input" / "RUN.fdf"),
                    relevant_output=str(Path(case) / "OUT" / "stdout.txt"),
                ).to_dict()
            )
        return jobs

    def _prepare_kpoint_cases(self, state):
        cases = self.science.prepare_kpoint_cases(
            state["parsed_request"]["dimensionality"]
        )
        return {
            "jobs": state.get("jobs", []) + self._jobs_for("kpoint", cases),
            "current_stage": "submit_kpoint_jobs",
        }

    def _submit_kpoint_jobs(self, state):
        return {
            "jobs": self.scheduler.submit_pending(state["jobs"]),
            "current_stage": "monitor_kpoint_jobs",
        }

    def _prepare_lattice_cases(self, state):
        cases = self.science.prepare_lattice_cases(
            state["parsed_request"]["dimensionality"]
        )
        return {
            "jobs": state["jobs"] + self._jobs_for("lattice", cases),
            "current_stage": "submit_lattice_jobs",
            "status": "running",
        }

    def _submit_lattice_jobs(self, state):
        return {
            "jobs": self.scheduler.submit_pending(state["jobs"]),
            "current_stage": "monitor_lattice_jobs",
        }

    def _prepare_geometry(self, state):
        cases = self.science.prepare_geometry_optimization(
            state["parsed_request"]["dimensionality"]
        )
        return {
            "jobs": state["jobs"] + self._jobs_for("geometry", cases),
            "current_stage": "submit_geometry_job",
            "status": "running",
        }

    def _submit_geometry(self, state):
        return {
            "jobs": self.scheduler.submit_pending(state["jobs"]),
            "current_stage": "monitor_geometry_job",
        }

    def _monitor(self, state, job_stage, completed_stage):
        jobs = self.scheduler.update(state["jobs"])
        relevant = [job for job in jobs if job["workflow_stage"] == job_stage]
        failed = [
            job for job in relevant if job["status"] in {"failed", "cancelled"}
        ]
        if failed:
            return {
                "jobs": jobs,
                "errors": state.get("errors", [])
                + ["%s scheduler job failed or was cancelled." % job_stage],
                "status": "failed",
                "current_stage": "write_summary",
            }
        if relevant and all(job["status"] == "completed" for job in relevant):
            return {
                "jobs": jobs,
                "status": "running",
                "current_stage": completed_stage,
            }
        return {
            "jobs": jobs,
            "status": "waiting",
            "current_stage": state["current_stage"],
        }

    def _monitor_kpoint_jobs(self, state):
        return self._monitor(
            state, "kpoint", "analyze_kpoint_convergence"
        )

    def _monitor_lattice_jobs(self, state):
        return self._monitor(state, "lattice", "analyze_lattice_results")

    def _monitor_geometry(self, state):
        return self._monitor(
            state, "geometry", "validate_geometry_optimization"
        )

    def _analyze_kpoints(self, state):
        try:
            result = _jsonable(self.science.analyze_kpoints())
            if result.get("converged_k") is None:
                raise ValueError(
                    "Configured k-point cases did not reach deterministic convergence."
                )
            return {
                "kpoint_result": result,
                "current_stage": "prepare_lattice_cases",
            }
        except Exception as exc:
            return self._analysis_failure(state, "k-point analysis", exc)

    def _analyze_lattice(self, state):
        try:
            result = _jsonable(
                self.science.analyze_lattice(
                    state["parsed_request"]["dimensionality"]
                )
            )
            return {
                "lattice_result": result,
                "current_stage": "prepare_geometry_optimization",
            }
        except Exception as exc:
            return self._analysis_failure(state, "lattice analysis", exc)

    def _validate_geometry(self, state):
        try:
            result = _jsonable(self.science.validate_geometry_optimization())
            return {
                "geometry_result": result,
                "current_stage": "generate_final_input",
            }
        except Exception as exc:
            return self._analysis_failure(
                state, "geometry optimization validation", exc
            )

    @staticmethod
    def _analysis_failure(state, stage, exc):
        return {
            "errors": state.get("errors", []) + ["%s failed: %s" % (stage, exc)],
            "status": "failed",
            "current_stage": "write_summary",
        }

    def _generate_final_input(self, state):
        try:
            location = self.science.generate_final_input(
                state["geometry_result"]
            )
            return {
                "final_input_location": str(location),
                "status": "completed",
                "current_stage": "write_summary",
            }
        except Exception as exc:
            return self._analysis_failure(state, "final input generation", exc)

    def _write_summary(self, state):
        summary_dir = self.root / ".py4siesta-agent"
        summary_dir.mkdir(parents=True, exist_ok=True)
        path = summary_dir / ("%s.md" % state["workflow_id"])
        selected = state.get("selected_structure") or {}
        lines = [
            "# py4siesta DFT setup workflow",
            "",
            "- Workflow ID: `%s`" % state["workflow_id"],
            "- Status: **%s**" % state.get("status", "unknown"),
            "- Selected agent: `DFTSetupAgent`",
            "- Original request: %s" % state["original_request"],
            "- Parsed request: `%s`"
            % json.dumps(state.get("parsed_request", {}), sort_keys=True),
            "- Structure: `%s` from %s"
            % (selected.get("identifier", "not selected"), selected.get("source", "unknown")),
            "- Structure selection reason: %s"
            % selected.get("selection_reason", "not available"),
            "- Structure query: `%s`"
            % json.dumps(state.get("structure_query", {}), sort_keys=True),
            "- K-point result: `%s`"
            % json.dumps(state.get("kpoint_result", {}), sort_keys=True),
            "- Lattice result: `%s`"
            % json.dumps(state.get("lattice_result", {}), sort_keys=True),
            "- Geometry result: `%s`"
            % json.dumps(state.get("geometry_result", {}), sort_keys=True),
            "- Final input: `%s`" % state.get("final_input_location", "not generated"),
            "",
            "## Scheduler jobs",
            "",
        ]
        for job in state.get("jobs", []):
            lines.append(
                "- %s `%s`: %s (%s nodes)"
                % (
                    job["workflow_stage"],
                    job.get("scheduler_job_id") or "pending",
                    job["status"],
                    job["requested_nodes"],
                )
            )
        lines.extend(["", "## Warnings and errors", ""])
        for warning in state.get("warnings", []):
            lines.append("- Warning: %s" % warning)
        for error in state.get("errors", []):
            lines.append("- Error: %s" % error)
        path.write_text("\n".join(lines) + "\n")
        return {
            "summary_location": str(path),
            "current_stage": "write_summary",
        }
