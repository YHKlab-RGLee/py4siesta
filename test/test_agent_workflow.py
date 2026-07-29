import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from py4siesta_agent.agents.dft_setup import DFTSetupAgent
from py4siesta_agent.agents.dft_setup.workflow_operations import (
    DFTWorkflowOperations,
)
from py4siesta_agent.router import AgentRouter
from py4siesta_agent.scheduler import SchedulerError, SchedulerManager
from py4siesta_agent.shared.model_client import create_model_client
from py4siesta_agent.shared.types import (
    AgentConfigurationError,
    DFTRequest,
    ModelOutputError,
    RouteDecision,
)


class FakeModel:
    def structured(self, instructions, user_input, schema):
        if schema is RouteDecision:
            if "dft" not in user_input.lower() and "res2" not in user_input.lower():
                return RouteDecision(agent="unsupported", reason="No matching capability.")
            return RouteDecision(agent="dft_setup", reason="DFT setup request.")
        if schema is DFTRequest:
            return DFTRequest(
                material="ReS2",
                dimensionality="2D",
                structure_type="monolayer",
                exchange_correlation_functional="LDA",
            )
        raise AssertionError("Unexpected model schema: %s" % schema)


class InvalidModel:
    def structured(self, instructions, user_input, schema):
        raise ModelOutputError("invalid structured output")


class FakeStructureSource:
    def query(self, material, dimensionality, structure_type):
        return [
            {
                "identifier": "c2db-ReS2",
                "source": "C2DB",
                "formula": material,
                "structure_path": "/mock/ReS2.fdf",
                "dimensionality": dimensionality,
                "structure_type": structure_type,
                "stable": True,
                "energy_above_hull": 0.0,
                "metadata": {},
            }
        ]


class FakeScientificWorkflow:
    def __init__(self, root):
        self.root = Path(root)
        self.script = self.root / "template.sh"
        self.script.write_text("#!/bin/bash\n#SBATCH --nodes=1\n")
        self.calls = []

    def resolve_scheduler_script(self, explicit=None):
        return str(self.script)

    def initial_kpoints(self, request):
        return [2, 2, 1]

    def initialize_origin(self, structure, functional, kpoints, scheduler):
        self.calls.append(("initialize_origin", structure, functional, kpoints))
        origin = self.root / "origin"
        origin.mkdir(exist_ok=True)
        return {"origin_dir": origin}

    def _cases(self, base, count):
        result = []
        for index in range(count):
            case = self.root / base / ("%02d" % index)
            case.mkdir(parents=True, exist_ok=True)
            (case / "slm_test").write_text("#SBATCH --nodes=1\n")
            result.append(str(case))
        return result

    def prepare_kpoint_cases(self, dimensionality):
        return self._cases("01.kpoint_sampling", 3)

    def prepare_lattice_cases(self, dimensionality):
        return self._cases("02.slab_eos", 2)

    def prepare_geometry_optimization(self, dimensionality):
        return self._cases("03.geometry", 1)

    @staticmethod
    def scheduler_script_for_case(case):
        return str(Path(case) / "slm_test")

    def analyze_kpoints(self):
        return {"converged_k": 6, "results": [[2, -1.0], [4, -1.1], [6, -1.1]]}

    def analyze_lattice(self, dimensionality):
        return {"optimized_lattice_A": 3.1}

    def validate_geometry_optimization(self):
        return {"converged": True, "structure_output": "/mock/ReS2.STRUCT_OUT"}

    def generate_final_input(self, geometry_result):
        destination = self.root / "final_input"
        destination.mkdir()
        return str(destination)


class FakeSchedulerBackend:
    name = "slurm"

    def __init__(self):
        self.next_id = 1
        self.states = {}
        self.submissions = []
        self.completion_cycles = 0

    def submit(self, case, script):
        job_id = str(self.next_id)
        self.next_id += 1
        self.states[job_id] = "queued"
        self.submissions.append((str(case), str(script)))
        return job_id

    def status(self, job_id):
        return self.states[job_id]

    def cancel(self, job_id):
        self.states[job_id] = "cancelled"

    def complete_all(self):
        self.completion_cycles += 1
        for job_id in self.states:
            self.states[job_id] = "completed"


class RouterAndModelTests(unittest.TestCase):
    def test_unsupported_request_is_not_forced_into_dft(self):
        agent = mock.Mock()
        agent.capability = "DFT setup"
        result = AgentRouter(FakeModel(), {"dft_setup": agent}).route(
            "plot an existing band structure"
        )
        self.assertEqual(result["status"], "unsupported")
        agent.run.assert_not_called()

    def test_invalid_structured_output_is_rejected(self):
        agent = mock.Mock()
        agent.capability = "DFT setup"
        with self.assertRaises(ModelOutputError):
            AgentRouter(InvalidModel(), {"dft_setup": agent}).route("setup dft")

    def test_missing_model_configuration_is_clear(self):
        with self.assertRaisesRegex(
            AgentConfigurationError, "PY4SIESTA_LLM_PROVIDER"
        ):
            create_model_client(environment={})


class SchedulerManagerTests(unittest.TestCase):
    def test_unparseable_nodes_are_rejected(self):
        with self.assertRaisesRegex(SchedulerError, "Cannot parse requested nodes"):
            SchedulerManager.parse_requested_nodes("#SBATCH --time=01:00:00", "run.sh")

    def test_aggregate_limit_leaves_jobs_pending(self):
        backend = FakeSchedulerBackend()
        manager = SchedulerManager(backend=backend, max_total_nodes=2)
        jobs = [
            {
                "calculation_directory": "/case/%d" % index,
                "scheduler_script": "/case/run.sh",
                "requested_nodes": 1,
                "status": "pending",
            }
            for index in range(3)
        ]
        manager.submit_pending(jobs)
        self.assertEqual([job["status"] for job in jobs], ["queued", "queued", "pending"])
        backend.complete_all()
        manager.update(jobs)
        self.assertEqual([job["status"] for job in jobs], ["completed", "completed", "queued"])

    def test_submission_uses_fixed_backend_arguments_not_model_commands(self):
        backend = FakeSchedulerBackend()
        manager = SchedulerManager(backend=backend, max_total_nodes=1)
        job = {
            "calculation_directory": "/safe/case",
            "scheduler_script": "/safe/case/run.sh",
            "requested_nodes": 1,
            "status": "pending",
            "model_command": "rm -rf /",
        }
        manager.submit_pending([job])
        self.assertEqual(backend.submissions, [("/safe/case", "/safe/case/run.sh")])

    def test_missing_scheduler_script_is_clear(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            workflow = DFTWorkflowOperations(tmpdir, environment={})
            with self.assertRaisesRegex(
                FileNotFoundError, "PY4SIESTA_SLURM_TEMPLATE"
            ):
                workflow.resolve_scheduler_script()

    def test_incomplete_geometry_output_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "03.geometry_optimization" / "OUT"
            output.mkdir(parents=True)
            workflow = DFTWorkflowOperations(tmpdir)
            with self.assertRaisesRegex(FileNotFoundError, "incomplete"):
                workflow.validate_geometry_optimization()


class PersistentWorkflowTests(unittest.TestCase):
    def test_single_command_monitors_and_finishes_full_workflow(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            backend = FakeSchedulerBackend()
            scheduler = SchedulerManager(backend=backend, max_total_nodes=2)
            science = FakeScientificWorkflow(tmpdir)
            checkpoint = Path(tmpdir) / "state.sqlite"
            agent = DFTSetupAgent(
                FakeModel(),
                calculation_root=tmpdir,
                scheduler_manager=scheduler,
                scientific_workflow=science,
                structure_source=FakeStructureSource(),
                checkpoint_path=checkpoint,
                poll_interval=0,
                sleeper=lambda _: backend.complete_all(),
            )
            result = AgentRouter(FakeModel(), {"dft_setup": agent}).route(
                "Calculate a ReS2 monolayer with LDA"
            )
            self.assertEqual(result["status"], "completed")
            self.assertEqual(
                result["state"]["parsed_request"]["exchange_correlation_functional"],
                "LDA",
            )
            self.assertEqual(result["state"]["selected_structure"]["source"], "C2DB")
            self.assertIn(("initialize_origin", "/mock/ReS2.fdf", "LDA", [2, 2, 1]), science.calls)
            self.assertEqual(result["state"]["kpoint_result"]["converged_k"], 6)
            self.assertEqual(result["state"]["lattice_result"]["optimized_lattice_A"], 3.1)
            self.assertTrue(result["state"]["geometry_result"]["converged"])
            self.assertTrue(Path(result["state"]["summary_location"]).is_file())
            self.assertTrue(Path(result["state"]["final_input_location"]).is_dir())
            self.assertTrue(checkpoint.is_file())
            self.assertGreaterEqual(backend.completion_cycles, 3)
            self.assertTrue(
                all(job["status"] == "completed" for job in result["state"]["jobs"])
            )
            agent.close()

    def test_failed_scheduler_job_stops_before_analysis(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            backend = FakeSchedulerBackend()
            science = FakeScientificWorkflow(tmpdir)
            science.analyze_kpoints = mock.Mock()

            def fail_all(_):
                for job_id in backend.states:
                    backend.states[job_id] = "failed"

            agent = DFTSetupAgent(
                FakeModel(),
                calculation_root=tmpdir,
                scheduler_manager=SchedulerManager(backend=backend, max_total_nodes=10),
                scientific_workflow=science,
                structure_source=FakeStructureSource(),
                checkpoint_path=Path(tmpdir) / "state.sqlite",
                poll_interval=0,
                sleeper=fail_all,
            )
            result = agent.run("setup dft for ReS2")
            self.assertEqual(result["status"], "failed")
            science.analyze_kpoints.assert_not_called()
            agent.close()


if __name__ == "__main__":
    unittest.main()
