"""Agent workflow resource accounting over shared Slurm operations."""

import os
from datetime import datetime, timezone

from py4siesta.scheduler import (
    SchedulerError, SlurmBackend, parse_requested_nodes, validate_scheduler_script,
)


TERMINAL_STATES = {"completed", "failed", "cancelled"}


class SchedulerManager:
    """Manage recorded agent jobs and enforce an aggregate workflow node budget."""

    def __init__(self, backend=None, max_total_nodes=None, environment=None):
        self.environment = os.environ if environment is None else environment
        backend_name = self.environment.get(
            "PY4SIESTA_SCHEDULER_BACKEND", "slurm"
        )
        if backend is None and backend_name != "slurm":
            raise SchedulerError("Unsupported scheduler backend: %s" % backend_name)
        self.backend = backend or SlurmBackend()
        self.max_total_nodes = int(
            max_total_nodes
            if max_total_nodes is not None
            else self.environment.get("PY4SIESTA_MAX_TOTAL_NODES", "10")
        )
        if self.max_total_nodes <= 0:
            raise SchedulerError("PY4SIESTA_MAX_TOTAL_NODES must be positive.")

    def validate_script(self, script_path):
        return validate_scheduler_script(script_path, self.backend)

    parse_requested_nodes = staticmethod(parse_requested_nodes)

    @staticmethod
    def active_nodes(jobs):
        return sum(
            int(job["requested_nodes"])
            for job in jobs
            if job.get("status") in {"queued", "running"}
        )

    def submit_pending(self, jobs):
        active = self.active_nodes(jobs)
        for job in jobs:
            if job.get("status") != "pending":
                continue
            requested = int(job["requested_nodes"])
            if active + requested > self.max_total_nodes:
                continue
            job["scheduler_job_id"] = self.backend.submit(
                job["calculation_directory"], job["scheduler_script"]
            )
            job["submitted_at"] = datetime.now(timezone.utc).isoformat()
            job["status"] = "queued"
            active += requested
        return jobs

    def update(self, jobs):
        for job in jobs:
            if job.get("status") not in {"queued", "running"}:
                continue
            try:
                job["status"] = self.backend.status(job["scheduler_job_id"])
                job.pop("query_error", None)
                if job["status"] in TERMINAL_STATES:
                    job["completed_at"] = datetime.now(timezone.utc).isoformat()
            except Exception as exc:
                job["query_error"] = "Scheduler status query failed: %s" % exc
        return self.submit_pending(jobs)

    def cancel(self, jobs, job_id):
        for job in jobs:
            if job.get("scheduler_job_id") == str(job_id):
                self.backend.cancel(job_id)
                job["status"] = "cancelled"
                return job
        raise KeyError("Managed scheduler job was not found: %s" % job_id)


__all__ = ["SchedulerError", "SchedulerManager", "SlurmBackend"]
