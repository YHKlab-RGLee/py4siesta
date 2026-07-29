"""Deterministic scheduler operations and agent workflow resource accounting."""

import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path


class SchedulerError(RuntimeError):
    pass


class SlurmBackend:
    name = "slurm"

    def submit(self, case_directory, script_path):
        result = subprocess.run(
            ["sbatch", str(script_path)],
            cwd=str(case_directory),
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        match = re.search(r"Submitted batch job\s+(\d+)", result.stdout)
        if not match:
            raise SchedulerError("Could not parse a Slurm job ID from sbatch output.")
        return match.group(1)

    def status(self, job_id):
        result = subprocess.run(
            ["squeue", "-h", "-j", str(job_id), "-o", "%T"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        state = result.stdout.strip().splitlines()
        if not state:
            result = subprocess.run(
                ["sacct", "-n", "-X", "-j", str(job_id), "--format=State"],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            state = result.stdout.strip().splitlines()
        return self.normalize_status(state[0] if state else "UNKNOWN")

    def cancel(self, job_id):
        subprocess.run(["scancel", str(job_id)], check=True)

    @staticmethod
    def normalize_status(value):
        state = value.strip().upper().split()[0].split("+")[0]
        if state in {"PENDING", "CONFIGURING", "REQUEUED"}:
            return "queued"
        if state in {"RUNNING", "COMPLETING"}:
            return "running"
        if state == "COMPLETED":
            return "completed"
        if state in {"CANCELLED", "PREEMPTED"}:
            return "cancelled"
        return "failed"


def parse_requested_nodes(text, script_path="<scheduler script>"):
    patterns = [
        r"^\s*#SBATCH\s+--nodes(?:=|\s+)(\d+)\s*(?:#.*)?$",
        r"^\s*#SBATCH\s+-N(?:=|\s+)(\d+)\s*(?:#.*)?$",
    ]
    for line in text.splitlines():
        for pattern in patterns:
            match = re.match(pattern, line)
            if match:
                value = int(match.group(1))
                if value > 0:
                    return value
    raise SchedulerError(
        "Cannot parse requested nodes from %s; add '#SBATCH --nodes=<count>'."
        % script_path
    )


def validate_scheduler_script(script_path, backend):
    path = Path(script_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError("Scheduler script does not exist: %s" % path)
    if not os.access(str(path), os.R_OK):
        raise SchedulerError("Scheduler script is not readable: %s" % path)
    text = path.read_text()
    if backend.name == "slurm" and "#SBATCH" not in text:
        raise SchedulerError(
            "Scheduler script %s is incompatible with the configured Slurm backend."
            % path
        )
    return {
        "path": str(path),
        "requested_nodes": parse_requested_nodes(text, path),
        "backend": backend.name,
    }


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
                if job["status"] in TERMINAL_STATES:
                    job["completed_at"] = datetime.now(timezone.utc).isoformat()
            except Exception as exc:
                job["status"] = "failed"
                job["completed_at"] = datetime.now(timezone.utc).isoformat()
                job["failure"] = "Scheduler status query failed: %s" % exc
        return self.submit_pending(jobs)

    def cancel(self, jobs, job_id):
        for job in jobs:
            if job.get("scheduler_job_id") == str(job_id):
                self.backend.cancel(job_id)
                job["status"] = "cancelled"
                return job
        raise KeyError("Managed scheduler job was not found: %s" % job_id)


__all__ = ["SchedulerError", "SchedulerManager", "SlurmBackend"]
