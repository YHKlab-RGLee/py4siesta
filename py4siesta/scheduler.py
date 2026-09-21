"""Shared Slurm operations for explicit jobs, legacy menus, and agent adapters."""

import os
import re
import subprocess
from pathlib import Path


class SchedulerError(RuntimeError):
    """Scheduler request failed; details preserve uncertain outcomes and evidence."""

    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = details or {}


def _job_id(value):
    value = str(value)
    if not re.fullmatch(r"[1-9][0-9]*(?:_[0-9]+)?", value):
        raise ValueError("job_id must be a positive Slurm job ID or an array task ID (123_4).")
    return value


def _cluster_args(cluster):
    if cluster is None:
        return []
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", cluster):
        raise ValueError("cluster must be a single Slurm cluster name.")
    return ["--clusters", cluster]


def _run(command, cwd=None, capture=True):
    options = dict(check=True)
    if cwd is not None:
        options["cwd"] = str(cwd)
    if capture:
        options.update(text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30)
    return subprocess.run(command, **options)


def _submit_script(script_path, case_directory=None, capture=False, parsable=False):
    # Legacy callers retain inherited output and their original script argument.
    command = ["sbatch"] + (["--parsable"] if parsable else []) + [str(script_path)]
    return _run(command, cwd=case_directory, capture=capture)


def _request(call, details):
    try:
        return call()
    except (subprocess.SubprocessError, OSError) as exc:
        def text(value):
            return value.decode(errors="replace") if isinstance(value, bytes) else (value or "")
        evidence = dict(details, stdout=text(getattr(exc, "stdout", None)),
                        stderr=text(getattr(exc, "stderr", None)))
        if isinstance(exc, subprocess.CalledProcessError):
            evidence["returncode"] = exc.returncode
        raise SchedulerError("Slurm request failed: %s" % exc, evidence) from exc


def submit_job(case_directory: str, script_path: str):
    """Submit exactly one script with sbatch --parsable; return job ID and paths.

    Resolve relative script_path against case_directory. No origin directory,
    filename prefix, or explicit node-count directive is required. The user owns
    script contents and cluster configuration. A 30-second timeout or unparseable
    response is uncertain: inspect scheduler records before any resubmission.
    """
    case = Path(case_directory).expanduser().resolve()
    if not case.is_dir():
        raise NotADirectoryError("Calculation directory does not exist: %s" % case)
    script = Path(script_path).expanduser()
    if not script.is_absolute():
        script = case / script
    script = script.resolve()
    if not script.is_file():
        raise FileNotFoundError("Scheduler script does not exist: %s" % script)
    if not os.access(str(script), os.R_OK):
        raise PermissionError("Scheduler script is not readable: %s" % script)
    paths = dict(case_directory=str(case), script_path=str(script))
    uncertain = dict(paths, status="submission_unknown", retry_safe=False)
    result = _request(lambda: _submit_script(script, case, capture=True, parsable=True), uncertain)
    match = re.fullmatch(r"([1-9][0-9]*)(?:;([A-Za-z0-9][A-Za-z0-9_.-]*))?", result.stdout.strip())
    if not match:
        raise SchedulerError("Slurm submission response has no unambiguous job ID; do not resubmit.",
                             dict(uncertain, stdout=result.stdout, stderr=result.stderr))
    return dict(paths, job_id=match.group(1), cluster=match.group(2), status="submitted",
                stdout=result.stdout, stderr=result.stderr)


def _normalize_status(value):
    words = value.strip().upper().split()
    state = words[0].split("+")[0] if words else "UNKNOWN"
    if state in {"PENDING", "CONFIGURING", "REQUEUED", "REQUEUE_FED", "REQUEUE_HOLD"}:
        return "queued"
    if state in {"RUNNING", "COMPLETING", "SUSPENDED", "RESIZING", "SIGNALING", "STAGE_OUT"}:
        return "running"
    if state == "COMPLETED":
        return "completed"
    if state in {"CANCELLED", "PREEMPTED"}:
        return "cancelled"
    if state in {"FAILED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "BOOT_FAIL", "DEADLINE", "REVOKED"}:
        return "failed"
    return "unknown"


def job_status(job_id: str, cluster: str = None):
    """Query squeue, then sacct for one allocation/task. Never submit or retry jobs.

    Return raw_state, normalized status, and exit_code (Slurm code:signal, or null).
    Missing/ambiguous records are unknown, not failed. Query errors raise
    SchedulerError with status=query_error. Array-parent aggregation is not provided.
    Scheduler completion does not establish SIESTA convergence.
    """
    job_id = _job_id(job_id)
    cluster_args = _cluster_args(cluster)
    base = dict(job_id=job_id, cluster=cluster)
    errors = dict(base, status="query_error")
    queue_error = None
    try:
        queue = _request(lambda: _run(["squeue", "-h", "-j", job_id, "-o", "%i|%T"] + cluster_args), errors)
    except SchedulerError as exc:
        # Some Slurm versions reject IDs already removed from the live queue.
        # Accounting can still provide their completed allocation record.
        queue_error = exc
        queue = subprocess.CompletedProcess([], 0, stdout="", stderr="")
    queue_rows = [line.strip().split("|") for line in queue.stdout.splitlines()
                  if line.strip() and not line.strip().startswith("CLUSTER:")]
    if len(queue_rows) == 1 and len(queue_rows[0]) == 2 and queue_rows[0][0].strip() == job_id:
        raw = queue_rows[0][1].strip()
        return dict(base, status=_normalize_status(raw), raw_state=raw, exit_code=None,
                    source="squeue")
    if queue_rows:
        return dict(base, status="unknown", raw_state=None, exit_code=None, source="squeue",
                    reason="No unique matching allocation/task; array parents are not aggregated.",
                    stdout=queue.stdout, stderr=queue.stderr)
    accounting = _request(lambda: _run([
        "sacct", "-n", "-X", "-P", "-j", job_id, "--format=JobID%64,State%40,ExitCode"
    ] + cluster_args), errors)
    rows = [line.strip().split("|") for line in accounting.stdout.splitlines() if line.strip()]
    matches = [row for row in rows if len(row) == 3 and row[0].strip() == job_id]
    if len(matches) != 1:
        if queue_error is not None:
            raise queue_error
        return dict(base, status="unknown", raw_state=None, exit_code=None, source="sacct",
                    reason="No unique matching accounting record.", stdout=accounting.stdout,
                    stderr=accounting.stderr)
    raw, code = matches[0][1].strip(), matches[0][2].strip()
    return dict(base, status=_normalize_status(raw), raw_state=raw, exit_code=code or None,
                source="sacct")


def cancel_job(job_id: str, cluster: str = None):
    """Request scancel for one job/task; cancellation must be confirmed by status."""
    job_id = _job_id(job_id)
    command = ["scancel"] + _cluster_args(cluster) + [job_id]
    base = dict(job_id=job_id, cluster=cluster)
    result = _request(lambda: _run(command), dict(base, status="cancellation_unknown"))
    return dict(base, status="cancel_requested", stdout=result.stdout, stderr=result.stderr)


class SlurmBackend:
    """Compatibility adapter for the existing agent's string-valued backend API."""

    name = "slurm"

    def submit(self, case_directory, script_path):
        # Preserve the legacy string result and standard sbatch output contract.
        result = _submit_script(script_path, case_directory, capture=True)
        match = re.search(r"Submitted batch job\s+(\d+)", result.stdout)
        if not match:
            raise SchedulerError("Could not parse a Slurm job ID from sbatch output.",
                                 dict(status="submission_unknown", retry_safe=False,
                                      stdout=result.stdout, stderr=result.stderr))
        return match.group(1)

    def status(self, job_id):
        result = job_status(job_id)
        if result["status"] == "unknown":
            raise SchedulerError("Slurm job status is unknown.", result)
        return result["status"]

    def cancel(self, job_id):
        cancel_job(job_id)

    normalize_status = staticmethod(_normalize_status)


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
