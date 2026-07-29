"""Validated, provider-neutral values used by the agent boundary."""

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator


class AgentError(RuntimeError):
    """Base class for user-facing agent failures."""


class AgentConfigurationError(AgentError):
    """Required agent configuration is missing or invalid."""


class ModelOutputError(AgentError, ValueError):
    """A structured language-model response failed local validation."""


class WorkflowError(AgentError):
    """A deterministic workflow operation could not continue."""


class RouteDecision(BaseModel):
    model_config = ConfigDict(extra="forbid")

    agent: str
    reason: str = Field(min_length=1)

    @field_validator("agent")
    @classmethod
    def supported_route(cls, value):
        if value not in {"dft_setup", "unsupported"}:
            raise ValueError("must be dft_setup or unsupported")
        return value


class DFTRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    material: str = Field(min_length=1)
    dimensionality: str
    structure_type: str
    exchange_correlation_functional: str
    local_structure: Optional[str] = None
    scheduler_script: Optional[str] = None
    initial_kpoints: Optional[List[int]] = None

    @field_validator("dimensionality")
    @classmethod
    def valid_dimensionality(cls, value):
        if value not in {"2D", "3D", "other"}:
            raise ValueError("must be 2D, 3D, or other")
        return value

    @field_validator("structure_type")
    @classmethod
    def valid_structure_type(cls, value):
        if value not in {"monolayer", "bulk", "local", "other"}:
            raise ValueError("must be monolayer, bulk, local, or other")
        return value

    @field_validator("exchange_correlation_functional")
    @classmethod
    def valid_functional(cls, value):
        normalized = value.upper()
        if normalized not in {"LDA", "GGA"}:
            raise ValueError("must be LDA or GGA")
        return normalized

    @field_validator("local_structure", "scheduler_script")
    @classmethod
    def normalize_path(cls, value):
        return None if value is None else str(Path(value).expanduser())

    @field_validator("initial_kpoints")
    @classmethod
    def valid_kpoints(cls, value):
        if value is not None and (len(value) != 3 or any(v <= 0 for v in value)):
            raise ValueError("must contain exactly three positive integers")
        return value


class StructureSelectionDecision(BaseModel):
    model_config = ConfigDict(extra="forbid")

    candidate_id: str = Field(min_length=1)
    reason: str = Field(min_length=1)


@dataclass
class JobRecord:
    workflow_stage: str
    calculation_directory: str
    scheduler_script: str
    requested_nodes: int
    relevant_input: Optional[str] = None
    relevant_output: Optional[str] = None
    scheduler_job_id: Optional[str] = None
    submitted_at: Optional[str] = None
    completed_at: Optional[str] = None
    status: str = "pending"
    failure: Optional[str] = None

    def to_dict(self):
        return asdict(self)


@dataclass
class WorkflowState:
    """Serializable state persisted by the LangGraph checkpointer."""

    original_request: str
    workflow_id: str
    current_stage: str = "parse_request"
    routing_decision: Dict[str, Any] = field(default_factory=dict)
    parsed_request: Dict[str, Any] = field(default_factory=dict)
    structure_query: Dict[str, Any] = field(default_factory=dict)
    structure_candidates: List[Dict[str, Any]] = field(default_factory=list)
    selected_structure: Dict[str, Any] = field(default_factory=dict)
    calculation_root: str = "."
    scheduler_script: Optional[str] = None
    jobs: List[Dict[str, Any]] = field(default_factory=list)
    kpoint_result: Dict[str, Any] = field(default_factory=dict)
    lattice_result: Dict[str, Any] = field(default_factory=dict)
    geometry_result: Dict[str, Any] = field(default_factory=dict)
    final_input_location: Optional[str] = None
    summary_location: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    status: str = "running"

    def to_dict(self):
        return asdict(self)
