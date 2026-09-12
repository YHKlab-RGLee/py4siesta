"""Deterministic structure-database queries used by the DFT setup agent."""

import json
import math
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


def reduced_composition(formula):
    """Return an order-independent reduced composition for simple formulas."""

    compact = str(formula).replace(" ", "")
    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", compact)
    if not tokens or "".join(symbol + count for symbol, count in tokens) != compact:
        return (compact.lower(),)
    counts = {}
    for symbol, count in tokens:
        counts[symbol] = counts.get(symbol, 0) + int(count or "1")
    divisor = 0
    for count in counts.values():
        divisor = math.gcd(divisor, count)
    return tuple(
        (symbol, count // divisor) for symbol, count in sorted(counts.items())
    )


def same_composition(first, second):
    return reduced_composition(first) == reduced_composition(second)


@dataclass
class StructureCandidate:
    identifier: str
    source: str
    formula: str
    structure_path: str
    dimensionality: str
    structure_type: str
    stable: Optional[bool] = None
    energy_above_hull: Optional[float] = None
    metadata: Dict[str, Any] = None

    def to_dict(self):
        value = asdict(self)
        value["metadata"] = value["metadata"] or {}
        return value


def _write_fdf(path, symbols, positions, cell):
    """Write database coordinates in the STRUCT.fdf format NanoCore consumes."""

    species = []
    for symbol in symbols:
        if symbol not in species:
            species.append(symbol)
    species_index = {symbol: index + 1 for index, symbol in enumerate(species)}

    try:
        from NanoCore.atomic_data import atomic_number

        numbers = [atomic_number(symbol) for symbol in species]
    except (KeyError, TypeError) as exc:
        raise ValueError("Database structure contains an unsupported element.") from exc

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write("NumberOfAtoms    %d\n" % len(symbols))
        stream.write("NumberOfSpecies  %d\n\n" % len(species))
        stream.write("%block ChemicalSpeciesLabel\n")
        for index, (number, symbol) in enumerate(zip(numbers, species), start=1):
            stream.write(" %d %d %s\n" % (index, number, symbol))
        stream.write("%endblock ChemicalSpeciesLabel\n\n")
        stream.write("LatticeConstant 1.0 Ang\n%block LatticeVectors\n")
        for vector in cell:
            stream.write(" %.12f %.12f %.12f\n" % tuple(vector))
        stream.write("%endblock LatticeVectors\n\n")
        stream.write("AtomicCoordinatesFormat Ang\n")
        stream.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        for serial, (symbol, position) in enumerate(zip(symbols, positions), start=1):
            stream.write(
                " %.12f %.12f %.12f %d %d\n"
                % (*position, species_index[symbol], serial)
            )
        stream.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
    return path


class StructureSource:
    """Query C2DB or Materials Project without model-generated database calls."""

    C2DB_OPTIMADE_URL = "https://c2db.fysik.dtu.dk/optimade/v1/structures"

    def __init__(
        self,
        output_directory,
        c2db_endpoint=None,
        materials_project_api_key=None,
        urlopen=None,
    ):
        self.structure_dir = Path(output_directory)
        self.c2db_endpoint = c2db_endpoint or self.C2DB_OPTIMADE_URL
        self.materials_project_api_key = materials_project_api_key
        self.urlopen = urlopen or urllib.request.urlopen

    def query(self, material, dimensionality, structure_type):
        if dimensionality == "2D" or structure_type == "monolayer":
            candidates = self._query_c2db(material)
            source = "C2DB"
        elif dimensionality == "3D" or structure_type == "bulk":
            candidates = self._query_materials_project(material)
            source = "Materials Project"
        else:
            raise ValueError(
                "A non-local structure request must identify a 2D/monolayer or 3D/bulk form."
            )
        filtered, decisions = self._filter(
            candidates, material, dimensionality, structure_type
        )
        self.last_query_record = {
            "source": source,
            "retrieved_candidate_ids": [
                candidate.identifier
                if isinstance(candidate, StructureCandidate)
                else candidate["identifier"]
                for candidate in candidates
            ],
            "deterministic_filtering": decisions,
        }
        return filtered

    def local(self, path, material, dimensionality, structure_type):
        structure = Path(path).expanduser().resolve()
        if not structure.is_file():
            raise FileNotFoundError("Local structure does not exist: %s" % structure)
        candidates = [
            StructureCandidate(
                identifier=structure.name,
                source="local",
                formula=material,
                structure_path=str(structure),
                dimensionality=dimensionality,
                structure_type=structure_type,
                metadata={"provided_by_user": True},
            ).to_dict()
        ]
        self.last_query_record = {
            "source": "local",
            "retrieved_candidate_ids": [structure.name],
            "deterministic_filtering": [
                {"identifier": structure.name, "accepted": True, "reasons": []}
            ],
        }
        return candidates

    @staticmethod
    def _filter(candidates, material, dimensionality, structure_type):
        filtered = []
        decisions = []
        for candidate in candidates:
            value = candidate.to_dict() if isinstance(candidate, StructureCandidate) else dict(candidate)
            reasons = []
            if not same_composition(value["formula"], material):
                reasons.append("formula mismatch")
            if dimensionality == "2D" and value.get("dimensionality") != "2D":
                reasons.append("dimensionality mismatch")
            if structure_type == "monolayer" and value.get("structure_type") != "monolayer":
                reasons.append("structure type mismatch")
            if value.get("stable") is False:
                reasons.append("marked unstable")
            accepted = not reasons
            decisions.append(
                {
                    "identifier": value["identifier"],
                    "accepted": accepted,
                    "reasons": reasons,
                }
            )
            if accepted:
                filtered.append(value)
        if not filtered:
            raise FileNotFoundError(
                "No compatible %s structure candidate was found for %s."
                % (structure_type, material)
            )
        return filtered, decisions

    def _query_c2db(self, material):
        endpoint = self.c2db_endpoint.rstrip("/")
        query = urllib.parse.urlencode(
            {"filter": 'chemical_formula_reduced="%s"' % material}
        )
        url = "%s?%s" % (endpoint, query)
        try:
            with self.urlopen(url, timeout=30) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except (OSError, urllib.error.URLError, json.JSONDecodeError) as exc:
            raise RuntimeError(
                "C2DB OPTIMADE request failed for %s: %s" % (url, exc)
            ) from exc

        records = payload.get("data")
        if not isinstance(records, list):
            raise RuntimeError("C2DB OPTIMADE response does not contain a data list.")

        candidates = []
        for record in records:
            attributes = record.get("attributes", {})
            formula = (
                attributes.get("chemical_formula_reduced")
                or attributes.get("chemical_formula_hill")
            )
            if not formula or not same_composition(formula, material):
                continue
            symbols = attributes.get("species_at_sites")
            positions = attributes.get("cartesian_site_positions")
            cell = attributes.get("lattice_vectors")
            dimensions = attributes.get("dimension_types")
            if (
                not isinstance(symbols, list)
                or not isinstance(positions, list)
                or len(symbols) != len(positions)
                or not isinstance(cell, list)
                or len(cell) != 3
                or dimensions != [1, 1, 0]
            ):
                continue
            identifier = str(record["id"])
            output = self.structure_dir / ("c2db-%s.fdf" % identifier)
            _write_fdf(output, symbols, positions, cell)
            candidates.append(
                StructureCandidate(
                    identifier=identifier,
                    source="C2DB",
                    formula=str(formula),
                    structure_path=str(output.resolve()),
                    dimensionality="2D",
                    structure_type="monolayer",
                    metadata={
                        "provider": payload.get("meta", {})
                        .get("provider", {})
                        .get("name", "C2DB"),
                        "optimade_id": identifier,
                        "optimade_url": endpoint,
                        "nperiodic_dimensions": attributes.get(
                            "nperiodic_dimensions"
                        ),
                        "nsites": attributes.get("nsites"),
                    },
                )
            )
        return candidates

    def _query_materials_project(self, material):
        try:
            from mp_api.client import MPRester
        except ImportError as exc:
            raise RuntimeError(
                "Materials Project access requires the optional mp-api package."
            ) from exc
        if not self.materials_project_api_key:
            raise RuntimeError("MP_API_KEY is required for Materials Project access.")

        candidates = []
        with MPRester(self.materials_project_api_key) as rester:
            documents = rester.materials.summary.search(
                formula=material,
                fields=[
                    "material_id",
                    "formula_pretty",
                    "structure",
                    "energy_above_hull",
                    "is_stable",
                ],
            )
        for document in documents:
            identifier = str(document.material_id)
            structure = document.structure
            output = self.structure_dir / ("%s.fdf" % identifier)
            _write_fdf(
                output,
                [str(site.specie.symbol) for site in structure.sites],
                [list(site.coords) for site in structure.sites],
                structure.lattice.matrix,
            )
            candidates.append(
                StructureCandidate(
                    identifier=identifier,
                    source="Materials Project",
                    formula=str(document.formula_pretty),
                    structure_path=str(output.resolve()),
                    dimensionality="3D",
                    structure_type="bulk",
                    stable=document.is_stable,
                    energy_above_hull=document.energy_above_hull,
                    metadata={},
                )
            )
        return candidates
