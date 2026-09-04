"""Target-independent leakage-safe splitting and CR/NCR benchmarks.

This module adds an exploratory API without changing the historical splitter.
Project-defined terminology is defined before its public keys are emitted:

``main_union_plus_criteria`` is the effective leakage-block policy.  Over the
complete release, it starts from ``main_union`` (the transitive connected
components formed by exact full CIF SHA-256, database-namespaced source
siblings, and release-authorized RAC5, MOFid-v2, and MOFid-v1 groups), adds an
edge for co-membership in every ordered user-selected explanatory criterion,
and takes connected-component closure again.  Missing criterion evidence is a
structure-specific singleton and therefore adds no shared edge.  This policy
is a partition guard, not proof of identity or common parentage.

``priority_main`` is the full-release conflict-aware explanatory hierarchy:
exact release-authorized RAC5 groups seed components, followed by exact
MOFid-v2 and then MOFid-v1 groups.  A lower group touching zero stronger
components creates one, one attaches only unresolved members, and two or more
records PARENT_METHOD_CONFLICT without merging stronger components and leaves
lower-only rows unresolved.  Remaining missing evidence becomes a
structure-specific singleton (or is explicitly excluded).  It excludes
Zeo++, CrystalNets, source ID, common name, CIF hash, provisional
source-ID/MOFid transitive identity groups, and StructureMatcher.

``RT`` aliases ``rac5_crystalnets``: exact equality of all 264 finite RAC5
binary64 values plus the complete successful current CrystalNets fingerprint.
That fingerprint includes network dimension; interpenetrated-subnet,
catenation, and subnet counts; top-level single-node/all-node nets and their
agreement; and, for every complete SingleNodes and AllNodes subnet, status,
dimension, topology key, topology name, topological genome, and per-subnet
agreement.  Incomplete input adds no edge.  ``M2T`` aliases
``mofid_v2_crystalnets``: exact equality of complete canonical
release-authorized MOFid-v2 plus that fingerprint.  In order, canonical text
converts the current non-null value to text, collapses Unicode-whitespace runs
to one ASCII space, trims, rejects an empty or declared whole-field
placeholder, applies Unicode NFKC, and applies default case-folding; this text
processing only does not change a CIF, coordinate, atom, occupancy, bond,
chemistry, topology, or unit cell.  Eligible statuses are exactly
SUCCESS, SUCCESS_TOPOLOGY_UNKNOWN, SUCCESS_TOPOLOGY_ERROR, and
SUCCESS_TOPOLOGY_TIMEOUT; every other status or incomplete input adds no edge.

``SM`` aliases ``structure_matcher_strict``: a convenience connected component
over authoritative direct pair edges whose forward and reverse
``fit(..., symmetric=True)`` both pass under Python 3.9, pymatgen 2024.2.8,
NumPy 1.26.4, and ElementComparator with ltol=stol=0.001,
angle_tol=0.01, primitive_cell=True, scale=False, attempt_supercell=True,
allow_subset=False, supercell_size=num_sites, and no ignored species.  The
direct pinned CifParser uses site/coordinate tolerances 0.0001, expands declared
symmetry, merges generated sites within tolerance, rounds coordinates near one
third/two thirds, checks occupancy, sorts the periodic Structure, preserves
disorder, and performs no manual repair, occupancy selection, atom deletion,
or chemistry edit.  Direct edges, not non-transitive components, are
authoritative; parser, timeout, OOM, matcher, asymmetric, or execution failures
are NOT_AVAILABLE rather than unmatched.  RT, M2T, and SM are reference-only
and enter neither priority_main nor the pre-existing main_union by themselves.

``representative`` is the versioned target-free diversity profile.  A complete
264-value finite RAC5 vector is used first.  Otherwise, a complete selected
Zeo++ vector uses the 13 intensive N2/He pore fields plus N2 channel and
bonded-framework dimensions.  A structure missing both remains eligible in an
explicit no-numeric-diversity tier.  Numeric tiers are median/IQR scaled with
zero-IQR fields retained using a unit divisor, never imputed; RAC5 is reduced to at
most 32 principal components; deterministic MiniBatchKMeans strata use sorted
IDs and a fixed profile seed.  The strata balance distributions only and can
never divide an effective leakage block.  The exact pinned numerical backend
is NumPy 1.26.4 plus scikit-learn 1.5.0; absence or version drift fails closed.

``fixed_pure_cr`` reserves one common approximately requested-size test set
from effective leakage blocks whose complete-release members are all strict
five-checker CR.  Strict CR means all five named checker results are available
and PASS; strict NCR means all five are available and FAIL.  The fixed test is
identical across every NCR-pool fraction and seed and has no effective-block
overlap with training or validation.  ``full_cr`` makes each cohort contain C
structures, where C is the eligible strict-CR count under the recorded cohort-
eligibility policy.  At requested NCR-pool fraction q, round-half-up(q*M) of
the M eligible strict-NCR structures replace the same number of eligible CR
structures.  Consequently a 100% endpoint uses every eligible NCR but is not
necessarily a 100%-NCR composition.  ``full_cr_diagnostic`` is the
supplementary prediction view over all CR structures; it reports exact-ID and
same-effective-block training overlap and is not the independent paper test.

Every result emitted here is exploratory and records ``official_split=false``.
Targets are neither inspected nor consumed by grouping, diversity, cohort
selection, or assignment; they may be attached only after a frozen assignment.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field, fields, is_dataclass
from decimal import Decimal, ROUND_HALF_UP
import hashlib
import json
import math
from numbers import Real
import os
from pathlib import Path
import platform
import re
import shutil
import tempfile
from types import MappingProxyType
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from CoREMOF import __version__
from CoREMOF._authority import (
    AuthorityStateError,
    IdentitySealRegistry,
    state_fingerprint,
)
from CoREMOF.labels import CHECKER_PRESETS, LABELS
from CoREMOF.parents import (
    SELECTABLE_PARENT_METHODS,
    ZEO_NUMERIC_FINGERPRINT_DEFINITION,
    ParentResolver,
    generated_output_terminology_definitions,
    parent_method_definition,
)


BENCHMARK_API_VERSION = "0.1.0"
BENCHMARK_NUMPY_VERSION = "1.26.4"
BENCHMARK_SKLEARN_VERSION = "1.5.0"
BENCHMARK_SCIPY_VERSION = "1.13.1"
BENCHMARK_JOBLIB_VERSION = "1.5.3"
BENCHMARK_THREADPOOLCTL_VERSION = "3.6.0"
DIVERSITY_PROFILE_VERSION = "coremof-representative-diversity/1.0"
DIVERSITY_PROFILE_SEED = 2602
PAIRED_LABEL_BALANCE_WEIGHT = 64.0
PAIRED_SWAP_CANDIDATE_COUNT = 96
EFFECTIVE_LEAKAGE_POLICY = "main_union_plus_criteria"
TEST_POLICY = "fixed_pure_cr"
LABEL_PURE_BLOCK_ELIGIBILITY = "complete_release_label_pure_effective_blocks"
_BENCHMARK_RESULT_FACTORY_TOKEN = object()
_BENCHMARK_RESULTS = IdentitySealRegistry("benchmark result generation")
TEST_POLICY_DEFINITION = (
    "fixed_pure_cr selects one common approximately requested-size test set "
    "by globally minimizing whole-block count deviation, preferring a count below "
    "the request on an equal-distance tie, then jointly assigning that test and a "
    "feasible exact CR-removal ladder with target-free diversity-balanced "
    "deterministic priority at profile seed 2602 from whole complete-release "
    "effective leakage blocks whose every member is strict five-checker CR. The same "
    "IDs are used for every q and seed; each selected block is consumed in full, no "
    "member of one of those blocks may enter train or validation, and requested-"
    "versus-achieved whole-block count deviation is reported."
)
EFFECTIVE_LEAKAGE_POLICY_DEFINITION = (
    "main_union_plus_criteria starts from full-release main_union, adds "
    "co-membership edges from every ordered selected criterion, and takes "
    "connected-component closure before any filter. Missing criterion "
    "evidence is a structure-specific singleton. The result is an "
    "indivisible partition guard, not identity or parent proof."
)

_PARTITIONS = ("train", "validation", "test")
_RAC_AVAILABILITY_COLUMN = "rac5_available"
_ZEO_AVAILABILITY_COLUMNS = ("n2_he_available", "periodicity_available")
_ZEO_NUMERIC_FIELDS = tuple(ZEO_NUMERIC_FINGERPRINT_DEFINITION["numeric_fields"]) + (
    "n2_channel_dimension",
    "structure_periodic_dimension",
)

_CRITERION_ALIASES = MappingProxyType(
    {
        "rt": "rac5_crystalnets",
        "rac5_topology": "rac5_crystalnets",
        "rac5_crystalnets": "rac5_crystalnets",
        "m2t": "mofid_v2_crystalnets",
        "mofid_v2_topology": "mofid_v2_crystalnets",
        "mofid_v2_crystalnets": "mofid_v2_crystalnets",
        "sm": "structure_matcher_strict",
        "structure_matcher": "structure_matcher_strict",
        "structure_matcher_strict": "structure_matcher_strict",
        **{method.casefold(): method for method in SELECTABLE_PARENT_METHODS},
    }
)

_ALIASES_BY_CRITERION = MappingProxyType(
    {
        method: tuple(
            alias
            for alias, canonical in (
                ("RT", "rac5_crystalnets"),
                ("rac5_topology", "rac5_crystalnets"),
                ("M2T", "mofid_v2_crystalnets"),
                ("mofid_v2_topology", "mofid_v2_crystalnets"),
                ("SM", "structure_matcher_strict"),
                ("structure_matcher", "structure_matcher_strict"),
            )
            if canonical == method
        )
        for method in SELECTABLE_PARENT_METHODS
    }
)


class BenchmarkError(ValueError):
    """Base error for the exploratory splitting and benchmark API."""


class BenchmarkDependencyError(BenchmarkError):
    """Raised when the pinned representative-diversity backend is absent."""


class BenchmarkFeasibilityError(BenchmarkError):
    """Raised when the fixed-size CR/NCR ladder cannot be constructed."""


def _freeze(value: object) -> object:
    if isinstance(value, Mapping):
        frozen = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise TypeError("mapping keys must be exact nonblank strings")
            frozen[key] = _freeze(item)
        return MappingProxyType(frozen)
    if isinstance(value, (list, tuple, set, frozenset)):
        return tuple(_freeze(item) for item in value)
    return value


def _jsonable(value: object) -> object:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(item) for item in value]
    return value


def _canonical_json(value: object) -> str:
    return json.dumps(
        _jsonable(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _digest(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _stable_key(*values: object) -> str:
    return hashlib.sha256(
        "\0".join(str(value) for value in values).encode("utf-8")
    ).hexdigest()


def _current_source_hashes() -> Dict[str, str]:
    package_root = Path(__file__).resolve().parent
    names = (
        "_authority.py",
        "benchmarks.py",
        "dataset.py",
        "labels.py",
        "parents.py",
        "targets.py",
        "_transactions.py",
    )
    return {
        name: hashlib.sha256((package_root / name).read_bytes()).hexdigest()
        for name in names
    }


_IMPORTED_SOURCE_HASHES = MappingProxyType(_current_source_hashes())


def _source_hashes() -> Mapping[str, str]:
    """Return import-bound source hashes and reject an on-disk code drift."""

    current = _current_source_hashes()
    if current != dict(_IMPORTED_SOURCE_HASHES):
        changed = sorted(
            name
            for name in set(current).union(_IMPORTED_SOURCE_HASHES)
            if current.get(name) != _IMPORTED_SOURCE_HASHES.get(name)
        )
        raise BenchmarkError(
            "CoREMOF implementation source changed after module import: {}".format(
                ", ".join(changed)
            )
        )
    return _IMPORTED_SOURCE_HASHES


def _seal_projection(value: object) -> object:
    """Convert frozen nested dataclass state to authority-fingerprint values."""

    if is_dataclass(value) and not isinstance(value, type):
        return (
            "dataclass",
            "{}.{}".format(type(value).__module__, type(value).__qualname__),
            tuple(
                (definition.name, _seal_projection(getattr(value, definition.name)))
                for definition in fields(value)
                if definition.name
                not in {
                    "_authority_factory_token",
                    "_authority_generation_marker",
                    "_dataset",
                }
            ),
        )
    if isinstance(value, Mapping):
        return (
            "mapping",
            tuple(
                (_seal_projection(key), _seal_projection(item))
                for key, item in value.items()
            ),
        )
    if type(value) in (tuple, list):
        return tuple(_seal_projection(item) for item in value)
    return value


def _benchmark_result_fingerprint(result: object) -> str:
    return state_fingerprint(
        ("coremof-benchmark-result-generation/1", _seal_projection(result))
    )


def _register_benchmark_result(result: object, dataset: object):
    if (
        getattr(result, "_authority_factory_token", None)
        is not _BENCHMARK_RESULT_FACTORY_TOKEN
    ):
        raise BenchmarkError("cannot register a non-factory benchmark result")
    _BENCHMARK_RESULTS.register(
        result,
        _benchmark_result_fingerprint(result),
        {
            "dataset": dataset,
            "release_binding": _jsonable(_release_binding(dataset)),
        },
    )
    return result


def _require_benchmark_result(result: object) -> None:
    try:
        context = _BENCHMARK_RESULTS.require(
            result, _benchmark_result_fingerprint(result)
        )
    except AuthorityStateError as error:
        raise BenchmarkError(str(error)) from error
    if not isinstance(context, Mapping):
        raise BenchmarkError("benchmark result authority context is malformed")
    dataset = context.get("dataset")
    if getattr(result, "_dataset", None) is not dataset:
        raise BenchmarkError("benchmark result dataset binding changed")
    if _jsonable(_release_binding(dataset)) != context.get("release_binding"):
        raise BenchmarkError("benchmark result release binding changed")
    _source_hashes()


def _release_binding(dataset: object) -> Mapping[str, object]:
    structure_ids = tuple(getattr(dataset, "structure_ids", ()))
    input_hashes = dict(getattr(dataset, "input_hashes", {}) or {})
    return MappingProxyType(
        {
            "dataset_version": getattr(dataset, "dataset_version", None),
            "structure_count": len(structure_ids),
            "ordered_structure_ids_sha256": hashlib.sha256(
                "\n".join(structure_ids).encode("utf-8")
            ).hexdigest(),
            "structure_ids_sha256": hashlib.sha256(
                "\n".join(sorted(structure_ids)).encode("utf-8")
            ).hexdigest(),
            "input_hashes": MappingProxyType(dict(sorted(input_hashes.items()))),
        }
    )


def _validate_exact_strings(values: object, name: str) -> Tuple[str, ...]:
    if isinstance(values, str):
        result = (values,)
    elif isinstance(values, (list, tuple)):
        result = tuple(values)
    else:
        raise TypeError("{} must be a string or ordered list/tuple".format(name))
    if not result:
        raise BenchmarkError("{} must not be empty".format(name))
    if any(
        not isinstance(value, str) or not value or value != value.strip()
        for value in result
    ):
        raise BenchmarkError("{} must contain exact nonblank strings".format(name))
    return result


def normalize_group_criteria(group_criteria: object) -> Tuple[str, ...]:
    """Normalize one criterion or an ordered criterion tuple.

    ``RT`` means exact complete 264-value finite RAC5 plus complete successful
    current CrystalNets equality.  ``M2T`` means exact complete canonical
    release-authorized MOFid-v2 plus that CrystalNets fingerprint and remains
    provisional with provisional MOFid input.  ``SM`` means the convenience
    connected-component view over authoritative direct symmetric matches from
    the pinned strict pymatgen protocol.  Incomplete evidence adds no match;
    all three are reference-only and enter neither ``priority_main`` nor the
    pre-existing ``main_union`` by themselves.
    """

    values = _validate_exact_strings(group_criteria, "group_criteria")
    result = []
    for value in values:
        canonical = _CRITERION_ALIASES.get(value.casefold())
        if canonical is None or canonical not in SELECTABLE_PARENT_METHODS:
            choices = sorted(set(SELECTABLE_PARENT_METHODS).union({"RT", "M2T", "SM"}))
            raise BenchmarkError(
                "unknown group criterion {!r}; choose one of {}".format(
                    value, ", ".join(choices)
                )
            )
        result.append(canonical)
    if len(set(result)) != len(result):
        raise BenchmarkError(
            "group_criteria normalizes to duplicate canonical criteria"
        )
    return tuple(result)


def _dataset_from(value: object) -> object:
    dataset = getattr(value, "dataset", value)
    # The legacy target-aware splitter intentionally accepts a
    # TargetMergedDataset.  This target-independent API instead performs every
    # grouping, feature lookup, and receipt hash lookup against its untouched
    # base release.  That keeps both target values and target input hashes out
    # of the frozen assignment by construction.
    target_columns = getattr(dataset, "target_columns", ())
    if target_columns:
        base_dataset = getattr(dataset, "base_dataset", None)
        if base_dataset is None or not hasattr(base_dataset, "metadata_rows"):
            raise BenchmarkError(
                "target-independent splitting cannot consume a dataset with "
                "target columns unless it exposes its validated base_dataset"
            )
        dataset = base_dataset
    if not hasattr(dataset, "metadata_rows"):
        raise TypeError("dataset must expose metadata_rows")
    return dataset


def _validate_classified_authority(
    classified: object, *, require_official: bool
) -> bool:
    """Revalidate package authority, requiring it for strict benchmarks.

    General exploratory splitting continues to accept the historical
    duck-typed classified input contract, but an object that claims package
    authority is always revalidated.  Strict CR/NCR pools accept only an
    authenticated, unmodified ``ClassifiedDataset`` generation produced by
    ``CoREMOFDataset.classify('5checker')``.
    """

    from CoREMOF.dataset import (
        ClassifiedDataset,
        _validate_classified_generation_if_present,
    )

    if not isinstance(classified, ClassifiedDataset):
        if require_official:
            raise BenchmarkError(
                "CR/NCR benchmarking requires an authenticated ClassifiedDataset "
                "produced from a validated release"
            )
        return False
    authenticated = _validate_classified_generation_if_present(classified)
    claimed = getattr(classified, "checker_view_official", None)
    if claimed is True and not authenticated:
        raise BenchmarkError(
            "classified checker authority could not be revalidated; rebuild the "
            "view from the validated release"
        )
    if require_official and not authenticated:
        raise BenchmarkError(
            "CR/NCR benchmarking requires the authenticated recomputed strict "
            "five-checker release view"
        )
    return authenticated


def available_group_criteria(dataset: object) -> Mapping[str, Mapping[str, object]]:
    """Report criterion requirements, aliases, and release availability.

    Availability means the release declares and validates the criterion.  It
    does not mean every structure has evidence: an unavailable row remains its
    own structure-specific singleton and never matches another null.
    """

    source = _dataset_from(dataset)
    resolver = ParentResolver(source, missing_parent="singleton")
    result: Dict[str, Mapping[str, object]] = {}
    for method in SELECTABLE_PARENT_METHODS:
        definition = parent_method_definition(method)
        reason = None
        try:
            resolver.resolve(method, missing_parent="singleton")
            available = True
        except ValueError as error:
            message = str(error)
            if (
                "not present in this release" not in message
                and "requires the exact validated" not in message
            ):
                raise
            available = False
            reason = message
        result[method] = _freeze(
            {
                "canonical_name": method,
                "display_name": definition["display_name"],
                "aliases": _ALIASES_BY_CRITERION.get(method, ()),
                "required_evidence": definition.get("input_fields", ()),
                "available": available,
                "unavailable_reason": reason,
                "missing_evidence_behavior": (
                    "structure-specific singleton; common nulls never match"
                ),
                "definition": definition,
            }
        )
    return MappingProxyType(result)


class _DisjointSet:
    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}
        self.rank = {value: 0 for value in self.parent}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        left_rank = self.rank[left_root]
        right_rank = self.rank[right_root]
        if left_rank < right_rank or (
            left_rank == right_rank and right_root < left_root
        ):
            left_root, right_root = right_root, left_root
            left_rank, right_rank = right_rank, left_rank
        self.parent[right_root] = left_root
        if left_rank == right_rank:
            self.rank[left_root] += 1


def _component_names(
    ids: Sequence[str], disjoint: _DisjointSet, prefix: str
) -> Mapping[str, str]:
    members: Dict[str, List[str]] = {}
    for structure_id in sorted(ids):
        members.setdefault(disjoint.find(structure_id), []).append(structure_id)
    names: Dict[str, str] = {}
    used: Dict[str, Tuple[str, ...]] = {}
    for root, values in sorted(members.items(), key=lambda item: tuple(item[1])):
        canonical = tuple(values)
        digest = hashlib.sha256("\0".join(canonical).encode("utf-8")).hexdigest()
        length = 16
        name = "{}-{}".format(prefix, digest[:length].upper())
        while name in used and used[name] != canonical:
            length += 1
            name = "{}-{}".format(prefix, digest[:length].upper())
        used[name] = canonical
        names[root] = name
    return MappingProxyType(
        {structure_id: names[disjoint.find(structure_id)] for structure_id in ids}
    )


@dataclass(frozen=True)
class _GroupingContext:
    criteria: Tuple[str, ...]
    criterion_groups: Mapping[str, Mapping[str, str]]
    effective_blocks: Mapping[str, str]
    criterion_diagnostics: Mapping[str, Mapping[str, str]]
    main_union_groups: Mapping[str, str]


def _grouping_context(dataset: object, criteria: Tuple[str, ...]) -> _GroupingContext:
    resolver = ParentResolver(dataset, missing_parent="singleton")
    main = resolver.resolve("main_union", missing_parent="singleton")
    ids = tuple(sorted(main.universe_ids))
    disjoint = _DisjointSet(ids)

    def add_groups(groups: Mapping[str, str]) -> None:
        members: Dict[str, List[str]] = {}
        for structure_id in ids:
            members.setdefault(groups[structure_id], []).append(structure_id)
        for group in sorted(members):
            values = members[group]
            for structure_id in values[1:]:
                disjoint.union(values[0], structure_id)

    main_groups = {structure_id: main.group_by_id[structure_id] for structure_id in ids}
    add_groups(main_groups)
    criterion_groups: Dict[str, Mapping[str, str]] = {}
    diagnostics: Dict[str, Mapping[str, str]] = {}
    for criterion in criteria:
        resolution = resolver.resolve(criterion, missing_parent="singleton")
        groups = {}
        criterion_diagnostics = {}
        for structure_id in ids:
            group = resolution.group_by_id.get(structure_id)
            if group is None:
                # A conflict in the explanatory hierarchy is not a license to
                # merge stronger parents.  For the criterion-specific view it
                # remains a singleton; main_union still supplies the broad guard.
                group = "SINGLETON:{}".format(structure_id)
                reason = resolution.exclusions.get(
                    structure_id, "MISSING_CRITERION_EVIDENCE"
                )
                criterion_diagnostics[structure_id] = reason
            groups[structure_id] = group
        add_groups(groups)
        criterion_groups[criterion] = MappingProxyType(groups)
        diagnostics[criterion] = MappingProxyType(criterion_diagnostics)
    return _GroupingContext(
        criteria=criteria,
        criterion_groups=MappingProxyType(criterion_groups),
        effective_blocks=_component_names(ids, disjoint, "LEAK"),
        criterion_diagnostics=MappingProxyType(diagnostics),
        main_union_groups=MappingProxyType(main_groups),
    )


def _load_backend():
    try:
        import joblib
        import numpy as np
        import scipy
        import sklearn
        import threadpoolctl
        from sklearn.cluster import MiniBatchKMeans
        from sklearn.decomposition import PCA
        from threadpoolctl import threadpool_info, threadpool_limits
    except (ImportError, ModuleNotFoundError) as error:
        raise BenchmarkDependencyError(
            "representative diversity requires the pinned benchmark extra: "
            "pip install 'CoREMOF-tools[benchmark]' (NumPy {}, "
            "scikit-learn {}, SciPy {}, joblib {}, and threadpoolctl {})".format(
                BENCHMARK_NUMPY_VERSION,
                BENCHMARK_SKLEARN_VERSION,
                BENCHMARK_SCIPY_VERSION,
                BENCHMARK_JOBLIB_VERSION,
                BENCHMARK_THREADPOOLCTL_VERSION,
            )
        ) from error
    expected = {
        "NumPy": BENCHMARK_NUMPY_VERSION,
        "scikit-learn": BENCHMARK_SKLEARN_VERSION,
        "SciPy": BENCHMARK_SCIPY_VERSION,
        "joblib": BENCHMARK_JOBLIB_VERSION,
        "threadpoolctl": BENCHMARK_THREADPOOLCTL_VERSION,
    }
    found = {
        "NumPy": np.__version__,
        "scikit-learn": sklearn.__version__,
        "SciPy": scipy.__version__,
        "joblib": joblib.__version__,
        "threadpoolctl": threadpoolctl.__version__,
    }
    if found != expected:
        raise BenchmarkDependencyError(
            "representative diversity requires exact backend versions {}; found {}. "
            "Install the pinned benchmark extra; version drift is not silently "
            "accepted.".format(
                ", ".join("{} {}".format(name, expected[name]) for name in expected),
                ", ".join("{} {}".format(name, found[name]) for name in found),
            )
        )
    backend = {
        "numpy": BENCHMARK_NUMPY_VERSION,
        "scikit_learn": BENCHMARK_SKLEARN_VERSION,
        "scipy": BENCHMARK_SCIPY_VERSION,
        "joblib": BENCHMARK_JOBLIB_VERSION,
        "threadpoolctl": BENCHMARK_THREADPOOLCTL_VERSION,
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "machine": platform.machine(),
        "thread_limit": 1,
        "extra": "benchmark",
        "version_drift_policy": "ERROR",
        "cross_architecture_bit_identity_guaranteed": False,
    }
    return np, PCA, MiniBatchKMeans, threadpool_limits, threadpool_info, backend


def _threadpool_runtime(threadpool_info) -> List[Mapping[str, object]]:
    """Return stable, non-path thread-library evidence for the receipt."""

    allowed = (
        "user_api",
        "internal_api",
        "prefix",
        "version",
        "threading_layer",
        "architecture",
        "num_threads",
    )
    records = []
    for raw in threadpool_info():
        records.append(
            {key: raw.get(key) for key in allowed if raw.get(key) is not None}
        )
    return sorted(records, key=lambda item: json.dumps(item, sort_keys=True))


def _truthy(value: object) -> bool:
    return isinstance(value, str) and value.strip().casefold() in {
        "true",
        "1",
        "yes",
        "success",
        "available",
    }


def _finite_vector(
    row: Mapping[str, object], fields: Sequence[str]
) -> Optional[Tuple[float, ...]]:
    values = []
    for name in fields:
        raw = row.get(name)
        if not isinstance(raw, str) or not raw.strip():
            return None
        try:
            value = float(raw)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(value):
            return None
        values.append(0.0 if value == 0.0 else value)
    return tuple(values)


def _read_release_feature_table(dataset: object, name: str):
    try:
        from CoREMOF.targets import _read_feature_table
    except ImportError as error:  # pragma: no cover - package corruption only
        raise BenchmarkError(
            "cannot import the release feature-table reader"
        ) from error
    try:
        return _read_feature_table(dataset, name)
    except Exception as error:
        raise BenchmarkError(
            "representative diversity requires the validated current {!r} "
            "feature table: {}".format(name, error)
        ) from error


@dataclass(frozen=True)
class DiversityIndex:
    """Immutable target-free feature-tier and cluster-stratum assignment."""

    strata_by_id: Mapping[str, str]
    tier_by_id: Mapping[str, str]
    topology_category_by_id: Mapping[str, str]
    profile: Mapping[str, object]
    input_receipts: Tuple[Mapping[str, object], ...]
    digest: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "strata_by_id", _freeze(self.strata_by_id))
        object.__setattr__(self, "tier_by_id", _freeze(self.tier_by_id))
        object.__setattr__(
            self, "topology_category_by_id", _freeze(self.topology_category_by_id)
        )
        object.__setattr__(self, "profile", _freeze(self.profile))
        object.__setattr__(
            self,
            "input_receipts",
            tuple(_freeze(receipt) for receipt in self.input_receipts),
        )


def _cluster_tier(
    ids: Sequence[str],
    vectors: Mapping[str, Tuple[float, ...]],
    tier: str,
    np,
    PCA,
    MiniBatchKMeans,
) -> Tuple[Mapping[str, str], Mapping[str, object]]:
    ordered = tuple(sorted(ids))
    if not ordered:
        return MappingProxyType({}), MappingProxyType(
            {"count": 0, "requested_clusters": 0, "effective_clusters": 0}
        )
    matrix = np.asarray(
        [vectors[structure_id] for structure_id in ordered], dtype=np.float64
    )
    medians = np.median(matrix, axis=0)
    q1 = np.percentile(matrix, 25.0, axis=0, method="linear")
    q3 = np.percentile(matrix, 75.0, axis=0, method="linear")
    iqr = q3 - q1
    zero_iqr = iqr == 0.0
    scales = iqr.copy()
    scales[zero_iqr] = 1.0
    scaled = (matrix - medians) / scales
    if not np.isfinite(scaled).all():
        raise BenchmarkError(
            "internal error: robust scaling produced a nonfinite value"
        )
    pca_components = None
    if tier == "rac5":
        component_count = min(32, scaled.shape[0], scaled.shape[1])
        if np.count_nonzero(scaled) == 0:
            transformed = np.zeros((scaled.shape[0], 1), dtype=np.float64)
            pca_components = 0
        else:
            transformed = PCA(
                n_components=component_count,
                svd_solver="full",
            ).fit_transform(scaled)
            pca_components = component_count
    else:
        transformed = scaled
    requested_k = min(
        len(ordered),
        256,
        max(16, int(math.ceil(math.sqrt(len(ordered))))),
    )
    if requested_k == 1 or np.count_nonzero(transformed) == 0:
        labels = np.zeros(len(ordered), dtype=int)
    else:
        model = MiniBatchKMeans(
            n_clusters=requested_k,
            random_state=DIVERSITY_PROFILE_SEED,
            n_init=10,
            batch_size=max(1024, requested_k * 3),
            reassignment_ratio=0.0,
        )
        labels = model.fit_predict(transformed)
    members: Dict[int, List[str]] = {}
    for structure_id, label in zip(ordered, labels.tolist()):
        members.setdefault(int(label), []).append(structure_id)
    canonical_labels = {
        old: index
        for index, (old, _) in enumerate(
            sorted(members.items(), key=lambda item: tuple(item[1]))
        )
    }
    strata = MappingProxyType(
        {
            structure_id: "{}:{:03d}".format(tier, canonical_labels[int(label)])
            for structure_id, label in zip(ordered, labels.tolist())
        }
    )
    details = MappingProxyType(
        {
            "count": len(ordered),
            "input_dimension": int(matrix.shape[1]),
            "zero_iqr_field_count": int(np.count_nonzero(zero_iqr)),
            "pca_components": pca_components,
            "requested_clusters": requested_k,
            "effective_clusters": len(members),
        }
    )
    return strata, details


def build_diversity_index(
    dataset: object, diversity: str = "representative"
) -> DiversityIndex:
    """Build the reusable target-free diversity index for a complete release."""

    source = _dataset_from(dataset)
    ids = tuple(sorted(str(row["structure_id"]) for row in source.metadata_rows))
    if diversity == "none":
        strata = {structure_id: "not_requested:000" for structure_id in ids}
        tiers = {structure_id: "not_requested" for structure_id in ids}
        topology = {structure_id: "not_requested" for structure_id in ids}
        profile = {
            "schema_version": DIVERSITY_PROFILE_VERSION,
            "name": "none",
            "target_columns_consumed": [],
            "scientific_feature_imputation": False,
        }
        digest = _digest(
            {"profile": profile, "strata": strata, "tiers": tiers, "topology": topology}
        )
        return DiversityIndex(strata, tiers, topology, profile, (), digest)
    if diversity != "representative":
        raise BenchmarkError("diversity must be 'representative' or 'none'")

    (
        np,
        PCA,
        MiniBatchKMeans,
        threadpool_limits,
        threadpool_info,
        backend,
    ) = _load_backend()
    rac_columns, rac_rows, rac_receipt = _read_release_feature_table(source, "rac5")
    zeo_columns, zeo_rows, zeo_receipt = _read_release_feature_table(source, "zeo")
    topology_columns, topology_rows, topology_receipt = _read_release_feature_table(
        source, "topology"
    )
    rac_fields = tuple(
        column for column in rac_columns if column != _RAC_AVAILABILITY_COLUMN
    )
    if len(rac_fields) != 264:
        raise BenchmarkError(
            "representative diversity requires exactly 264 RAC5 descriptor columns; "
            "found {}".format(len(rac_fields))
        )
    missing_zeo = set(_ZEO_AVAILABILITY_COLUMNS + _ZEO_NUMERIC_FIELDS).difference(
        zeo_columns
    )
    if missing_zeo:
        raise BenchmarkError(
            "representative diversity Zeo++ table is missing: {}".format(
                ", ".join(sorted(missing_zeo))
            )
        )
    required_topology = {
        "topology_available",
        "network_dimension",
        "single_node_net",
        "all_node_net",
        "single_all_agree",
    }
    missing_topology = required_topology.difference(topology_columns)
    if missing_topology:
        raise BenchmarkError(
            "representative diversity topology table is missing: {}".format(
                ", ".join(sorted(missing_topology))
            )
        )

    rac_vectors: Dict[str, Tuple[float, ...]] = {}
    zeo_vectors: Dict[str, Tuple[float, ...]] = {}
    tiers: Dict[str, str] = {}
    topology_categories: Dict[str, str] = {}
    for structure_id in ids:
        rac_vector = None
        if _truthy(rac_rows[structure_id].get(_RAC_AVAILABILITY_COLUMN)):
            rac_vector = _finite_vector(rac_rows[structure_id], rac_fields)
        if rac_vector is not None:
            rac_vectors[structure_id] = rac_vector
            tiers[structure_id] = "rac5"
        else:
            zeo_vector = None
            if all(
                _truthy(zeo_rows[structure_id].get(column))
                for column in _ZEO_AVAILABILITY_COLUMNS
            ):
                zeo_vector = _finite_vector(zeo_rows[structure_id], _ZEO_NUMERIC_FIELDS)
            if zeo_vector is not None:
                zeo_vectors[structure_id] = zeo_vector
                tiers[structure_id] = "zeo"
            else:
                tiers[structure_id] = "no_numeric"
        top = topology_rows[structure_id]
        if _truthy(top.get("topology_available")):
            values = tuple(
                str(top.get(name, "")).strip() or "<MISSING>"
                for name in (
                    "network_dimension",
                    "single_node_net",
                    "all_node_net",
                    "single_all_agree",
                )
            )
            topology_categories[structure_id] = "|".join(values)
        else:
            topology_categories[structure_id] = "<NO_CURRENT_SUCCESS>"

    with threadpool_limits(limits=1):
        backend["threadpool_libraries"] = _threadpool_runtime(threadpool_info)
        rac_strata, rac_details = _cluster_tier(
            tuple(rac_vectors), rac_vectors, "rac5", np, PCA, MiniBatchKMeans
        )
        zeo_strata, zeo_details = _cluster_tier(
            tuple(zeo_vectors), zeo_vectors, "zeo", np, PCA, MiniBatchKMeans
        )
    strata = {}
    for structure_id in ids:
        if structure_id in rac_strata:
            strata[structure_id] = rac_strata[structure_id]
        elif structure_id in zeo_strata:
            strata[structure_id] = zeo_strata[structure_id]
        else:
            strata[structure_id] = "no_numeric:000"
    profile = {
        "schema_version": DIVERSITY_PROFILE_VERSION,
        "name": "representative",
        "purpose": "target-free distribution-balancing strata, never identity or leakage groups",
        "profile_seed": DIVERSITY_PROFILE_SEED,
        "backend": backend,
        "tier_precedence": ["complete_rac5", "complete_selected_zeo", "no_numeric"],
        "rac5_fields": list(rac_fields),
        "zeo_fields": list(_ZEO_NUMERIC_FIELDS),
        "robust_scaling": {
            "center": "median",
            "scale": "interquartile_range",
            "percentile_method": "linear",
            "zero_iqr_behavior": (
                "subtract_the_median_use_a_unit_divisor_and_retain_centered_"
                "deviations_without_forcing_them_to_zero"
            ),
            "imputation": False,
        },
        "rac5_reduction": {
            "algorithm": "PCA",
            "maximum_components": 32,
            "svd_solver": "full",
            "all_zero_scaled_matrix_behavior": "one_zero_coordinate_without_fit",
        },
        "clustering": {
            "algorithm": "MiniBatchKMeans",
            "k_formula": "min(n,256,max(16,ceil(sqrt(n))))",
            "sorted_structure_ids": True,
            "n_init": 10,
            "batch_size_formula": "max(1024,3*k)",
            "reassignment_ratio": 0.0,
        },
        "balance_alongside": [
            "source_database",
            "structure_variant",
            "current_topology_category",
            "feature_availability_tier",
        ],
        "target_columns_consumed": [],
        "scientific_feature_imputation": False,
        "tier_diagnostics": {
            "rac5": dict(rac_details),
            "zeo": dict(zeo_details),
            "no_numeric": {
                "count": sum(value == "no_numeric" for value in tiers.values())
            },
        },
    }
    payload = {
        "profile": profile,
        "inputs": [rac_receipt, zeo_receipt, topology_receipt],
        "strata": strata,
        "tiers": tiers,
        "topology_categories": topology_categories,
    }
    return DiversityIndex(
        strata,
        tiers,
        topology_categories,
        profile,
        (rac_receipt, zeo_receipt, topology_receipt),
        _digest(payload),
    )


def _validate_fractions(
    train: object, val: object, test: object
) -> Tuple[float, float, float]:
    raw = (train, val, test)
    if any(isinstance(value, bool) or not isinstance(value, Real) for value in raw):
        raise TypeError("train, val, and test must be finite non-boolean numbers")
    result = tuple(float(value) for value in raw)
    if any(not math.isfinite(value) or value < 0.0 for value in result):
        raise BenchmarkError("train, val, and test must be finite and non-negative")
    if not math.isclose(sum(result), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise BenchmarkError("train, val, and test must sum to 1.0")
    if not any(value > 0.0 for value in result):
        raise BenchmarkError("at least one partition fraction must be positive")
    return result  # type: ignore[return-value]


def _optional_filter(values: object, name: str) -> Optional[Tuple[str, ...]]:
    if values is None:
        return None
    return _validate_exact_strings(values, name)


def _row_index(
    dataset: object,
) -> Tuple[Mapping[str, Mapping[str, object]], Mapping[str, int]]:
    rows = tuple(dataset.metadata_rows)
    by_id: Dict[str, Mapping[str, object]] = {}
    index: Dict[str, int] = {}
    for position, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise TypeError("metadata rows must be mappings")
        structure_id = row.get("structure_id")
        if (
            not isinstance(structure_id, str)
            or not structure_id
            or structure_id != structure_id.strip()
        ):
            raise BenchmarkError("metadata rows require exact nonblank structure_id")
        if structure_id in by_id:
            raise BenchmarkError("duplicate structure_id: {}".format(structure_id))
        by_id[structure_id] = row
        index[structure_id] = position
    return MappingProxyType(by_id), MappingProxyType(index)


def _label_mapping(classified: object) -> Mapping[str, str]:
    value = getattr(classified, "label_by_id", None)
    if not isinstance(value, Mapping):
        raise TypeError("classified dataset must expose label_by_id")
    result = {}
    for structure_id, label in value.items():
        if not isinstance(structure_id, str) or not isinstance(label, str):
            raise TypeError("classification labels must map string IDs to strings")
        canonical = label.upper()
        if canonical not in LABELS:
            raise BenchmarkError(
                "unknown checker-consensus label {!r} for {}".format(
                    label, structure_id
                )
            )
        result[structure_id] = canonical
    return MappingProxyType(result)


def _balance_values(
    ids: Sequence[str],
    rows: Mapping[str, Mapping[str, object]],
    labels: Mapping[str, str],
    diversity: DiversityIndex,
) -> Mapping[str, Tuple[str, ...]]:
    result = {}
    for structure_id in ids:
        row = rows[structure_id]
        result[structure_id] = (
            labels.get(structure_id, "<MISSING>"),
            str(row.get("source_database", "") or "<MISSING>"),
            str(row.get("structure_variant", "") or "<MISSING>"),
            diversity.topology_category_by_id[structure_id],
            diversity.tier_by_id[structure_id],
            diversity.strata_by_id[structure_id],
        )
    return MappingProxyType(result)


def _assign_blocks(
    ids: Sequence[str],
    block_by_id: Mapping[str, str],
    strata_by_id: Mapping[str, Tuple[str, ...]],
    partition_fractions: Mapping[str, float],
    seed: object,
    namespace: str,
) -> Mapping[str, str]:
    active = tuple(
        partition
        for partition, fraction in partition_fractions.items()
        if fraction > 0.0
    )
    if not active:
        return MappingProxyType({})
    members_by_block: Dict[str, List[str]] = {}
    for structure_id in sorted(ids):
        members_by_block.setdefault(block_by_id[structure_id], []).append(structure_id)
    total = len(ids)
    target_total = {
        partition: partition_fractions[partition] * total for partition in active
    }
    stratum_total = Counter(strata_by_id[structure_id] for structure_id in ids)
    target_stratum = {
        partition: {
            stratum: partition_fractions[partition] * count
            for stratum, count in stratum_total.items()
        }
        for partition in active
    }
    block_strata = {
        block: Counter(strata_by_id[structure_id] for structure_id in members)
        for block, members in members_by_block.items()
    }
    blocks = sorted(
        members_by_block,
        key=lambda block: (
            -len(members_by_block[block]),
            -max(block_strata[block].values()),
            _stable_key(seed, namespace, "block", block),
        ),
    )
    current_total = {partition: 0 for partition in active}
    current_stratum = {partition: Counter() for partition in active}
    assignment = {}
    for block_index, block in enumerate(blocks):
        empty = [partition for partition in active if current_total[partition] == 0]
        candidates = (
            empty if empty and len(blocks) - block_index <= len(empty) else active
        )
        scores = []
        for partition in candidates:
            old_total = current_total[partition]
            new_total = old_total + len(members_by_block[block])
            target = target_total[partition]
            scale = max(target, 1.0)
            score = ((new_total - target) / scale) ** 2 - (
                (old_total - target) / scale
            ) ** 2
            score += 2.0 * (
                (max(0.0, new_total - target) / scale) ** 2
                - (max(0.0, old_total - target) / scale) ** 2
            )
            marginal = 0.0
            for stratum, addition in block_strata[block].items():
                old = current_stratum[partition].get(stratum, 0)
                new = old + addition
                stratum_target = target_stratum[partition][stratum]
                stratum_scale = max(stratum_target, 1.0)
                marginal += ((new - stratum_target) / stratum_scale) ** 2 - (
                    (old - stratum_target) / stratum_scale
                ) ** 2
            if stratum_total:
                score += marginal / len(stratum_total)
            scores.append(
                (
                    score,
                    _stable_key(seed, namespace, "partition", block, partition),
                    partition,
                )
            )
        chosen = min(scores)[2]
        assignment[block] = chosen
        current_total[chosen] += len(members_by_block[block])
        current_stratum[chosen].update(block_strata[block])
    return MappingProxyType(assignment)


def _assign_paired_blocks(
    pool_ids: Sequence[str],
    seed_cohorts: Sequence[CRNCRCohort],
    block_by_id: Mapping[str, str],
    strata_by_id: Mapping[str, Tuple[str, ...]],
    partition_fractions: Mapping[str, float],
    seed: int,
) -> Mapping[str, str]:
    """Assign blocks once while minimizing deviations over the whole ladder."""

    partitions = tuple(
        partition
        for partition, fraction in partition_fractions.items()
        if fraction > 0.0
    )
    if not partitions:
        return MappingProxyType({})
    selected_sets = tuple(
        set(cohort.structure_ids).intersection(pool_ids) for cohort in seed_cohorts
    )
    members_by_block: Dict[str, List[str]] = {}
    for structure_id in sorted(pool_ids):
        members_by_block.setdefault(block_by_id[structure_id], []).append(structure_id)
    exposure = {
        block: tuple(
            sum(structure_id in selected for structure_id in members)
            for selected in selected_sets
        )
        for block, members in members_by_block.items()
    }
    targets = {
        partition: tuple(
            partition_fractions[partition] * len(selected) for selected in selected_sets
        )
        for partition in partitions
    }
    level_strata = tuple(
        Counter(strata_by_id[value] for value in selected) for selected in selected_sets
    )
    level_labels = tuple(
        Counter(strata_by_id[value][0] for value in selected)
        for selected in selected_sets
    )
    block_level_strata = {
        block: tuple(
            Counter(strata_by_id[value] for value in members if value in selected)
            for selected in selected_sets
        )
        for block, members in members_by_block.items()
    }
    block_level_labels = {
        block: tuple(
            Counter(strata_by_id[value][0] for value in members if value in selected)
            for selected in selected_sets
        )
        for block, members in members_by_block.items()
    }
    if len(partitions) == 1:
        return MappingProxyType({block: partitions[0] for block in members_by_block})
    if set(partitions) != {"train", "validation"}:
        raise RuntimeError(
            "internal error: paired assignment expects train/validation partitions"
        )

    # Since train is the exact complement of validation at every ladder level,
    # start with a deterministic globally stratified block assignment, then
    # improve the validation objective jointly across all requested levels.
    # Single-block flips correct size drift; bounded deterministic swaps can
    # replace a persistent CR block with a later NCR prefix even when neither
    # flip helps alone.  This avoids the ordering trap where CR fills
    # validation before NCR is considered, while retaining one immutable block
    # assignment per seed.
    validation = "validation"
    initial = dict(
        _assign_blocks(
            pool_ids,
            block_by_id,
            strata_by_id,
            partition_fractions,
            seed,
            "paired-ladder-initial",
        )
    )
    selected_validation: Set[str] = {
        block for block, partition in initial.items() if partition == validation
    }
    current = [
        sum(exposure[block][level] for block in selected_validation)
        for level in range(len(selected_sets))
    ]
    current_strata = [Counter() for _ in selected_sets]
    current_labels = [Counter() for _ in selected_sets]
    for block in selected_validation:
        for level in range(len(selected_sets)):
            current_strata[level].update(block_level_strata[block][level])
            current_labels[level].update(block_level_labels[block][level])
    total_stratum_terms = sum(len(values) for values in level_strata)
    total_label_terms = sum(len(values) for values in level_labels)

    def change_cost(changes: Sequence[Tuple[str, int]]) -> float:
        count_changes = [0 for _ in selected_sets]
        stratum_changes = [Counter() for _ in selected_sets]
        label_changes = [Counter() for _ in selected_sets]
        for block, direction in changes:
            for level, addition in enumerate(exposure[block]):
                count_changes[level] += direction * addition
                for stratum, count in block_level_strata[block][level].items():
                    stratum_changes[level][stratum] += direction * count
                for label, count in block_level_labels[block][level].items():
                    label_changes[level][label] += direction * count
        score = 0.0
        for level, addition in enumerate(count_changes):
            old = current[level]
            new = old + addition
            target = targets[validation][level]
            scale = max(target, 1.0)
            score += ((new - target) / scale) ** 2 - ((old - target) / scale) ** 2
            score += 2.0 * (
                (max(0.0, new - target) / scale) ** 2
                - (max(0.0, old - target) / scale) ** 2
            )
        if selected_sets:
            score /= len(selected_sets)
        # Label composition is an explicit benchmark objective: validation is
        # meant to mimic the mixed training population at every ladder level.
        # Keep this term separate from the fine-grained composite strata so it
        # cannot be diluted when source/variant/topology strata are numerous.
        label_score = 0.0
        for level, additions in enumerate(label_changes):
            for label, addition in additions.items():
                if addition == 0:
                    continue
                old = current_labels[level].get(label, 0)
                new = old + addition
                target = partition_fractions[validation] * level_labels[level][label]
                scale = max(target, 1.0)
                label_score += ((new - target) / scale) ** 2 - (
                    (old - target) / scale
                ) ** 2
        if total_label_terms:
            score += PAIRED_LABEL_BALANCE_WEIGHT * label_score / total_label_terms
        stratum_score = 0.0
        for level, additions in enumerate(stratum_changes):
            for stratum, addition in additions.items():
                if addition == 0:
                    continue
                old = current_strata[level].get(stratum, 0)
                new = old + addition
                target = partition_fractions[validation] * level_strata[level][stratum]
                scale = max(target, 1.0)
                stratum_score += ((new - target) / scale) ** 2 - (
                    (old - target) / scale
                ) ** 2
        if total_stratum_terms:
            score += stratum_score / total_stratum_terms
        return score

    def apply_changes(changes: Sequence[Tuple[str, int]]) -> None:
        for block, direction in changes:
            if direction > 0:
                selected_validation.add(block)
            else:
                selected_validation.remove(block)
            for level, addition in enumerate(exposure[block]):
                current[level] += direction * addition
                for stratum, count in block_level_strata[block][level].items():
                    updated = current_strata[level][stratum] + direction * count
                    if updated:
                        current_strata[level][stratum] = updated
                    else:
                        del current_strata[level][stratum]
                for label, count in block_level_labels[block][level].items():
                    updated = current_labels[level][label] + direction * count
                    if updated:
                        current_labels[level][label] = updated
                    else:
                        del current_labels[level][label]

    active_blocks = tuple(block for block in members_by_block if any(exposure[block]))
    maximum_improvements = max(16, 8 * len(selected_sets))
    pair_candidate_count = PAIRED_SWAP_CANDIDATE_COUNT
    for _ in range(maximum_improvements):
        single_candidates = []
        for block in active_blocks:
            direction = -1 if block in selected_validation else 1
            single_candidates.append(
                (
                    change_cost(((block, direction),)),
                    _stable_key(seed, "paired-ladder-single", direction, block),
                    block,
                    direction,
                )
            )
        best_single = min(single_candidates, default=None)
        if best_single is not None and best_single[0] < -1e-15:
            apply_changes(((best_single[2], best_single[3]),))
            continue

        removals = sorted(
            (
                change_cost(((block, -1),)),
                _stable_key(seed, "paired-ladder-removal", block),
                block,
            )
            for block in selected_validation
            if any(exposure[block])
        )[:pair_candidate_count]
        additions = sorted(
            (
                change_cost(((block, 1),)),
                _stable_key(seed, "paired-ladder-addition", block),
                block,
            )
            for block in active_blocks
            if block not in selected_validation
        )[:pair_candidate_count]
        pair_candidates = []
        for _, removal_key, removal in removals:
            for _, addition_key, addition in additions:
                pair_candidates.append(
                    (
                        change_cost(((removal, -1), (addition, 1))),
                        removal_key,
                        addition_key,
                        removal,
                        addition,
                    )
                )
        best_pair = min(pair_candidates, default=None)
        if best_pair is None or best_pair[0] >= -1e-15:
            break
        apply_changes(((best_pair[3], -1), (best_pair[4], 1)))
    return MappingProxyType(
        {
            block: validation if block in selected_validation else "train"
            for block in members_by_block
        }
    )


def _leakage_audit(
    assignments: Mapping[str, str], blocks: Mapping[str, str]
) -> Mapping[str, object]:
    partitions_by_block: Dict[str, Set[str]] = {}
    sizes = Counter()
    for structure_id, partition in assignments.items():
        block = blocks[structure_id]
        partitions_by_block.setdefault(block, set()).add(partition)
        sizes[block] += 1
    crossed = tuple(
        sorted(
            block for block, values in partitions_by_block.items() if len(values) > 1
        )
    )
    return MappingProxyType(
        {
            "block_count": len(partitions_by_block),
            "cross_split_block_count": len(crossed),
            "cross_split_blocks": crossed,
            "max_selected_block_size": max(sizes.values()) if sizes else 0,
            "passed": not crossed,
        }
    )


def _assignment_digest(
    assignments: Mapping[str, str], exclusions: Mapping[str, str]
) -> str:
    lines = []
    for structure_id in sorted(set(assignments).union(exclusions)):
        lines.append(
            "\t".join(
                (
                    structure_id,
                    assignments.get(structure_id, "EXCLUDED"),
                    exclusions.get(structure_id, ""),
                )
            )
        )
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def _safe_stem(stem: str) -> str:
    if (
        not isinstance(stem, str)
        or not stem
        or stem in {".", ".."}
        or Path(stem).name != stem
    ):
        raise BenchmarkError(
            "stem must be a non-empty filename stem without directories"
        )
    return stem


def _publish_directory(staging: Path, target: Path, overwrite: bool) -> Path:
    from CoREMOF._transactions import publish_directory

    return publish_directory(staging, target, overwrite=overwrite)


@dataclass(frozen=True)
class DataSplitResult:
    """Immutable general split with criterion and effective-block evidence."""

    train_ids: Tuple[str, ...]
    validation_ids: Tuple[str, ...]
    test_ids: Tuple[str, ...]
    assignments: Mapping[str, str]
    exclusions: Mapping[str, str]
    labels: Mapping[str, str]
    criterion_groups: Mapping[str, Mapping[str, str]]
    effective_leakage_blocks: Mapping[str, str]
    diversity_strata: Mapping[str, str]
    diversity_tiers: Mapping[str, str]
    diagnostics: Mapping[str, object]
    group_criteria: Tuple[str, ...]
    leakage_guard: str
    diversity: str
    fractions: Tuple[float, float, float]
    random_state: str
    dataset_version: Optional[str]
    checker_view: Optional[str]
    input_hashes: Mapping[str, str]
    diversity_index_hash: str
    assignment_digest: str
    official_split: bool
    _receipt: Mapping[str, object] = field(repr=False, compare=False)
    _dataset: object = field(repr=False, compare=False)
    _authority_factory_token: object = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        for name in (
            "train_ids",
            "validation_ids",
            "test_ids",
            "group_criteria",
            "fractions",
        ):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        for name in (
            "assignments",
            "exclusions",
            "labels",
            "criterion_groups",
            "effective_leakage_blocks",
            "diversity_strata",
            "diversity_tiers",
            "diagnostics",
            "input_hashes",
            "_receipt",
        ):
            object.__setattr__(self, name, _freeze(getattr(self, name)))
        if self.official_split is not False:
            raise BenchmarkError(
                "new exploratory splits must record official_split=false"
            )

    @property
    def counts(self) -> Mapping[str, int]:
        _require_benchmark_result(self)
        return MappingProxyType(
            {
                "train": len(self.train_ids),
                "validation": len(self.validation_ids),
                "test": len(self.test_ids),
                "excluded": len(self.exclusions),
            }
        )

    @property
    def leakage_audit(self) -> Mapping[str, object]:
        _require_benchmark_result(self)
        return _leakage_audit(self.assignments, self.effective_leakage_blocks)

    def receipt(self) -> Mapping[str, object]:
        _require_benchmark_result(self)
        return _jsonable(self._receipt)  # type: ignore[return-value]

    def assignment_rows(self) -> Tuple[Mapping[str, object], ...]:
        _require_benchmark_result(self)
        ids = sorted(set(self.assignments).union(self.exclusions))
        return tuple(
            {
                "structure_id": structure_id,
                "partition": self.assignments.get(structure_id, "EXCLUDED"),
                "label": self.labels.get(structure_id, ""),
                "effective_leakage_block": self.effective_leakage_blocks.get(
                    structure_id, ""
                ),
                "diversity_tier": self.diversity_tiers.get(structure_id, ""),
                "diversity_stratum": self.diversity_strata.get(structure_id, ""),
                "exclusion_reason": self.exclusions.get(structure_id, ""),
            }
            for structure_id in ids
        )

    def write(
        self,
        output_directory: object,
        stem: str = "coremof_data_split",
        overwrite: bool = False,
    ) -> Tuple[Path, Path]:
        """Transactionally publish assignment CSV and receipt JSON."""

        _require_benchmark_result(self)
        if type(overwrite) is not bool:
            raise TypeError("overwrite must be a boolean")
        stem = _safe_stem(stem)
        directory = Path(output_directory)
        directory.mkdir(parents=True, exist_ok=True)
        csv_path = directory / (stem + ".csv")
        json_path = directory / (stem + ".json")
        if not overwrite and (csv_path.exists() or json_path.exists()):
            raise FileExistsError("Refusing to overwrite existing data-split output")
        staging = Path(tempfile.mkdtemp(prefix=".{}-".format(stem), dir=str(directory)))
        staged_csv = staging / csv_path.name
        staged_json = staging / json_path.name
        preserve_staging = False
        try:
            rows = self.assignment_rows()
            fields = (
                tuple(rows[0])
                if rows
                else (
                    "structure_id",
                    "partition",
                    "label",
                    "effective_leakage_block",
                    "diversity_tier",
                    "diversity_stratum",
                    "exclusion_reason",
                )
            )
            with staged_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
                handle.flush()
                os.fsync(handle.fileno())
            with staged_json.open("w", encoding="utf-8") as handle:
                json.dump(
                    self.receipt(), handle, indent=2, sort_keys=True, allow_nan=False
                )
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            from CoREMOF._transactions import publish_file_bundle

            publish_file_bundle(
                (staged_csv, staged_json),
                (csv_path, json_path),
                overwrite=overwrite,
            )
        except BaseException as error:
            preserve_staging = hasattr(error, "coremof_preserved_staging_directory")
            raise
        finally:
            if not preserve_staging:
                shutil.rmtree(staging, ignore_errors=True)
        return csv_path, json_path

    def attach_targets(
        self, target_data: object, missing: str = "keep", **kwargs: object
    ):
        """Attach targets to the frozen assignment without resplitting."""

        _require_benchmark_result(self)
        from CoREMOF.attachments import attach_targets

        return attach_targets(
            self, target_data, dataset=self._dataset, missing=missing, **kwargs
        )


def data_split(
    classified: object,
    group_criteria: object = "priority_main",
    train: object = 0.8,
    val: object = 0.1,
    test: object = 0.1,
    leakage_guard: str = EFFECTIVE_LEAKAGE_POLICY,
    diversity: str = "representative",
    random_state: object = 42,
    labels: object = ("CR", "NCR"),
    sources: object = None,
    variants: object = None,
    metals: object = None,
    structure_ids: object = None,
    official: bool = False,
) -> DataSplitResult:
    """Build a target-independent exploratory split over effective blocks."""

    if leakage_guard != EFFECTIVE_LEAKAGE_POLICY:
        raise BenchmarkError(
            "leakage_guard must be {!r}; this API always unions main_union with "
            "every selected criterion".format(EFFECTIVE_LEAKAGE_POLICY)
        )
    if type(official) is not bool:
        raise TypeError("official must be a boolean")
    if official:
        raise BenchmarkError(
            "No audited official assignment manifest is available; use official=False"
        )
    _validate_classified_authority(classified, require_official=False)
    if isinstance(random_state, bool) or not isinstance(random_state, (int, str)):
        raise TypeError("random_state must be a non-boolean integer or exact string")
    if isinstance(random_state, str) and (
        not random_state or random_state != random_state.strip()
    ):
        raise BenchmarkError("random_state string must be exact and nonblank")
    seed = str(random_state)
    fractions = _validate_fractions(train, val, test)
    criteria = normalize_group_criteria(group_criteria)
    dataset = _dataset_from(classified)
    rows, _ = _row_index(dataset)
    label_by_id = _label_mapping(classified)
    view_ids = tuple(getattr(classified, "structure_ids", tuple(label_by_id)))
    if (
        len(set(view_ids)) != len(view_ids)
        or set(view_ids).difference(rows)
        or set(label_by_id) != set(view_ids)
    ):
        raise BenchmarkError(
            "classified view IDs and labels must be the same unique subset of the release"
        )
    grouping = _grouping_context(dataset, criteria)
    diversity_index = build_diversity_index(dataset, diversity=diversity)

    label_filter = _optional_filter(labels, "labels")
    source_filter = _optional_filter(sources, "sources")
    variant_filter = _optional_filter(variants, "variants")
    metal_filter = _optional_filter(metals, "metals")
    id_filter = _optional_filter(structure_ids, "structure_ids")
    view_set = set(view_ids)
    allowed_labels = {value.upper() for value in label_filter} if label_filter else None
    if allowed_labels is not None:
        unknown_labels = allowed_labels.difference(LABELS)
        if unknown_labels:
            raise BenchmarkError(
                "unknown labels: {}; expected only {}".format(
                    ", ".join(sorted(unknown_labels)),
                    ", ".join(sorted(LABELS)),
                )
            )
    allowed_sources = (
        {value.casefold() for value in source_filter} if source_filter else None
    )
    allowed_variants = (
        {value.casefold() for value in variant_filter} if variant_filter else None
    )
    allowed_metals = (
        {value.casefold() for value in metal_filter} if metal_filter else None
    )
    allowed_ids = set(id_filter) if id_filter else None
    if allowed_ids is not None:
        unknown = allowed_ids.difference(rows)
        if unknown:
            raise KeyError(
                "unknown structure IDs: {}".format(", ".join(sorted(unknown)[:5]))
            )
    exclusions: Dict[str, str] = {}
    eligible = []
    for structure_id in sorted(rows):
        row = rows[structure_id]
        reason = None
        if structure_id not in view_set:
            reason = "PRESELECTION_FILTER"
        elif structure_id not in label_by_id:
            reason = "LABEL_NOT_AVAILABLE"
        elif (
            allowed_labels is not None
            and label_by_id[structure_id] not in allowed_labels
        ):
            reason = "LABEL_FILTER"
        elif (
            allowed_sources is not None
            and str(row.get("source_database", "")).casefold() not in allowed_sources
        ):
            reason = "SOURCE_FILTER"
        elif (
            allowed_variants is not None
            and str(row.get("structure_variant", "")).casefold() not in allowed_variants
        ):
            reason = "VARIANT_FILTER"
        elif allowed_ids is not None and structure_id not in allowed_ids:
            reason = "STRUCTURE_ID_FILTER"
        elif allowed_metals is not None:
            metals_value = row.get("metal_elements", "")
            present = {
                value.casefold()
                for value in re.split(r"[,;|\s]+", str(metals_value))
                if value
            }
            if not present.intersection(allowed_metals):
                reason = "METAL_FILTER"
        if reason is None:
            eligible.append(structure_id)
        else:
            exclusions[structure_id] = reason
    balance = _balance_values(eligible, rows, label_by_id, diversity_index)
    fractions_by_partition = dict(zip(_PARTITIONS, fractions))
    by_block = _assign_blocks(
        eligible,
        grouping.effective_blocks,
        balance,
        fractions_by_partition,
        seed,
        "general-data-split",
    )
    assignments = {
        structure_id: by_block[grouping.effective_blocks[structure_id]]
        for structure_id in eligible
    }
    audit = _leakage_audit(assignments, grouping.effective_blocks)
    if not audit["passed"]:
        raise RuntimeError(
            "internal error: an effective leakage block crossed partitions"
        )
    ids_by_partition = {
        partition: tuple(
            sorted(
                structure_id
                for structure_id, value in assignments.items()
                if value == partition
            )
        )
        for partition in _PARTITIONS
    }
    assignment_digest = _assignment_digest(assignments, exclusions)
    counts = {partition: len(ids_by_partition[partition]) for partition in _PARTITIONS}
    counts["excluded"] = len(exclusions)
    target_counts = {
        partition: fractions_by_partition[partition] * len(eligible)
        for partition in _PARTITIONS
    }
    deviations = {
        partition: counts[partition] - target_counts[partition]
        for partition in _PARTITIONS
    }
    input_hashes = dict(getattr(dataset, "input_hashes", {}) or {})
    for item in diversity_index.input_receipts:
        input_hashes[str(item["release_path"])] = str(item["sha256"])
    complete_block_sizes = Counter(grouping.effective_blocks.values())
    selected_labels_by_block: Dict[str, Set[str]] = {}
    for structure_id in eligible:
        selected_labels_by_block.setdefault(
            grouping.effective_blocks[structure_id], set()
        ).add(label_by_id[structure_id])
    mixed_label_blocks = tuple(
        sorted(
            block
            for block, block_labels in selected_labels_by_block.items()
            if len(block_labels) > 1
        )
    )
    receipt = {
        "schema_version": "coremof-data-split-receipt/1.0",
        "benchmark_api_version": BENCHMARK_API_VERSION,
        "implementation": {
            "package": "CoREMOF-tools",
            "package_version": __version__,
            "source_sha256": dict(_source_hashes()),
        },
        "contract_definitions": _jsonable(generated_output_terminology_definitions()),
        "effective_leakage_policy_definition": EFFECTIVE_LEAKAGE_POLICY_DEFINITION,
        "dataset_version": getattr(dataset, "dataset_version", None),
        "checker_view": getattr(classified, "checker_view", None),
        "official_split": False,
        "group_criteria": list(criteria),
        "group_criterion_definitions": {
            criterion: _jsonable(parent_method_definition(criterion))
            for criterion in criteria
        },
        "criterion_groups_sha256": {
            criterion: _digest(sorted(groups.items()))
            for criterion, groups in grouping.criterion_groups.items()
        },
        "leakage_guard": leakage_guard,
        "effective_leakage_blocks": {
            "sha256": _digest(sorted(grouping.effective_blocks.items())),
            "complete_release_block_count": len(complete_block_sizes),
            "maximum_complete_release_block_size": max(complete_block_sizes.values())
            if complete_block_sizes
            else 0,
            "selected_mixed_label_block_count": len(mixed_label_blocks),
            "selected_mixed_label_blocks": list(mixed_label_blocks),
        },
        "diversity": diversity,
        "diversity_profile": _jsonable(diversity_index.profile),
        "diversity_index_sha256": diversity_index.digest,
        "target_columns_consumed": [],
        "fractions": dict(zip(_PARTITIONS, fractions)),
        "random_state": seed,
        "filters": {
            "labels": list(label_filter) if label_filter else None,
            "sources": list(source_filter) if source_filter else None,
            "variants": list(variant_filter) if variant_filter else None,
            "metals": list(metal_filter) if metal_filter else None,
            "structure_ids": list(id_filter) if id_filter else None,
        },
        "counts": counts,
        "requested_counts": target_counts,
        "partition_count_deviations": deviations,
        "criterion_diagnostics": _jsonable(grouping.criterion_diagnostics),
        "leakage_audit": dict(audit),
        "assignment_sha256": assignment_digest,
        "input_hashes": dict(sorted(input_hashes.items())),
        "release_binding": _jsonable(_release_binding(dataset)),
        "partitions": {
            partition: list(ids_by_partition[partition]) for partition in _PARTITIONS
        },
        "exclusions": dict(sorted(exclusions.items())),
    }
    return _register_benchmark_result(
        DataSplitResult(
            train_ids=ids_by_partition["train"],
            validation_ids=ids_by_partition["validation"],
            test_ids=ids_by_partition["test"],
            assignments=assignments,
            exclusions=exclusions,
            labels={
                structure_id: label_by_id[structure_id] for structure_id in label_by_id
            },
            criterion_groups=grouping.criterion_groups,
            effective_leakage_blocks=grouping.effective_blocks,
            diversity_strata=diversity_index.strata_by_id,
            diversity_tiers=diversity_index.tier_by_id,
            diagnostics={
                "criterion": grouping.criterion_diagnostics,
                "partition_count_deviations": deviations,
                "mixed_label_blocks": mixed_label_blocks,
                "leakage_audit": audit,
            },
            group_criteria=criteria,
            leakage_guard=leakage_guard,
            diversity=diversity,
            fractions=fractions,
            random_state=seed,
            dataset_version=getattr(dataset, "dataset_version", None),
            checker_view=getattr(classified, "checker_view", None),
            input_hashes=input_hashes,
            diversity_index_hash=diversity_index.digest,
            assignment_digest=assignment_digest,
            official_split=False,
            _receipt=receipt,
            _dataset=dataset,
            _authority_factory_token=_BENCHMARK_RESULT_FACTORY_TOKEN,
        ),
        dataset,
    )


def _round_half_up(value: Decimal) -> int:
    return int(value.quantize(Decimal("1"), rounding=ROUND_HALF_UP))


def _fraction_values(values: object) -> Tuple[Decimal, ...]:
    if not isinstance(values, (list, tuple)) or not values:
        raise TypeError("ncr_pool_fractions must be a non-empty ordered list/tuple")
    result = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError("NCR-pool fractions must be finite non-boolean numbers")
        decimal = Decimal(str(value))
        if not decimal.is_finite() or decimal < 0 or decimal > 1:
            raise BenchmarkError("NCR-pool fractions must lie in [0, 1]")
        result.append(decimal)
    if len(set(result)) != len(result):
        raise BenchmarkError("ncr_pool_fractions must not contain duplicates")
    if tuple(result) != tuple(sorted(result)):
        raise BenchmarkError("ncr_pool_fractions must be strictly increasing")
    return tuple(result)


def _seed_values(values: object) -> Tuple[int, ...]:
    if not isinstance(values, (list, tuple)) or not values:
        raise TypeError("seeds must be a non-empty ordered list/tuple")
    result = tuple(values)
    if any(isinstance(value, bool) or not isinstance(value, int) for value in result):
        raise TypeError("seeds must contain non-boolean integers")
    if len(set(result)) != len(result):
        raise BenchmarkError("seeds must not contain duplicates")
    return result


def _balanced_order(
    ids: Iterable[str],
    strata: Mapping[str, Tuple[str, ...]],
    seed: object,
    namespace: str,
) -> Tuple[str, ...]:
    buckets: Dict[Tuple[str, ...], List[str]] = {}
    for structure_id in sorted(ids):
        buckets.setdefault(strata[structure_id], []).append(structure_id)
    ranked = []
    for stratum, members in buckets.items():
        members.sort(key=lambda value: _stable_key(seed, namespace, stratum, value))
        for index, structure_id in enumerate(members):
            ranked.append(
                (
                    Decimal(2 * index + 1) / Decimal(2 * len(members)),
                    _stable_key(seed, namespace, "interleave", stratum, structure_id),
                    structure_id,
                )
            )
    return tuple(item[2] for item in sorted(ranked))


def _fixed_test_blocks(
    cr_ids: Sequence[str],
    label_by_id: Mapping[str, str],
    blocks: Mapping[str, str],
    strata: Mapping[str, Tuple[str, ...]],
    desired_count: int,
) -> Tuple[str, ...]:
    if desired_count < 0:
        raise BenchmarkError("fixed pure-CR test count must be non-negative")
    if desired_count == 0:
        return ()
    members_by_block: Dict[str, List[str]] = {}
    for structure_id, block in blocks.items():
        members_by_block.setdefault(block, []).append(structure_id)
    candidates = {
        block: tuple(sorted(members))
        for block, members in members_by_block.items()
        if members and all(label_by_id.get(value) == "CR" for value in members)
    }
    cr_set = set(cr_ids)
    candidates = {
        block: tuple(value for value in members if value in cr_set)
        for block, members in candidates.items()
    }
    candidates = {block: members for block, members in candidates.items() if members}
    if desired_count > 0 and not candidates:
        raise BenchmarkFeasibilityError(
            "fixed_pure_cr requested a nonempty test set, but no complete-release "
            "effective leakage block contains only strict five-checker CR rows"
        )
    # First solve the whole-block cardinality problem globally.  A greedy
    # nearest-addition rule can miss an exact combination (for example block
    # sizes 8, 7, and 3 for a requested count of 10), and that can turn a
    # feasible NCR-replacement gate into a false rejection.  Python integers
    # provide a compact subset-sum bitset, so the complete release remains
    # inexpensive at this scale.
    reachable_counts = 1
    for members in candidates.values():
        reachable_counts |= reachable_counts << len(members)
    total_candidate_count = sum(len(members) for members in candidates.values())
    achieved_count = min(
        (
            count
            for count in range(1, total_candidate_count + 1)
            if reachable_counts & (1 << count)
        ),
        key=lambda count: (
            abs(count - desired_count),
            count > desired_count,
            count,
        ),
    )

    # Use the target-free balanced structure order as the deterministic
    # secondary priority among subsets with the globally optimal cardinality.
    # The exact count is never traded away for this distributional preference.
    balanced_ids = _balanced_order(
        (value for members in candidates.values() for value in members),
        strata,
        DIVERSITY_PROFILE_SEED,
        "fixed-pure-cr-block-priority",
    )
    rank_by_id = {
        structure_id: index for index, structure_id in enumerate(balanced_ids)
    }
    block_order = sorted(
        candidates,
        key=lambda block: (
            tuple(sorted(rank_by_id[value] for value in candidates[block])),
            _stable_key(DIVERSITY_PROFILE_SEED, "fixed-pure-cr", block),
            block,
        ),
    )

    # Reconstruct one exact whole-block subset in the balanced priority order.
    # Each attainable sum is recorded only on first discovery, so predecessor
    # storage is O(achieved_count) and no leakage block can be reused.
    mask = (1 << (achieved_count + 1)) - 1
    reachable = 1
    parent_sum = [-1] * (achieved_count + 1)
    parent_block: List[Optional[str]] = [None] * (achieved_count + 1)
    for block in block_order:
        size = len(candidates[block])
        newly_reachable = ((reachable << size) & mask) & ~reachable
        pending = newly_reachable
        while pending:
            bit = pending & -pending
            count = bit.bit_length() - 1
            parent_sum[count] = count - size
            parent_block[count] = block
            pending ^= bit
        reachable |= newly_reachable
        if reachable & (1 << achieved_count):
            break
    if not reachable & (1 << achieved_count):
        raise RuntimeError(
            "internal error: failed to reconstruct fixed-test block subset"
        )
    selected_blocks = []
    count = achieved_count
    while count:
        block = parent_block[count]
        if block is None or parent_sum[count] < 0:
            raise RuntimeError(
                "internal error: incomplete fixed-test subset predecessor"
            )
        selected_blocks.append(block)
        count = parent_sum[count]
    return tuple(
        sorted(
            structure_id
            for block in selected_blocks
            for structure_id in candidates[block]
        )
    )


def _projected_block_members(
    ids: Iterable[str], blocks: Mapping[str, str]
) -> Mapping[str, Tuple[str, ...]]:
    members: Dict[str, List[str]] = {}
    for structure_id in sorted(ids):
        members.setdefault(blocks[structure_id], []).append(structure_id)
    return MappingProxyType(
        {block: tuple(values) for block, values in sorted(members.items())}
    )


def _balanced_block_order(
    members_by_block: Mapping[str, Sequence[str]],
    strata: Mapping[str, Tuple[str, ...]],
    seed: object,
    namespace: str,
) -> Tuple[str, ...]:
    balanced_ids = _balanced_order(
        (
            structure_id
            for members in members_by_block.values()
            for structure_id in members
        ),
        strata,
        seed,
        namespace,
    )
    rank = {structure_id: index for index, structure_id in enumerate(balanced_ids)}
    return tuple(
        sorted(
            members_by_block,
            key=lambda block: (
                tuple(sorted(rank[value] for value in members_by_block[block])),
                _stable_key(seed, namespace, "block", block),
                block,
            ),
        )
    )


def _exact_block_subset(
    members_by_block: Mapping[str, Sequence[str]],
    block_order: Sequence[str],
    target_count: int,
    *,
    label: str,
    seed: int,
    transition: str,
) -> Tuple[str, ...]:
    """Choose one deterministic exact whole-block subset or fail closed."""

    if target_count < 0:
        raise BenchmarkError("whole-block target count must be non-negative")
    if target_count == 0:
        return ()
    available_count = sum(len(members_by_block[block]) for block in block_order)
    if target_count > available_count:
        raise BenchmarkFeasibilityError(
            "whole-block {} transition {} for seed {} requests {} rows but only "
            "{} remain".format(label, transition, seed, target_count, available_count)
        )

    mask = (1 << (target_count + 1)) - 1
    reachable = 1
    parent_sum = [-1] * (target_count + 1)
    parent_block: List[Optional[str]] = [None] * (target_count + 1)
    for block in block_order:
        size = len(members_by_block[block])
        newly_reachable = ((reachable << size) & mask) & ~reachable
        pending = newly_reachable
        while pending:
            bit = pending & -pending
            count = bit.bit_length() - 1
            parent_sum[count] = count - size
            parent_block[count] = block
            pending ^= bit
        reachable |= newly_reachable
        if reachable & (1 << target_count):
            break
    if not reachable & (1 << target_count):
        full_reachable = 1
        for block in block_order:
            full_reachable |= full_reachable << len(members_by_block[block])
        lower = max(
            count for count in range(0, target_count) if full_reachable & (1 << count)
        )
        upper = next(
            (
                count
                for count in range(target_count + 1, available_count + 1)
                if full_reachable & (1 << count)
            ),
            None,
        )
        raise BenchmarkFeasibilityError(
            "exact whole-block {} transition {} for seed {} is infeasible: "
            "requested {}, nearest reachable below {}, nearest reachable above {}; "
            "no block was split and no exploratory deviation was emitted".format(
                label, transition, seed, target_count, lower, upper
            )
        )
    selected = []
    count = target_count
    while count:
        block = parent_block[count]
        if block is None or parent_sum[count] < 0:
            raise RuntimeError("internal error: incomplete whole-block predecessor")
        selected.append(block)
        count = parent_sum[count]
    return tuple(selected)


def _nested_exact_block_memberships(
    ids: Sequence[str],
    blocks: Mapping[str, str],
    strata: Mapping[str, Tuple[str, ...]],
    target_counts: Sequence[int],
    seed: int,
    namespace: str,
    label: str,
) -> Mapping[int, Tuple[str, ...]]:
    """Return globally feasible nested exact whole-block memberships.

    Every positive transition is solved jointly as one exact multi-bin
    assignment.  Blocks not needed by the largest requested membership enter
    an explicit unused bin.  Size-one blocks are interchangeable exact fillers;
    larger blocks are placed first in deterministic target-free priority order.
    A bounded exhaustive fallback is used only when the greedy placement cannot
    be certified by the available singleton slack.
    """

    if tuple(target_counts) != tuple(sorted(target_counts)):
        raise RuntimeError("internal error: nested whole-block targets are not ordered")
    if not target_counts:
        raise RuntimeError("internal error: nested whole-block targets are empty")
    unique_targets = tuple(dict.fromkeys(target_counts))
    if unique_targets[0] < 0:
        raise BenchmarkError("whole-block target count must be non-negative")
    members_by_block = _projected_block_members(ids, blocks)
    priority = _balanced_block_order(members_by_block, strata, seed, namespace)
    available_count = sum(len(values) for values in members_by_block.values())
    final_target = unique_targets[-1]
    if final_target > available_count:
        raise BenchmarkFeasibilityError(
            "joint whole-block {} assignment for seed {} requests {} rows but only "
            "{} are available".format(label, seed, final_target, available_count)
        )

    transitions = []
    previous = 0
    for target in unique_targets:
        delta = target - previous
        if delta:
            # This cheap one-dimensional proof preserves precise nearest-count
            # diagnostics when even one transition is individually impossible.
            _exact_block_subset(
                members_by_block,
                priority,
                delta,
                label=label,
                seed=seed,
                transition="{}->{}".format(previous, target),
            )
            transitions.append((target, delta))
        previous = target
    bin_names = [("transition", target) for target, _ in transitions]
    capacities = [delta for _, delta in transitions]
    unused_count = available_count - final_target
    if unused_count:
        bin_names.append(("unused", None))
        capacities.append(unused_count)
    if sum(capacities) != available_count:
        raise RuntimeError(
            "internal error: joint whole-block capacities are incomplete"
        )

    singleton_blocks = [
        block for block in priority if len(members_by_block[block]) == 1
    ]
    priority_rank = {block: index for index, block in enumerate(priority)}
    multi_blocks = sorted(
        (block for block in priority if len(members_by_block[block]) > 1),
        key=lambda block: (-len(members_by_block[block]), priority_rank[block]),
    )

    def candidate_bins(
        block: str, size: int, residual: Sequence[int]
    ) -> Tuple[int, ...]:
        return tuple(
            sorted(
                (index for index, value in enumerate(residual) if value >= size),
                key=lambda index: (
                    1 if bin_names[index][0] == "unused" else 0,
                    index,
                    residual[index] - size,
                    _stable_key(seed, namespace, block, bin_names[index]),
                ),
            )
        )

    residual = list(capacities)
    assigned_bins: List[int] = []
    greedy_complete = True
    for block in multi_blocks:
        size = len(members_by_block[block])
        choices = candidate_bins(block, size, residual)
        if not choices:
            greedy_complete = False
            break
        chosen = choices[0]
        residual[chosen] -= size
        assigned_bins.append(chosen)

    if not greedy_complete:
        # Exact fallback.  Equal residual capacities are symmetric for
        # feasibility, while the first reconstructed labeled path remains
        # deterministic.  A search-limit failure is reported as such and never
        # mislabeled as mathematical infeasibility.
        if len(multi_blocks) >= 900:
            raise BenchmarkFeasibilityError(
                "joint exact whole-block {} assignment for seed {} exceeded the "
                "bounded-search input limit; no suite was emitted and mathematical "
                "infeasibility was not inferred".format(label, seed)
            )
        dead: Set[Tuple[int, Tuple[int, ...]]] = set()
        explored = [0]
        search_limit = 2_000_000

        def search(index: int, state: Tuple[int, ...]) -> Optional[Tuple[int, ...]]:
            if index == len(multi_blocks):
                return ()
            key = (index, tuple(sorted(state)))
            if key in dead:
                return None
            explored[0] += 1
            if explored[0] > search_limit:
                raise BenchmarkFeasibilityError(
                    "joint exact whole-block {} assignment for seed {} exceeded "
                    "the {}-state search limit; no suite was emitted and mathematical "
                    "infeasibility was not inferred".format(label, seed, search_limit)
                )
            block = multi_blocks[index]
            size = len(members_by_block[block])
            seen_residuals: Set[int] = set()
            for chosen in candidate_bins(block, size, state):
                if state[chosen] in seen_residuals:
                    continue
                seen_residuals.add(state[chosen])
                updated = list(state)
                updated[chosen] -= size
                suffix = search(index + 1, tuple(updated))
                if suffix is not None:
                    return (chosen,) + suffix
            dead.add(key)
            return None

        solution = search(0, tuple(capacities))
        if solution is None:
            raise BenchmarkFeasibilityError(
                "joint exact whole-block {} assignment for seed {} is infeasible; "
                "no block was split and no exploratory deviation was emitted".format(
                    label, seed
                )
            )
        assigned_bins = list(solution)
        residual = list(capacities)
        for block, chosen in zip(multi_blocks, assigned_bins):
            residual[chosen] -= len(members_by_block[block])

    if sum(residual) != len(singleton_blocks):
        raise RuntimeError("internal error: singleton filler accounting differs")
    blocks_by_bin: List[List[str]] = [[] for _ in capacities]
    for block, chosen in zip(multi_blocks, assigned_bins):
        blocks_by_bin[chosen].append(block)
    cursor = 0
    for index, count in enumerate(residual):
        blocks_by_bin[index].extend(singleton_blocks[cursor : cursor + count])
        cursor += count
    if cursor != len(singleton_blocks):
        raise RuntimeError(
            "internal error: singleton fillers were not consumed exactly"
        )

    selected_blocks: Set[str] = set()
    result: Dict[int, Tuple[str, ...]] = {}
    transition_bin = 0
    for target in unique_targets:
        if target:
            while (
                transition_bin < len(bin_names)
                and bin_names[transition_bin][0] == "transition"
                and bin_names[transition_bin][1] <= target
            ):
                selected_blocks.update(blocks_by_bin[transition_bin])
                transition_bin += 1
        members = tuple(
            sorted(
                structure_id
                for block in selected_blocks
                for structure_id in members_by_block[block]
            )
        )
        if len(members) != target:
            raise RuntimeError(
                "internal error: joint whole-block target was not reached"
            )
        result[target] = members
    return MappingProxyType(result)


def _whole_block_membership_audit(
    selected_ids: Iterable[str], blocks: Mapping[str, str]
) -> Mapping[str, object]:
    selected = set(selected_ids)
    complete: Dict[str, Set[str]] = {}
    for structure_id, block in blocks.items():
        complete.setdefault(block, set()).add(structure_id)
    touched = {
        block for structure_id, block in blocks.items() if structure_id in selected
    }
    partial = tuple(
        sorted(block for block in touched if not complete[block].issubset(selected))
    )
    return MappingProxyType(
        {
            "selected_effective_block_count": len(touched),
            "partially_selected_block_count": len(partial),
            "partially_selected_blocks": partial,
            "passed": not partial,
        }
    )


@dataclass(frozen=True)
class CRNCRCohort:
    """One fixed-size, target-independent CR/NCR cohort membership."""

    seed: int
    requested_ncr_pool_fraction: str
    selected_ncr_count: int
    selected_cr_count: int
    actual_ncr_ratio: float
    structure_ids: Tuple[str, ...]
    cr_ids: Tuple[str, ...]
    ncr_ids: Tuple[str, ...]
    whole_block_membership_audit: Mapping[str, object]

    def __post_init__(self) -> None:
        for name in ("structure_ids", "cr_ids", "ncr_ids"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        object.__setattr__(
            self,
            "whole_block_membership_audit",
            _freeze(self.whole_block_membership_audit),
        )


@dataclass(frozen=True)
class CRNCRCohortSuite:
    """Immutable cohort memberships before train/validation assignment."""

    cohorts: Tuple[CRNCRCohort, ...]
    fixed_test_ids: Tuple[str, ...]
    full_cr_ids: Tuple[str, ...]
    full_ncr_ids: Tuple[str, ...]
    raw_strict_cr_ids: Tuple[str, ...]
    raw_strict_ncr_ids: Tuple[str, ...]
    label_accounting: Mapping[str, str]
    cohort_eligibility_policy: str
    cohort_eligibility_exclusions: Mapping[str, str]
    cr_priority_orders: Mapping[str, Tuple[str, ...]]
    ncr_priority_orders: Mapping[str, Tuple[str, ...]]
    effective_leakage_blocks: Mapping[str, str]
    criterion_groups: Mapping[str, Mapping[str, str]]
    diversity_strata: Mapping[str, str]
    diversity_tiers: Mapping[str, str]
    balance_strata: Mapping[str, Tuple[str, ...]]
    group_criteria: Tuple[str, ...]
    ncr_pool_fractions: Tuple[str, ...]
    seeds: Tuple[int, ...]
    fractions: Tuple[float, float, float]
    diversity: str
    diversity_index_hash: str
    dataset_version: Optional[str]
    checker_view: Optional[str]
    input_hashes: Mapping[str, str]
    _receipt: Mapping[str, object] = field(repr=False, compare=False)
    _dataset: object = field(repr=False, compare=False)
    _authority_factory_token: object = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        for name in (
            "cohorts",
            "fixed_test_ids",
            "full_cr_ids",
            "full_ncr_ids",
            "raw_strict_cr_ids",
            "raw_strict_ncr_ids",
            "group_criteria",
            "ncr_pool_fractions",
            "seeds",
            "fractions",
        ):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        for name in (
            "cr_priority_orders",
            "ncr_priority_orders",
            "label_accounting",
            "cohort_eligibility_exclusions",
            "effective_leakage_blocks",
            "criterion_groups",
            "diversity_strata",
            "diversity_tiers",
            "balance_strata",
            "input_hashes",
            "_receipt",
        ):
            object.__setattr__(self, name, _freeze(getattr(self, name)))

    @property
    def pool_counts(self) -> Mapping[str, int]:
        """Return CR/NCR counts after the declared cohort-eligibility policy."""

        _require_benchmark_result(self)
        return MappingProxyType(
            {"CR": len(self.full_cr_ids), "NCR": len(self.full_ncr_ids)}
        )

    @property
    def raw_pool_counts(self) -> Mapping[str, int]:
        """Return complete-release strict CR/NCR counts before eligibility."""

        _require_benchmark_result(self)
        return MappingProxyType(
            {"CR": len(self.raw_strict_cr_ids), "NCR": len(self.raw_strict_ncr_ids)}
        )

    def receipt(self) -> Mapping[str, object]:
        _require_benchmark_result(self)
        return _jsonable(self._receipt)  # type: ignore[return-value]

    def data_split(
        self,
        train: object = None,
        val: object = None,
        test: object = None,
        include_full_cr_diagnostic: bool = True,
    ) -> "CRNCRBenchmarkSuite":
        """Assign fixed cohort members without changing any membership."""

        _require_benchmark_result(self)
        stored = self.fractions
        requested = (
            stored
            if train is None and val is None and test is None
            else _validate_fractions(
                stored[0] if train is None else train,
                stored[1] if val is None else val,
                stored[2] if test is None else test,
            )
        )
        if requested != stored:
            raise BenchmarkError(
                "cohort fixed-test membership was frozen for fractions {}; rebuild "
                "the cohorts before changing train/validation/test fractions".format(
                    stored
                )
            )
        return _split_cr_ncr_cohorts(self, include_full_cr_diagnostic)

    def write(
        self,
        output_directory: object,
        stem: str = "coremof_cr_ncr_cohorts",
        overwrite: bool = False,
    ) -> Path:
        _require_benchmark_result(self)
        stem = _safe_stem(stem)
        parent = Path(output_directory)
        parent.mkdir(parents=True, exist_ok=True)
        target = parent / stem
        staging = Path(tempfile.mkdtemp(prefix=".{}-".format(stem), dir=str(parent)))
        try:
            with (staging / "membership_manifest.csv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=(
                        "seed",
                        "requested_ncr_pool_fraction",
                        "structure_id",
                        "label",
                    ),
                )
                writer.writeheader()
                for cohort in self.cohorts:
                    cr = set(cohort.cr_ids)
                    for structure_id in cohort.structure_ids:
                        writer.writerow(
                            {
                                "seed": cohort.seed,
                                "requested_ncr_pool_fraction": cohort.requested_ncr_pool_fraction,
                                "structure_id": structure_id,
                                "label": "CR" if structure_id in cr else "NCR",
                            }
                        )
            with (staging / "label_accounting_manifest.csv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                eligible_ids = set(self.full_cr_ids).union(self.full_ncr_ids)
                fields = (
                    "structure_id",
                    "label",
                    "strict_pool_member",
                    "cohort_eligible",
                    "exclusion_reason",
                    "effective_leakage_block",
                )
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for structure_id in sorted(self.label_accounting):
                    label = self.label_accounting[structure_id]
                    writer.writerow(
                        {
                            "structure_id": structure_id,
                            "label": label,
                            "strict_pool_member": (
                                "true" if label in {"CR", "NCR"} else "false"
                            ),
                            "cohort_eligible": (
                                "true" if structure_id in eligible_ids else "false"
                            ),
                            "exclusion_reason": (
                                ""
                                if structure_id in eligible_ids
                                else (
                                    "BLOCK_NOT_SINGLE_STRICT_LABEL"
                                    if label in {"CR", "NCR"}
                                    else (
                                        "PASS_FAIL_MIXTURE"
                                        if label == "AMBIGUOUS"
                                        else "REQUIRED_CHECKER_NOT_AVAILABLE"
                                    )
                                )
                            ),
                            "effective_leakage_block": self.effective_leakage_blocks[
                                structure_id
                            ],
                        }
                    )
            with (staging / "receipt.json").open("w", encoding="utf-8") as handle:
                json.dump(
                    self.receipt(), handle, indent=2, sort_keys=True, allow_nan=False
                )
                handle.write("\n")
            files = sorted(
                path
                for path in staging.iterdir()
                if path.is_file() and path.name != "SHA256SUMS"
            )
            with (staging / "SHA256SUMS").open("w", encoding="utf-8") as handle:
                for path in files:
                    handle.write(
                        "{}  {}\n".format(
                            hashlib.sha256(path.read_bytes()).hexdigest(), path.name
                        )
                    )
            return _publish_directory(staging, target, overwrite)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise


def build_cr_ncr_cohorts(
    classified: object,
    ncr_pool_fractions: object = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    seeds: object = (42, 43, 44, 45, 46),
    total_size: str = "full_cr",
    train: object = 0.8,
    val: object = 0.1,
    test: object = 0.1,
    group_criteria: object = "priority_main",
    cohort_eligibility: object = None,
    diversity: str = "representative",
    test_policy: str = TEST_POLICY,
) -> CRNCRCohortSuite:
    """Freeze paired fixed-size strict five-checker CR/NCR memberships."""

    fractions = _validate_fractions(train, val, test)
    q_values = _fraction_values(ncr_pool_fractions)
    seed_values = _seed_values(seeds)
    if total_size != "full_cr":
        raise BenchmarkError("total_size must be 'full_cr'")
    if test_policy != TEST_POLICY:
        raise BenchmarkError("test_policy must be 'fixed_pure_cr'")
    if getattr(classified, "checker_view", None) != "5checker":
        raise BenchmarkError(
            "CR/NCR benchmarking requires the recomputed strict 5checker view"
        )
    _validate_classified_authority(classified, require_official=True)
    dataset = _dataset_from(classified)
    rows, _ = _row_index(dataset)
    labels = _label_mapping(classified)
    view_ids = tuple(getattr(classified, "structure_ids", ()))
    if set(view_ids) != set(rows) or set(labels) != set(rows):
        raise BenchmarkError(
            "CR/NCR benchmarking requires the unfiltered complete-release checker view"
        )
    raw_cr_ids = tuple(
        sorted(structure_id for structure_id in rows if labels[structure_id] == "CR")
    )
    raw_ncr_ids = tuple(
        sorted(structure_id for structure_id in rows if labels[structure_id] == "NCR")
    )
    complete_label_counts = Counter(labels.values())
    criteria = normalize_group_criteria(group_criteria)
    grouping = _grouping_context(dataset, criteria)
    complete_members_by_block: Dict[str, List[str]] = {}
    for structure_id, block in grouping.effective_blocks.items():
        complete_members_by_block.setdefault(block, []).append(structure_id)
    label_pure_blocks = {
        label: {
            block
            for block, members in complete_members_by_block.items()
            if members
            and all(labels[structure_id] == label for structure_id in members)
        }
        for label in ("CR", "NCR")
    }
    pure_cr_ids = tuple(
        value
        for value in raw_cr_ids
        if grouping.effective_blocks[value] in label_pure_blocks["CR"]
    )
    pure_ncr_ids = tuple(
        value
        for value in raw_ncr_ids
        if grouping.effective_blocks[value] in label_pure_blocks["NCR"]
    )
    pure_ids = set(pure_cr_ids).union(pure_ncr_ids)
    cohort_eligibility_exclusions = {
        structure_id: "BLOCK_NOT_SINGLE_STRICT_LABEL"
        for structure_id in tuple(raw_cr_ids) + tuple(raw_ncr_ids)
        if structure_id not in pure_ids
    }
    if cohort_eligibility is None:
        if cohort_eligibility_exclusions:
            excluded_counts = Counter(
                labels[structure_id] for structure_id in cohort_eligibility_exclusions
            )
            raise BenchmarkFeasibilityError(
                "the complete strict-row CR/NCR pools cannot be selected as whole "
                "effective leakage blocks: {} CR and {} NCR structures share a block "
                "with another label. Explicitly request cohort_eligibility={!r} for "
                "the label-pure-block sensitivity cohort; no partial block was "
                "selected.".format(
                    excluded_counts["CR"],
                    excluded_counts["NCR"],
                    LABEL_PURE_BLOCK_ELIGIBILITY,
                )
            )
        eligibility_policy = "complete_strict_pool_whole_effective_blocks"
        cr_ids = raw_cr_ids
        ncr_ids = raw_ncr_ids
    elif cohort_eligibility == LABEL_PURE_BLOCK_ELIGIBILITY:
        eligibility_policy = LABEL_PURE_BLOCK_ELIGIBILITY
        cr_ids = pure_cr_ids
        ncr_ids = pure_ncr_ids
    else:
        raise BenchmarkError(
            "cohort_eligibility must be null or {!r}".format(
                LABEL_PURE_BLOCK_ELIGIBILITY
            )
        )
    c_count = len(cr_ids)
    m_count = len(ncr_ids)
    if c_count == 0:
        raise BenchmarkFeasibilityError(
            "fixed-size full_cr cohorts require at least one strict five-checker CR structure"
        )
    if m_count > c_count:
        raise BenchmarkFeasibilityError(
            "fixed-size CR/NCR ladder infeasible: strict five-checker pool counts "
            "are C={} CR and M={} NCR, so M>C. The maximum feasible complete-NCR-"
            "pool fraction before reserving a clean test is {:.12g}; no count was "
            "capped, duplicated, or resized.".format(
                c_count, m_count, c_count / m_count if m_count else 1.0
            )
        )

    diversity_index = build_diversity_index(dataset, diversity=diversity)
    balance = _balance_values(tuple(rows), rows, labels, diversity_index)
    desired_test_count = _round_half_up(Decimal(str(fractions[2])) * Decimal(c_count))
    preliminary_fixed_test = _fixed_test_blocks(
        cr_ids,
        labels,
        grouping.effective_blocks,
        balance,
        desired_test_count,
    )
    test_count = len(preliminary_fixed_test)
    target_ncr_counts = tuple(_round_half_up(q * Decimal(m_count)) for q in q_values)
    if m_count > c_count - test_count:
        maximum = (
            Decimal(c_count - test_count) / Decimal(m_count) if m_count else Decimal(1)
        )
        raise BenchmarkFeasibilityError(
            "fixed-size CR/NCR ladder infeasible after reserving the common pure-CR "
            "test: C={} CR, M={} NCR, requested test count={}, achieved whole-block "
            "test count={}, leaving C-test_count={} replaceable CR. M exceeds that "
            "capacity. The maximum feasible NCR-pool fraction is {}; no count was "
            "capped, duplicated, or resized.".format(
                c_count,
                m_count,
                desired_test_count,
                test_count,
                c_count - test_count,
                format(maximum, "f"),
            )
        )
    # Select the fixed test jointly with one feasible CR-removal ladder.  Using
    # an independently chosen same-cardinality test can consume the only block
    # sizes that make later exact q transitions feasible.
    fixed_test_and_removal = _nested_exact_block_memberships(
        cr_ids,
        grouping.effective_blocks,
        balance,
        (test_count,) + tuple(test_count + value for value in target_ncr_counts),
        DIVERSITY_PROFILE_SEED,
        "fixed-test-plus-cr-removal-feasibility",
        "fixed test plus CR removal",
    )
    fixed_test = fixed_test_and_removal[test_count]
    fixed_test_set = set(fixed_test)
    cohorts = []
    cr_orders = {}
    ncr_orders = {}
    for seed in seed_values:
        remaining_cr_ids = tuple(
            value for value in cr_ids if value not in fixed_test_set
        )
        selected_ncr_by_count = _nested_exact_block_memberships(
            ncr_ids,
            grouping.effective_blocks,
            balance,
            target_ncr_counts,
            seed,
            "ncr-whole-block-priority",
            "NCR addition",
        )
        removed_cr_by_count = _nested_exact_block_memberships(
            remaining_cr_ids,
            grouping.effective_blocks,
            balance,
            target_ncr_counts,
            seed,
            "cr-whole-block-removal-priority",
            "CR removal",
        )
        cr_order = tuple(fixed_test) + _balanced_order(
            remaining_cr_ids, balance, seed, "cr-priority-reporting"
        )
        ncr_order = _balanced_order(ncr_ids, balance, seed, "ncr-priority-reporting")
        cr_orders[str(seed)] = cr_order
        ncr_orders[str(seed)] = ncr_order
        for q, selected_ncr_count in zip(q_values, target_ncr_counts):
            selected_cr_count = c_count - selected_ncr_count
            selected_ncr = selected_ncr_by_count[selected_ncr_count]
            removed_cr = set(removed_cr_by_count[selected_ncr_count])
            selected_cr = tuple(
                sorted(value for value in cr_ids if value not in removed_cr)
            )
            if (
                len(selected_ncr) != selected_ncr_count
                or len(selected_cr) != selected_cr_count
            ):
                raise RuntimeError(
                    "internal error: whole-block cohort counts differ from formula"
                )
            selected = tuple(sorted(selected_cr + selected_ncr))
            if len(selected) != c_count:
                raise RuntimeError("internal error: cohort size differs from C")
            whole_block_audit = _whole_block_membership_audit(
                selected, grouping.effective_blocks
            )
            if not whole_block_audit["passed"]:
                raise RuntimeError(
                    "internal error: cohort membership partially selected an "
                    "effective leakage block"
                )
            cohorts.append(
                CRNCRCohort(
                    seed=seed,
                    requested_ncr_pool_fraction=format(q, "f"),
                    selected_ncr_count=selected_ncr_count,
                    selected_cr_count=selected_cr_count,
                    actual_ncr_ratio=(selected_ncr_count / c_count),
                    structure_ids=selected,
                    cr_ids=selected_cr,
                    ncr_ids=selected_ncr,
                    whole_block_membership_audit=whole_block_audit,
                )
            )
    nesting_audit = []
    by_seed = {
        seed: [cohort for cohort in cohorts if cohort.seed == seed]
        for seed in seed_values
    }
    for seed, values in by_seed.items():
        values.sort(key=lambda value: Decimal(value.requested_ncr_pool_fraction))
        for previous, current in zip(values, values[1:]):
            nesting_audit.append(
                {
                    "seed": seed,
                    "from": previous.requested_ncr_pool_fraction,
                    "to": current.requested_ncr_pool_fraction,
                    "ncr_nested_addition": set(previous.ncr_ids).issubset(
                        current.ncr_ids
                    ),
                    "cr_reverse_nested_removal": set(current.cr_ids).issubset(
                        previous.cr_ids
                    ),
                }
            )
    if any(
        not item["ncr_nested_addition"] or not item["cr_reverse_nested_removal"]
        for item in nesting_audit
    ):
        raise RuntimeError("internal error: paired cohort nesting failed")
    full_endpoint = [
        cohort
        for cohort in cohorts
        if Decimal(cohort.requested_ncr_pool_fraction) == Decimal(1)
    ]
    if full_endpoint and any(
        set(cohort.ncr_ids) != set(ncr_ids) for cohort in full_endpoint
    ):
        raise RuntimeError("internal error: 100% endpoint omitted NCR structures")
    input_hashes = dict(getattr(dataset, "input_hashes", {}) or {})
    for item in diversity_index.input_receipts:
        input_hashes[str(item["release_path"])] = str(item["sha256"])
    block_sizes = Counter(grouping.effective_blocks.values())
    maximum_feasible_fraction = (
        min(Decimal(1), Decimal(c_count - test_count) / Decimal(m_count))
        if m_count
        else Decimal(1)
    )
    receipt = {
        "schema_version": "coremof-cr-ncr-cohort-suite-receipt/1.0",
        "benchmark_api_version": BENCHMARK_API_VERSION,
        "implementation": {
            "package": "CoREMOF-tools",
            "package_version": __version__,
            "source_sha256": dict(_source_hashes()),
        },
        "contract_definitions": _jsonable(generated_output_terminology_definitions()),
        "dataset_version": getattr(dataset, "dataset_version", None),
        "checker_view": "5checker",
        "checker_names": list(CHECKER_PRESETS["5checker"]),
        "strict_pool_definition": {
            "CR": "all five named checkers available and PASS",
            "NCR": "all five named checkers available and FAIL",
            "AMBIGUOUS": "all five checker votes available with both PASS and FAIL",
            "UNCHECKED": "at least one of the five required checkers is NOT_AVAILABLE",
            "execution_problem_behavior": "NOT_AVAILABLE non-vote; never FAIL or NCR",
        },
        "complete_release_label_accounting": {
            "counts": {
                label: complete_label_counts.get(label, 0)
                for label in ("CR", "NCR", "AMBIGUOUS", "UNCHECKED")
            },
            "membership_sha256": {
                label: _digest(
                    tuple(
                        sorted(
                            structure_id
                            for structure_id, value in labels.items()
                            if value == label
                        )
                    )
                )
                for label in ("CR", "NCR", "AMBIGUOUS", "UNCHECKED")
            },
            "companion_manifest": "label_accounting_manifest.csv",
            "non_strict_rows_enter_cohort": False,
        },
        "raw_strict_pool_counts": {
            "CR": len(raw_cr_ids),
            "NCR": len(raw_ncr_ids),
        },
        "cohort_eligibility_policy": {
            "name": eligibility_policy,
            "definition": (
                "A structure is cohort-eligible only when every complete-release "
                "member of its effective leakage block has the same strict "
                "five-checker CR or NCR label."
                if eligibility_policy == LABEL_PURE_BLOCK_ELIGIBILITY
                else "Every strict CR/NCR row is eligible because the complete "
                "release contains no strict row in a label-impure effective block."
            ),
            "scope": (
                "explicit label-pure-block sensitivity cohort; not the complete "
                "strict-row population"
                if eligibility_policy == LABEL_PURE_BLOCK_ELIGIBILITY
                else "complete strict-row population"
            ),
        },
        "eligible_pool_counts": {"C_CR": c_count, "M_NCR": m_count},
        "excluded_strict_membership": {
            "CR": sum(
                labels[structure_id] == "CR"
                for structure_id in cohort_eligibility_exclusions
            ),
            "NCR": sum(
                labels[structure_id] == "NCR"
                for structure_id in cohort_eligibility_exclusions
            ),
            "reason": "BLOCK_NOT_SINGLE_STRICT_LABEL",
            "membership_sha256": _digest(tuple(sorted(cohort_eligibility_exclusions))),
            "companion_manifest": "label_accounting_manifest.csv",
        },
        "strict_pool_membership_sha256": {
            "CR": _digest(cr_ids),
            "NCR": _digest(ncr_ids),
        },
        "raw_strict_pool_membership_sha256": {
            "CR": _digest(raw_cr_ids),
            "NCR": _digest(raw_ncr_ids),
        },
        "formula": {
            "selected_NCR": "round_half_up(q*M)",
            "selected_CR": "C-selected_NCR",
            "total_size": "C",
            "actual_NCR_ratio": "selected_NCR/C",
            "meaning_of_q_1": "use every NCR structure, not 100 percent NCR composition",
            "formula_pool_basis": "eligible_pool_counts",
        },
        "nested_whole_block_solver": {
            "name": "joint_transition_bin_assignment_v1",
            "definition": (
                "All requested q-transition deltas and one optional unused bin "
                "are solved together. Every effective leakage block is assigned "
                "once and size-one blocks fill exact residual capacities. A "
                "bounded exhaustive fallback distinguishes a search limit from "
                "mathematical infeasibility."
            ),
            "silent_count_deviation": False,
        },
        "ncr_pool_fractions": [format(value, "f") for value in q_values],
        "seeds": list(seed_values),
        "fractions": dict(zip(_PARTITIONS, fractions)),
        "total_size": total_size,
        "group_criteria": list(criteria),
        "group_criterion_definitions": {
            criterion: _jsonable(parent_method_definition(criterion))
            for criterion in criteria
        },
        "criterion_groups_sha256": {
            criterion: _digest(sorted(groups.items()))
            for criterion, groups in grouping.criterion_groups.items()
        },
        "leakage_guard": EFFECTIVE_LEAKAGE_POLICY,
        "effective_leakage_policy_definition": EFFECTIVE_LEAKAGE_POLICY_DEFINITION,
        "effective_leakage_blocks": {
            "sha256": _digest(sorted(grouping.effective_blocks.items())),
            "block_count": len(block_sizes),
            "maximum_complete_release_block_size": max(block_sizes.values())
            if block_sizes
            else 0,
        },
        "test_policy": test_policy,
        "test_policy_definition": TEST_POLICY_DEFINITION,
        "feasibility_audit": {
            "M_le_C": m_count <= c_count,
            "M_le_C_minus_test_count": m_count <= c_count - test_count,
            "replaceable_CR_count": c_count - test_count,
            "maximum_feasible_NCR_pool_fraction": format(
                maximum_feasible_fraction, "f"
            ),
            "silent_cap_duplicate_or_resize": False,
        },
        "fixed_test_selection": {
            "profile_seed": DIVERSITY_PROFILE_SEED,
            "cardinality_optimizer": "global nearest positive subset sum; equal-distance tie prefers below requested",
            "membership_solver": (
                "joint_transition_bin_assignment_v1 with the CR-removal ladder; "
                "the test is not selected independently"
            ),
            "requested_count": desired_test_count,
            "achieved_whole_block_count": test_count,
            "pure_CR_complete_release_blocks": True,
            "ids": list(fixed_test),
        },
        "diversity": diversity,
        "diversity_profile": _jsonable(diversity_index.profile),
        "diversity_index_sha256": diversity_index.digest,
        "target_columns_consumed": [],
        "input_hashes": dict(sorted(input_hashes.items())),
        "release_binding": _jsonable(_release_binding(dataset)),
        "official_split": False,
        "publication_authorized": False,
        "benchmark_contract_compliance": "PASS",
        "cohort_count": len(cohorts),
        "cohorts": [
            {
                "seed": cohort.seed,
                "requested_ncr_pool_fraction": cohort.requested_ncr_pool_fraction,
                "selected_CR": cohort.selected_cr_count,
                "selected_NCR": cohort.selected_ncr_count,
                "total": len(cohort.structure_ids),
                "actual_NCR_ratio": cohort.actual_ncr_ratio,
                "membership_sha256": _digest(cohort.structure_ids),
                "whole_block_membership_audit": _jsonable(
                    cohort.whole_block_membership_audit
                ),
            }
            for cohort in cohorts
        ],
        "nesting_audit": nesting_audit,
        "all_ncr_present_at_full_endpoint": all(
            set(cohort.ncr_ids) == set(ncr_ids) for cohort in full_endpoint
        )
        if full_endpoint
        else None,
        "all_zero_partially_selected_effective_blocks": all(
            cohort.whole_block_membership_audit["passed"] for cohort in cohorts
        ),
    }
    return _register_benchmark_result(
        CRNCRCohortSuite(
            cohorts=tuple(cohorts),
            fixed_test_ids=fixed_test,
            full_cr_ids=cr_ids,
            full_ncr_ids=ncr_ids,
            raw_strict_cr_ids=raw_cr_ids,
            raw_strict_ncr_ids=raw_ncr_ids,
            label_accounting=labels,
            cohort_eligibility_policy=eligibility_policy,
            cohort_eligibility_exclusions=cohort_eligibility_exclusions,
            cr_priority_orders=cr_orders,
            ncr_priority_orders=ncr_orders,
            effective_leakage_blocks=grouping.effective_blocks,
            criterion_groups=grouping.criterion_groups,
            diversity_strata=diversity_index.strata_by_id,
            diversity_tiers=diversity_index.tier_by_id,
            balance_strata=balance,
            group_criteria=criteria,
            ncr_pool_fractions=tuple(format(value, "f") for value in q_values),
            seeds=seed_values,
            fractions=fractions,
            diversity=diversity,
            diversity_index_hash=diversity_index.digest,
            dataset_version=getattr(dataset, "dataset_version", None),
            checker_view="5checker",
            input_hashes=input_hashes,
            _receipt=receipt,
            _dataset=dataset,
            _authority_factory_token=_BENCHMARK_RESULT_FACTORY_TOKEN,
        ),
        dataset,
    )


@dataclass(frozen=True)
class CRNCRBenchmarkRun:
    """One paired cohort with a fixed pure-CR test assignment."""

    seed: int
    requested_ncr_pool_fraction: str
    actual_ncr_ratio: float
    train_ids: Tuple[str, ...]
    validation_ids: Tuple[str, ...]
    test_ids: Tuple[str, ...]
    cr_ids: Tuple[str, ...]
    ncr_ids: Tuple[str, ...]
    assignments: Mapping[str, str]
    requested_counts: Mapping[str, int]
    achieved_counts: Mapping[str, int]
    count_deviations: Mapping[str, int]
    label_counts_by_partition: Mapping[str, Mapping[str, int]]
    requested_label_counts_by_partition: Mapping[str, Mapping[str, float]]
    label_count_deviations_by_partition: Mapping[str, Mapping[str, float]]
    leakage_audit: Mapping[str, object]
    assignment_digest: str
    full_cr_diagnostic: Optional[Mapping[str, object]]

    def __post_init__(self) -> None:
        for name in ("train_ids", "validation_ids", "test_ids", "cr_ids", "ncr_ids"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        for name in (
            "assignments",
            "requested_counts",
            "achieved_counts",
            "count_deviations",
            "label_counts_by_partition",
            "requested_label_counts_by_partition",
            "label_count_deviations_by_partition",
            "leakage_audit",
        ):
            object.__setattr__(self, name, _freeze(getattr(self, name)))
        if self.full_cr_diagnostic is not None:
            object.__setattr__(
                self, "full_cr_diagnostic", _freeze(self.full_cr_diagnostic)
            )

    @property
    def run_key(self) -> str:
        return "seed{}_q{}".format(
            self.seed, self.requested_ncr_pool_fraction.replace(".", "p")
        )


@dataclass(frozen=True)
class CRNCRBenchmarkSuite:
    """Immutable paired benchmark suite with one common clean test manifest."""

    runs: Tuple[CRNCRBenchmarkRun, ...]
    fixed_test_ids: Tuple[str, ...]
    full_cr_diagnostic_ids: Tuple[str, ...]
    eligible_cr_ids: Tuple[str, ...]
    eligible_ncr_ids: Tuple[str, ...]
    label_accounting: Mapping[str, str]
    cohort_eligibility_exclusions: Mapping[str, str]
    effective_leakage_blocks: Mapping[str, str]
    criterion_groups: Mapping[str, Mapping[str, str]]
    diversity_strata: Mapping[str, str]
    diversity_tiers: Mapping[str, str]
    group_criteria: Tuple[str, ...]
    seeds: Tuple[int, ...]
    ncr_pool_fractions: Tuple[str, ...]
    fractions: Tuple[float, float, float]
    diversity_index_hash: str
    dataset_version: Optional[str]
    checker_view: str
    input_hashes: Mapping[str, str]
    official_split: bool
    assignment_digest: str
    _receipt: Mapping[str, object] = field(repr=False, compare=False)
    _dataset: object = field(repr=False, compare=False)
    _authority_factory_token: object = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        for name in (
            "runs",
            "fixed_test_ids",
            "full_cr_diagnostic_ids",
            "eligible_cr_ids",
            "eligible_ncr_ids",
            "group_criteria",
            "seeds",
            "ncr_pool_fractions",
            "fractions",
        ):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        for name in (
            "effective_leakage_blocks",
            "label_accounting",
            "cohort_eligibility_exclusions",
            "criterion_groups",
            "diversity_strata",
            "diversity_tiers",
            "input_hashes",
            "_receipt",
        ):
            object.__setattr__(self, name, _freeze(getattr(self, name)))
        if self.official_split is not False:
            raise BenchmarkError("benchmark suites must remain exploratory")

    def receipt(self) -> Mapping[str, object]:
        _require_benchmark_result(self)
        return _jsonable(self._receipt)  # type: ignore[return-value]

    def attach_targets(
        self, target_data: object, missing: str = "keep", **kwargs: object
    ):
        """Attach targets to every frozen run without selection or resplitting."""

        _require_benchmark_result(self)
        from CoREMOF.attachments import attach_targets

        return attach_targets(
            self, target_data, dataset=self._dataset, missing=missing, **kwargs
        )

    def write(
        self,
        output_directory: object,
        stem: str = "coremof_cr_ncr_benchmark",
        overwrite: bool = False,
    ) -> Path:
        """Transactionally write manifests, per-run views, receipt, and checksums."""

        _require_benchmark_result(self)
        stem = _safe_stem(stem)
        parent = Path(output_directory)
        parent.mkdir(parents=True, exist_ok=True)
        target = parent / stem
        staging = Path(tempfile.mkdtemp(prefix=".{}-".format(stem), dir=str(parent)))
        try:
            runs_directory = staging / "runs"
            runs_directory.mkdir()
            index_rows = []
            membership_path = staging / "membership_manifest.csv"
            with membership_path.open("w", encoding="utf-8", newline="") as handle:
                fields = (
                    "run_key",
                    "seed",
                    "requested_ncr_pool_fraction",
                    "actual_ncr_ratio",
                    "structure_id",
                    "label",
                    "partition",
                    "effective_leakage_block",
                    "diversity_tier",
                    "diversity_stratum",
                )
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for run in self.runs:
                    cr = set(run.cr_ids)
                    run_rows = []
                    for structure_id in sorted(run.assignments):
                        row = {
                            "run_key": run.run_key,
                            "seed": run.seed,
                            "requested_ncr_pool_fraction": run.requested_ncr_pool_fraction,
                            "actual_ncr_ratio": format(run.actual_ncr_ratio, ".17g"),
                            "structure_id": structure_id,
                            "label": "CR" if structure_id in cr else "NCR",
                            "partition": run.assignments[structure_id],
                            "effective_leakage_block": self.effective_leakage_blocks[
                                structure_id
                            ],
                            "diversity_tier": self.diversity_tiers[structure_id],
                            "diversity_stratum": self.diversity_strata[structure_id],
                        }
                        writer.writerow(row)
                        run_rows.append(row)
                    run_path = runs_directory / (run.run_key + ".csv")
                    with run_path.open("w", encoding="utf-8", newline="") as run_handle:
                        run_writer = csv.DictWriter(run_handle, fieldnames=fields)
                        run_writer.writeheader()
                        run_writer.writerows(run_rows)
                    index_rows.append(
                        {
                            "run_key": run.run_key,
                            "seed": run.seed,
                            "requested_ncr_pool_fraction": run.requested_ncr_pool_fraction,
                            "actual_ncr_ratio": run.actual_ncr_ratio,
                            "selected_CR": len(run.cr_ids),
                            "selected_NCR": len(run.ncr_ids),
                            "train": len(run.train_ids),
                            "validation": len(run.validation_ids),
                            "test": len(run.test_ids),
                            "assignment_sha256": run.assignment_digest,
                            "run_manifest": "runs/{}.csv".format(run.run_key),
                        }
                    )
            with (staging / "fixed_test_manifest.csv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("structure_id", "label", "effective_leakage_block"),
                )
                writer.writeheader()
                for structure_id in self.fixed_test_ids:
                    writer.writerow(
                        {
                            "structure_id": structure_id,
                            "label": "CR",
                            "effective_leakage_block": self.effective_leakage_blocks[
                                structure_id
                            ],
                        }
                    )
            if self.full_cr_diagnostic_ids:
                with (staging / "full_cr_diagnostic_manifest.csv").open(
                    "w", encoding="utf-8", newline=""
                ) as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=("structure_id", "label", "effective_leakage_block"),
                    )
                    writer.writeheader()
                    for structure_id in self.full_cr_diagnostic_ids:
                        writer.writerow(
                            {
                                "structure_id": structure_id,
                                "label": "CR",
                                "effective_leakage_block": self.effective_leakage_blocks[
                                    structure_id
                                ],
                            }
                        )
            with (staging / "label_accounting_manifest.csv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                eligible_ids = set(self.eligible_cr_ids).union(self.eligible_ncr_ids)
                fields = (
                    "structure_id",
                    "label",
                    "strict_pool_member",
                    "cohort_eligible",
                    "exclusion_reason",
                    "effective_leakage_block",
                )
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for structure_id in sorted(self.label_accounting):
                    label = self.label_accounting[structure_id]
                    writer.writerow(
                        {
                            "structure_id": structure_id,
                            "label": label,
                            "strict_pool_member": (
                                "true" if label in {"CR", "NCR"} else "false"
                            ),
                            "cohort_eligible": (
                                "true" if structure_id in eligible_ids else "false"
                            ),
                            "exclusion_reason": (
                                ""
                                if structure_id in eligible_ids
                                else (
                                    self.cohort_eligibility_exclusions.get(
                                        structure_id, "BLOCK_NOT_SINGLE_STRICT_LABEL"
                                    )
                                    if label in {"CR", "NCR"}
                                    else (
                                        "PASS_FAIL_MIXTURE"
                                        if label == "AMBIGUOUS"
                                        else "REQUIRED_CHECKER_NOT_AVAILABLE"
                                    )
                                )
                            ),
                            "effective_leakage_block": self.effective_leakage_blocks[
                                structure_id
                            ],
                        }
                    )
            with (staging / "suite_index.json").open("w", encoding="utf-8") as handle:
                json.dump(index_rows, handle, indent=2, sort_keys=True, allow_nan=False)
                handle.write("\n")
            with (staging / "receipt.json").open("w", encoding="utf-8") as handle:
                json.dump(
                    self.receipt(), handle, indent=2, sort_keys=True, allow_nan=False
                )
                handle.write("\n")
            files = sorted(
                path
                for path in staging.rglob("*")
                if path.is_file() and path.name != "SHA256SUMS"
            )
            with (staging / "SHA256SUMS").open("w", encoding="utf-8") as handle:
                for path in files:
                    digest = hashlib.sha256(path.read_bytes()).hexdigest()
                    handle.write(
                        "{}  {}\n".format(digest, path.relative_to(staging).as_posix())
                    )
            return _publish_directory(staging, target, overwrite)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise


def _split_cr_ncr_cohorts(
    cohorts: CRNCRCohortSuite,
    include_full_cr_diagnostic: bool,
) -> CRNCRBenchmarkSuite:
    if type(include_full_cr_diagnostic) is not bool:
        raise TypeError("include_full_cr_diagnostic must be a boolean")
    train_fraction, val_fraction, _ = cohorts.fractions
    non_test_total = train_fraction + val_fraction
    if non_test_total <= 0.0:
        raise BenchmarkError(
            "fixed_pure_cr requires a nonzero train or validation fraction"
        )
    normalized = {
        "train": train_fraction / non_test_total,
        "validation": val_fraction / non_test_total,
    }
    test_set = set(cohorts.fixed_test_ids)
    all_pool = tuple(
        sorted(
            set(cohorts.full_cr_ids).union(cohorts.full_ncr_ids).difference(test_set)
        )
    )
    blocks_with_test = {
        cohorts.effective_leakage_blocks[structure_id] for structure_id in test_set
    }
    if any(
        cohorts.effective_leakage_blocks[structure_id] in blocks_with_test
        for structure_id in all_pool
    ):
        raise RuntimeError(
            "internal error: fixed test did not consume whole leakage blocks"
        )
    partition_by_seed = {}
    for seed in cohorts.seeds:
        seed_cohorts = tuple(
            sorted(
                (cohort for cohort in cohorts.cohorts if cohort.seed == seed),
                key=lambda cohort: Decimal(cohort.requested_ncr_pool_fraction),
            )
        )
        by_block = _assign_paired_blocks(
            all_pool,
            seed_cohorts,
            cohorts.effective_leakage_blocks,
            cohorts.balance_strata,
            normalized,
            seed,
        )
        partition_by_seed[seed] = {
            structure_id: by_block[cohorts.effective_leakage_blocks[structure_id]]
            for structure_id in all_pool
        }
    runs = []
    target_train = _round_half_up(
        Decimal(str(cohorts.fractions[0])) * Decimal(len(cohorts.full_cr_ids))
    )
    target_test = _round_half_up(
        Decimal(str(cohorts.fractions[2])) * Decimal(len(cohorts.full_cr_ids))
    )
    target_validation = len(cohorts.full_cr_ids) - target_train - target_test
    if target_validation < 0:
        target_validation = _round_half_up(
            Decimal(str(cohorts.fractions[1])) * Decimal(len(cohorts.full_cr_ids))
        )
        target_train = len(cohorts.full_cr_ids) - target_validation - target_test
    requested_counts = {
        "train": target_train,
        "validation": target_validation,
        "test": target_test,
        "total": len(cohorts.full_cr_ids),
    }
    for cohort in cohorts.cohorts:
        assignments = {
            structure_id: (
                "test"
                if structure_id in test_set
                else partition_by_seed[cohort.seed][structure_id]
            )
            for structure_id in cohort.structure_ids
        }
        ids_by_partition = {
            partition: tuple(
                sorted(
                    value
                    for value, assigned in assignments.items()
                    if assigned == partition
                )
            )
            for partition in _PARTITIONS
        }
        audit = _leakage_audit(assignments, cohorts.effective_leakage_blocks)
        if not audit["passed"]:
            raise RuntimeError(
                "internal error: benchmark effective block crossed partitions"
            )
        achieved = {
            partition: len(ids_by_partition[partition]) for partition in _PARTITIONS
        }
        achieved["total"] = len(assignments)
        deviations = {
            partition: achieved[partition] - requested_counts[partition]
            for partition in ("train", "validation", "test", "total")
        }
        cr_set = set(cohort.cr_ids)
        label_counts_by_partition = {
            partition: {
                "CR": sum(value in cr_set for value in ids_by_partition[partition]),
                "NCR": sum(
                    value not in cr_set for value in ids_by_partition[partition]
                ),
            }
            for partition in _PARTITIONS
        }
        remaining_label_counts = {
            "CR": len(cohort.cr_ids) - len(cohorts.fixed_test_ids),
            "NCR": len(cohort.ncr_ids),
        }
        requested_label_counts_by_partition = {
            "train": {
                label: normalized["train"] * count
                for label, count in remaining_label_counts.items()
            },
            "validation": {
                label: normalized["validation"] * count
                for label, count in remaining_label_counts.items()
            },
            "test": {
                "CR": float(len(cohorts.fixed_test_ids)),
                "NCR": 0.0,
            },
        }
        label_count_deviations_by_partition = {
            partition: {
                label: (
                    label_counts_by_partition[partition][label]
                    - requested_label_counts_by_partition[partition][label]
                )
                for label in ("CR", "NCR")
            }
            for partition in _PARTITIONS
        }
        diagnostic = None
        if include_full_cr_diagnostic:
            train_ids = set(ids_by_partition["train"])
            train_blocks = {
                cohorts.effective_leakage_blocks[value] for value in train_ids
            }
            exact_overlap = tuple(
                sorted(set(cohorts.raw_strict_cr_ids).intersection(train_ids))
            )
            block_overlap = tuple(
                sorted(
                    value
                    for value in cohorts.raw_strict_cr_ids
                    if cohorts.effective_leakage_blocks[value] in train_blocks
                )
            )
            diagnostic = {
                "definition": (
                    "supplementary prediction view over the complete raw strict-CR "
                    "pool, including CR rows excluded by the label-pure-block "
                    "cohort policy; "
                    "it includes seen structures and is not the independent fixed test"
                ),
                "structure_count": len(cohorts.raw_strict_cr_ids),
                "exact_training_overlap_count": len(exact_overlap),
                "exact_training_overlap_ids": exact_overlap,
                "same_effective_block_training_overlap_count": len(block_overlap),
                "same_effective_block_training_overlap_ids": block_overlap,
            }
        assignment_digest = _assignment_digest(assignments, {})
        runs.append(
            CRNCRBenchmarkRun(
                seed=cohort.seed,
                requested_ncr_pool_fraction=cohort.requested_ncr_pool_fraction,
                actual_ncr_ratio=cohort.actual_ncr_ratio,
                train_ids=ids_by_partition["train"],
                validation_ids=ids_by_partition["validation"],
                test_ids=ids_by_partition["test"],
                cr_ids=cohort.cr_ids,
                ncr_ids=cohort.ncr_ids,
                assignments=assignments,
                requested_counts=requested_counts,
                achieved_counts=achieved,
                count_deviations=deviations,
                label_counts_by_partition=label_counts_by_partition,
                requested_label_counts_by_partition=requested_label_counts_by_partition,
                label_count_deviations_by_partition=label_count_deviations_by_partition,
                leakage_audit=audit,
                assignment_digest=assignment_digest,
                full_cr_diagnostic=diagnostic,
            )
        )
    # Persistent selected structures must retain their preassigned partition
    # within one seed as q changes.
    stability_audit = []
    for seed in cohorts.seeds:
        seed_runs = [run for run in runs if run.seed == seed]
        stable = True
        for left_index, left in enumerate(seed_runs):
            for right in seed_runs[left_index + 1 :]:
                common = set(left.assignments).intersection(right.assignments)
                if any(
                    left.assignments[value] != right.assignments[value]
                    for value in common
                ):
                    stable = False
                    break
        stability_audit.append(
            {"seed": seed, "persistent_partition_assignments_stable": stable}
        )
    if not all(
        item["persistent_partition_assignments_stable"] for item in stability_audit
    ):
        raise RuntimeError("internal error: paired partition stability failed")
    fixed_test_digest = _digest(cohorts.fixed_test_ids)
    suite_assignment_digest = _digest(
        [(run.run_key, run.assignment_digest) for run in runs]
    )
    receipt = {
        "schema_version": "coremof-cr-ncr-benchmark-suite-receipt/1.0",
        "benchmark_api_version": BENCHMARK_API_VERSION,
        "implementation": {
            "package": "CoREMOF-tools",
            "package_version": __version__,
            "source_sha256": dict(_source_hashes()),
        },
        "cohort_receipt": cohorts.receipt(),
        "complete_release_label_accounting": cohorts.receipt()[
            "complete_release_label_accounting"
        ],
        "dataset_version": cohorts.dataset_version,
        "checker_view": cohorts.checker_view,
        "input_hashes": dict(cohorts.input_hashes),
        "official_split": False,
        "publication_authorized": False,
        "benchmark_contract_compliance": "PASS",
        "test_policy": TEST_POLICY,
        "test_policy_definition": TEST_POLICY_DEFINITION,
        "fixed_test_sha256": fixed_test_digest,
        "fixed_test_ids": list(cohorts.fixed_test_ids),
        "fractions": dict(zip(_PARTITIONS, cohorts.fractions)),
        "requested_counts": requested_counts,
        "paired_partition_assignment_profile": {
            "definition": (
                "assign each non-test effective leakage block once per seed; "
                "jointly minimize count, strict-label composition, and target-free "
                "balance-stratum deviations over every NCR ladder level"
            ),
            "strict_label_composition_weight": PAIRED_LABEL_BALANCE_WEIGHT,
            "maximum_improvement_rounds": max(16, 8 * len(cohorts.ncr_pool_fractions)),
            "swap_candidate_limit_per_side": PAIRED_SWAP_CANDIDATE_COUNT,
        },
        "run_count": len(runs),
        "runs": [
            {
                "run_key": run.run_key,
                "seed": run.seed,
                "requested_ncr_pool_fraction": run.requested_ncr_pool_fraction,
                "actual_ncr_ratio": run.actual_ncr_ratio,
                "selected_CR": len(run.cr_ids),
                "selected_NCR": len(run.ncr_ids),
                "achieved_counts": dict(run.achieved_counts),
                "count_deviations": dict(run.count_deviations),
                "label_counts_by_partition": _jsonable(run.label_counts_by_partition),
                "requested_label_counts_by_partition": _jsonable(
                    run.requested_label_counts_by_partition
                ),
                "label_count_deviations_by_partition": _jsonable(
                    run.label_count_deviations_by_partition
                ),
                "leakage_audit": dict(run.leakage_audit),
                "assignment_sha256": run.assignment_digest,
                "full_cr_diagnostic": _jsonable(run.full_cr_diagnostic),
            }
            for run in runs
        ],
        "persistent_partition_stability_audit": stability_audit,
        "all_zero_crossed_effective_blocks": all(
            run.leakage_audit["passed"] for run in runs
        ),
        "common_test_identical_across_runs": all(
            run.test_ids == cohorts.fixed_test_ids for run in runs
        ),
        "suite_assignment_sha256": suite_assignment_digest,
    }
    return _register_benchmark_result(
        CRNCRBenchmarkSuite(
            runs=tuple(runs),
            fixed_test_ids=cohorts.fixed_test_ids,
            full_cr_diagnostic_ids=(
                cohorts.raw_strict_cr_ids if include_full_cr_diagnostic else ()
            ),
            eligible_cr_ids=cohorts.full_cr_ids,
            eligible_ncr_ids=cohorts.full_ncr_ids,
            label_accounting=cohorts.label_accounting,
            cohort_eligibility_exclusions=cohorts.cohort_eligibility_exclusions,
            effective_leakage_blocks=cohorts.effective_leakage_blocks,
            criterion_groups=cohorts.criterion_groups,
            diversity_strata=cohorts.diversity_strata,
            diversity_tiers=cohorts.diversity_tiers,
            group_criteria=cohorts.group_criteria,
            seeds=cohorts.seeds,
            ncr_pool_fractions=cohorts.ncr_pool_fractions,
            fractions=cohorts.fractions,
            diversity_index_hash=cohorts.diversity_index_hash,
            dataset_version=cohorts.dataset_version,
            checker_view=str(cohorts.checker_view),
            input_hashes=cohorts.input_hashes,
            official_split=False,
            assignment_digest=suite_assignment_digest,
            _receipt=receipt,
            _dataset=cohorts._dataset,
            _authority_factory_token=_BENCHMARK_RESULT_FACTORY_TOKEN,
        ),
        cohorts._dataset,
    )


def build_cr_ncr_benchmark(
    classified: object,
    ncr_pool_fractions: object = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    seeds: object = (42, 43, 44, 45, 46),
    total_size: str = "full_cr",
    train: object = 0.8,
    val: object = 0.1,
    test: object = 0.1,
    group_criteria: object = "priority_main",
    cohort_eligibility: object = None,
    diversity: str = "representative",
    test_policy: str = TEST_POLICY,
    include_full_cr_diagnostic: bool = True,
) -> CRNCRBenchmarkSuite:
    """Build cohorts and paired partitions in one convenience call."""

    cohort_suite = build_cr_ncr_cohorts(
        classified,
        ncr_pool_fractions=ncr_pool_fractions,
        seeds=seeds,
        total_size=total_size,
        train=train,
        val=val,
        test=test,
        group_criteria=group_criteria,
        cohort_eligibility=cohort_eligibility,
        diversity=diversity,
        test_policy=test_policy,
    )
    return cohort_suite.data_split(
        include_full_cr_diagnostic=include_full_cr_diagnostic
    )


def _classified_data_split(
    self,
    group_criteria: object = "priority_main",
    train: object = 0.8,
    val: object = 0.1,
    test: object = 0.1,
    leakage_guard: str = EFFECTIVE_LEAKAGE_POLICY,
    diversity: str = "representative",
    random_state: object = 42,
    labels: object = ("CR", "NCR"),
    sources: object = None,
    variants: object = None,
    metals: object = None,
    structure_ids: object = None,
    official: bool = False,
) -> DataSplitResult:
    """Bound convenience method for ``ClassifiedDataset``."""

    return data_split(
        self,
        group_criteria=group_criteria,
        train=train,
        val=val,
        test=test,
        leakage_guard=leakage_guard,
        diversity=diversity,
        random_state=random_state,
        labels=labels,
        sources=sources,
        variants=variants,
        metals=metals,
        structure_ids=structure_ids,
        official=official,
    )


def _classified_build_cr_ncr_cohorts(
    self,
    ncr_pool_fractions: object = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    seeds: object = (42, 43, 44, 45, 46),
    total_size: str = "full_cr",
    train: object = 0.8,
    val: object = 0.1,
    test: object = 0.1,
    group_criteria: object = "priority_main",
    cohort_eligibility: object = None,
    diversity: str = "representative",
    test_policy: str = TEST_POLICY,
) -> CRNCRCohortSuite:
    """Bound convenience method freezing strict CR/NCR memberships."""

    return build_cr_ncr_cohorts(
        self,
        ncr_pool_fractions=ncr_pool_fractions,
        seeds=seeds,
        total_size=total_size,
        train=train,
        val=val,
        test=test,
        group_criteria=group_criteria,
        cohort_eligibility=cohort_eligibility,
        diversity=diversity,
        test_policy=test_policy,
    )


def _classified_build_cr_ncr_benchmark(
    self,
    ncr_pool_fractions: object = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    seeds: object = (42, 43, 44, 45, 46),
    total_size: str = "full_cr",
    train: object = 0.8,
    val: object = 0.1,
    test: object = 0.1,
    group_criteria: object = "priority_main",
    cohort_eligibility: object = None,
    diversity: str = "representative",
    test_policy: str = TEST_POLICY,
    include_full_cr_diagnostic: bool = True,
) -> CRNCRBenchmarkSuite:
    """Bound convenience method building the paired benchmark suite."""

    return build_cr_ncr_benchmark(
        self,
        ncr_pool_fractions=ncr_pool_fractions,
        seeds=seeds,
        total_size=total_size,
        train=train,
        val=val,
        test=test,
        group_criteria=group_criteria,
        cohort_eligibility=cohort_eligibility,
        diversity=diversity,
        test_policy=test_policy,
        include_full_cr_diagnostic=include_full_cr_diagnostic,
    )


def _install_classified_dataset_extensions() -> None:
    """Install additive APIs without changing legacy hash-bound source files.

    Historical split receipts bind the exact bytes of ``dataset.py`` and
    ``splitters.py``.  Keeping these conveniences in this separately versioned
    module preserves those receipt bytes while normal package imports still
    expose the documented methods on ``ClassifiedDataset``.
    """

    from CoREMOF.dataset import ClassifiedDataset

    additions = {
        "data_split": _classified_data_split,
        "build_cr_ncr_cohorts": _classified_build_cr_ncr_cohorts,
        "build_cr_ncr_benchmark": _classified_build_cr_ncr_benchmark,
    }
    for name, method in additions.items():
        existing = getattr(ClassifiedDataset, name, None)
        if existing is not None and existing is not method:
            raise RuntimeError(
                "ClassifiedDataset already defines an incompatible {} method".format(
                    name
                )
            )
        setattr(ClassifiedDataset, name, method)


_install_classified_dataset_extensions()


__all__ = [
    "BENCHMARK_API_VERSION",
    "BENCHMARK_NUMPY_VERSION",
    "BENCHMARK_SKLEARN_VERSION",
    "BENCHMARK_SCIPY_VERSION",
    "BENCHMARK_JOBLIB_VERSION",
    "BENCHMARK_THREADPOOLCTL_VERSION",
    "LABEL_PURE_BLOCK_ELIGIBILITY",
    "BenchmarkDependencyError",
    "BenchmarkError",
    "BenchmarkFeasibilityError",
    "CRNCRBenchmarkRun",
    "CRNCRBenchmarkSuite",
    "CRNCRCohort",
    "CRNCRCohortSuite",
    "DataSplitResult",
    "DiversityIndex",
    "available_group_criteria",
    "build_cr_ncr_benchmark",
    "build_cr_ncr_cohorts",
    "build_diversity_index",
    "data_split",
    "normalize_group_criteria",
]
