"""Hash-bound source subsets with grouping computed on the complete release.

The contract is a delivery artifact, not a digital signature or publication
permission. A receiver supplies its expected checksum from a trusted handoff.
Only the private factory can authenticate a loaded in-memory projection.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from types import MappingProxyType
from typing import Mapping

from .dataset import (
    CoREMOFDataset, ReleaseValidationError, StructureRecord, _DATASET_GENERATIONS,
    _ParentGroupView, _capture_release_file, _deep_freeze, _index_unique,
    _parent_prefixes, _read_csv, _read_json_object, _register_dataset_generation,
    _reject_retired_or_reserved_keys, _require_columns, _require_dataset_generation,
    _validate_cif_manifest, _validate_dataset_info, _validate_metadata_identity,
    _validate_parent_methods, _validate_parent_rows, _validate_published_labels,
    _LEGACY_TARGET_COLUMNS, _validate_target_free_layout,
    _without_legacy_target_information,
)
from .labels import CHECKER_COLUMNS, CHECKER_PRESETS, classify_checker_row

SCHEMA = "coremof-source-projection/1.0"
_LABELS = ("CR", "NCR", "AMBIGUOUS", "UNCHECKED")
_SHA = re.compile(r"[0-9a-f]{64}")
_CORE_FILES = ("metadata/metadata.csv", "parent_groups/parent_groups.csv",
               "parent_groups/parent_group_methods.json", "manifests/cif_manifest.csv")


def _plain(value):
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    return value


def _digest(value):
    return hashlib.sha256(json.dumps(_plain(value), sort_keys=True, separators=(",", ":"),
                                      ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def _keys(value, keys, description):
    if not isinstance(value, Mapping) or set(value) != set(keys):
        raise ReleaseValidationError("invalid {} fields".format(description))


def _hash(value):
    if not isinstance(value, str) or not _SHA.fullmatch(value):
        raise ReleaseValidationError("expected an exact lowercase SHA-256")
    return value


def _count(value):
    if type(value) is not int or value < 0:
        raise ReleaseValidationError("counts must be nonnegative integers")
    return value


def _exact_map(value, ids, description):
    if not isinstance(value, Mapping) or set(value) != set(ids):
        raise ReleaseValidationError("{} membership differs".format(description))
    if any(not isinstance(v, str) or not v or v != v.strip() for v in value.values()):
        raise ReleaseValidationError("{} has invalid values".format(description))


def _tables(root):
    snapshots = {name: _capture_release_file(root / name, name) for name in _CORE_FILES}
    tables = {}
    fields = {}
    for name in (_CORE_FILES[0], _CORE_FILES[1], _CORE_FILES[3]):
        columns, rows = _read_csv(snapshots[name], name)
        indexed, order = _index_unique(rows, root / name)
        tables[name] = (indexed, order)
        fields[name] = columns
    ids = set(tables[_CORE_FILES[0]][0])
    if not ids or any(set(indexed) != ids for indexed, _ in tables.values()):
        raise ReleaseValidationError("source tables require one identical nonempty ID set")
    _require_columns(fields[_CORE_FILES[0]], (
        "structure_id", "cif_file", "source_database", "source_id", "structure_variant",
        "metal_elements", *CHECKER_COLUMNS.values(), "label_3checker", "label_4checker", "label_5checker",
    ), root / _CORE_FILES[0])
    _require_columns(fields[_CORE_FILES[3]], ("structure_id", "cif_file", "size_bytes", "sha256"), root / _CORE_FILES[3])
    return snapshots, tables, fields


def projection_context(dataset):
    """Return only factory-authenticated projection state, including target views."""
    entry = _DATASET_GENERATIONS.entry(dataset)
    if entry is None:
        return None
    context = _require_dataset_generation(dataset)
    if context["kind"] == "target_merged":
        return projection_context(context["base_dataset"])
    if context["kind"] == "source_projection":
        return dataset._authority_extra_state
    return None


def _profile(contract, criteria):
    for profile in contract["profiles"]:
        if tuple(profile["criteria"]) == tuple(criteria):
            return profile
    raise ReleaseValidationError(
        "source projection does not contain this ordered grouping profile: {}. "
        "Export it from the complete release; separately projected criteria "
        "must not be recombined.".format(tuple(criteria))
    )


def projected_grouping(dataset, criteria):
    contract = projection_context(dataset)
    if contract is None:
        return None
    from .benchmarks import _GroupingContext
    profile = _profile(contract, criteria)
    return _GroupingContext(tuple(criteria), profile["criterion_groups"],
                            profile["effective_groups"], profile["criterion_diagnostics"],
                            profile["main_union_groups"])


def projected_diversity(dataset, diversity):
    contract = projection_context(dataset)
    if contract is None or diversity == "none":
        return None
    if diversity != "representative" or contract["diversity"] is None:
        raise ReleaseValidationError("source projection has no verified representative diversity index")
    from .benchmarks import DiversityIndex
    cached = contract["diversity"]
    profile = dict(cached["profile"])
    profile["scope"] = "selected-source projection of a saved complete-release diversity index"
    profile["complete_release_index_sha256"] = cached["complete_release_index_sha256"]
    return DiversityIndex(cached["strata"], cached["tiers"], cached["topology_categories"],
                          profile, (), _digest(cached))


def projected_pure_groups(dataset, criteria):
    contract = projection_context(dataset)
    if contract is None:
        return None
    counts = _profile(contract, criteria)["complete_group_label_counts"]
    return {label: {group for group, row in counts.items()
                    if row[label] > 0 and sum(row.values()) == row[label]}
            for label in ("CR", "NCR")}


def projection_receipt(dataset, criteria, receipt):
    """Add explicit selected/full-universe scopes only for projection outputs."""
    contract = projection_context(dataset)
    if contract is None:
        return
    profile = _profile(contract, criteria)
    receipt["source_projection"] = {
        "schema_version": SCHEMA, "sources": list(contract["sources"]),
        "selected_structure_count": contract["dataset_info"]["structure_count"],
        "complete_release_structure_count": contract["complete_release"]["structure_count"],
        "contract_sha256": dataset.input_hashes["source_projection.json"],
        "grouping_and_label_purity_scope": "complete release before source selection",
        "cohort_pool_and_prediction_diagnostic_scope": "selected sources only",
        "omitted_source_IDs_included": False,
        "publication_authorized": False,
    }
    receipt["implementation"]["source_sha256"]["projections.py"] = _implementation_hash()
    groups = receipt.get("effective_leakage_blocks", {})
    if "maximum_complete_release_block_size" in groups:
        groups["maximum_complete_release_block_size"] = profile["complete_maximum_group_size"]
    if "complete_release_block_count" in groups:
        groups["complete_release_block_count"] = profile["complete_group_count"]
    groups["selected_source_group_count"] = len(set(profile["effective_groups"].values()))
    if "complete_release_label_accounting" in receipt:
        receipt["selected_source_label_accounting"] = receipt["complete_release_label_accounting"]
        receipt["complete_release_label_accounting"] = {
            **_plain(contract["complete_release"]["label_accounting"]),
            "companion_manifest": None,
            "scope": "complete-release summary, omitted structure IDs are not exported",
            "non_strict_rows_enter_cohort": False,
        }
        receipt["raw_strict_pool_counts_scope"] = "selected sources before optional eligibility filtering"
        receipt["test_policy_definition"] += (
            " For a source projection, every selected eligible member of a chosen "
            "full-release group is reserved; omitted-source members are not model inputs."
        )


def _implementation_hash():
    value = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if value != _IMPORTED_HASH:
        raise ReleaseValidationError("projection implementation changed after import")
    return value


_IMPORTED_HASH = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def export_source_projection(dataset, source_root, contract_path, *, sources,
                             group_profiles=(("priority_main",), ("RT", "M2T")),
                             diversity="representative"):
    """Write one new contract for existing, unchanged source-subset tables.

    The producer needs the complete authorized release. The receiver does not.
    ``diversity='none'`` explicitly omits cached numeric strata. The output must
    not exist. Source CIFs and metadata are never edited or copied by this call.
    """
    from .benchmarks import _grouping_context, build_diversity_index, normalize_group_criteria
    from .targets import CURRENT_FEATURE_TABLES
    from ._transactions import publish_file_bundle

    if _require_dataset_generation(dataset)["kind"] != "validated_release":
        raise ReleaseValidationError("projection export requires an authenticated complete release")
    if isinstance(sources, (str, bytes)) or not isinstance(sources, (tuple, list)):
        raise TypeError("sources must be an explicit list or tuple")
    if not sources or len(set(sources)) != len(sources) or set(sources).difference({"COD", "CSD", "SI"}):
        raise ValueError("sources must be distinct COD, CSD or SI values")
    if not isinstance(group_profiles, (tuple, list)) or not group_profiles:
        raise ValueError("group_profiles must be a nonempty sequence of ordered criteria")
    profiles = tuple(normalize_group_criteria(value) for value in group_profiles)
    if len(set(profiles)) != len(profiles):
        raise ValueError("duplicate grouping profiles")
    root = Path(source_root).resolve()
    destination = Path(contract_path).absolute()
    checksum_path = Path(str(destination) + ".sha256")
    if destination.exists() or checksum_path.exists() or destination.is_symlink() or checksum_path.is_symlink():
        raise FileExistsError("projection contract output already exists")
    snapshots, tables, fields = _tables(root)
    metadata, order = tables[_CORE_FILES[0]]
    selected = {record.structure_id for record in dataset.records if record.source_database in sources}
    if set(metadata) != selected:
        raise ReleaseValidationError("source subset does not equal the requested complete-release selection")
    for sid in order:
        if (metadata[sid] != dict(dataset[sid].metadata)
                or tables[_CORE_FILES[1]][0][sid] != dict(dataset.parent_by_id[sid])
                or tables[_CORE_FILES[3]][0][sid] != dict(dataset[sid].cif_manifest or {})):
            raise ReleaseValidationError("source record differs from complete-release evidence: " + sid)
    if snapshots[_CORE_FILES[2]].sha256 != dataset.input_hashes[_CORE_FILES[2]]:
        raise ReleaseValidationError("source parent methods differ from complete release")
    for name in CURRENT_FEATURE_TABLES.values():
        if not (root / name).is_file():
            continue
        snapshot = _capture_release_file(root / name, name)
        columns, selected_rows = _read_csv(snapshot, name)
        full_columns, full_rows = _read_csv(_capture_release_file(dataset.release_root / name, name), name)
        selected_index, _ = _index_unique(selected_rows, root / name)
        full_index, _ = _index_unique(full_rows, dataset.release_root / name)
        if columns != full_columns or set(selected_index) != selected or set(full_index) != set(dataset.structure_ids):
            raise ReleaseValidationError("source feature membership/schema differs: " + name)
        if any(selected_index[sid] != full_index[sid] for sid in selected):
            raise ReleaseValidationError("source feature values differ: " + name)
        snapshots[name] = snapshot
    info = {key: _plain(dataset.dataset_info[key]) for key in (
        "dataset_version", "classification_definitions", "release_status", "definitions", "parent_grouping",
        "metadata_layout",
    ) if key in dataset.dataset_info}
    info["structure_count"] = len(selected)
    labels = {record.structure_id: classify_checker_row(record.metadata, CHECKER_PRESETS["5checker"])
              for record in dataset.records}
    counts = Counter(labels.values())
    parent_sizes = {}
    for prefix in _parent_prefixes(fields[_CORE_FILES[1]]):
        parent_sizes[prefix] = {row[prefix + "_group"]: int(row[prefix + "_size"])
                                for row in tables[_CORE_FILES[1]][0].values()}
    payload = {"schema_version": SCHEMA, "sources": sorted(sources), "dataset_info": info,
        "selected_files": {name: {"sha256": snap.sha256, "size_bytes": snap.size_bytes}
                           for name, snap in sorted(snapshots.items())},
        "complete_release": {"dataset_version": dataset.dataset_version,
            "structure_count": len(dataset), "input_sha256": dict(dataset.input_hashes),
            "structure_ids_sha256": _digest(sorted(dataset.structure_ids)),
            "label_accounting": {"counts": {label: counts[label] for label in _LABELS},
                "membership_sha256": {label: _digest(sorted(sid for sid, value in labels.items() if value == label))
                                       for label in _LABELS}}},
        "parent_group_sizes": parent_sizes, "profiles": [], "diversity": None,
        "official_split": False, "publication_authorized": False,
        "producer_sha256": _implementation_hash()}
    for criteria in profiles:
        grouping = _grouping_context(dataset, criteria)
        label_counts = defaultdict(Counter)
        for sid, group in grouping.effective_blocks.items():
            label_counts[group][labels[sid]] += 1
        represented = {grouping.effective_blocks[sid] for sid in selected}
        payload["profiles"].append({"criteria": list(criteria),
            "effective_groups": {sid: grouping.effective_blocks[sid] for sid in sorted(selected)},
            "main_union_groups": {sid: grouping.main_union_groups[sid] for sid in sorted(selected)},
            "criterion_groups": {key: {sid: value[sid] for sid in sorted(selected)}
                                 for key, value in grouping.criterion_groups.items()},
            "criterion_diagnostics": {key: {sid: value[sid] for sid in sorted(selected) if sid in value}
                                      for key, value in grouping.criterion_diagnostics.items()},
            "complete_group_label_counts": {group: {label: label_counts[group][label] for label in _LABELS}
                                            for group in sorted(represented)},
            "complete_group_count": len(label_counts),
            "complete_maximum_group_size": max(sum(row.values()) for row in label_counts.values())})
    if diversity not in ("none", "representative"):
        raise ValueError("diversity must be none or representative")
    if diversity == "representative":
        index = build_diversity_index(dataset, diversity=diversity)
        payload["diversity"] = {
            "profile": _plain(index.profile), "complete_release_index_sha256": index.digest,
            "complete_release_input_receipts": _plain(index.input_receipts),
            "strata": {sid: index.strata_by_id[sid] for sid in sorted(selected)},
            "tiers": {sid: index.tier_by_id[sid] for sid in sorted(selected)},
            "topology_categories": {sid: index.topology_category_by_id[sid] for sid in sorted(selected)}}
    # Validate the emitted meaning before making the contract visible.
    _validate_contract(payload, tables, fields)
    _require_dataset_generation(dataset)
    for name, snapshot in snapshots.items():
        if _capture_release_file(root / name, name).sha256 != snapshot.sha256:
            raise ReleaseValidationError("source changed during projection export")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".coremof-projection-", dir=str(destination.parent)))
    preserve_staging = False
    try:
        staged = staging / destination.name
        staged.write_text(json.dumps(payload, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
        digest = hashlib.sha256(staged.read_bytes()).hexdigest()
        staged_checksum = staging / checksum_path.name
        staged_checksum.write_text(digest + "  " + destination.name + "\n", encoding="utf-8")
        publish_file_bundle((staged, staged_checksum), (destination, checksum_path), overwrite=False)
    except BaseException as error:
        preserve_staging = hasattr(error, "coremof_preserved_staging_directory")
        raise
    finally:
        if not preserve_staging:
            shutil.rmtree(staging, ignore_errors=True)
    return MappingProxyType({"path": str(destination), "sha256": digest,
                             "structure_count": len(selected), "publication_authorized": False})


def _validate_contract(contract, tables, fields):
    from .benchmarks import normalize_group_criteria
    from .targets import CURRENT_FEATURE_TABLES
    _keys(contract, ("schema_version", "sources", "dataset_info", "selected_files", "complete_release",
                     "parent_group_sizes", "profiles", "diversity", "official_split", "publication_authorized",
                     "producer_sha256"), "source projection")
    if contract["schema_version"] != SCHEMA or contract["official_split"] is not False or contract["publication_authorized"] is not False:
        raise ReleaseValidationError("source projection is exploratory and not publication authorization")
    _hash(contract["producer_sha256"])
    sources = contract["sources"]
    if not isinstance(sources, list) or not sources or any(type(v) is not str for v in sources) or sources != sorted(set(sources)) or set(sources).difference({"COD", "CSD", "SI"}):
        raise ReleaseValidationError("invalid source selection")
    metadata, order = tables[_CORE_FILES[0]]
    ids = set(order)
    if {row["source_database"] for row in metadata.values()} != set(sources):
        raise ReleaseValidationError("source metadata differs from contract selection")
    info = contract["dataset_info"]
    _reject_retired_or_reserved_keys(info, "dataset_info")
    _validate_dataset_info(info, len(order))
    _validate_metadata_identity(tuple(metadata.values()))
    _validate_published_labels(tuple(metadata.values()), info)
    files = contract["selected_files"]
    if not isinstance(files, dict) or not set(_CORE_FILES).issubset(files) or set(files).difference(set(_CORE_FILES) | set(CURRENT_FEATURE_TABLES.values())):
        raise ReleaseValidationError("unsupported source projection files")
    for binding in files.values():
        _keys(binding, ("sha256", "size_bytes"), "selected file binding")
        _hash(binding["sha256"])
        _count(binding["size_bytes"])
    full = contract["complete_release"]
    _keys(full, ("dataset_version", "structure_count", "input_sha256", "structure_ids_sha256", "label_accounting"), "complete release")
    if full["dataset_version"] != info["dataset_version"] or _count(full["structure_count"]) < len(ids):
        raise ReleaseValidationError("complete-release identity/count differs")
    _hash(full["structure_ids_sha256"])
    if not isinstance(full["input_sha256"], dict) or not set(_CORE_FILES).issubset(full["input_sha256"]):
        raise ReleaseValidationError("complete-release input bindings are missing")
    for digest in full["input_sha256"].values():
        _hash(digest)
    if full["input_sha256"][_CORE_FILES[2]] != files[_CORE_FILES[2]]["sha256"]:
        raise ReleaseValidationError("source and complete-release parent methods differ")
    _keys(full["label_accounting"], ("counts", "membership_sha256"), "complete label accounting")
    _keys(full["label_accounting"]["counts"], _LABELS, "complete labels")
    _keys(full["label_accounting"]["membership_sha256"], _LABELS, "complete label memberships")
    if sum(_count(n) for n in full["label_accounting"]["counts"].values()) != full["structure_count"]:
        raise ReleaseValidationError("complete label counts do not reconcile")
    for value in full["label_accounting"]["membership_sha256"].values():
        _hash(value)
    sizes = contract["parent_group_sizes"]
    prefixes = _parent_prefixes(fields[_CORE_FILES[1]])
    _keys(sizes, prefixes, "complete parent sizes")
    for prefix in prefixes:
        expected_groups = {row[prefix + "_group"] for row in tables[_CORE_FILES[1]][0].values()}
        _keys(sizes[prefix], expected_groups, "complete parent size groups")
        for value in sizes[prefix].values():
            if _count(value) < 1 or value > full["structure_count"]:
                raise ReleaseValidationError("invalid complete parent size")
    _validate_parent_rows(fields[_CORE_FILES[1]], tuple(tables[_CORE_FILES[1]][0].values()),
                          metadata_by_id=metadata, _complete_group_sizes=sizes)
    profiles = contract["profiles"]
    if not isinstance(profiles, list) or not profiles:
        raise ReleaseValidationError("source projection needs grouping profiles")
    seen = set()
    main_union = None
    for profile in profiles:
        _keys(profile, ("criteria", "effective_groups", "main_union_groups", "criterion_groups", "criterion_diagnostics",
                        "complete_group_label_counts", "complete_group_count", "complete_maximum_group_size"), "grouping profile")
        criteria = normalize_group_criteria(profile["criteria"])
        if list(criteria) != profile["criteria"] or criteria in seen:
            raise ReleaseValidationError("noncanonical or duplicate projection criteria")
        seen.add(criteria)
        _exact_map(profile["effective_groups"], ids, "effective groups")
        _exact_map(profile["main_union_groups"], ids, "main-union groups")
        if main_union is not None and main_union != profile["main_union_groups"]:
            raise ReleaseValidationError("main-union groups differ between projection profiles")
        main_union = profile["main_union_groups"]
        _keys(profile["criterion_groups"], criteria, "criterion groups")
        _keys(profile["criterion_diagnostics"], criteria, "criterion diagnostics")
        for criterion in criteria:
            _exact_map(profile["criterion_groups"][criterion], ids, "criterion groups")
            diagnostic = profile["criterion_diagnostics"][criterion]
            if not isinstance(diagnostic, dict) or set(diagnostic).difference(ids) or any(not isinstance(v, str) or not v for v in diagnostic.values()):
                raise ReleaseValidationError("invalid criterion diagnostics")
        for mapping in [profile["main_union_groups"], *profile["criterion_groups"].values()]:
            destinations = {}
            for sid, group in mapping.items():
                if destinations.setdefault(group, profile["effective_groups"][sid]) != profile["effective_groups"][sid]:
                    raise ReleaseValidationError("effective projection splits a selected criterion group")
        groups = set(profile["effective_groups"].values())
        _keys(profile["complete_group_label_counts"], groups, "complete group labels")
        observed = defaultdict(Counter)
        for sid, group in profile["effective_groups"].items():
            observed[group][metadata[sid]["label_5checker"]] += 1
        total = Counter()
        for group, counts in profile["complete_group_label_counts"].items():
            _keys(counts, _LABELS, "group label counts")
            for label, count in counts.items():
                if _count(count) < observed[group][label]:
                    raise ReleaseValidationError("complete group label counts omit selected members")
                total[label] += count
            if sum(counts.values()) > _count(profile["complete_maximum_group_size"]):
                raise ReleaseValidationError("group size exceeds complete-release maximum")
        if any(total[label] > full["label_accounting"]["counts"][label] for label in _LABELS):
            raise ReleaseValidationError("represented groups exceed complete label totals")
        if not len(groups) <= _count(profile["complete_group_count"]) <= full["structure_count"] or not 1 <= profile["complete_maximum_group_size"] <= full["structure_count"]:
            raise ReleaseValidationError("invalid complete group summary")
    cached = contract["diversity"]
    if cached is not None:
        _keys(cached, ("profile", "complete_release_index_sha256", "complete_release_input_receipts",
                       "strata", "tiers", "topology_categories"), "diversity projection")
        _hash(cached["complete_release_index_sha256"])
        for field in ("strata", "tiers", "topology_categories"):
            _exact_map(cached[field], ids, field)
        if set(cached["tiers"].values()).difference({"rac5", "zeo", "no_numeric"}):
            raise ReleaseValidationError("invalid diversity feature-availability tier")
        profile = cached["profile"]
        if not isinstance(profile, dict) or profile.get("name") != "representative" or profile.get("target_columns_consumed") != [] or profile.get("scientific_feature_imputation") is not False:
            raise ReleaseValidationError("invalid target-free diversity profile")
        if not isinstance(cached["complete_release_input_receipts"], list) or len(cached["complete_release_input_receipts"]) != 3:
            raise ReleaseValidationError("representative index requires three input bindings")
        receipts = cached["complete_release_input_receipts"]
        if {item.get("name") for item in receipts if isinstance(item, dict)} != {"rac5", "zeo", "topology"}:
            raise ReleaseValidationError("invalid representative feature bindings")
        for item in receipts:
            if item.get("release_path") != CURRENT_FEATURE_TABLES[item["name"]]:
                raise ReleaseValidationError("invalid representative input path")
            _hash(item["sha256"])
            _count(item["size_bytes"])
            if _count(item["row_count"]) != full["structure_count"]:
                raise ReleaseValidationError("representative inputs do not cover the complete release")


def load_source_projection(source_root, contract_path, *, expected_sha256,
                           verify_cif_files=False, dataset_class=CoREMOFDataset,
                           include_legacy_targets=False):
    """Validate selected rows and a trusted contract without reading omitted rows.

    Original file and contract hashes remain bound even when a historical
    combined table is exposed as a target-free view. The explicit compatibility
    option is not allowed to bypass a declared target-free layout.
    """
    _hash(expected_sha256)
    if type(verify_cif_files) is not bool:
        raise TypeError("verify_cif_files must be a boolean")
    if type(include_legacy_targets) is not bool:
        raise TypeError("include_legacy_targets must be a boolean")
    root = Path(source_root).expanduser().resolve()
    path = Path(contract_path).expanduser().resolve()
    if path.stat().st_size > 128 * 1024 * 1024:
        raise ReleaseValidationError("source projection contract exceeds the 128-MiB limit")
    snapshot = _capture_release_file(path, "source projection contract")
    if snapshot.sha256 != expected_sha256:
        raise ReleaseValidationError("source projection contract checksum mismatch")

    def pairs(values):
        result = {}
        for key, value in values:
            if key in result:
                raise ReleaseValidationError("duplicate source projection JSON key")
            result[key] = value
        return result

    def nonfinite(value):
        raise ReleaseValidationError("non-finite source projection JSON value: " + value)

    contract = json.loads(snapshot.data.decode("utf-8"), object_pairs_hook=pairs, parse_constant=nonfinite)
    snapshots, tables, fields = _tables(root)
    _validate_contract(contract, tables, fields)
    for name, binding in contract["selected_files"].items():
        current = snapshots.get(name)
        if current is None:
            current = _capture_release_file(root / name, name)
            snapshots[name] = current
        if current.sha256 != binding["sha256"] or current.size_bytes != binding["size_bytes"]:
            raise ReleaseValidationError("source projection input checksum mismatch: " + name)
    metadata, order = tables[_CORE_FILES[0]]
    parents = tables[_CORE_FILES[1]][0]
    manifests = tables[_CORE_FILES[3]][0]
    methods = _read_json_object(snapshots[_CORE_FILES[2]], "parent methods")
    _reject_retired_or_reserved_keys(metadata, "metadata_rows")
    _reject_retired_or_reserved_keys(methods, "parent_group_methods")
    _validate_parent_methods(methods, contract["dataset_info"], fields[_CORE_FILES[1]])
    _validate_cif_manifest(metadata, manifests, root=root, verify_files=verify_cif_files)
    _validate_target_free_layout(contract["dataset_info"], fields[_CORE_FILES[0]])
    if not include_legacy_targets:
        metadata = {
            sid: {name: value for name, value in row.items() if name not in _LEGACY_TARGET_COLUMNS}
            for sid, row in metadata.items()
        }
    info = (
        contract["dataset_info"] if include_legacy_targets
        else _without_legacy_target_information(contract["dataset_info"])
    )
    prefixes = _parent_prefixes(fields[_CORE_FILES[1]])
    immutable_parents = {sid: MappingProxyType(parents[sid]) for sid in order}
    records = [StructureRecord(sid, MappingProxyType(metadata[sid]),
                 _ParentGroupView(immutable_parents[sid], prefixes), MappingProxyType(manifests[sid]))
               for sid in order]
    hashes = {name: value.sha256 for name, value in sorted(snapshots.items())}
    hashes["source_projection.json"] = expected_sha256
    result = dataset_class(root, records, _deep_freeze(info),
        _deep_freeze(methods), MappingProxyType(immutable_parents), MappingProxyType(hashes), verify_cif_files)
    result._authority_extra_state = _deep_freeze(contract)
    _register_dataset_generation(result, kind="source_projection", official_release_source=True)
    return result
