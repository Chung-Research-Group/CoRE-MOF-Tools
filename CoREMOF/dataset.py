"""Lightweight access to a validated CoRE-MOF release.

The loader in this module uses only the Python standard library.  Importing it
does not import pandas, NumPy, crystallographic software, or network clients.
It treats release metadata as an audited contract: identifiers must agree
across tables, parent-group declarations must be internally consistent, and
all published checker labels are independently recomputed before use.
"""

import csv
import hashlib
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import (
    Dict,
    Iterable,
    Iterator,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)

from .labels import (
    CHECKER_COLUMNS,
    CHECKER_PRESETS,
    LABELS,
    classify_checker_row,
    resolve_checker_view,
)


PARENT_STATUSES = frozenset({"MATCHED", "UNMATCHED", "NOT_AVAILABLE"})
PUBLIC_CHECKER_STATUSES = frozenset({"PASS", "FAIL", "NOT_AVAILABLE"})
PUBLIC_PARENT_METHOD_PREFIXES = MappingProxyType(
    {
        "rac5_zeo": "rac_zeo",
        "rac5": "rac",
        "zeo": "zeo",
        "source_id": "source",
        "mofid_v2": "mofid2",
        "mofid_v1": "mofid1",
        "common_name": "name",
        "identity_union": "identity",
    }
)
PARENT_METHOD_ALIASES = MappingProxyType(
    {
        **dict(PUBLIC_PARENT_METHOD_PREFIXES),
        "rac_zeo": "rac_zeo",
        "rac": "rac",
        "source": "source",
        "mofid2": "mofid2",
        "mofid1": "mofid1",
        "name": "name",
        "identity": "identity",
    }
)

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_STRUCTURE_ID_RE = re.compile(
    r"^(ASR|FSR|ION)-(COD|CSD|SI)-(\d{4}|UNKN)-(\d{4})$"
)
_PARENT_GROUP_PREFIXES = MappingProxyType(
    {
        "rac_zeo": "RZ-",
        "rac": "R-",
        "zeo": "Z-",
        "source": "S-",
        "mofid2": "M2-",
        "mofid1": "M1-",
        "name": "N-",
        "identity": "I-",
    }
)


class ReleaseValidationError(ValueError):
    """Raised when release tables do not satisfy the public data contract."""


def _deep_freeze(value: object) -> object:
    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _deep_freeze(item) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple, set, frozenset)):
        return tuple(_deep_freeze(item) for item in value)
    return value


def normalize_parent_method(method: str) -> str:
    """Translate a documented criterion name to its CSV column prefix."""

    if not isinstance(method, str):
        raise TypeError("parent method must be a string")
    try:
        return PARENT_METHOD_ALIASES[method]
    except KeyError:
        raise ValueError(
            "unknown parent method {!r}; choose one of {}".format(
                method, ", ".join(sorted(PARENT_METHOD_ALIASES))
            )
        )


@dataclass(frozen=True)
class ParentGroup:
    """One structure's membership under one parent criterion."""

    method: str
    status: str
    group_id: str
    size: int

    @property
    def available(self) -> bool:
        return self.status != "NOT_AVAILABLE"

    @property
    def matched(self) -> bool:
        return self.status == "MATCHED"


@dataclass(frozen=True)
class StructureRecord:
    """One row of release metadata joined to its parent and CIF records."""

    structure_id: str
    metadata: Mapping[str, str]
    parent_groups: Mapping[str, ParentGroup]
    cif_manifest: Optional[Mapping[str, str]] = None

    def __getitem__(self, key: str) -> str:
        return self.metadata[key]

    def get(self, key: str, default: Optional[str] = None) -> Optional[str]:
        return self.metadata.get(key, default)

    @property
    def source_database(self) -> str:
        return self.metadata.get("source_database", "")

    @property
    def structure_variant(self) -> str:
        return self.metadata.get("structure_variant", "")

    @property
    def metal_elements(self) -> Tuple[str, ...]:
        value = self.metadata.get("metal_elements", "")
        return tuple(element for element in value.split(";") if element)

    def parent_group(self, method: str) -> ParentGroup:
        prefix = normalize_parent_method(method)
        try:
            return self.parent_groups[prefix]
        except KeyError:
            raise KeyError(
                "parent method {!r} is not present in this release".format(method)
            )


class _ParentGroupView(Mapping[str, ParentGroup]):
    """Lazy immutable view over one raw parent-group CSV row."""

    __slots__ = ("_row", "_prefixes")

    def __init__(self, row: Mapping[str, str], prefixes: Sequence[str]) -> None:
        self._row = row
        self._prefixes = tuple(prefixes)

    def __getitem__(self, prefix: str) -> ParentGroup:
        if prefix not in self._prefixes:
            raise KeyError(prefix)
        return ParentGroup(
            method=prefix,
            status=self._row["{}_status".format(prefix)],
            group_id=self._row["{}_group".format(prefix)],
            size=int(self._row["{}_size".format(prefix)]),
        )

    def __iter__(self) -> Iterator[str]:
        return iter(self._prefixes)

    def __len__(self) -> int:
        return len(self._prefixes)


@dataclass(frozen=True)
class ClassifiedRecord:
    """A release record with a recomputed checker-consensus label."""

    record: StructureRecord
    checker_view: str
    label: str
    checker_statuses: Mapping[str, str]

    @property
    def checker_preset(self) -> str:
        """Backward-compatible alias for :attr:`checker_view`."""

        return self.checker_view

    @property
    def structure_id(self) -> str:
        return self.record.structure_id

    @property
    def metadata(self) -> Mapping[str, str]:
        return self.record.metadata

    @property
    def parent_groups(self) -> Mapping[str, ParentGroup]:
        return self.record.parent_groups

    @property
    def source_database(self) -> str:
        return self.record.source_database

    @property
    def structure_variant(self) -> str:
        return self.record.structure_variant

    @property
    def metal_elements(self) -> Tuple[str, ...]:
        return self.record.metal_elements

    def parent_group(self, method: str) -> ParentGroup:
        return self.record.parent_group(method)


class CoREMOFDataset:
    """A strict, read-only view of one on-disk CoRE-MOF release."""

    def __init__(
        self,
        release_root: Path,
        records: Sequence[StructureRecord],
        dataset_info: Mapping[str, object],
        parent_group_methods: Mapping[str, object],
        parent_by_id: Mapping[str, Mapping[str, str]],
        input_hashes: Mapping[str, str],
        cif_files_verified: bool,
    ) -> None:
        self.release_root = release_root
        self.records = tuple(records)
        self.metadata_rows = tuple(record.metadata for record in self.records)
        self.dataset_info = dataset_info
        self.parent_group_methods = parent_group_methods
        self.parent_by_id = parent_by_id
        self.cif_hashes = MappingProxyType(
            {
                record.structure_id: record.cif_manifest["sha256"]
                for record in self.records
                if record.cif_manifest is not None
            }
        )
        self.input_hashes = input_hashes
        self.cif_files_verified = cif_files_verified
        self._by_id = {record.structure_id: record for record in self.records}

    @classmethod
    def from_release(
        cls,
        release_root: Union[str, Path],
        verify_cif_files: bool = False,
    ) -> "CoREMOFDataset":
        """Load and validate a release directory.

        Required files are ``metadata/metadata.csv``,
        ``parent_groups/parent_groups.csv``,
        ``parent_groups/parent_group_methods.json``, and ``dataset_info.json``.
        If ``manifests/cif_manifest.csv`` is present, it is joined and validated
        against the same exact identifier set.
        """

        root = Path(release_root).expanduser()
        if not root.is_dir():
            raise FileNotFoundError("release directory does not exist: {}".format(root))
        root = root.resolve()

        metadata_path = root / "metadata" / "metadata.csv"
        parents_path = root / "parent_groups" / "parent_groups.csv"
        methods_path = root / "parent_groups" / "parent_group_methods.json"
        info_path = root / "dataset_info.json"
        manifest_path = root / "manifests" / "cif_manifest.csv"

        metadata_fields, metadata_rows = _read_csv(
            metadata_path, "release metadata"
        )
        _require_columns(
            metadata_fields,
            (
                "structure_id",
                "cif_file",
                "source_database",
                "source_id",
                "structure_variant",
                "metal_elements",
                *tuple(CHECKER_COLUMNS.values()),
                "label_3checker",
                "label_4checker",
                "label_5checker",
            ),
            metadata_path,
        )
        metadata_by_id, metadata_order = _index_unique(
            metadata_rows, metadata_path
        )

        parent_fields, parent_rows = _read_csv(parents_path, "parent groups")
        required_parent_columns = ["structure_id"]
        for prefix in sorted(set(PUBLIC_PARENT_METHOD_PREFIXES.values())):
            required_parent_columns.extend(
                (
                    "{}_status".format(prefix),
                    "{}_group".format(prefix),
                    "{}_size".format(prefix),
                )
            )
        _require_columns(parent_fields, tuple(required_parent_columns), parents_path)
        parent_by_id, _ = _index_unique(parent_rows, parents_path)

        dataset_info = _read_json_object(info_path, "dataset information")
        parent_methods = _read_json_object(methods_path, "parent-group methods")

        metadata_ids = set(metadata_by_id)
        _require_exact_id_set(
            metadata_ids, set(parent_by_id), "metadata.csv", "parent_groups.csv"
        )

        manifest_by_id = None
        if manifest_path.is_file():
            manifest_fields, manifest_rows = _read_csv(
                manifest_path, "CIF manifest"
            )
            _require_columns(
                manifest_fields,
                ("structure_id", "cif_file", "size_bytes", "sha256"),
                manifest_path,
            )
            manifest_by_id, _ = _index_unique(manifest_rows, manifest_path)
            _require_exact_id_set(
                metadata_ids,
                set(manifest_by_id),
                "metadata.csv",
                "cif_manifest.csv",
            )

        _validate_dataset_info(dataset_info, len(metadata_rows))
        _validate_metadata_identity(metadata_rows)
        _validate_parent_methods(parent_methods, dataset_info, parent_fields)
        parent_prefixes = _validate_parent_rows(parent_fields, parent_rows)
        _validate_published_labels(metadata_rows, dataset_info)
        if manifest_by_id is not None:
            _validate_cif_manifest(
                metadata_by_id,
                manifest_by_id,
                root=root,
                verify_files=verify_cif_files,
            )
        elif verify_cif_files:
            raise ReleaseValidationError(
                "verify_cif_files=True requires manifests/cif_manifest.csv"
            )

        records = []
        for structure_id in metadata_order:
            metadata = MappingProxyType(metadata_by_id[structure_id])
            raw_parent = MappingProxyType(parent_by_id[structure_id])
            manifest = None
            if manifest_by_id is not None:
                manifest = MappingProxyType(manifest_by_id[structure_id])
            records.append(
                StructureRecord(
                    structure_id=structure_id,
                    metadata=metadata,
                    parent_groups=_ParentGroupView(raw_parent, parent_prefixes),
                    cif_manifest=manifest,
                )
            )

        input_paths = {
            "metadata/metadata.csv": metadata_path,
            "parent_groups/parent_groups.csv": parents_path,
            "parent_groups/parent_group_methods.json": methods_path,
            "dataset_info.json": info_path,
        }
        if manifest_path.is_file():
            input_paths["manifests/cif_manifest.csv"] = manifest_path
        input_hashes = MappingProxyType(
            {
                name: _sha256_file(path)
                for name, path in sorted(input_paths.items())
            }
        )
        immutable_parent_by_id = MappingProxyType(
            {
                structure_id: MappingProxyType(parent_by_id[structure_id])
                for structure_id in metadata_order
            }
        )

        return cls(
            release_root=root,
            records=records,
            dataset_info=_deep_freeze(dataset_info),
            parent_group_methods=_deep_freeze(parent_methods),
            parent_by_id=immutable_parent_by_id,
            input_hashes=input_hashes,
            cif_files_verified=verify_cif_files,
        )

    @property
    def dataset_version(self) -> str:
        value = self.dataset_info.get("dataset_version", "")
        return str(value)

    @property
    def structure_ids(self) -> Tuple[str, ...]:
        return tuple(record.structure_id for record in self.records)

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self) -> Iterator[StructureRecord]:
        return iter(self.records)

    def __getitem__(
        self, key: Union[int, slice, str]
    ) -> Union[StructureRecord, Tuple[StructureRecord, ...]]:
        if isinstance(key, str):
            return self._by_id[key]
        return self.records[key]

    def classify(
        self,
        checkers: Union[int, str, Sequence[str]] = "5checker",
        labels: Optional[Iterable[str]] = None,
        sources: Optional[Iterable[str]] = None,
        variants: Optional[Iterable[str]] = None,
        metals: Optional[Iterable[str]] = None,
        structure_ids: Optional[Iterable[str]] = None,
    ) -> "ClassifiedDataset":
        """Recompute one checker view and optionally filter its records."""

        view_name, checker_names, view_is_official = resolve_checker_view(checkers)
        records = []
        for record in self.records:
            label = classify_checker_row(record.metadata, checker_names)
            statuses = {
                checker: record.metadata[column]
                for checker, column in CHECKER_COLUMNS.items()
            }
            records.append(
                ClassifiedRecord(
                    record=record,
                    checker_view=view_name,
                    label=label,
                    checker_statuses=MappingProxyType(statuses),
                )
            )
        classified = ClassifiedDataset(
            self, view_name, records, checker_view_official=view_is_official
        )
        return classified.filter(
            labels=labels,
            sources=sources,
            variants=variants,
            metals=metals,
            structure_ids=structure_ids,
        )


class ClassifiedDataset:
    """A checker-labelled, filterable dataset ready for parent-safe splitting."""

    def __init__(
        self,
        dataset: CoREMOFDataset,
        checker_view: str,
        records: Sequence[ClassifiedRecord],
        selection_steps: Sequence[Mapping[str, object]] = (),
        checker_view_official: bool = True,
    ) -> None:
        self.dataset = dataset
        self.checker_view = checker_view
        self.checker_preset = checker_view
        self.checker_view_official = bool(checker_view_official)
        self.records = tuple(records)
        self._selection_steps = tuple(dict(step) for step in selection_steps)
        self._by_id = {record.structure_id: record for record in self.records}
        self.label_by_id = MappingProxyType(
            {record.structure_id: record.label for record in self.records}
        )
        selected_digest = hashlib.sha256(
            "\n".join(sorted(self._by_id)).encode("utf-8")
        ).hexdigest()
        self.selection_filters = _deep_freeze(
            {
                "applied": bool(self._selection_steps),
                "steps": tuple(dict(step) for step in self._selection_steps),
                "selected_count": len(self.records),
                "selected_ids_sha256": selected_digest,
            }
        )

    @property
    def dataset_version(self) -> str:
        return self.dataset.dataset_version

    @property
    def release_root(self) -> Path:
        return self.dataset.release_root

    @property
    def structure_ids(self) -> Tuple[str, ...]:
        return tuple(record.structure_id for record in self.records)

    @property
    def labels(self) -> Tuple[str, ...]:
        return tuple(record.label for record in self.records)

    def ids_for_label(self, label: str) -> Tuple[str, ...]:
        """Return IDs with one canonical consensus label in view order."""

        if not isinstance(label, str) or label.upper() not in LABELS:
            raise ValueError(
                "label must be one of {}".format(", ".join(sorted(LABELS)))
            )
        canonical = label.upper()
        return tuple(
            record.structure_id
            for record in self.records
            if record.label == canonical
        )

    @property
    def cr_ids(self) -> Tuple[str, ...]:
        return self.ids_for_label("CR")

    @property
    def ncr_ids(self) -> Tuple[str, ...]:
        return self.ids_for_label("NCR")

    @property
    def ambiguous_ids(self) -> Tuple[str, ...]:
        return self.ids_for_label("AMBIGUOUS")

    @property
    def unchecked_ids(self) -> Tuple[str, ...]:
        return self.ids_for_label("UNCHECKED")

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self) -> Iterator[ClassifiedRecord]:
        return iter(self.records)

    def __getitem__(
        self, key: Union[int, slice, str]
    ) -> Union[ClassifiedRecord, Tuple[ClassifiedRecord, ...]]:
        if isinstance(key, str):
            return self._by_id[key]
        return self.records[key]

    def filter(
        self,
        labels: Optional[Iterable[str]] = None,
        sources: Optional[Iterable[str]] = None,
        variants: Optional[Iterable[str]] = None,
        metals: Optional[Iterable[str]] = None,
        structure_ids: Optional[Iterable[str]] = None,
    ) -> "ClassifiedDataset":
        """Return a new view matching all supplied filters.

        A structure passes ``metals`` when it contains at least one requested
        metal.  Other filters use exact values from public metadata.
        """

        label_set = _filter_set(labels, "labels")
        if label_set is not None:
            label_set = {value.upper() for value in label_set}
            invalid = label_set.difference(LABELS)
            if invalid:
                raise ValueError(
                    "unknown label(s): {}".format(", ".join(sorted(invalid)))
                )
        source_set = _filter_set(sources, "sources")
        variant_set = _filter_set(variants, "variants")
        metal_set = _filter_set(metals, "metals")
        structure_id_set = _filter_set(structure_ids, "structure_ids")
        if structure_id_set is not None:
            unknown = structure_id_set.difference(self.dataset.structure_ids)
            if unknown:
                raise KeyError(
                    "unknown structure_id filter value(s): {}".format(
                        ", ".join(sorted(unknown)[:10])
                    )
                )

        if all(
            value is None
            for value in (
                label_set,
                source_set,
                variant_set,
                metal_set,
                structure_id_set,
            )
        ):
            return self
        folded_sources = (
            {value.casefold() for value in source_set}
            if source_set is not None
            else None
        )
        folded_variants = (
            {value.casefold() for value in variant_set}
            if variant_set is not None
            else None
        )
        folded_metals = (
            {value.casefold() for value in metal_set}
            if metal_set is not None
            else None
        )

        selected = []
        for record in self.records:
            if structure_id_set is not None and record.structure_id not in structure_id_set:
                continue
            if label_set is not None and record.label not in label_set:
                continue
            if (
                folded_sources is not None
                and record.source_database.casefold()
                not in folded_sources
            ):
                continue
            if (
                folded_variants is not None
                and record.structure_variant.casefold()
                not in folded_variants
            ):
                continue
            if folded_metals is not None and not folded_metals.intersection(
                element.casefold() for element in record.metal_elements
            ):
                continue
            selected.append(record)
        selection_step = {
            "labels": tuple(sorted(label_set)) if label_set is not None else None,
            "sources": tuple(sorted(source_set)) if source_set is not None else None,
            "variants": tuple(sorted(variant_set)) if variant_set is not None else None,
            "metals": tuple(sorted(metal_set)) if metal_set is not None else None,
            "structure_ids": (
                tuple(sorted(structure_id_set)) if structure_id_set is not None else None
            ),
        }
        return ClassifiedDataset(
            self.dataset,
            self.checker_view,
            selected,
            selection_steps=self._selection_steps + (selection_step,),
            checker_view_official=self.checker_view_official,
        )

    select = filter

    def label_counts(self) -> Mapping[str, int]:
        return MappingProxyType(dict(Counter(self.labels)))

    def train_valid_test_split(
        self,
        parent_method: str = "priority_main",
        fractions: Sequence[float] = (0.8, 0.1, 0.1),
        random_state: int = 42,
        labels: Optional[Iterable[str]] = ("CR", "NCR"),
        sources: Optional[Iterable[str]] = None,
        variants: Optional[Iterable[str]] = None,
        metals: Optional[Iterable[str]] = None,
        structure_ids: Optional[Iterable[str]] = None,
        missing_parent: str = "singleton",
        leakage_guard: str = "auto",
        stratify_by: Sequence[str] = ("label",),
        **splitter_options: object
    ) -> object:
        """Create a parent-safe train/validation/test split.

        The splitter is imported only when requested, keeping metadata loading
        independent from optional training dependencies and future splitter
        implementations.
        """

        from .splitters import ParentGroupSplitter

        splitter = ParentGroupSplitter(
            self,
            parent_method=parent_method,
            random_state=random_state,
            missing_parent=missing_parent,
            leakage_guard=leakage_guard,
            stratify_by=stratify_by,
            labels=labels,
            sources=sources,
            variants=variants,
            metals=metals,
            structure_ids=structure_ids,
            **splitter_options
        )
        return splitter.train_valid_test_split(fractions=fractions)

    def split(self, **kwargs: object) -> object:
        """Alias for :meth:`train_valid_test_split`."""

        return self.train_valid_test_split(**kwargs)


def _read_csv(path: Path, description: str) -> Tuple[Tuple[str, ...], List[Dict[str, str]]]:
    if not path.is_file():
        raise FileNotFoundError("{} file does not exist: {}".format(description, path))
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ReleaseValidationError("CSV has no header: {}".format(path))
        fields = tuple(reader.fieldnames)
        if len(set(fields)) != len(fields):
            raise ReleaseValidationError(
                "CSV has duplicate header names: {}".format(path)
            )
        rows = []
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ReleaseValidationError(
                    "{}:{} contains more values than its header".format(
                        path, row_number
                    )
                )
            missing = [field for field, value in row.items() if value is None]
            if missing:
                raise ReleaseValidationError(
                    "{}:{} is missing value(s) for {}".format(
                        path, row_number, ", ".join(missing)
                    )
                )
            rows.append(dict(row))
    return fields, rows


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json_object(path: Path, description: str) -> Dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError("{} file does not exist: {}".format(description, path))
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseValidationError("cannot read {} {}: {}".format(description, path, error))
    if not isinstance(value, dict):
        raise ReleaseValidationError("{} must contain a JSON object: {}".format(description, path))
    return value


def _require_columns(
    fields: Sequence[str], required: Sequence[str], path: Path
) -> None:
    missing = [column for column in required if column not in fields]
    if missing:
        raise ReleaseValidationError(
            "{} is missing required column(s): {}".format(path, ", ".join(missing))
        )


def _index_unique(
    rows: Sequence[Mapping[str, str]], path: Path
) -> Tuple[Dict[str, Mapping[str, str]], Tuple[str, ...]]:
    indexed = {}
    order = []
    for row_number, row in enumerate(rows, start=2):
        structure_id = row.get("structure_id", "")
        if not structure_id:
            raise ReleaseValidationError(
                "{}:{} has an empty structure_id".format(path, row_number)
            )
        if structure_id in indexed:
            raise ReleaseValidationError(
                "{} contains duplicate structure_id {!r}".format(path, structure_id)
            )
        indexed[structure_id] = row
        order.append(structure_id)
    return indexed, tuple(order)


def _require_exact_id_set(
    expected: set, observed: set, expected_name: str, observed_name: str
) -> None:
    if expected == observed:
        return
    missing = sorted(expected.difference(observed))[:5]
    extra = sorted(observed.difference(expected))[:5]
    raise ReleaseValidationError(
        "structure_id sets differ between {} and {}: missing={} extra={}".format(
            expected_name, observed_name, missing, extra
        )
    )


def _validate_dataset_info(info: Mapping[str, object], row_count: int) -> None:
    version = info.get("dataset_version")
    if not isinstance(version, str) or not version:
        raise ReleaseValidationError("dataset_info.json has no dataset_version")
    declared_count = info.get("structure_count")
    if not isinstance(declared_count, int) or isinstance(declared_count, bool):
        raise ReleaseValidationError(
            "dataset_info.json structure_count must be an integer"
        )
    if declared_count != row_count:
        raise ReleaseValidationError(
            "dataset_info.json structure_count {} does not match {} metadata rows".format(
                declared_count, row_count
            )
        )

    definitions = info.get("classification_definitions")
    if not isinstance(definitions, dict):
        raise ReleaseValidationError("classification_definitions must be an object")
    for name, expected in CHECKER_PRESETS.items():
        if definitions.get(name) != list(expected):
            raise ReleaseValidationError(
                (
                    "dataset_info.json {} checker definition does not match "
                    "the package contract"
                ).format(name)
            )


def _validate_metadata_identity(rows: Sequence[Mapping[str, str]]) -> None:
    """Validate public naming and its denormalized identity fields."""

    for row_number, row in enumerate(rows, start=2):
        structure_id = row["structure_id"]
        match = _STRUCTURE_ID_RE.fullmatch(structure_id)
        if match is None:
            raise ReleaseValidationError(
                "metadata.csv:{} has invalid public structure_id {!r}".format(
                    row_number, structure_id
                )
            )
        variant, source, _, _ = match.groups()
        if row["structure_variant"] != variant:
            raise ReleaseValidationError(
                "{} structure_variant {!r} disagrees with its public ID".format(
                    structure_id, row["structure_variant"]
                )
            )
        if row["source_database"] != source:
            raise ReleaseValidationError(
                "{} source_database {!r} disagrees with its public ID".format(
                    structure_id, row["source_database"]
                )
            )
        if not row["source_id"].strip():
            raise ReleaseValidationError(
                "{} has an empty source_id".format(structure_id)
            )
        expected_cif = "cifs/{}.cif".format(structure_id)
        if row["cif_file"] != expected_cif:
            raise ReleaseValidationError(
                "{} must use CIF path {!r}, not {!r}".format(
                    structure_id, expected_cif, row["cif_file"]
                )
            )


def _parent_prefixes(fields: Sequence[str]) -> Tuple[str, ...]:
    prefixes = []
    for field in fields:
        if field.endswith("_status"):
            prefix = field[: -len("_status")]
            if (
                "{}_group".format(prefix) not in fields
                or "{}_size".format(prefix) not in fields
            ):
                raise ReleaseValidationError(
                    "parent criterion {!r} does not have status/group/size columns".format(
                        prefix
                    )
                )
            prefixes.append(prefix)
    if not prefixes:
        raise ReleaseValidationError("parent_groups.csv contains no parent criteria")
    return tuple(prefixes)


def _validate_parent_methods(
    methods: Mapping[str, object],
    info: Mapping[str, object],
    fields: Sequence[str],
) -> None:
    method_version = methods.get("dataset_version")
    if method_version != info.get("dataset_version"):
        raise ReleaseValidationError(
            "parent_group_methods.json dataset_version does not match dataset_info.json"
        )
    prefixes = set(_parent_prefixes(fields))
    declared = methods.get("csv_column_prefixes")
    if not isinstance(declared, dict):
        raise ReleaseValidationError("csv_column_prefixes must be an object")
    if declared != dict(PUBLIC_PARENT_METHOD_PREFIXES):
        raise ReleaseValidationError(
            "parent_group_methods.json csv_column_prefixes does not match the "
            "public method contract"
        )
    if set(declared.values()) != prefixes:
        raise ReleaseValidationError(
            "parent-group CSV prefixes do not match parent_group_methods.json"
        )


def _validate_parent_rows(
    fields: Sequence[str], rows: Sequence[Mapping[str, str]]
) -> Tuple[str, ...]:
    prefixes = _parent_prefixes(fields)
    counts = {prefix: Counter() for prefix in prefixes}
    parsed_sizes = {}

    for row_number, row in enumerate(rows, start=2):
        structure_id = row["structure_id"]
        for prefix in prefixes:
            status = row["{}_status".format(prefix)]
            group = row["{}_group".format(prefix)]
            raw_size = row["{}_size".format(prefix)]
            if status not in PARENT_STATUSES:
                raise ReleaseValidationError(
                    "parent_groups.csv:{} {} has unknown {} status {!r}".format(
                        row_number, structure_id, prefix, status
                    )
                )
            if not group:
                raise ReleaseValidationError(
                    "parent_groups.csv:{} {} has an empty {} group".format(
                        row_number, structure_id, prefix
                    )
                )
            expected_prefix = _PARENT_GROUP_PREFIXES.get(prefix)
            if expected_prefix is not None and re.fullmatch(
                re.escape(expected_prefix) + r"[0-9A-F]{8,64}", group
            ) is None:
                raise ReleaseValidationError(
                    "parent_groups.csv:{} {} has invalid {} group ID {!r}".format(
                        row_number, structure_id, prefix, group
                    )
                )
            try:
                size = int(raw_size)
            except (TypeError, ValueError):
                raise ReleaseValidationError(
                    "parent_groups.csv:{} {} has invalid {} size {!r}".format(
                        row_number, structure_id, prefix, raw_size
                    )
                )
            if size < 1:
                raise ReleaseValidationError(
                    "parent_groups.csv:{} {} has non-positive {} size".format(
                        row_number, structure_id, prefix
                    )
                )
            counts[prefix][group] += 1
            parsed_sizes[(structure_id, prefix)] = size

    for row in rows:
        structure_id = row["structure_id"]
        for prefix in prefixes:
            group = row["{}_group".format(prefix)]
            status = row["{}_status".format(prefix)]
            observed_size = counts[prefix][group]
            declared_size = parsed_sizes[(structure_id, prefix)]
            if declared_size != observed_size:
                raise ReleaseValidationError(
                    "{} {} declares size {} but group {} contains {} rows".format(
                        structure_id,
                        prefix,
                        declared_size,
                        group,
                        observed_size,
                    )
                )
            if observed_size > 1 and status != "MATCHED":
                raise ReleaseValidationError(
                    "{} {} group {} has {} rows but status is {}".format(
                        structure_id, prefix, group, observed_size, status
                    )
                )
            if observed_size == 1 and status == "MATCHED":
                raise ReleaseValidationError(
                    "{} {} singleton group cannot have MATCHED status".format(
                        structure_id, prefix
                    )
                )
    return prefixes


def _validate_published_labels(
    rows: Sequence[Mapping[str, str]], info: Mapping[str, object]
) -> None:
    observed_counts = {
        "label_3checker": Counter(),
        "label_4checker": Counter(),
        "label_5checker": Counter(),
    }
    for row in rows:
        structure_id = row["structure_id"]
        for checker, status_column in CHECKER_COLUMNS.items():
            status = row[status_column]
            if status not in PUBLIC_CHECKER_STATUSES:
                raise ReleaseValidationError(
                    "{} has non-public {} status {!r}; release metadata permits only "
                    "PASS, FAIL, and NOT_AVAILABLE".format(
                        structure_id, checker, status
                    )
                )
        for preset_name in CHECKER_PRESETS:
            label_column = "label_{}".format(preset_name)
            published = row[label_column]
            if published not in LABELS:
                raise ReleaseValidationError(
                    "{} has unknown published {} value {!r}".format(
                        structure_id, label_column, published
                    )
                )
            try:
                computed = classify_checker_row(row, preset_name)
            except (KeyError, ValueError, TypeError) as error:
                raise ReleaseValidationError(
                    "cannot classify {} for {}: {}".format(
                        structure_id, preset_name, error
                    )
                )
            if published != computed:
                raise ReleaseValidationError(
                    "{} {} is {!r}, but checker statuses recompute to {!r}".format(
                        structure_id, label_column, published, computed
                    )
                )
            observed_counts[label_column][published] += 1

    declared_counts = info.get("label_counts")
    if declared_counts is not None:
        if not isinstance(declared_counts, dict):
            raise ReleaseValidationError("dataset_info.json label_counts must be an object")
        for column, counts in observed_counts.items():
            declared = declared_counts.get(column)
            if declared is not None and declared != dict(counts):
                raise ReleaseValidationError(
                    "dataset_info.json {} does not match recomputed metadata counts".format(
                        column
                    )
                )


def _validate_cif_manifest(
    metadata: Mapping[str, Mapping[str, str]],
    manifest: Mapping[str, Mapping[str, str]],
    root: Path,
    verify_files: bool,
) -> None:
    for structure_id, row in manifest.items():
        cif_file = row["cif_file"]
        expected_cif = "cifs/{}.cif".format(structure_id)
        if cif_file != expected_cif:
            raise ReleaseValidationError(
                "{} must use manifest CIF path {!r}, not {!r}".format(
                    structure_id, expected_cif, cif_file
                )
            )
        metadata_cif = metadata[structure_id].get("cif_file")
        if metadata_cif != cif_file:
            raise ReleaseValidationError(
                "{} CIF path differs between metadata.csv and cif_manifest.csv".format(
                    structure_id
                )
            )
        try:
            size = int(row["size_bytes"])
        except (TypeError, ValueError):
            raise ReleaseValidationError(
                "{} has invalid CIF size {!r}".format(structure_id, row["size_bytes"])
            )
        if size < 1:
            raise ReleaseValidationError(
                "{} has a non-positive CIF size".format(structure_id)
            )
        if not _SHA256_RE.fullmatch(row["sha256"]):
            raise ReleaseValidationError(
                "{} has an invalid lowercase SHA-256".format(structure_id)
            )
        if verify_files:
            cif_path = root / cif_file
            if not cif_path.is_file():
                raise ReleaseValidationError(
                    "{} CIF file does not exist: {}".format(structure_id, cif_path)
                )
            try:
                cif_path.resolve().relative_to(root)
            except ValueError:
                raise ReleaseValidationError(
                    "{} CIF path resolves outside the release root".format(
                        structure_id
                    )
                )
            observed_size = cif_path.stat().st_size
            if observed_size != size:
                raise ReleaseValidationError(
                    "{} CIF size is {}, but the manifest declares {}".format(
                        structure_id, observed_size, size
                    )
                )
            observed_hash = _sha256_file(cif_path)
            if observed_hash != row["sha256"]:
                raise ReleaseValidationError(
                    "{} CIF SHA-256 does not match the manifest".format(structure_id)
                )


def _filter_set(
    values: Optional[Iterable[str]], description: str
) -> Optional[set]:
    if values is None:
        return None
    if isinstance(values, str):
        result = {values}
    else:
        try:
            result = set(values)
        except TypeError:
            raise TypeError("{} filter must be a string or iterable".format(description))
    if any(not isinstance(value, str) or not value for value in result):
        raise ValueError("{} filter values must be non-empty strings".format(description))
    return result


__all__ = [
    "ClassifiedDataset",
    "ClassifiedRecord",
    "CoREMOFDataset",
    "PARENT_METHOD_ALIASES",
    "PUBLIC_PARENT_METHOD_PREFIXES",
    "PUBLIC_CHECKER_STATUSES",
    "ParentGroup",
    "ReleaseValidationError",
    "StructureRecord",
    "normalize_parent_method",
]
