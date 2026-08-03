"""Parent-structure resolution for leakage-safe CoRE-MOF splits.

The release metadata contains several independent notions of a parent
structure.  This module deliberately keeps those notions separate.  The
``priority_main`` resolver is the opinionated public default: an available
RAC5 group is an anchor, MOFid v2 may attach otherwise unresolved structures,
and MOFid v1 is the final fallback.  Lower-priority evidence is never allowed
to merge two stronger components silently.

The implementation uses only the Python standard library and operates on a
small dataset protocol (``metadata_rows``, ``parent_by_id`` and
``cif_hashes``).  It therefore remains importable without any of the optional
scientific dependencies used elsewhere in :mod:`CoREMOF`.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from types import MappingProxyType
from typing import Dict, Iterable, List, Mapping, Optional, Set, Tuple


DIRECT_PARENT_METHODS = ("rac5", "mofid_v2", "mofid_v1")
REFERENCE_PARENT_METHODS = (
    "rac5_zeo",
    "zeo",
    "source_id",
    "common_name",
    "identity_union",
    "none",
)
COMPUTED_PARENT_METHODS = ("priority_main", "main_union")
PARENT_METHODS = DIRECT_PARENT_METHODS + REFERENCE_PARENT_METHODS + COMPUTED_PARENT_METHODS


_METHOD_KEYS = {
    "rac5": ("rac5", "rac", "rac5_group", "rac_group"),
    "mofid_v2": ("mofid_v2", "mofid2", "mofid_v2_group", "mofid2_group"),
    "mofid_v1": ("mofid_v1", "mofid1", "mofid_v1_group", "mofid1_group"),
    "rac5_zeo": ("rac5_zeo", "rac_zeo", "rac5_zeo_group", "rac_zeo_group"),
    "zeo": ("zeo", "zeo_group"),
    "source_id": ("source_id", "source", "source_id_group", "source_group"),
    "common_name": ("common_name", "name", "common_name_group", "name_group"),
    "identity_union": (
        "identity_union",
        "identity",
        "identity_union_group",
        "identity_group",
    ),
}

_RELEASE_BASE = {
    "rac5": "rac",
    "mofid_v2": "mofid2",
    "mofid_v1": "mofid1",
    "rac5_zeo": "rac_zeo",
    "zeo": "zeo",
    "source_id": "source",
    "common_name": "name",
    "identity_union": "identity",
}


_MISSING_TEXT = {"", "NA", "N/A", "NONE", "NULL", "NOT_AVAILABLE", "UNAVAILABLE"}


class _DisjointSet:
    """Small deterministic union-find implementation."""

    def __init__(self, values: Iterable[str]):
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
        if left_rank < right_rank or (left_rank == right_rank and left_root > right_root):
            left_root, right_root = right_root, left_root
            left_rank, right_rank = right_rank, left_rank
        self.parent[right_root] = left_root
        if left_rank == right_rank:
            self.rank[left_root] += 1


@dataclass(frozen=True)
class ParentConflict:
    """One lower-priority relation that spans stronger parent components.

    The conflict is diagnostic: it never authorizes merging the stronger
    components.  ``component_members`` retains the stronger component and the
    members of the lower-priority group already anchored to it.  Rows in
    ``unresolved_ids`` have only the conflicting lower-priority evidence and
    therefore cannot be placed by ``priority_main`` itself.
    """

    lower_method: str
    lower_group: str
    stronger_components: Tuple[str, ...]
    member_ids: Tuple[str, ...]
    unresolved_ids: Tuple[str, ...]
    component_members: Tuple[Tuple[str, Tuple[str, ...]], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "stronger_components", tuple(self.stronger_components))
        object.__setattr__(self, "member_ids", tuple(self.member_ids))
        object.__setattr__(self, "unresolved_ids", tuple(self.unresolved_ids))
        object.__setattr__(
            self,
            "component_members",
            tuple(
                (str(component), tuple(members))
                for component, members in self.component_members
            ),
        )


@dataclass(frozen=True)
class ParentResolution:
    """Resolved parent component for every usable structure.

    ``group_by_id`` and ``exclusions`` are disjoint mappings.  Exclusion
    reasons are stable machine-readable strings.  ``evidence_by_id`` records
    which evidence channels were available; it is diagnostic and is not used
    as an implicit classification label.
    """

    method: str
    group_by_id: Mapping[str, str]
    exclusions: Mapping[str, str]
    evidence_by_id: Mapping[str, Tuple[str, ...]]
    universe_ids: Tuple[str, ...]
    missing_parent: str
    conflicts: Tuple[ParentConflict, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "group_by_id", MappingProxyType(dict(self.group_by_id)))
        object.__setattr__(self, "exclusions", MappingProxyType(dict(self.exclusions)))
        object.__setattr__(
            self,
            "evidence_by_id",
            MappingProxyType(
                {
                    structure_id: tuple(channels)
                    for structure_id, channels in self.evidence_by_id.items()
                }
            ),
        )
        object.__setattr__(self, "conflicts", tuple(self.conflicts))

    @property
    def groups(self) -> Mapping[str, str]:
        """Alias for ``group_by_id`` for concise interactive use."""

        return self.group_by_id

    def members_by_group(self) -> Dict[str, Tuple[str, ...]]:
        """Return stable group-to-members mapping."""

        members: Dict[str, List[str]] = {}
        for structure_id in self.universe_ids:
            group = self.group_by_id.get(structure_id)
            if group is not None:
                members.setdefault(group, []).append(structure_id)
        return {group: tuple(values) for group, values in sorted(members.items())}


class ParentResolver:
    """Resolve release parent groups without mutating the dataset.

    Parameters
    ----------
    dataset
        Dataset-like object with ``metadata_rows`` and ``parent_by_id``.
        ``cif_hashes`` is additionally used by ``main_union``.
    missing_parent
        ``"singleton"`` (default) keeps each missing structure in its own
        non-matching group, as required by the release parent contract.
        ``"exclude"`` is an explicit stricter option.
    """

    def __init__(self, dataset, missing_parent: str = "singleton"):
        self.dataset = dataset
        self.missing_parent = self._validate_missing_parent(missing_parent)
        rows = tuple(getattr(dataset, "metadata_rows"))
        self._rows = rows
        self._row_by_id: Dict[str, Mapping[str, object]] = {}
        self._ids: List[str] = []
        for row in rows:
            structure_id = str(row.get("structure_id", "")).strip()
            if not structure_id:
                raise ValueError("Every metadata row must have a non-empty structure_id")
            if structure_id in self._row_by_id:
                raise ValueError("Duplicate structure_id in metadata_rows: %s" % structure_id)
            self._ids.append(structure_id)
            self._row_by_id[structure_id] = row
        self._parent_by_id = getattr(dataset, "parent_by_id", {}) or {}
        self._cif_hashes = getattr(dataset, "cif_hashes", {}) or {}
        self._validate_release_parent_table()

    @staticmethod
    def _validate_missing_parent(value: str) -> str:
        if value not in ("exclude", "singleton"):
            raise ValueError("missing_parent must be 'exclude' or 'singleton'")
        return value

    @property
    def universe_ids(self) -> Tuple[str, ...]:
        return tuple(self._ids)

    def _validate_release_parent_table(self) -> None:
        """Fail closed on inconsistent release-style status/group/size triples."""

        if not isinstance(self._parent_by_id, Mapping):
            return
        for method, keys in _METHOD_KEYS.items():
            group_keys = tuple(key for key in keys if key.endswith("_group"))
            declared: List[Tuple[str, Mapping[str, object], str]] = []
            for structure_id in self._ids:
                entry = self._parent_by_id.get(structure_id)
                if not isinstance(entry, Mapping):
                    continue
                for key in group_keys:
                    if key in entry:
                        group = str(entry[key]).strip()
                        if not group:
                            raise ValueError("%s release parent group may not be empty" % method)
                        declared.append((structure_id, entry, group))
                        break
            if not declared:
                continue
            counts: Dict[str, int] = {}
            for _, _, group in declared:
                counts[group] = counts.get(group, 0) + 1
            for structure_id, entry, group in declared:
                status = self._release_status(entry, method)
                self._validate_release_size(entry, method, status)
                base = _RELEASE_BASE[method]
                advertised = entry.get("%s_size" % base, entry.get("%s_size" % method))
                if int(advertised) != counts[group]:
                    raise ValueError(
                        "%s release parent size for %s does not match observed group %s"
                        % (method, structure_id, group)
                    )

    def resolve(
        self,
        method: str = "priority_main",
        structure_ids: Optional[Iterable[str]] = None,
        missing_parent: Optional[str] = None,
    ) -> ParentResolution:
        """Resolve *method* over the full dataset, then optionally subset it.

        Computing first and subsetting second is intentional: a structure
        filtered out of a particular experiment can still be a parent bridge
        and must not allow related structures to leak across split partitions.
        """

        if not isinstance(method, str):
            raise TypeError("parent method must be a string")
        method = method.strip().lower()
        if method not in PARENT_METHODS:
            raise ValueError(
                "Unknown parent method %r; choose one of %s"
                % (method, ", ".join(PARENT_METHODS))
            )
        missing = self._validate_missing_parent(missing_parent or self.missing_parent)
        if method == "none":
            resolution = self._resolve_none(missing)
        elif method == "priority_main":
            resolution = self._resolve_priority(missing)
        elif method == "main_union":
            resolution = self._resolve_main_union(missing)
        else:
            resolution = self._resolve_direct(method, missing)
        if structure_ids is None:
            return resolution
        requested = set(structure_ids)
        unknown = requested.difference(self._row_by_id)
        if unknown:
            raise KeyError("Unknown structure_id(s): %s" % ", ".join(sorted(unknown)))
        ordered = tuple(structure_id for structure_id in self._ids if structure_id in requested)
        return ParentResolution(
            method=resolution.method,
            group_by_id={
                structure_id: resolution.group_by_id[structure_id]
                for structure_id in ordered
                if structure_id in resolution.group_by_id
            },
            exclusions={
                structure_id: resolution.exclusions[structure_id]
                for structure_id in ordered
                if structure_id in resolution.exclusions
            },
            evidence_by_id={
                structure_id: resolution.evidence_by_id.get(structure_id, ())
                for structure_id in ordered
            },
            universe_ids=ordered,
            missing_parent=resolution.missing_parent,
            conflicts=tuple(
                conflict
                for conflict in resolution.conflicts
                if requested.intersection(conflict.member_ids)
            ),
        )

    def _parent_value(self, structure_id: str, method: str) -> Optional[str]:
        """Read either ID-major or method-major parent metadata."""

        keys = _METHOD_KEYS[method]
        source = self._parent_by_id
        candidate = None
        candidate_key = None
        release_entry = None
        if isinstance(source, Mapping):
            entry = source.get(structure_id)
            if isinstance(entry, Mapping):
                release_entry = entry
                for key in keys:
                    if key in entry:
                        candidate = entry[key]
                        candidate_key = key
                        break
            elif entry is not None and method == "identity_union":
                candidate = entry
            if candidate is None:
                for key in keys:
                    method_map = source.get(key)
                    if isinstance(method_map, Mapping) and structure_id in method_map:
                        candidate = method_map[structure_id]
                        candidate_key = key
                        break
        if (
            release_entry is not None
            and candidate_key is not None
            and candidate_key.endswith("_group")
        ):
            status = self._release_status(release_entry, method)
            self._validate_release_size(release_entry, method, status)
            if status == "NOT_AVAILABLE":
                return None
        if isinstance(candidate, Mapping):
            for key in ("group", "group_id", "parent_group", "value"):
                if key in candidate:
                    candidate = candidate[key]
                    break
            else:
                candidate = None
        if candidate is None:
            return None
        text = str(candidate).strip()
        if text.upper() in _MISSING_TEXT:
            return None
        return text

    @staticmethod
    def _release_status(entry: Mapping[str, object], method: str) -> str:
        base = _RELEASE_BASE[method]
        status = None
        for key in ("%s_status" % base, "%s_status" % method):
            if key in entry:
                status = str(entry[key]).strip().upper()
                break
        if status not in {"MATCHED", "UNMATCHED", "NOT_AVAILABLE"}:
            raise ValueError(
                "%s release parent status must be MATCHED, UNMATCHED, or "
                "NOT_AVAILABLE (got %r)" % (method, status)
            )
        return status

    @staticmethod
    def _validate_release_size(
        entry: Mapping[str, object], method: str, status: str
    ) -> None:
        base = _RELEASE_BASE[method]
        value = None
        for key in ("%s_size" % base, "%s_size" % method):
            if key in entry:
                value = entry[key]
                break
        try:
            size = int(value)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            raise ValueError("%s release parent size must be a positive integer" % method)
        if size < 1:
            raise ValueError("%s release parent size must be a positive integer" % method)
        if status == "MATCHED" and size < 2:
            raise ValueError("%s MATCHED parent group must have size >= 2" % method)
        if status in {"UNMATCHED", "NOT_AVAILABLE"} and size != 1:
            raise ValueError("%s %s parent group must have size 1" % (method, status))

    def _source_evidence(self, structure_id: str) -> Optional[str]:
        """Return a database-namespaced source group."""

        row = self._row_by_id[structure_id]
        database = str(row.get("source_database", "")).strip().upper()
        source_group = self._parent_value(structure_id, "source_id")
        if source_group is None:
            entry = (
                self._parent_by_id.get(structure_id)
                if isinstance(self._parent_by_id, Mapping)
                else None
            )
            if isinstance(entry, Mapping) and any(
                key in entry
                for key in ("source_group", "source_id_group", "source_status", "source_id_status")
            ):
                # A release row explicitly declared this criterion
                # unavailable.  Do not resurrect it from display metadata.
                return None
            source_id = str(row.get("source_id", "")).strip()
            if not source_id:
                return None
            source_group = " ".join(source_id.casefold().split())
        if not database:
            database = "UNKNOWN"
        return "%s:%s" % (database, source_group)

    def _full_cif_hash(self, structure_id: str) -> Optional[str]:
        candidate = None
        if isinstance(self._cif_hashes, Mapping):
            candidate = self._cif_hashes.get(structure_id)
        if candidate is None:
            row = self._row_by_id[structure_id]
            for key in ("cif_sha256", "sha256", "cif_hash_sha256"):
                if row.get(key):
                    candidate = row[key]
                    break
        if candidate is None:
            return None
        text = str(candidate).strip().lower()
        if len(text) != 64 or any(character not in "0123456789abcdef" for character in text):
            raise ValueError(
                "CIF integrity hashes used for main_union must be full lowercase SHA-256"
            )
        return text

    def _resolve_none(self, missing: str) -> ParentResolution:
        groups = {structure_id: "SINGLETON:%s" % structure_id for structure_id in self._ids}
        return ParentResolution(
            method="none",
            group_by_id=groups,
            exclusions={},
            evidence_by_id={structure_id: ("singleton",) for structure_id in self._ids},
            universe_ids=tuple(self._ids),
            missing_parent=missing,
        )

    def _resolve_direct(self, method: str, missing: str) -> ParentResolution:
        groups: Dict[str, str] = {}
        exclusions: Dict[str, str] = {}
        evidence: Dict[str, Tuple[str, ...]] = {}
        for structure_id in self._ids:
            group = self._parent_value(structure_id, method)
            if group is not None:
                groups[structure_id] = group
                evidence[structure_id] = (method,)
            elif missing == "singleton":
                groups[structure_id] = "SINGLETON:%s" % structure_id
                evidence[structure_id] = ()
            else:
                exclusions[structure_id] = "MISSING_PARENT_EVIDENCE"
                evidence[structure_id] = ()
        return ParentResolution(
            method=method,
            group_by_id=groups,
            exclusions=exclusions,
            evidence_by_id=evidence,
            universe_ids=tuple(self._ids),
            missing_parent=missing,
        )

    def _resolve_priority(self, missing: str) -> ParentResolution:
        raw = {
            method: {
                structure_id: self._parent_value(structure_id, method)
                for structure_id in self._ids
            }
            for method in DIRECT_PARENT_METHODS
        }
        assigned: Dict[str, str] = {}
        exclusions: Dict[str, str] = {}
        evidence: Dict[str, Tuple[str, ...]] = {}
        conflicts: List[ParentConflict] = []
        for structure_id in self._ids:
            evidence[structure_id] = tuple(
                method for method in DIRECT_PARENT_METHODS if raw[method][structure_id] is not None
            )

        rac_members = self._members_for_raw_groups(raw["rac5"])
        for group, members in sorted(rac_members.items()):
            component = "RAC5:%s" % group
            for structure_id in members:
                assigned[structure_id] = component

        for method, prefix in (("mofid_v2", "MOFID_V2"), ("mofid_v1", "MOFID_V1")):
            members_by_group = self._members_for_raw_groups(raw[method])
            for group, members in sorted(members_by_group.items()):
                active_members = [value for value in members if value not in exclusions]
                stronger = {
                    assigned[value] for value in active_members if value in assigned
                }
                unresolved = [value for value in active_members if value not in assigned]
                if len(stronger) > 1:
                    # Stronger anchors/components remain intact.  A row that
                    # only lower-priority evidence could place is unsafe.
                    stronger_components = tuple(sorted(stronger))
                    conflicts.append(
                        ParentConflict(
                            lower_method=method,
                            lower_group=group,
                            stronger_components=stronger_components,
                            member_ids=tuple(sorted(active_members)),
                            unresolved_ids=tuple(sorted(unresolved)),
                            component_members=tuple(
                                (
                                    component,
                                    tuple(
                                        sorted(
                                            value
                                            for value in active_members
                                            if assigned.get(value) == component
                                        )
                                    ),
                                )
                                for component in stronger_components
                            ),
                        )
                    )
                    for structure_id in unresolved:
                        exclusions[structure_id] = "PARENT_METHOD_CONFLICT"
                    continue
                if len(stronger) == 1:
                    component = next(iter(stronger))
                else:
                    component = "%s:%s" % (prefix, group)
                for structure_id in unresolved:
                    assigned[structure_id] = component

        for structure_id in self._ids:
            if structure_id in assigned or structure_id in exclusions:
                continue
            if missing == "singleton":
                assigned[structure_id] = "SINGLETON:%s" % structure_id
            else:
                exclusions[structure_id] = "MISSING_PARENT_EVIDENCE"

        return ParentResolution(
            method="priority_main",
            group_by_id={
                structure_id: assigned[structure_id]
                for structure_id in self._ids
                if structure_id in assigned
            },
            exclusions={
                structure_id: exclusions[structure_id]
                for structure_id in self._ids
                if structure_id in exclusions
            },
            evidence_by_id=evidence,
            universe_ids=tuple(self._ids),
            missing_parent=missing,
            conflicts=tuple(conflicts),
        )

    @staticmethod
    def _members_for_raw_groups(
        group_by_id: Mapping[str, Optional[str]]
    ) -> Dict[str, List[str]]:
        members: Dict[str, List[str]] = {}
        for structure_id, group in group_by_id.items():
            if group is not None:
                members.setdefault(group, []).append(structure_id)
        for group in members:
            members[group].sort()
        return members

    def _resolve_main_union(self, missing: str) -> ParentResolution:
        cif_hash_by_id = {
            structure_id: self._full_cif_hash(structure_id)
            for structure_id in self._ids
        }
        missing_hashes = [
            structure_id
            for structure_id, cif_hash in cif_hash_by_id.items()
            if cif_hash is None
        ]
        if missing_hashes:
            raise ValueError(
                "main_union requires a full CIF SHA-256 for every structure; "
                "missing {} of {} (examples: {})".format(
                    len(missing_hashes),
                    len(self._ids),
                    ", ".join(sorted(missing_hashes)[:5]),
                )
            )
        disjoint = _DisjointSet(self._ids)
        tokens_by_id: Dict[str, List[str]] = {structure_id: [] for structure_id in self._ids}
        members_by_token: Dict[str, List[str]] = {}

        for structure_id in self._ids:
            cif_hash = cif_hash_by_id[structure_id]
            token = "cif_sha256:%s" % cif_hash
            tokens_by_id[structure_id].append("cif_sha256")
            members_by_token.setdefault(token, []).append(structure_id)
            source = self._source_evidence(structure_id)
            if source is not None:
                token = "source_id:%s" % source
                tokens_by_id[structure_id].append("source_id")
                members_by_token.setdefault(token, []).append(structure_id)
            for method in DIRECT_PARENT_METHODS:
                group = self._parent_value(structure_id, method)
                if group is not None:
                    token = "%s:%s" % (method, group)
                    tokens_by_id[structure_id].append(method)
                    members_by_token.setdefault(token, []).append(structure_id)

        for token in sorted(members_by_token):
            members = sorted(members_by_token[token])
            for structure_id in members[1:]:
                disjoint.union(members[0], structure_id)

        components: Dict[str, List[str]] = {}
        available_ids: Set[str] = set()
        for structure_id in self._ids:
            if tokens_by_id[structure_id]:
                available_ids.add(structure_id)
                components.setdefault(disjoint.find(structure_id), []).append(structure_id)

        component_name: Dict[str, str] = {}
        used_names: Dict[str, Tuple[str, ...]] = {}
        for root, members in sorted(components.items(), key=lambda item: tuple(sorted(item[1]))):
            canonical_members = tuple(sorted(members))
            digest = hashlib.sha256("\0".join(canonical_members).encode("utf-8")).hexdigest()
            length = 16
            name = "MAIN-%s" % digest[:length]
            while name in used_names and used_names[name] != canonical_members:
                length += 1
                name = "MAIN-%s" % digest[:length]
            used_names[name] = canonical_members
            component_name[root] = name

        groups: Dict[str, str] = {}
        exclusions: Dict[str, str] = {}
        for structure_id in self._ids:
            if structure_id in available_ids:
                groups[structure_id] = component_name[disjoint.find(structure_id)]
            elif missing == "singleton":
                groups[structure_id] = "SINGLETON:%s" % structure_id
            else:
                exclusions[structure_id] = "MISSING_PARENT_EVIDENCE"

        return ParentResolution(
            method="main_union",
            group_by_id=groups,
            exclusions=exclusions,
            evidence_by_id={
                structure_id: tuple(tokens_by_id[structure_id]) for structure_id in self._ids
            },
            universe_ids=tuple(self._ids),
            missing_parent=missing,
        )


__all__ = [
    "COMPUTED_PARENT_METHODS",
    "DIRECT_PARENT_METHODS",
    "PARENT_METHODS",
    "REFERENCE_PARENT_METHODS",
    "ParentConflict",
    "ParentResolution",
    "ParentResolver",
]
