"""Lightweight access to a validated CoRE-MOF release.

The loader in this module uses only the Python standard library.  Importing it
does not import pandas, NumPy, crystallographic software, or network clients.
It treats release metadata as an audited contract: identifiers must agree
across tables, parent-group declarations must be internally consistent, and
all published checker labels are independently recomputed before use.

Project terminology in the embedded machine declarations is closed here.
Canonical identifier text is text processing only: convert a non-null value
to text, collapse Unicode whitespace runs to one ASCII space, trim, reject an
empty or declared whole-field placeholder, apply Unicode NFKC, then casefold.
It does not change CIF bytes, coordinates, atoms, occupancies, bonds, chemistry,
topology, or unit cells. ``identity_union`` means provisional source-ID/MOFid
transitive groups freshly computed per release: connect equal
database-namespaced ``(source_database, source_id)`` keys or equal complete
MOFid-v2/MOFid-v1 text and take connected components. Missing values add no
edge; it excludes CIF coordinates/hashes, RAC5, Zeo++, CrystalNets, common
names, DOI, and the optional strict StructureMatcher method, whose pinned
symmetric direct edges are authoritative, whose connected components are
convenience views rather than duplicate proof, and whose parser, timeout, or
execution errors are NOT_AVAILABLE rather than unmatched. The identity
group is not proof of identity or parentage.

``priority_main`` is the conflict-aware explanatory hierarchy over the full
release: seed exact release-authorized RAC5 groups, then process MOFid-v2 and
MOFid-v1 groups. Zero stronger components creates a component, exactly one
attaches unresolved rows, and two or more never merge, record
``PARENT_METHOD_CONFLICT``, and leave lower-only rows unresolved; remaining
missing evidence is a structure-specific singleton or explicit exclusion. It
does not use Zeo++, topology, source ID, common name, CIF hash, or
StructureMatcher. ``main_union`` is the separate leakage guard, not a parent
or explanatory method and not proof: before filtering, it takes full-release
transitive connected components over exact full lowercase 64-hex CIF SHA-256,
database-namespaced source siblings, and available release-authorized RAC5,
MOFid-v2, and MOFid-v1 edges. A CIF hash is mandatory and malformed/missing
input fails closed; missing optional evidence adds no edge. ``parent_only``
means only the selected explanatory parent relation supplies split blocks.
``auto`` selects ``main_union`` for ``priority_main`` and ``parent_only`` for
every explicit direct/reference parent method.

``rac5_crystalnets`` / ``RT-`` is an optional non-decisive reference requiring
exact equality of all 264 available finite RAC5 binary64 values plus the
complete current-success CrystalNets fingerprint; incomplete/error/timeout
input is ``NOT_AVAILABLE`` and adds no match. ``mofid_v2_crystalnets`` /
``M2T-`` replaces RAC5 with the complete canonical MOFid-v2 text. Its eligible
statuses are exactly ``SUCCESS``, ``SUCCESS_TOPOLOGY_UNKNOWN``,
``SUCCESS_TOPOLOGY_ERROR``, and ``SUCCESS_TOPOLOGY_TIMEOUT``; the latter two
are successful calculated identifiers with embedded topology qualifiers, not
execution failures. Every other status or incomplete input is
``NOT_AVAILABLE`` and adds no edge. It is provisional, non-decisive, excluded
from the default parent/leakage methods, and must be rebuilt if authorized
MOFid-v2 values change. ``structure_matcher_strict`` / ``SM-`` is the optional
connected-component convenience view over symmetric direct matches from the
pinned strict pymatgen protocol; direct edges are authoritative, components
are not duplicate proof, and parse/timeout/error cases are ``NOT_AVAILABLE``
rather than unmatched. It enters neither default parent nor leakage method.
"""

import csv
import hashlib
import io
import json
import os
import re
import stat
import unicodedata
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

from ._authority import (
    AuthorityStateError,
    IdentitySealRegistry,
    reject_sealed_copy,
    state_fingerprint,
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
M2T_ELIGIBLE_MOFID_V2_STATUSES = frozenset(
    {
        "SUCCESS",
        "SUCCESS_TOPOLOGY_UNKNOWN",
        "SUCCESS_TOPOLOGY_ERROR",
        "SUCCESS_TOPOLOGY_TIMEOUT",
    }
)
BASE_PUBLIC_PARENT_METHOD_PREFIXES = MappingProxyType(
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
OPTIONAL_PUBLIC_PARENT_METHOD_PREFIXES = MappingProxyType(
    {
        "rac5_crystalnets": "rac_crystalnets",
        "mofid_v2_crystalnets": "mofid2_crystalnets",
        "structure_matcher_strict": "sm",
    }
)
PUBLIC_PARENT_METHOD_PREFIXES = MappingProxyType(
    {
        **dict(BASE_PUBLIC_PARENT_METHOD_PREFIXES),
        **dict(OPTIONAL_PUBLIC_PARENT_METHOD_PREFIXES),
    }
)
PARENT_METHOD_ALIASES = MappingProxyType(
    {
        **dict(PUBLIC_PARENT_METHOD_PREFIXES),
        "rac_zeo": "rac_zeo",
        "rac_crystalnets": "rac_crystalnets",
        "mofid2_crystalnets": "mofid2_crystalnets",
        "structure_matcher": "sm",
        "structure_matcher_strict": "sm",
        "sm": "sm",
        "rac": "rac",
        "source": "source",
        "mofid2": "mofid2",
        "mofid1": "mofid1",
        "name": "name",
        "identity": "identity",
    }
)

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_RETIRED_PARENT_TOPOLOGY_KEY_RE = re.compile(
    r"^(?:rac5|rac|mofid_v2|mofid2)_topology(?:_(?:status|group|size))?$"
)
_RESERVED_AUTHORITY_KEYS = frozenset(
    {
        "_validated_release_token",
        "_official_checker_view_token",
        "_authority_generation_marker",
        "_authority_locked",
        "_authority_extra_state",
        "_authority_extra_bindings",
        "_authority_factory_token",
        "_authority_source",
        "_factory_token",
    }
)
_NORMALIZED_RESERVED_AUTHORITY_KEYS = frozenset(
    re.sub(
        r"[^a-z0-9]+",
        "",
        unicodedata.normalize("NFKC", key).casefold(),
    )
    for key in _RESERVED_AUTHORITY_KEYS
)
_NORMALIZED_PUBLIC_AUTHORITY_CLAIM_KEYS = frozenset(
    {
        "publicationauthorized",
        "publicationstatus",
        "publicationauthorization",
        "publicationready",
        "publishable",
        "authoritative",
        "officialsplit",
        "releaseauthority",
        "releasestatus",
        "releasestate",
        "candidatestatus",
        "candidatestate",
        "evidencestate",
        "prioritymainchanged",
        "mainunionchanged",
        "checkerviewofficial",
    }
)
_STAGED_RELEASE_STATUS = "PROVISIONAL_LATEST_AUDITED_SNAPSHOT"
_STAGED_REFERENCE_CANDIDATE_STATE = "STAGED_CANDIDATE"
_STAGED_REFERENCE_EVIDENCE_STATE = (
    "AUDITED_STAGED_CANDIDATE_REFERENCE; full input hashes and local "
    "evidence paths are retained only in the external integration receipt"
)
_RT_RELEASE_STATE = "STAGED_CANDIDATE_REFERENCE_NOT_PUBLICATION_AUTHORIZATION"
_M2T_RELEASE_STATE = (
    "PROVISIONAL_STAGED_CANDIDATE_REFERENCE_NOT_PUBLICATION_AUTHORIZATION; "
    "candidate-only while the null-unresolved MOFid bundle is stage-only"
)
_ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS = {
    ("dataset_info", "release_status"): _STAGED_RELEASE_STATUS,
}
_ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS.update(
    {
        (
            "parent_group_methods",
            "criteria",
            method,
            "decision_contract",
            "publication_authorized",
        ): False
        for method in ("rac5_crystalnets", "mofid_v2_crystalnets")
    }
)
_ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS.update(
    {
        (
            "parent_group_methods",
            "definitions",
            method,
            "decision_contract",
            "publication_authorized",
        ): False
        for method in ("rac5_crystalnets", "mofid_v2_crystalnets")
    }
)
_ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS.update(
    {
        (
            "dataset_info",
            "definitions",
            method,
            "decision_contract",
            "publication_authorized",
        ): False
        for method in ("rac5_crystalnets", "mofid_v2_crystalnets")
    }
)
_ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS.update(
    {
        (
            "dataset_info",
            "parent_grouping",
            "criterion_contracts",
            method,
            "decision_contract",
            "publication_authorized",
        ): False
        for method in ("rac5_crystalnets", "mofid_v2_crystalnets")
    }
)
for _method, _release_state in (
    ("rac5_crystalnets", _RT_RELEASE_STATE),
    ("mofid_v2_crystalnets", _M2T_RELEASE_STATE),
):
    for _contract_path in (
        ("parent_group_methods", "criteria", _method),
        ("parent_group_methods", "definitions", _method),
        ("dataset_info", "definitions", _method),
        (
            "dataset_info",
            "parent_grouping",
            "criterion_contracts",
            _method,
        ),
    ):
        _ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS[
            _contract_path + ("release_state",)
        ] = _release_state
for _integration_path in (
    ("parent_group_methods", "crystalnets_reference_integration"),
    (
        "dataset_info",
        "parent_grouping",
        "crystalnets_reference_integration",
    ),
):
    _ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS.update(
        {
            _integration_path + ("candidate_state",): (
                _STAGED_REFERENCE_CANDIDATE_STATE
            ),
            _integration_path + ("evidence_state",): (
                _STAGED_REFERENCE_EVIDENCE_STATE
            ),
            _integration_path + ("main_union_changed",): False,
            _integration_path + ("priority_main_changed",): False,
            _integration_path + ("publication_authorized",): False,
        }
    )
_ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS = MappingProxyType(
    _ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS
)
del _contract_path, _integration_path, _method, _release_state
_M2T_MOFID_V2_PLACEHOLDERS = frozenset(
    {
        "",
        "-",
        "nan",
        "none",
        "null",
        "n/a",
        "na",
        "unknown",
        "missing",
        "timeout",
        "timed out",
        "error",
        "failed",
        "fail",
        "fail process",
        "failed process",
        "process failed",
        "not_available",
        "not available",
        "unavailable",
    }
)
_STRUCTURE_ID_RE = re.compile(
    r"^(ASR|FSR|ION)-(COD|CSD|SI)-(\d{4}|UNKN)-(\d{4})$"
)
_PARENT_GROUP_PREFIXES = MappingProxyType(
    {
        "rac_zeo": "RZ-",
        "rac_crystalnets": "RT-",
        "mofid2_crystalnets": "M2T-",
        "sm": "SM-",
        "rac": "R-",
        "zeo": "Z-",
        "source": "S-",
        "mofid2": "M2-",
        "mofid1": "M1-",
        "name": "N-",
        "identity": "I-",
    }
)

_STRUCTURE_MATCHER_METHOD_ID = "pymatgen_structure_matcher_strict_v2"
_STRUCTURE_MATCHER_METHOD_SCHEMA = "coremof-structure-matcher-method/2.0"
_STRUCTURE_MATCHER_RECEIPT_SCHEMA = (
    "coremof-structure-matcher-release-adapter-receipt/1.0"
)
_STRUCTURE_MATCHER_RECEIPT_FILE = (
    "parent_groups/structure_matcher_strict_evidence_receipt.json"
)
_STRUCTURE_MATCHER_METHOD_CONTRACT = MappingProxyType(
    {
        "role": "OPTIONAL_REFERENCE",
        "method_id": _STRUCTURE_MATCHER_METHOD_ID,
        "method_schema_version": _STRUCTURE_MATCHER_METHOD_SCHEMA,
        "authoritative_evidence": "DIRECT_SYMMETRIC_STRICT_MATCH_EDGES",
        "fit_policy": "FIT_SYMMETRIC_TRUE_REQUIRED",
        "component_semantics": (
            "CONNECTED_COMPONENT_CONVENIENCE_DIRECT_EDGES_AUTHORITATIVE"
        ),
        "component_completeness_policy": (
            "INCOMPLETE_COMPONENTS_NOT_AVAILABLE_UNIQUE_SINGLETON"
        ),
        "public_status_policy": ["MATCHED", "UNMATCHED", "NOT_AVAILABLE"],
        "included_in_priority_main": False,
        "included_in_main_union": False,
        "historical_relaxed_executed": False,
        "historical_relaxed_exposed": False,
    }
)
_CRYSTALNETS_REFERENCE_METHODS = MappingProxyType(
    {
        "rac5_crystalnets": ("rac_crystalnets", "RT-"),
        "mofid_v2_crystalnets": ("mofid2_crystalnets", "M2T-"),
    }
)
_CRYSTALNETS_REFERENCE_CONTRACTS = MappingProxyType(
    json.loads(
        r'''
{
  "mofid_v2_crystalnets": {
    "comparison": "exact equality of the complete combined key; zero tolerance",
    "decision_contract": {
      "classification_role": "NON_DECISIVE_REFERENCE_ONLY",
      "explicit_direct_or_reference_selection": {
        "leakage_guard_auto": "RESOLVE_TO_PARENT_ONLY",
        "split_block_source": "SELECTED_PARENT_RELATION_ONLY"
      },
      "field_meanings": {
        "authorized_criterion_order": "the release-authorized explanatory priority sequence",
        "classification_role": "how this criterion may be used for classification",
        "closure": "how direct leakage edges become leakage blocks",
        "construction_scope": "which release rows exist when leakage blocks are constructed",
        "decision_contract": "the closed machine-readable non-decisive use rules",
        "edge_relation_logic": "whether one or every listed relation is needed for a leakage edge",
        "full_cif_content_digest_requirement": "the mandatory complete lowercase 64-hex CIF-content digest rule; missing or malformed input fails construction closed; this method contract does not expose the digest value",
        "mofid_v2_crystalnets_rebuild_trigger": "the upstream MOFid-v2 change that invalidates M2T groups",
        "semantic_role": "what a leakage block means and does not prove"
      },
      "main_union": {
        "changed": false,
        "closure": "TRANSITIVE_CONNECTED_COMPONENT",
        "conflict_rows": "RETAIN_WITH_PARENT_METHOD_CONFLICT_DIAGNOSTIC",
        "construction_scope": "FULL_RELEASE_BEFORE_FILTERING",
        "consumed": false,
        "edge_relation_logic": "ANY_ONE",
        "edge_relations": [
          "exact_full_cif_sha256",
          "database_namespaced_source_sibling",
          "available_release_authorized_rac5",
          "available_release_authorized_mofid_v2",
          "available_release_authorized_mofid_v1"
        ],
        "full_cif_content_digest_requirement": {
          "algorithm": "SHA-256",
          "encoding": "LOWERCASE_HEX",
          "length": 64,
          "missing_or_malformed": "FAIL_CLOSED",
          "required": true
        },
        "missing_optional_relation": "NO_EDGE",
        "parent_conflict_ledger": "RECORD_EVERY_SPANNING_LOWER_GROUP",
        "parent_only_conflict_rows": "EXCLUDE_UNRESOLVED_ROW",
        "semantic_role": "SPLIT_LEAKAGE_GUARD_NOT_PARENT_PROOF"
      },
      "mofid_v2_crystalnets_rebuild_trigger": "REBUILD_M2T_IF_AUTHORIZED_MOFID_V2_VALUES_CHANGE",
      "priority_main": {
        "authorized_criterion_order": [
          "rac5",
          "mofid_v2",
          "mofid_v1"
        ],
        "changed": false,
        "consumed": false,
        "excluded_evidence": [
          "zeo",
          "crystalnets",
          "source_id",
          "common_name",
          "cif_sha256",
          "structure_matcher"
        ],
        "leakage_guard_auto": "RESOLVE_TO_MAIN_UNION",
        "missing_parent_exclude": "EXCLUDE_ROW",
        "multiple_stronger_components": "DO_NOT_MERGE_RECORD_PARENT_METHOD_CONFLICT",
        "one_stronger_component": "ATTACH_UNRESOLVED_ROWS",
        "remaining_unresolved_default": "STRUCTURE_SPECIFIC_SINGLETON",
        "zero_stronger_components": "CREATE_COMPONENT_FROM_UNRESOLVED_ROWS"
      },
      "publication_authorized": false
    },
    "display_name": "exact MOFid-v2 plus current-success CrystalNets groups",
    "excludes": [
      "RAC5",
      "Zeo++",
      "source ID",
      "DOI",
      "common name",
      "CIF hash",
      "StructureMatcher"
    ],
    "group_prefix": "M2T-",
    "inputs": {
      "crystalnets": {
        "availability_gate": "effective execution_status=SUCCESS, topology_available=true, error=null, and a nonempty list of complete SingleNodes/AllNodes subnet results",
        "canonicalization": "project only the included scientific fields; sort each projected subnet by canonical JSON with sorted object keys; retain duplicate subnet multiplicity; serialize the complete projection as canonical JSON and hash with SHA-256",
        "excluded_fields": [
          "runtime_seconds",
          "paths",
          "CIF hashes",
          "errors and diagnostics",
          "software and method boilerplate",
          "original subnet index or order"
        ],
        "included_fields": [
          "network_dimension",
          "interpenetrated_subnet_count",
          "catenation_degree",
          "subnet_count",
          "single_node_net",
          "all_node_net",
          "single_all_agree",
          "every canonically sorted subnet's single_all_agree and both node status, dimension, topology_key, topology_name, and topological_genome"
        ]
      },
      "mofid_v2": {
        "canonical_text": {
          "comparison": "compare the complete resulting MOFid-v2 string exactly",
          "operations_in_order": [
            "convert a non-null value to text",
            "collapse every run of Unicode whitespace to one ASCII space",
            "strip leading and trailing whitespace",
            "reject an empty value or a declared whole-field placeholder",
            "apply Unicode NFKC",
            "apply Unicode casefold()"
          ],
          "scientific_effect": "No CIF byte, coordinate, atom, occupancy, bond, chemistry, topology, or unit-cell value is changed."
        },
        "required_statuses": [
          "SUCCESS",
          "SUCCESS_TOPOLOGY_ERROR",
          "SUCCESS_TOPOLOGY_TIMEOUT",
          "SUCCESS_TOPOLOGY_UNKNOWN"
        ],
        "status_gate": "SUCCESS, SUCCESS_TOPOLOGY_UNKNOWN, SUCCESS_TOPOLOGY_ERROR, and SUCCESS_TOPOLOGY_TIMEOUT are the only eligible statuses; SUCCESS_TOPOLOGY_ERROR and SUCCESS_TOPOLOGY_TIMEOUT are successful calculated identifiers whose embedded topology qualifier is ERROR or TIMEOUT, not execution failures; every other MOFid-v2 status is ineligible and adds no match; the projected parent result is a structure-specific NOT_AVAILABLE singleton"
      }
    },
    "key": "mofid_v2_crystalnets",
    "missing_error_behavior": "any incomplete, missing, error, timeout, partial, or non-success input adds no match and receives a structure-specific NOT_AVAILABLE singleton; nulls never match one another",
    "purpose": "non-decisive combined identifier/topology reference screening",
    "relation_to_main_union": "not consumed; main_union is unchanged. Before label, source, variant, metal, ID, target, or other row filtering across the full release, it forms a split-leakage graph: add a direct edge for any one of an equal exact full CIF SHA-256, an equal database-namespaced source-sibling key, an equal available release-authorized RAC5 group, an equal available release-authorized MOFid-v2 group, or an equal available release-authorized MOFid-v1 group; then take transitive connected-component closure. Full CIF SHA-256 is mandatory for every row; if it is missing or unavailable, main_union construction fails closed instead of emitting a block. Missing or unavailable optional source, RAC5, MOFid-v2, or MOFid-v1 evidence adds no corresponding edge. It is not an explanatory parent claim or proof of identity",
    "relation_to_priority_main": "not consumed; priority_main is unchanged. Across the complete release, available release-authorized exact RAC5 groups seed components. For available release-authorized exact MOFid-v2 groups and then MOFid-v1 groups: when a group touches zero stronger components, its unresolved rows form a new component; when it touches exactly one, its unresolved rows attach to that component; when it touches two or more, it never merges the stronger components, records PARENT_METHOD_CONFLICT, and its lower-only unresolved rows remain unresolved. Under the default splitting behavior, any row still unresolved after priority processing is a structure-specific singleton; when a caller explicitly selects missing_parent=exclude, that row is excluded instead. main_union keeps conflict cases inside leakage blocks, while parent_only uses only the selected parent relation for split blocks and excludes unresolved rows. It excludes Zeo++, CrystalNets, source ID, common name, CIF hash, and StructureMatcher, including this CrystalNets-combined criterion",
    "release_state": "PROVISIONAL_STAGED_CANDIDATE_REFERENCE_NOT_PUBLICATION_AUTHORIZATION; candidate-only while the null-unresolved MOFid bundle is stage-only",
    "result_mapping": {
      "MATCHED": "the complete combined key belongs to an equivalence class with more than one structure; size is the complete class size",
      "NOT_AVAILABLE": "the combined key is incomplete or an input is non-success; the result is a structure-specific singleton with size 1",
      "UNMATCHED": "the combined key is complete but its equivalence class contains only this structure; size is 1"
    },
    "scope": "non-decisive reference screening evidence, not proof of crystallographic identity or common parentage",
    "transitivity": "one equal-key equivalence class; no graph chaining"
  },
  "rac5_crystalnets": {
    "comparison": "exact equality of the complete combined key; zero tolerance",
    "decision_contract": {
      "classification_role": "NON_DECISIVE_REFERENCE_ONLY",
      "explicit_direct_or_reference_selection": {
        "leakage_guard_auto": "RESOLVE_TO_PARENT_ONLY",
        "split_block_source": "SELECTED_PARENT_RELATION_ONLY"
      },
      "field_meanings": {
        "authorized_criterion_order": "the release-authorized explanatory priority sequence",
        "classification_role": "how this criterion may be used for classification",
        "closure": "how direct leakage edges become leakage blocks",
        "construction_scope": "which release rows exist when leakage blocks are constructed",
        "decision_contract": "the closed machine-readable non-decisive use rules",
        "edge_relation_logic": "whether one or every listed relation is needed for a leakage edge",
        "full_cif_content_digest_requirement": "the mandatory complete lowercase 64-hex CIF-content digest rule; missing or malformed input fails construction closed; this method contract does not expose the digest value",
        "semantic_role": "what a leakage block means and does not prove"
      },
      "main_union": {
        "changed": false,
        "closure": "TRANSITIVE_CONNECTED_COMPONENT",
        "conflict_rows": "RETAIN_WITH_PARENT_METHOD_CONFLICT_DIAGNOSTIC",
        "construction_scope": "FULL_RELEASE_BEFORE_FILTERING",
        "consumed": false,
        "edge_relation_logic": "ANY_ONE",
        "edge_relations": [
          "exact_full_cif_sha256",
          "database_namespaced_source_sibling",
          "available_release_authorized_rac5",
          "available_release_authorized_mofid_v2",
          "available_release_authorized_mofid_v1"
        ],
        "full_cif_content_digest_requirement": {
          "algorithm": "SHA-256",
          "encoding": "LOWERCASE_HEX",
          "length": 64,
          "missing_or_malformed": "FAIL_CLOSED",
          "required": true
        },
        "missing_optional_relation": "NO_EDGE",
        "parent_conflict_ledger": "RECORD_EVERY_SPANNING_LOWER_GROUP",
        "parent_only_conflict_rows": "EXCLUDE_UNRESOLVED_ROW",
        "semantic_role": "SPLIT_LEAKAGE_GUARD_NOT_PARENT_PROOF"
      },
      "priority_main": {
        "authorized_criterion_order": [
          "rac5",
          "mofid_v2",
          "mofid_v1"
        ],
        "changed": false,
        "consumed": false,
        "excluded_evidence": [
          "zeo",
          "crystalnets",
          "source_id",
          "common_name",
          "cif_sha256",
          "structure_matcher"
        ],
        "leakage_guard_auto": "RESOLVE_TO_MAIN_UNION",
        "missing_parent_exclude": "EXCLUDE_ROW",
        "multiple_stronger_components": "DO_NOT_MERGE_RECORD_PARENT_METHOD_CONFLICT",
        "one_stronger_component": "ATTACH_UNRESOLVED_ROWS",
        "remaining_unresolved_default": "STRUCTURE_SPECIFIC_SINGLETON",
        "zero_stronger_components": "CREATE_COMPONENT_FROM_UNRESOLVED_ROWS"
      },
      "publication_authorized": false
    },
    "display_name": "exact RAC5 plus current-success CrystalNets groups",
    "excludes": [
      "MOFid",
      "Zeo++",
      "source ID",
      "DOI",
      "common name",
      "CIF hash",
      "StructureMatcher"
    ],
    "group_prefix": "RT-",
    "inputs": {
      "crystalnets": {
        "availability_gate": "effective execution_status=SUCCESS, topology_available=true, error=null, and a nonempty list of complete SingleNodes/AllNodes subnet results",
        "canonicalization": "project only the included scientific fields; sort each projected subnet by canonical JSON with sorted object keys; retain duplicate subnet multiplicity; serialize the complete projection as canonical JSON and hash with SHA-256",
        "excluded_fields": [
          "runtime_seconds",
          "paths",
          "CIF hashes",
          "errors and diagnostics",
          "software and method boilerplate",
          "original subnet index or order"
        ],
        "included_fields": [
          "network_dimension",
          "interpenetrated_subnet_count",
          "catenation_degree",
          "subnet_count",
          "single_node_net",
          "all_node_net",
          "single_all_agree",
          "every canonically sorted subnet's single_all_agree and both node status, dimension, topology_key, topology_name, and topological_genome"
        ]
      },
      "rac5": "all 264 ordered finite depth-5 RAC descriptors parsed as IEEE-754 binary64, with negative zero mapped to positive zero and serialized with float.hex()"
    },
    "key": "rac5_crystalnets",
    "missing_error_behavior": "any incomplete, missing, error, timeout, partial, or non-success input adds no match and receives a structure-specific NOT_AVAILABLE singleton; nulls never match one another",
    "numeric_tolerance": {
      "absolute": 0.0,
      "relative": 0.0,
      "scaling": "none"
    },
    "purpose": "high-specificity non-decisive reference screening",
    "relation_to_main_union": "not consumed; main_union is unchanged. Before label, source, variant, metal, ID, target, or other row filtering across the full release, it forms a split-leakage graph: add a direct edge for any one of an equal exact full CIF SHA-256, an equal database-namespaced source-sibling key, an equal available release-authorized RAC5 group, an equal available release-authorized MOFid-v2 group, or an equal available release-authorized MOFid-v1 group; then take transitive connected-component closure. Full CIF SHA-256 is mandatory for every row; if it is missing or unavailable, main_union construction fails closed instead of emitting a block. Missing or unavailable optional source, RAC5, MOFid-v2, or MOFid-v1 evidence adds no corresponding edge. It is not an explanatory parent claim or proof of identity",
    "relation_to_priority_main": "not consumed; priority_main is unchanged. Across the complete release, available release-authorized exact RAC5 groups seed components. For available release-authorized exact MOFid-v2 groups and then MOFid-v1 groups: when a group touches zero stronger components, its unresolved rows form a new component; when it touches exactly one, its unresolved rows attach to that component; when it touches two or more, it never merges the stronger components, records PARENT_METHOD_CONFLICT, and its lower-only unresolved rows remain unresolved. Under the default splitting behavior, any row still unresolved after priority processing is a structure-specific singleton; when a caller explicitly selects missing_parent=exclude, that row is excluded instead. main_union keeps conflict cases inside leakage blocks, while parent_only uses only the selected parent relation for split blocks and excludes unresolved rows. It excludes Zeo++, CrystalNets, source ID, common name, CIF hash, and StructureMatcher, including this CrystalNets-combined criterion",
    "release_state": "STAGED_CANDIDATE_REFERENCE_NOT_PUBLICATION_AUTHORIZATION",
    "result_mapping": {
      "MATCHED": "the complete combined key belongs to an equivalence class with more than one structure; size is the complete class size",
      "NOT_AVAILABLE": "the combined key is incomplete or an input is non-success; the result is a structure-specific singleton with size 1",
      "UNMATCHED": "the combined key is complete but its equivalence class contains only this structure; size is 1"
    },
    "scope": "non-decisive reference screening evidence, not proof of crystallographic identity or common parentage",
    "transitivity": "one equal-key equivalence class; no graph chaining"
  }
}
'''
    )
)
_CRYSTALNETS_REFERENCE_INTEGRATION = MappingProxyType(
    {
        "candidate_state": _STAGED_REFERENCE_CANDIDATE_STATE,
        "evidence_state": _STAGED_REFERENCE_EVIDENCE_STATE,
        "main_union_changed": False,
        "priority_main_changed": False,
        "publication_authorized": False,
    }
)

# These identities are intentionally module-private.  Public constructors can
# build useful dataset-like and classified objects, but they cannot mint a
# claim that release files were validated or that a canonical checker preset
# was recomputed by this package.
_VALIDATED_RELEASE_TOKEN = object()
_OFFICIAL_CHECKER_VIEW_TOKEN = object()
_DATASET_FACTORY_TOKEN = object()
_CLASSIFIED_FACTORY_TOKEN = object()
_DATASET_GENERATIONS = IdentitySealRegistry("validated dataset generation")
_CLASSIFIED_GENERATIONS = IdentitySealRegistry("official checker-view generation")


class ReleaseValidationError(ValueError):
    """Raised when release tables do not satisfy the public data contract."""


def _normalized_declaration_key(value: str) -> str:
    """Normalize Unicode compatibility, case, and punctuation for key checks."""

    return re.sub(
        r"[^a-z0-9]+",
        "",
        unicodedata.normalize("NFKC", value).casefold(),
    )


def _reject_retired_or_reserved_keys(
    value: object,
    location: str,
    *,
    _root: Optional[str] = None,
    _path: Tuple[str, ...] = (),
) -> None:
    """Close retired aliases and public/private authority claims recursively.

    Public JSON/CSV keys are compared after Unicode NFKC, case folding, and
    punctuation removal so compatibility characters or spelling variants
    cannot bypass the authority boundary.  Authority/state-shaped fields are
    accepted only at exact trusted paths with exact type-sensitive staged or
    non-decisive values.
    """

    root = location if _root is None else _root

    if isinstance(value, Mapping):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ReleaseValidationError(
                    "{} contains a non-string declaration key".format(location)
                )
            normalized = _normalized_declaration_key(key)
            path = (root,) + _path + (key,)
            retired = re.fullmatch(
                r"(?:rac5|rac|mofidv2|mofid2)topology(?:status|group|size)?",
                normalized,
            ) is not None
            private_reserved = normalized in _NORMALIZED_RESERVED_AUTHORITY_KEYS
            public_authority = normalized in _NORMALIZED_PUBLIC_AUTHORITY_CLAIM_KEYS
            allowed_public_state = (
                path in _ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS
                and type(item)
                is type(_ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS[path])
                and item == _ALLOWED_PUBLIC_AUTHORITY_STATE_DECLARATIONS[path]
            )
            if retired or private_reserved or (
                public_authority and not allowed_public_state
            ):
                raise ReleaseValidationError(
                    "{} contains retired or reserved declaration key {!r}; an "
                    "exact closed terminal reference contract and exact validated "
                    "release are required".format(location, key)
                )
            _reject_retired_or_reserved_keys(
                item,
                "{}.{}".format(location, key),
                _root=root,
                _path=_path + (key,),
            )
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _reject_retired_or_reserved_keys(
                item,
                "{}[{}]".format(location, index),
                _root=root,
                _path=_path + ("[{}]".format(index),),
            )


def _canonical_complete_mofid_v2(value: object, structure_id: str) -> str:
    """Return complete M2T MOFid-v2 text under the release canonicalization."""

    if not isinstance(value, str):
        raise ReleaseValidationError(
            "{} M2T evidence requires mofid_v2 to be an exact string".format(
                structure_id
            )
        )
    collapsed = " ".join(value.split())
    if collapsed.casefold() in _M2T_MOFID_V2_PLACEHOLDERS:
        raise ReleaseValidationError(
            "{} M2T evidence requires a complete nonplaceholder mofid_v2".format(
                structure_id
            )
        )
    canonical = unicodedata.normalize("NFKC", collapsed).casefold()
    if not canonical or canonical in _M2T_MOFID_V2_PLACEHOLDERS:
        raise ReleaseValidationError(
            "{} M2T evidence requires a complete nonplaceholder mofid_v2".format(
                structure_id
            )
        )
    return canonical


@dataclass(frozen=True)
class _ReleaseFileSnapshot:
    """One stable release-file byte generation and its integrity values."""

    path: Path
    data: bytes
    sha256: str
    size_bytes: int


def _release_stat_signature(
    value: os.stat_result,
) -> Tuple[int, int, int, int, int]:
    return (
        int(value.st_dev),
        int(value.st_ino),
        int(value.st_size),
        int(getattr(value, "st_mtime_ns", int(value.st_mtime * 1_000_000_000))),
        int(getattr(value, "st_ctime_ns", int(value.st_ctime * 1_000_000_000))),
    )


def _read_release_snapshot_bytes(handle) -> bytes:
    return handle.read()


def _capture_release_file(
    path: Path, description: str
) -> _ReleaseFileSnapshot:
    """Capture one stable regular file and verify it through the same open fd."""

    if not path.is_file():
        raise FileNotFoundError(
            "{} file does not exist: {}".format(description, path)
        )
    path_before = path.stat()
    if not stat.S_ISREG(path_before.st_mode):
        raise ReleaseValidationError(
            "{} is not a regular file: {}".format(description, path)
        )
    with path.open("rb") as handle:
        descriptor_before = os.fstat(handle.fileno())
        if _release_stat_signature(path_before) != _release_stat_signature(
            descriptor_before
        ):
            raise ReleaseValidationError(
                "{} changed or was replaced while being opened: {}".format(
                    description, path
                )
            )
        data = _read_release_snapshot_bytes(handle)
        snapshot_sha256 = hashlib.sha256(data).hexdigest()
        handle.seek(0)
        verification_digest = hashlib.sha256()
        verification_size = 0
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            verification_digest.update(block)
            verification_size += len(block)
        descriptor_after = os.fstat(handle.fileno())
    path_after = path.stat()
    signatures = {
        _release_stat_signature(path_before),
        _release_stat_signature(descriptor_before),
        _release_stat_signature(descriptor_after),
        _release_stat_signature(path_after),
    }
    if (
        len(signatures) != 1
        or len(data) != descriptor_after.st_size
        or verification_size != len(data)
        or verification_digest.hexdigest() != snapshot_sha256
    ):
        raise ReleaseValidationError(
            "{} changed or was replaced during byte capture: {}".format(
                description, path
            )
        )
    return _ReleaseFileSnapshot(
        path=path,
        data=data,
        sha256=snapshot_sha256,
        size_bytes=len(data),
    )


def _deep_freeze(value: object) -> object:
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or key != key.strip():
                raise TypeError("mapping keys must be exact nonblank strings")
            result[key] = _deep_freeze(item)
        return MappingProxyType(result)
    if isinstance(value, (list, tuple)):
        return tuple(_deep_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        raise TypeError("unordered sets are not valid release contract values")
    return value


def _strict_json_equal(actual: object, expected: object) -> bool:
    """Compare decoded JSON without Python's bool/int equality shortcuts."""

    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return set(actual) == set(expected) and all(
            _strict_json_equal(actual[key], expected[key]) for key in expected
        )
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            _strict_json_equal(actual_item, expected_item)
            for actual_item, expected_item in zip(actual, expected)
        )
    return actual == expected


def _strict_frozen_json_equal(actual: object, expected: object) -> bool:
    """Compare a frozen in-memory JSON value to one decoded JSON contract."""

    if isinstance(expected, dict):
        return isinstance(actual, Mapping) and set(actual) == set(expected) and all(
            _strict_frozen_json_equal(actual[key], expected[key])
            for key in expected
        )
    if isinstance(expected, list):
        return isinstance(actual, (list, tuple)) and len(actual) == len(expected) and all(
            _strict_frozen_json_equal(actual_item, expected_item)
            for actual_item, expected_item in zip(actual, expected)
        )
    return type(actual) is type(expected) and actual == expected


def _has_authenticated_crystalnets_reference_contract(
    dataset: object, method: str
) -> bool:
    """Return whether every public RT/M2T declaration copy is exact.

    This check is deliberately reusable by the generic dataset-like parent
    resolver.  A syntactically valid ``RT-`` or ``M2T-`` triad alone does not
    authorize use of a release-defined scientific method name.
    """

    if method not in _CRYSTALNETS_REFERENCE_METHODS:
        return False
    methods = getattr(dataset, "parent_group_methods", None)
    info = getattr(dataset, "dataset_info", None)
    if not isinstance(methods, Mapping) or not isinstance(info, Mapping):
        return False
    try:
        _reject_retired_or_reserved_keys(methods, "parent_group_methods")
        _reject_retired_or_reserved_keys(info, "dataset_info")
    except ReleaseValidationError:
        return False
    prefix, _ = _CRYSTALNETS_REFERENCE_METHODS[method]
    declared = methods.get("csv_column_prefixes")
    criteria = methods.get("criteria")
    definitions = methods.get("definitions")
    methods_integration = methods.get("crystalnets_reference_integration")
    info_definitions = info.get("definitions")
    parent_grouping = info.get("parent_grouping")
    if not all(
        isinstance(value, Mapping)
        for value in (
            declared,
            criteria,
            definitions,
            info_definitions,
            parent_grouping,
        )
    ):
        return False
    info_contracts = parent_grouping.get("criterion_contracts")
    info_integration = parent_grouping.get("crystalnets_reference_integration")
    if not isinstance(info_contracts, Mapping):
        return False
    expected = _CRYSTALNETS_REFERENCE_CONTRACTS[method]
    integration = dict(_CRYSTALNETS_REFERENCE_INTEGRATION)
    return (
        declared.get(method) == prefix
        and _strict_frozen_json_equal(criteria.get(method), expected)
        and _strict_frozen_json_equal(definitions.get(method), expected)
        and _strict_frozen_json_equal(methods_integration, integration)
        and _strict_frozen_json_equal(info_definitions.get(method), expected)
        and _strict_frozen_json_equal(info_contracts.get(method), expected)
        and _strict_frozen_json_equal(info_integration, integration)
    )


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


def _exact_nonblank_text(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value != value.strip()
    ):
        raise AuthorityStateError("{} must be an exact nonblank string".format(name))
    return value


def _record_authority_state(record: object) -> Tuple[object, ...]:
    """Return one exact, type-sensitive record payload for generation sealing."""

    if not isinstance(record, StructureRecord):
        raise AuthorityStateError(
            "authenticated datasets must contain StructureRecord rows"
        )
    structure_id = _exact_nonblank_text(record.structure_id, "structure_id")
    if not isinstance(record.metadata, Mapping):
        raise AuthorityStateError("record metadata must be a mapping")
    if record.metadata.get("structure_id") != structure_id:
        raise AuthorityStateError(
            "record structure_id differs from its metadata structure_id"
        )
    if not isinstance(record.parent_groups, Mapping):
        raise AuthorityStateError("record parent_groups must be a mapping")
    if isinstance(record.parent_groups, _ParentGroupView):
        parent_state = (
            "release_parent_view",
            record.parent_groups._prefixes,
            state_fingerprint(record.parent_groups._row),
        )
    else:
        groups = []
        for prefix, group in record.parent_groups.items():
            if (
                not isinstance(prefix, str)
                or not prefix
                or prefix != prefix.strip()
                or not isinstance(group, ParentGroup)
            ):
                raise AuthorityStateError(
                    "record parent_groups must map exact strings to ParentGroup values"
                )
            groups.append(
                (prefix, group.method, group.status, group.group_id, group.size)
            )
        parent_state = ("materialized_parent_groups", tuple(groups))
    manifest = record.cif_manifest
    if manifest is not None and not isinstance(manifest, Mapping):
        raise AuthorityStateError("record CIF manifest must be a mapping or None")
    return (
        structure_id,
        state_fingerprint(record.metadata),
        parent_state,
        state_fingerprint(manifest),
    )


def _dataset_authority_payload(dataset: object) -> Tuple[object, ...]:
    """Validate and expose all state that can support an authoritative claim."""

    records = getattr(dataset, "records", None)
    metadata_rows = getattr(dataset, "metadata_rows", None)
    by_id = getattr(dataset, "_by_id", None)
    if type(records) is not tuple or type(metadata_rows) is not tuple:
        raise AuthorityStateError(
            "authenticated dataset records and metadata_rows must be tuples"
        )
    if not isinstance(by_id, Mapping):
        raise AuthorityStateError("authenticated dataset _by_id must be a mapping")
    if len(records) != len(metadata_rows) or len(records) != len(by_id):
        raise AuthorityStateError(
            "authenticated dataset record, metadata, and index sizes differ"
        )

    record_digests = []
    structure_ids = []
    seen_structure_ids = set()
    for index, record in enumerate(records):
        state = _record_authority_state(record)
        structure_id = state[0]
        if structure_id in seen_structure_ids:
            raise AuthorityStateError(
                "authenticated dataset contains duplicate structure IDs"
            )
        if metadata_rows[index] is not record.metadata:
            raise AuthorityStateError(
                "metadata_rows must reference each record's exact metadata mapping"
            )
        if by_id.get(structure_id) is not record:
            raise AuthorityStateError(
                "_by_id must reference each exact authenticated record"
            )
        structure_ids.append(structure_id)
        seen_structure_ids.add(structure_id)
        record_digests.append(state_fingerprint(state))
    if tuple(by_id) != tuple(structure_ids):
        raise AuthorityStateError(
            "authenticated dataset index order differs from record order"
        )

    for name in (
        "dataset_info",
        "parent_group_methods",
        "parent_by_id",
        "cif_hashes",
        "input_hashes",
    ):
        if not isinstance(getattr(dataset, name, None), Mapping):
            raise AuthorityStateError("{} must be a mapping".format(name))
    parent_by_id = dataset.parent_by_id
    cif_hashes = dataset.cif_hashes
    input_hashes = dataset.input_hashes
    unknown_parent_ids = set(parent_by_id).difference(structure_ids)
    unknown_cif_ids = set(cif_hashes).difference(structure_ids)
    if unknown_parent_ids or unknown_cif_ids:
        raise AuthorityStateError(
            "parent/CIF mappings contain IDs outside the dataset generation"
        )
    if any(not isinstance(row, Mapping) for row in parent_by_id.values()):
        raise AuthorityStateError("parent_by_id values must be mappings")
    for structure_id in structure_ids:
        raw_parent = parent_by_id.get(structure_id)
        if raw_parent is None:
            raise AuthorityStateError(
                "authenticated release parent_by_id must cover every structure"
            )
        record = by_id[structure_id]
        parent_view = record.parent_groups
        if isinstance(parent_view, _ParentGroupView) and parent_view._row != raw_parent:
            raise AuthorityStateError(
                "record parent view is not bound to its exact release parent row"
            )
    for structure_id, digest in cif_hashes.items():
        if (
            not isinstance(digest, str)
            or _SHA256_RE.fullmatch(digest) is None
        ):
            raise AuthorityStateError(
                "cif_hashes must contain exact lowercase SHA-256 strings"
            )
        record = by_id[structure_id]
        if (
            record.cif_manifest is None
            or record.cif_manifest.get("sha256") != digest
        ):
            raise AuthorityStateError(
                "cif_hashes differs from the joined record manifest"
            )
    if any(
        not isinstance(key, str)
        or not key
        or key != key.strip()
        or not isinstance(value, str)
        or _SHA256_RE.fullmatch(value) is None
        for key, value in input_hashes.items()
    ):
        raise AuthorityStateError(
            "input_hashes must map exact nonblank names to lowercase SHA-256 strings"
        )
    if type(getattr(dataset, "cif_files_verified", None)) is not bool:
        raise AuthorityStateError("cif_files_verified must be an exact boolean")
    release_root = getattr(dataset, "release_root", None)
    if not isinstance(release_root, Path):
        raise AuthorityStateError("release_root must be a pathlib.Path")

    return (
        "coremof-dataset-generation/1",
        "{}.{}".format(type(dataset).__module__, type(dataset).__qualname__),
        release_root,
        tuple(record_digests),
        state_fingerprint(dataset.dataset_info),
        state_fingerprint(dataset.parent_group_methods),
        state_fingerprint(cif_hashes),
        state_fingerprint(input_hashes),
        dataset.cif_files_verified,
        getattr(dataset, "_authority_extra_state", None),
    )


def _dataset_generation_fingerprint(dataset: object) -> str:
    return state_fingerprint(_dataset_authority_payload(dataset))


def _dataset_identity_snapshot(dataset: object) -> Tuple[object, ...]:
    """Capture identities of every recursively immutable generation member."""

    records = dataset.records
    record_snapshot = []
    for record in records:
        parent_view = record.parent_groups
        if isinstance(parent_view, _ParentGroupView):
            parent_identity = (
                id(parent_view._row),
                id(parent_view._prefixes),
                parent_view._prefixes,
            )
        else:
            parent_identity = (
                id(parent_view),
                tuple(
                    (key, id(group), group.method, group.status, group.group_id, group.size)
                    for key, group in parent_view.items()
                ),
            )
        record_snapshot.append(
            (
                id(record),
                record.structure_id,
                id(record.metadata),
                id(parent_view),
                parent_identity,
                id(record.cif_manifest),
            )
        )
    bindings = getattr(dataset, "_authority_extra_bindings", None)
    if bindings is None:
        binding_snapshot = None
    elif isinstance(bindings, Mapping):
        binding_items = []
        for name, expected in bindings.items():
            marker = object()
            if (
                not isinstance(name, str)
                or not name
                or name != name.strip()
                or getattr(dataset, name, marker) is not expected
            ):
                raise AuthorityStateError(
                    "target dataset authority binding changed after creation"
                )
            binding_items.append((name, id(expected)))
        binding_snapshot = (id(bindings), tuple(binding_items))
    else:
        raise AuthorityStateError("dataset authority bindings must be a mapping or None")
    return (
        id(dataset.records),
        id(dataset.metadata_rows),
        id(dataset.dataset_info),
        id(dataset.parent_group_methods),
        id(dataset.parent_by_id),
        id(dataset.cif_hashes),
        id(dataset.input_hashes),
        id(dataset._by_id),
        id(dataset._authority_extra_state),
        binding_snapshot,
        dataset.release_root,
        dataset.cif_files_verified,
        tuple(record_snapshot),
    )


def _register_dataset_generation(
    dataset: object,
    *,
    kind: str,
    official_release_source: bool,
    base_dataset: object = None,
) -> str:
    """Seal one exact factory generation and return its state fingerprint."""

    if kind not in {"validated_release", "target_merged"}:
        raise AuthorityStateError("unknown dataset generation kind")
    if type(official_release_source) is not bool:
        raise TypeError("official_release_source must be a boolean")
    fingerprint = _dataset_seal_fingerprint(dataset, kind)
    _DATASET_GENERATIONS.register(
        dataset,
        fingerprint,
        {
            "kind": kind,
            "official_release_source": official_release_source,
            "base_dataset": base_dataset,
            "identity_snapshot": _dataset_identity_snapshot(dataset),
        },
    )
    object.__setattr__(dataset, "_authority_locked", True)
    return fingerprint


def _dataset_seal_fingerprint(dataset: object, kind: str) -> str:
    """Recompute the compact, source-bound fingerprint for one generation."""

    # A release generation is already bound to exact source bytes by the
    # loader's input SHA-256 map.  Hash those byte commitments plus every
    # independent public contract/CIF commitment and the ordered universe;
    # the recursively immutable object-identity snapshot below then prevents
    # any parsed member from being replaced while retaining that authority.
    if kind == "validated_release":
        return state_fingerprint(
            (
                "coremof-validated-release-generation/1",
                dataset.release_root,
                dataset.input_hashes,
                dataset.cif_hashes,
                dataset.dataset_info,
                dataset.parent_group_methods,
                tuple(record.structure_id for record in dataset.records),
                dataset.cif_files_verified,
            )
        )
    else:
        # TargetMergedDataset validates the exact transformation before this
        # call.  Its receipt commits source/feature files and target values;
        # the immutable target extra-state commits provenance and all joined
        # containers.  Avoid reserializing tens of thousands of wide rows.
        return state_fingerprint(
            (
                "coremof-target-merged-generation/1",
                dataset.input_hashes,
                dataset.cif_hashes,
                dataset.dataset_info,
                dataset.parent_group_methods,
                tuple(record.structure_id for record in dataset.records),
                dataset._authority_extra_state,
                dataset.cif_files_verified,
            )
        )
def _validate_dataset_generation_if_present(dataset: object) -> bool:
    """Revalidate a sealed generation; return whether it has release authority."""

    entry = _DATASET_GENERATIONS.entry(dataset)
    if entry is None:
        return False
    expected_fingerprint, context = entry
    if not isinstance(context, Mapping):  # pragma: no cover - private invariant
        raise AuthorityStateError("dataset authority context is malformed")
    if _dataset_identity_snapshot(dataset) != context.get("identity_snapshot"):
        raise AuthorityStateError(
            "validated dataset generation changed after authentication; rebuild it"
        )
    if _dataset_seal_fingerprint(dataset, context.get("kind")) != expected_fingerprint:
        raise AuthorityStateError(
            "validated dataset generation fingerprint changed; rebuild it"
        )
    if context.get("kind") == "target_merged" and context.get(
        "official_release_source"
    ) is True:
        if not _validate_dataset_generation_if_present(context.get("base_dataset")):
            raise AuthorityStateError(
                "target generation lost its authenticated base release"
            )
    return context.get("official_release_source") is True


def _release_receipt_state(
    dataset: object,
) -> Tuple[Optional[str], Optional[str], bool]:
    """Return closed parent/dataset status text and provisional state.

    Public ``release_status`` text is descriptive, not publication authority.
    The only accepted status on current public declaration surfaces is the
    exact staged-snapshot marker, while the validated RT/M2T integration
    contract is explicitly non-decisive and publication-authorized false.
    Consequently no caller-provided status pair can mint a non-provisional
    receipt; a future final state requires a separate authenticated package
    contract that does not yet exist.
    """

    authenticated = _validate_dataset_generation_if_present(dataset)
    parent_methods = getattr(dataset, "parent_group_methods", None)
    dataset_info = getattr(dataset, "dataset_info", None)
    if parent_methods is None:
        parent_methods = {}
    if dataset_info is None:
        dataset_info = {}
    if not isinstance(parent_methods, Mapping):
        raise TypeError("parent_group_methods must be a mapping when present")
    if not isinstance(dataset_info, Mapping):
        raise TypeError("dataset_info must be a mapping when present")
    _reject_retired_or_reserved_keys(
        parent_methods, "parent_group_methods"
    )
    _reject_retired_or_reserved_keys(dataset_info, "dataset_info")

    # Factory-loaded releases were scanned before their immutable generation
    # was sealed.  Generic protocol objects carry no such authentication, so
    # scan their metadata declarations at every receipt boundary as well.
    if not authenticated:
        metadata_rows = getattr(dataset, "metadata_rows", ())
        if not isinstance(metadata_rows, (list, tuple)):
            raise TypeError("metadata_rows must be an ordered list or tuple")
        _reject_retired_or_reserved_keys(metadata_rows, "metadata_rows")

    parent_status = parent_methods.get("release_status")
    dataset_status = dataset_info.get("release_status")
    for name, status in (
        ("parent release_status", parent_status),
        ("dataset release_status", dataset_status),
    ):
        if status is not None and (
            not isinstance(status, str)
            or not status
            or status != status.strip()
        ):
            raise ReleaseValidationError(
                "{} must be an exact nonblank string or None".format(name)
            )

    # There is intentionally no public FINAL-status shortcut.  The current
    # closed contract authenticates only staged/non-publication state.
    return parent_status, dataset_status, True


def _require_dataset_generation(dataset: object) -> Mapping[str, object]:
    """Require any package-factory dataset generation and return its context."""

    current = _DATASET_GENERATIONS.entry(dataset)
    if current is None:
        raise AuthorityStateError(
            "dataset was not produced by an authenticated package factory"
        )
    expected_fingerprint, context = current
    if not isinstance(context, Mapping):  # pragma: no cover - private invariant
        raise AuthorityStateError("dataset authority context is malformed")
    if _dataset_identity_snapshot(dataset) != context.get("identity_snapshot"):
        raise AuthorityStateError(
            "authenticated dataset generation changed; rebuild it"
        )
    if _dataset_seal_fingerprint(dataset, context.get("kind")) != expected_fingerprint:
        raise AuthorityStateError(
            "authenticated dataset generation fingerprint changed; rebuild it"
        )
    if context.get("kind") == "target_merged" and context.get(
        "official_release_source"
    ) is True:
        if not _validate_dataset_generation_if_present(context.get("base_dataset")):
            raise AuthorityStateError(
                "target generation lost its authenticated base release"
            )
    return context


def _registered_dataset_fingerprint(dataset: object) -> str:
    """Return a generation fingerprint only after its live state revalidates."""

    _require_dataset_generation(dataset)
    current = _DATASET_GENERATIONS.entry(dataset)
    if current is None:  # pragma: no cover - checked above
        raise AuthorityStateError("dataset generation is not registered")
    return current[0]


def _classified_authority_payload(value: object) -> Tuple[object, ...]:
    dataset = getattr(value, "dataset", None)
    if not _validate_dataset_generation_if_present(dataset):
        raise AuthorityStateError(
            "official checker views require an unchanged validated release generation"
        )
    checker_view = _exact_nonblank_text(
        getattr(value, "checker_view", None), "checker_view"
    )
    if checker_view not in CHECKER_PRESETS:
        raise AuthorityStateError(
            "official checker views require a canonical checker preset"
        )
    if getattr(value, "checker_view_official", None) is not True:
        raise AuthorityStateError("official checker-view flag changed after creation")
    if getattr(value, "checker_preset", None) != checker_view:
        raise AuthorityStateError("checker_preset differs from checker_view")
    records = getattr(value, "records", None)
    by_id = getattr(value, "_by_id", None)
    label_by_id = getattr(value, "label_by_id", None)
    steps = getattr(value, "_selection_steps", None)
    filters = getattr(value, "selection_filters", None)
    if type(records) is not tuple or type(steps) is not tuple:
        raise AuthorityStateError(
            "official checker records and selection steps must be tuples"
        )
    if not isinstance(by_id, Mapping) or not isinstance(label_by_id, Mapping):
        raise AuthorityStateError("official checker indexes must be mappings")
    if not isinstance(filters, Mapping):
        raise AuthorityStateError("official selection_filters must be a mapping")

    checker_names = CHECKER_PRESETS[checker_view]
    ids = []
    seen_ids = set()
    ordered_ids_digest = hashlib.sha256()
    for record in records:
        if not isinstance(record, ClassifiedRecord):
            raise AuthorityStateError(
                "official checker views require ClassifiedRecord rows"
            )
        structure_id = _exact_nonblank_text(
            record.structure_id, "classified structure_id"
        )
        if structure_id in seen_ids:
            raise AuthorityStateError(
                "official checker view contains duplicate structure IDs"
            )
        try:
            source_record = dataset._by_id[structure_id]
        except (AttributeError, KeyError) as error:
            raise AuthorityStateError(
                "official checker row is absent from its source generation"
            ) from error
        if record.record is not source_record:
            raise AuthorityStateError(
                "official checker row is not bound to the exact source record"
            )
        expected_statuses = {
            checker: source_record.metadata[column]
            for checker, column in CHECKER_COLUMNS.items()
        }
        if (
            record.checker_view != checker_view
            or not isinstance(record.checker_statuses, Mapping)
            or dict(record.checker_statuses) != expected_statuses
            or record.label
            != classify_checker_row(source_record.metadata, checker_names)
        ):
            raise AuthorityStateError(
                "official checker row differs from source-bound recomputation"
            )
        if by_id.get(structure_id) is not record:
            raise AuthorityStateError(
                "official checker _by_id does not reference the exact row"
            )
        if label_by_id.get(structure_id) != record.label:
            raise AuthorityStateError(
                "official checker label index differs from the exact row"
            )
        ids.append(structure_id)
        seen_ids.add(structure_id)
        ordered_ids_digest.update(structure_id.encode("utf-8"))
        ordered_ids_digest.update(b"\0")
    if tuple(by_id) != tuple(ids) or tuple(label_by_id) != tuple(ids):
        raise AuthorityStateError(
            "official checker index order differs from record order"
        )
    selected_digest = hashlib.sha256(
        "\n".join(sorted(ids)).encode("utf-8")
    ).hexdigest()
    if (
        filters.get("applied") is not bool(steps)
        or filters.get("steps") != steps
        or filters.get("selected_count") != len(records)
        or filters.get("selected_ids_sha256") != selected_digest
    ):
        raise AuthorityStateError(
            "official checker selection contract differs from its records"
        )
    return (
        "coremof-classified-generation/1",
        _registered_dataset_fingerprint(dataset),
        checker_view,
        True,
        ordered_ids_digest.hexdigest(),
        len(ids),
        steps,
        filters,
    )


def _classified_identity_snapshot(value: object) -> Tuple[object, ...]:
    records = value.records
    return (
        id(value.dataset),
        value.checker_view,
        value.checker_preset,
        value.checker_view_official,
        id(value.records),
        id(value._selection_steps),
        id(value._by_id),
        id(value.label_by_id),
        id(value.selection_filters),
        tuple(
            (
                id(record),
                id(record.record),
                record.structure_id,
                record.checker_view,
                record.label,
                id(record.checker_statuses),
            )
            for record in records
        ),
    )


def _register_classified_generation(value: object) -> str:
    payload = _classified_authority_payload(value)
    fingerprint = state_fingerprint(payload)
    _CLASSIFIED_GENERATIONS.register(
        value,
        fingerprint,
        {
            "dataset": value.dataset,
            "identity_snapshot": _classified_identity_snapshot(value),
        },
    )
    object.__setattr__(value, "_authority_locked", True)
    return fingerprint


def _require_classified_generation(value: object) -> Mapping[str, object]:
    current = _CLASSIFIED_GENERATIONS.entry(value)
    if current is None:
        raise AuthorityStateError(
            "checker view was not produced by the authenticated package factory"
        )
    expected_fingerprint, context = current
    if not isinstance(context, Mapping):  # pragma: no cover - private invariant
        raise AuthorityStateError("classified authority context is malformed")
    if _classified_identity_snapshot(value) != context.get("identity_snapshot"):
        raise AuthorityStateError(
            "official checker-view generation changed; rebuild it"
        )
    # Re-derive the exact statuses and labels from the still-authenticated
    # source generation.  An identity snapshot alone is insufficient because
    # ``object.__setattr__`` can bypass a frozen dataclass implementation.
    current_fingerprint = state_fingerprint(_classified_authority_payload(value))
    if current_fingerprint != expected_fingerprint:
        raise AuthorityStateError(
            "official checker-view fingerprint changed; rebuild it"
        )
    return context


def _validate_classified_generation_if_present(value: object) -> bool:
    entry = _CLASSIFIED_GENERATIONS.entry(value)
    if entry is None:
        return False
    _require_classified_generation(value)
    return True


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
        *,
        _validated_release_token: object = None,
    ) -> None:
        if type(cif_files_verified) is not bool:
            raise TypeError("cif_files_verified must be a boolean")
        object.__setattr__(self, "_authority_locked", False)
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
        # Kept as a non-authoritative compatibility attribute.  Authority is
        # identity-bound in ``_DATASET_GENERATIONS`` and never follows a token
        # copied onto another object.
        self._validated_release_token = None
        self._authority_extra_state = None
        self._authority_extra_bindings = None
        self._by_id = MappingProxyType(
            {record.structure_id: record for record in self.records}
        )

    def __setattr__(self, name: str, value: object) -> None:
        if getattr(self, "_authority_locked", False):
            raise AttributeError(
                "validated dataset generations are immutable; rebuild from the "
                "package factory"
            )
        object.__setattr__(self, name, value)

    def __copy__(self):
        reject_sealed_copy(self, "copied")
        raise TypeError("generic dataset copying is not part of the public contract")

    def __deepcopy__(self, memo):
        reject_sealed_copy(self, "deep-copied")
        raise TypeError("generic dataset copying is not part of the public contract")

    def __reduce_ex__(self, protocol):
        reject_sealed_copy(self, "pickled")
        raise TypeError("generic dataset pickling is not part of the public contract")

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
        against the same exact identifier set. A release declaring optional
        ``sm_*`` parent columns must also provide the exact, hash-bound
        StructureMatcher method contract and release-adapter receipt.
        """

        if type(verify_cif_files) is not bool:
            raise TypeError("verify_cif_files must be a boolean")

        root = Path(release_root).expanduser()
        if not root.is_dir():
            raise FileNotFoundError("release directory does not exist: {}".format(root))
        root = root.resolve()

        metadata_path = root / "metadata" / "metadata.csv"
        parents_path = root / "parent_groups" / "parent_groups.csv"
        methods_path = root / "parent_groups" / "parent_group_methods.json"
        info_path = root / "dataset_info.json"
        manifest_path = root / "manifests" / "cif_manifest.csv"

        metadata_snapshot = _capture_release_file(
            metadata_path, "release metadata"
        )
        parents_snapshot = _capture_release_file(parents_path, "parent groups")
        methods_snapshot = _capture_release_file(
            methods_path, "parent-group methods"
        )
        info_snapshot = _capture_release_file(info_path, "dataset information")

        metadata_fields, metadata_rows = _read_csv(
            metadata_snapshot, "release metadata"
        )
        _reject_retired_or_reserved_keys(metadata_rows, "metadata_rows")
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

        parent_fields, parent_rows = _read_csv(parents_snapshot, "parent groups")
        required_parent_columns = ["structure_id"]
        # The eight original criteria remain mandatory.  New audited
        # sensitivity criteria are optional so this package stays compatible
        # with already published releases while accepting an extended parent
        # table once it is available.
        for prefix in sorted(set(BASE_PUBLIC_PARENT_METHOD_PREFIXES.values())):
            required_parent_columns.extend(
                (
                    "{}_status".format(prefix),
                    "{}_group".format(prefix),
                    "{}_size".format(prefix),
                )
            )
        _require_columns(parent_fields, tuple(required_parent_columns), parents_path)
        parent_by_id, _ = _index_unique(parent_rows, parents_path)

        dataset_info = _read_json_object(info_snapshot, "dataset information")
        parent_methods = _read_json_object(
            methods_snapshot, "parent-group methods"
        )
        _reject_retired_or_reserved_keys(dataset_info, "dataset_info")
        _reject_retired_or_reserved_keys(
            parent_methods, "parent_group_methods"
        )

        metadata_ids = set(metadata_by_id)
        _require_exact_id_set(
            metadata_ids, set(parent_by_id), "metadata.csv", "parent_groups.csv"
        )

        manifest_by_id = None
        manifest_snapshot = None
        if manifest_path.is_file():
            manifest_snapshot = _capture_release_file(manifest_path, "CIF manifest")
            manifest_fields, manifest_rows = _read_csv(
                manifest_snapshot, "CIF manifest"
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
        parent_prefixes = _validate_parent_rows(
            parent_fields, parent_rows, metadata_by_id=metadata_by_id
        )
        structure_matcher_inputs = _validate_structure_matcher_evidence(
            parent_methods,
            dataset_info,
            parent_fields,
            parent_rows,
            root=root,
            parents_snapshot=parents_snapshot,
        )
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

        input_snapshots = {
            "metadata/metadata.csv": metadata_snapshot,
            "parent_groups/parent_groups.csv": parents_snapshot,
            "parent_groups/parent_group_methods.json": methods_snapshot,
            "dataset_info.json": info_snapshot,
        }
        if manifest_snapshot is not None:
            input_snapshots["manifests/cif_manifest.csv"] = manifest_snapshot
        input_snapshots.update(structure_matcher_inputs)
        input_hashes = MappingProxyType(
            {
                name: snapshot.sha256
                for name, snapshot in sorted(input_snapshots.items())
            }
        )
        immutable_parent_by_id = MappingProxyType(
            {
                structure_id: MappingProxyType(parent_by_id[structure_id])
                for structure_id in metadata_order
            }
        )

        result = cls(
            release_root=root,
            records=records,
            dataset_info=_deep_freeze(dataset_info),
            parent_group_methods=_deep_freeze(parent_methods),
            parent_by_id=immutable_parent_by_id,
            input_hashes=input_hashes,
            cif_files_verified=verify_cif_files,
        )
        _register_dataset_generation(
            result,
            kind="validated_release",
            official_release_source=True,
        )
        return result

    @property
    def dataset_version(self) -> str:
        if not isinstance(self.dataset_info, Mapping):
            raise TypeError("dataset_info must be a mapping")
        value = self.dataset_info.get("dataset_version")
        if (
            not isinstance(value, str)
            or not value
            or value != value.strip()
        ):
            raise ValueError("dataset_version must be an exact nonblank string")
        return value

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

        validated_release_source = _validate_dataset_generation_if_present(self)
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
        authenticated_official = view_is_official and validated_release_source
        classified = ClassifiedDataset(
            self,
            view_name,
            records,
            checker_view_official=authenticated_official,
            _official_checker_view_token=(
                _OFFICIAL_CHECKER_VIEW_TOKEN if authenticated_official else None
            ),
            _factory_token=(
                _CLASSIFIED_FACTORY_TOKEN if authenticated_official else None
            ),
        )
        return classified.filter(
            labels=labels,
            sources=sources,
            variants=variants,
            metals=metals,
            structure_ids=structure_ids,
        )

    def merge_targets(
        self,
        sources: Sequence[object],
        alias_registry: Optional[object] = None,
        feature_tables: Sequence[str] = (),
    ) -> object:
        """Return release rows joined to user targets and optional features.

        Importing the target layer lazily keeps ordinary release loading
        lightweight.  Earlier identifiers are accepted only when an explicit
        :class:`CoREMOF.targets.AliasRegistry` is supplied.
        """

        from .targets import merge_targets

        return merge_targets(
            self,
            sources,
            alias_registry=alias_registry,
            feature_tables=feature_tables,
        )


class ClassifiedDataset:
    """A checker-labelled, filterable dataset ready for parent-safe splitting."""

    def __init__(
        self,
        dataset: CoREMOFDataset,
        checker_view: str,
        records: Sequence[ClassifiedRecord],
        selection_steps: Sequence[Mapping[str, object]] = (),
        checker_view_official: bool = False,
        *,
        _official_checker_view_token: object = None,
        _factory_token: object = None,
    ) -> None:
        if not isinstance(checker_view, str) or not checker_view or checker_view != checker_view.strip():
            raise ValueError("checker_view must be an exact nonblank string")
        if type(checker_view_official) is not bool:
            raise TypeError("checker_view_official must be a boolean")
        if checker_view_official:
            if checker_view not in CHECKER_PRESETS:
                raise ValueError(
                    "checker_view_official=True requires a canonical official checker preset"
                )
            if (
                _factory_token is not _CLASSIFIED_FACTORY_TOKEN
                or
                _official_checker_view_token is not _OFFICIAL_CHECKER_VIEW_TOKEN
                or not _validate_dataset_generation_if_present(dataset)
            ):
                raise ValueError(
                    "checker_view_official=True requires the internal validated-release "
                    "checker-recomputation path"
                )
        object.__setattr__(self, "_authority_locked", False)
        self.dataset = dataset
        self.checker_view = checker_view
        self.checker_preset = checker_view
        self.checker_view_official = checker_view_official
        self._official_checker_view_token = (
            _OFFICIAL_CHECKER_VIEW_TOKEN if checker_view_official else None
        )
        self.records = tuple(records)
        for record in self.records:
            structure_id = getattr(record, "structure_id", None)
            if (
                not isinstance(structure_id, str)
                or not structure_id
                or structure_id != structure_id.strip()
            ):
                raise ValueError(
                    "classified records must expose exact nonblank string structure IDs"
                )
        self._selection_steps = tuple(
            _deep_freeze(dict(step)) for step in selection_steps
        )
        self._by_id = MappingProxyType(
            {record.structure_id: record for record in self.records}
        )
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
        if checker_view_official:
            _register_classified_generation(self)

    def __setattr__(self, name: str, value: object) -> None:
        if getattr(self, "_authority_locked", False):
            raise AttributeError(
                "official checker-view generations are immutable; rebuild them "
                "from the validated dataset"
            )
        object.__setattr__(self, name, value)

    def __copy__(self):
        reject_sealed_copy(self, "copied")
        raise TypeError(
            "generic checker-view copying is not part of the public contract"
        )

    def __deepcopy__(self, memo):
        reject_sealed_copy(self, "deep-copied")
        raise TypeError(
            "generic checker-view copying is not part of the public contract"
        )

    def __reduce_ex__(self, protocol):
        reject_sealed_copy(self, "pickled")
        raise TypeError(
            "generic checker-view pickling is not part of the public contract"
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

    def _validate_official_records(self) -> None:
        """Recompute every official label and bind its complete status evidence."""

        _require_classified_generation(self)

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

        authenticated_official = _validate_classified_generation_if_present(self)

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
            checker_view_official=authenticated_official,
            _official_checker_view_token=(
                _OFFICIAL_CHECKER_VIEW_TOKEN
                if authenticated_official
                else None
            ),
            _factory_token=(
                _CLASSIFIED_FACTORY_TOKEN if authenticated_official else None
            ),
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
        required_targets: Optional[Iterable[str]] = None,
        required_target_mode: str = "all",
        missing_parent: str = "singleton",
        leakage_guard: str = "auto",
        stratify_by: Sequence[str] = ("label",),
        **splitter_options: object
    ) -> object:
        """Create a parent-safe train/validation/test split.

        The splitter is imported only when requested, keeping metadata loading
        independent from optional training dependencies and future splitter
        implementations.  ``priority_main`` and ``auto`` are project-defined
        policies, not standard statistical terms.  Inspect their exact
        machine-readable contracts with
        :func:`CoREMOF.parents.parent_method_definition` and
        :func:`CoREMOF.parents.leakage_guard_definition`.
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
            required_targets=required_targets,
            required_target_mode=required_target_mode,
            **splitter_options
        )
        return splitter.train_valid_test_split(fractions=fractions)

    def split(self, **kwargs: object) -> object:
        """Alias for :meth:`train_valid_test_split`."""

        return self.train_valid_test_split(**kwargs)


def _is_authenticated_official_checker_view(value: object) -> bool:
    """Return true only for an internally recomputed canonical release view."""

    if not isinstance(value, ClassifiedDataset):
        return False
    try:
        _require_classified_generation(value)
    except (AttributeError, KeyError, TypeError, ValueError):
        return False
    return True


def _read_csv(
    snapshot: _ReleaseFileSnapshot, description: str
) -> Tuple[Tuple[str, ...], List[Dict[str, str]]]:
    path = snapshot.path
    try:
        text = snapshot.data.decode("utf-8-sig")
    except UnicodeDecodeError as error:
        raise ReleaseValidationError(
            "cannot read {} {}: {}".format(description, path, error)
        )
    with io.StringIO(text, newline="") as handle:
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


def _read_json_object(
    snapshot: _ReleaseFileSnapshot, description: str
) -> Dict[str, object]:
    path = snapshot.path
    try:
        value = json.loads(snapshot.data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
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
    if (
        not isinstance(version, str)
        or not version
        or version != version.strip()
    ):
        raise ReleaseValidationError(
            "dataset_info.json dataset_version must be an exact nonblank string"
        )
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
        if (
            not row["source_id"]
            or row["source_id"] != row["source_id"].strip()
        ):
            raise ReleaseValidationError(
                "{} source_id must be exact nonblank text".format(structure_id)
            )
        expected_cif = "cifs/{}.cif".format(structure_id)
        if row["cif_file"] != expected_cif:
            raise ReleaseValidationError(
                "{} must use CIF path {!r}, not {!r}".format(
                    structure_id, expected_cif, row["cif_file"]
                )
            )


def _parent_prefixes(fields: Sequence[str]) -> Tuple[str, ...]:
    if any("relaxed" in field.casefold() for field in fields):
        raise ReleaseValidationError(
            "historical relaxed StructureMatcher evidence may not be exposed "
            "in parent_groups.csv"
        )

    suffixes = ("_status", "_group", "_size")
    expected_suffixes = set(suffixes)
    prefixes = []
    observed = {}
    for field in fields:
        for suffix in suffixes:
            if not field.endswith(suffix):
                continue
            prefix = field[: -len(suffix)]
            if not prefix:
                raise ReleaseValidationError(
                    "parent_groups.csv contains an empty parent criterion prefix"
                )
            if prefix not in observed:
                observed[prefix] = set()
                prefixes.append(prefix)
            observed[prefix].add(suffix)
            break
    if not prefixes:
        raise ReleaseValidationError("parent_groups.csv contains no parent criteria")
    for prefix in prefixes:
        if observed[prefix] != expected_suffixes:
            missing = expected_suffixes.difference(observed[prefix])
            raise ReleaseValidationError(
                "parent criterion {!r} must have a complete status/group/size "
                "triad; missing {}".format(
                    prefix,
                    ", ".join(sorted(value.lstrip("_") for value in missing)),
                )
            )
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
    if any(
        not isinstance(method, str) or not isinstance(prefix, str)
        for method, prefix in declared.items()
    ):
        raise ReleaseValidationError(
            "csv_column_prefixes keys and values must all be strings"
        )
    expected_base = dict(BASE_PUBLIC_PARENT_METHOD_PREFIXES)
    supported = dict(PUBLIC_PARENT_METHOD_PREFIXES)
    if any(declared.get(method) != prefix for method, prefix in expected_base.items()):
        raise ReleaseValidationError(
            "parent_group_methods.json public method contract is missing or "
            "changes a required parent method"
        )
    unsupported = {
        method: prefix
        for method, prefix in declared.items()
        if supported.get(method) != prefix
    }
    if unsupported:
        raise ReleaseValidationError(
            "parent_group_methods.json contains unsupported parent methods: {}".format(
                ", ".join(sorted(unsupported))
            )
        )
    if set(declared.values()) != prefixes:
        raise ReleaseValidationError(
            "parent-group CSV prefixes do not match parent_group_methods.json"
        )
    expected_fields = {"structure_id"}
    for prefix in declared.values():
        expected_fields.update(
            {
                "{}_status".format(prefix),
                "{}_group".format(prefix),
                "{}_size".format(prefix),
            }
        )
    if set(fields) != expected_fields:
        unexpected = sorted(set(fields).difference(expected_fields))
        missing = sorted(expected_fields.difference(fields))
        detail = []
        if unexpected:
            detail.append("unexpected={}".format(", ".join(unexpected)))
        if missing:
            detail.append("missing={}".format(", ".join(missing)))
        raise ReleaseValidationError(
            "parent_groups.csv must contain exactly structure_id and the "
            "declared status/group/size triads ({})".format("; ".join(detail))
        )

    criteria = methods.get("criteria")
    declared_crystalnets = {
        method
        for method in _CRYSTALNETS_REFERENCE_METHODS
        if method in declared
    }
    for method, (prefix, _) in _CRYSTALNETS_REFERENCE_METHODS.items():
        if method not in declared:
            if isinstance(criteria, dict) and method in criteria:
                raise ReleaseValidationError(
                    "criteria.{} is declared without {}_* columns".format(
                        method, prefix
                    )
                )
            continue
        if not isinstance(criteria, dict):
            raise ReleaseValidationError(
                "criteria.{} requires a criteria object".format(method)
            )
        criterion = criteria.get(method)
        if not isinstance(criterion, dict):
            raise ReleaseValidationError(
                "criteria.{} must be an object".format(method)
            )
        expected_criterion = _CRYSTALNETS_REFERENCE_CONTRACTS[method]
        if not _strict_json_equal(criterion, expected_criterion):
            raise ReleaseValidationError(
                "criteria.{} must match the exact closed terminal reference "
                "contract, including its inputs, comparison, missing/error, "
                "release-state, result, non-decisive-use, and rebuild rules".format(
                    method
                )
            )

    method_definitions = methods.get("definitions")
    if declared_crystalnets:
        if not isinstance(method_definitions, dict):
            raise ReleaseValidationError(
                "parent_group_methods.definitions must contain the declared "
                "RT/M2T contracts"
            )
        method_integration = methods.get("crystalnets_reference_integration")
        if not _strict_json_equal(
            method_integration, dict(_CRYSTALNETS_REFERENCE_INTEGRATION)
        ):
            raise ReleaseValidationError(
                "parent_group_methods.crystalnets_reference_integration must "
                "match the exact closed staged non-decisive contract"
            )
    for method in _CRYSTALNETS_REFERENCE_METHODS:
        method_definition_present = (
            isinstance(method_definitions, dict) and method in method_definitions
        )
        if method in declared_crystalnets:
            if not _strict_json_equal(
                method_definitions.get(method),
                _CRYSTALNETS_REFERENCE_CONTRACTS[method],
            ):
                raise ReleaseValidationError(
                    "parent_group_methods.definitions.{} must exactly match the "
                    "trusted criterion contract".format(method)
                )
        elif method_definition_present:
            raise ReleaseValidationError(
                "parent_group_methods.definitions declares {} without matching "
                "parent-group columns".format(method)
            )

    info_definitions = info.get("definitions")
    parent_grouping = info.get("parent_grouping")
    info_contracts = (
        parent_grouping.get("criterion_contracts")
        if isinstance(parent_grouping, dict)
        else None
    )
    if declared_crystalnets:
        if not isinstance(info_definitions, dict):
            raise ReleaseValidationError(
                "dataset_info.definitions must contain the declared RT/M2T contracts"
            )
        if not isinstance(parent_grouping, dict) or not isinstance(info_contracts, dict):
            raise ReleaseValidationError(
                "dataset_info.parent_grouping.criterion_contracts must contain "
                "the declared RT/M2T contracts"
            )
        integration = parent_grouping.get("crystalnets_reference_integration")
        if not _strict_json_equal(
            integration, dict(_CRYSTALNETS_REFERENCE_INTEGRATION)
        ):
            raise ReleaseValidationError(
                "dataset_info.parent_grouping.crystalnets_reference_integration "
                "must match the exact closed staged non-decisive contract"
            )
    for method in _CRYSTALNETS_REFERENCE_METHODS:
        expected = _CRYSTALNETS_REFERENCE_CONTRACTS[method]
        definition_present = (
            isinstance(info_definitions, dict) and method in info_definitions
        )
        parent_contract_present = (
            isinstance(info_contracts, dict) and method in info_contracts
        )
        if method in declared_crystalnets:
            if not _strict_json_equal(info_definitions.get(method), expected):
                raise ReleaseValidationError(
                    "dataset_info.definitions.{} must exactly match the trusted "
                    "parent-method criterion contract".format(method)
                )
            if not _strict_json_equal(info_contracts.get(method), expected):
                raise ReleaseValidationError(
                    "dataset_info.parent_grouping.criterion_contracts.{} must "
                    "exactly match the trusted parent-method criterion contract".format(
                        method
                    )
                )
        elif definition_present or parent_contract_present:
            raise ReleaseValidationError(
                "dataset_info declares {} without matching parent-group columns".format(
                    method
                )
            )

    if "structure_matcher_strict" not in declared:
        if isinstance(criteria, dict) and "structure_matcher_strict" in criteria:
            raise ReleaseValidationError(
                "structure_matcher_strict criteria are declared without sm_* columns"
            )
        return

    if not isinstance(criteria, dict):
        raise ReleaseValidationError(
            "criteria.structure_matcher_strict requires a criteria object"
        )
    criterion = criteria.get("structure_matcher_strict")
    if not isinstance(criterion, dict):
        raise ReleaseValidationError(
            "criteria.structure_matcher_strict must be an object"
        )
    expected_keys = set(_STRUCTURE_MATCHER_METHOD_CONTRACT).union(
        {"evidence_receipt"}
    )
    if set(criterion) != expected_keys:
        raise ReleaseValidationError(
            "criteria.structure_matcher_strict does not match the exact "
            "optional-reference contract"
        )
    for key, expected in _STRUCTURE_MATCHER_METHOD_CONTRACT.items():
        if not _strict_json_equal(criterion.get(key), expected):
            raise ReleaseValidationError(
                "criteria.structure_matcher_strict {!r} must be {!r}".format(
                    key, expected
                )
            )
    evidence = criterion.get("evidence_receipt")
    if not isinstance(evidence, dict) or set(evidence) != {"file", "sha256"}:
        raise ReleaseValidationError(
            "criteria.structure_matcher_strict evidence_receipt must contain "
            "exactly file and sha256"
        )
    if evidence.get("file") != _STRUCTURE_MATCHER_RECEIPT_FILE:
        raise ReleaseValidationError(
            "criteria.structure_matcher_strict evidence receipt must use {!r}".format(
                _STRUCTURE_MATCHER_RECEIPT_FILE
            )
        )
    evidence_sha256 = evidence.get("sha256")
    if (
        not isinstance(evidence_sha256, str)
        or _SHA256_RE.fullmatch(evidence_sha256) is None
    ):
        raise ReleaseValidationError(
            "criteria.structure_matcher_strict evidence receipt SHA-256 must "
            "be 64 lowercase hexadecimal characters"
        )


def _validate_structure_matcher_evidence(
    methods: Mapping[str, object],
    info: Mapping[str, object],
    fields: Sequence[str],
    rows: Sequence[Mapping[str, str]],
    *,
    root: Path,
    parents_snapshot: _ReleaseFileSnapshot,
) -> Dict[str, _ReleaseFileSnapshot]:
    """Validate the hash-bound optional StructureMatcher release receipt."""

    if "sm" not in _parent_prefixes(fields):
        return {}
    criteria = methods["criteria"]
    criterion = criteria["structure_matcher_strict"]  # type: ignore[index]
    evidence = criterion["evidence_receipt"]  # type: ignore[index]
    receipt_path = root / _STRUCTURE_MATCHER_RECEIPT_FILE
    try:
        receipt_snapshot = _capture_release_file(
            receipt_path, "StructureMatcher release evidence receipt"
        )
    except FileNotFoundError:
        raise ReleaseValidationError(
            "structure_matcher_strict evidence receipt is missing: {}".format(
                receipt_path
            )
        )
    observed_receipt_hash = receipt_snapshot.sha256
    if observed_receipt_hash != evidence["sha256"]:  # type: ignore[index]
        raise ReleaseValidationError(
            "structure_matcher_strict evidence receipt SHA-256 does not match"
        )
    receipt = _read_json_object(
        receipt_snapshot, "StructureMatcher release evidence receipt"
    )
    expected_receipt_keys = {
        "schema_version",
        "status",
        "dataset_version",
        "method_id",
        "method_schema_version",
        "parent_groups_sha256",
        "strict_pair_ledger_sha256",
        "structure_count",
        "candidate_pair_count",
        "successful_pair_count",
        "unresolved_pair_count",
        "strict_direct_match_edge_count",
        "not_available_structure_count",
        "historical_relaxed_executed",
        "historical_relaxed_exposed",
    }
    if set(receipt) != expected_receipt_keys:
        raise ReleaseValidationError(
            "StructureMatcher evidence receipt does not match the exact "
            "release-adapter schema"
        )
    expected_values = {
        "schema_version": _STRUCTURE_MATCHER_RECEIPT_SCHEMA,
        "status": "PASS",
        "dataset_version": info.get("dataset_version"),
        "method_id": _STRUCTURE_MATCHER_METHOD_ID,
        "method_schema_version": _STRUCTURE_MATCHER_METHOD_SCHEMA,
        "historical_relaxed_executed": False,
        "historical_relaxed_exposed": False,
    }
    for key, expected in expected_values.items():
        if not _strict_json_equal(receipt.get(key), expected):
            raise ReleaseValidationError(
                "StructureMatcher evidence receipt {!r} must be {!r}".format(
                    key, expected
                )
            )
    for key in ("parent_groups_sha256", "strict_pair_ledger_sha256"):
        value = receipt.get(key)
        if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
            raise ReleaseValidationError(
                "StructureMatcher evidence receipt {} must be a full lowercase "
                "SHA-256".format(key)
            )
    if receipt["parent_groups_sha256"] != parents_snapshot.sha256:
        raise ReleaseValidationError(
            "StructureMatcher evidence receipt is not bound to parent_groups.csv"
        )
    expected_structure_count = info.get("structure_count")
    integer_fields = (
        "structure_count",
        "candidate_pair_count",
        "successful_pair_count",
        "unresolved_pair_count",
        "strict_direct_match_edge_count",
        "not_available_structure_count",
    )
    for key in integer_fields:
        value = receipt.get(key)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ReleaseValidationError(
                "StructureMatcher evidence receipt {} must be a non-negative "
                "integer".format(key)
            )
    if receipt["structure_count"] != expected_structure_count or len(rows) != expected_structure_count:
        raise ReleaseValidationError(
            "StructureMatcher evidence receipt structure_count does not match "
            "the release"
        )
    if (
        receipt["successful_pair_count"] + receipt["unresolved_pair_count"]
        != receipt["candidate_pair_count"]
    ):
        raise ReleaseValidationError(
            "StructureMatcher successful and unresolved pair counts do not "
            "sum to candidate_pair_count"
        )
    if receipt["strict_direct_match_edge_count"] > receipt["successful_pair_count"]:
        raise ReleaseValidationError(
            "StructureMatcher direct match edge count exceeds successful pairs"
        )
    observed_not_available = sum(row["sm_status"] == "NOT_AVAILABLE" for row in rows)
    if receipt["not_available_structure_count"] != observed_not_available:
        raise ReleaseValidationError(
            "StructureMatcher evidence receipt NOT_AVAILABLE count does not "
            "match parent_groups.csv"
        )
    if receipt["unresolved_pair_count"] and not observed_not_available:
        raise ReleaseValidationError(
            "unresolved StructureMatcher pairs require NOT_AVAILABLE public rows"
        )
    return {_STRUCTURE_MATCHER_RECEIPT_FILE: receipt_snapshot}


def _validate_parent_rows(
    fields: Sequence[str],
    rows: Sequence[Mapping[str, str]],
    *,
    metadata_by_id: Optional[Mapping[str, Mapping[str, str]]] = None,
) -> Tuple[str, ...]:
    prefixes = _parent_prefixes(fields)
    counts = {prefix: Counter() for prefix in prefixes}
    parsed_sizes = {}
    m2t_mofid_by_group: Dict[str, str] = {}

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
            if (
                prefix == "mofid2_crystalnets"
                and status in {"MATCHED", "UNMATCHED"}
            ):
                if metadata_by_id is None:
                    raise ReleaseValidationError(
                        "mofid_v2_crystalnets parent validation requires release metadata"
                    )
                metadata = metadata_by_id.get(structure_id)
                mofid_status = (
                    metadata.get("mofid_v2_status")
                    if isinstance(metadata, Mapping)
                    else None
                )
                if mofid_status not in M2T_ELIGIBLE_MOFID_V2_STATUSES:
                    raise ReleaseValidationError(
                        "{} has {} status {} but ineligible mofid_v2_status {!r}".format(
                            structure_id, prefix, status, mofid_status
                        )
                    )
                canonical_mofid = _canonical_complete_mofid_v2(
                    metadata.get("mofid_v2"), structure_id
                )
                previous_mofid = m2t_mofid_by_group.setdefault(
                    group, canonical_mofid
                )
                if previous_mofid != canonical_mofid:
                    raise ReleaseValidationError(
                        "{} group {} spans conflicting canonical mofid_v2 values".format(
                            prefix, group
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
            size_text = raw_size
            if (
                not size_text
                or not size_text.isascii()
                or not size_text.isdigit()
                or (len(size_text) > 1 and size_text.startswith("0"))
            ):
                raise ReleaseValidationError(
                    "parent_groups.csv:{} {} has invalid {} size {!r}".format(
                        row_number, structure_id, prefix, raw_size
                    )
                )
            size = int(size_text)
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
        if not isinstance(values, (list, tuple)):
            raise TypeError(
                "{} filter must be an exact string or ordered list/tuple of "
                "strings".format(description)
            )
        result = set(values)
    if any(
        not isinstance(value, str) or not value or value != value.strip()
        for value in result
    ):
        raise ValueError(
            "{} filter values must be exact nonblank strings".format(description)
        )
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
