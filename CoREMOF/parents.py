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

Project-defined terms used by this module are deliberately defined here before
their machine keys.  ``identity_union`` is the provisional source-ID/MOFid
transitive group: it forms connected components from exact database-namespaced
``source_database``/``source_id`` equality plus eligible complete MOFid-v2 or
MOFid-v1 equality, treats missing values as no edge, and does not use RAC5,
Zeo++, CrystalNets, CIF hashes, or StructureMatcher.  It is not proof of
identity and is excluded from ``main_union``.  ``main_union`` is the split-
leakage guard, not a parent or explanatory method and not proof of identity or
parentage: before filtering it forms transitive connected components over the
complete release from exact full lowercase 64-hex CIF SHA-256, database-
namespaced source siblings, and available RAC5, MOFid-v2, and MOFid-v1 edges.
Missing or malformed CIF hashes fail closed; missing optional evidence adds no
edge.  The strict StructureMatcher reference is the pinned strict pymatgen
protocol: symmetric forward-and-reverse direct edges are authoritative, while
its connected components are convenience views rather than duplicate proof.
Parser, timeout, OOM, matcher, asymmetric, or execution errors are
NOT_AVAILABLE rather than unmatched, and this reference enters neither
``priority_main`` nor ``main_union``.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from types import MappingProxyType
from typing import Dict, Iterable, List, Mapping, Optional, Set, Tuple


DIRECT_PARENT_METHODS = ("rac5", "mofid_v2", "mofid_v1")
OPTIONAL_CRYSTALNETS_PARENT_METHODS = (
    "rac5_crystalnets",
    "mofid_v2_crystalnets",
)
OPTIONAL_REFERENCE_PARENT_METHODS = OPTIONAL_CRYSTALNETS_PARENT_METHODS + (
    "structure_matcher_strict",
)
REFERENCE_PARENT_METHODS = (
    "rac5_zeo",
    "rac5_crystalnets",
    "mofid_v2_crystalnets",
    "structure_matcher_strict",
    "zeo",
    "source_id",
    "common_name",
    "identity_union",
    "none",
)
COMPUTED_PARENT_METHODS = ("priority_main", "main_union")
PARENT_METHODS = DIRECT_PARENT_METHODS + REFERENCE_PARENT_METHODS + COMPUTED_PARENT_METHODS
SELECTABLE_PARENT_METHODS = (
    ("priority_main",) + DIRECT_PARENT_METHODS + REFERENCE_PARENT_METHODS
)


_METHOD_KEYS = {
    "rac5": ("rac5", "rac", "rac5_group", "rac_group"),
    "mofid_v2": ("mofid_v2", "mofid2", "mofid_v2_group", "mofid2_group"),
    "mofid_v1": ("mofid_v1", "mofid1", "mofid_v1_group", "mofid1_group"),
    "rac5_zeo": ("rac5_zeo", "rac_zeo", "rac5_zeo_group", "rac_zeo_group"),
    "rac5_crystalnets": (
        "rac5_crystalnets",
        "rac_crystalnets",
        "rac5_crystalnets_group",
        "rac_crystalnets_group",
    ),
    "mofid_v2_crystalnets": (
        "mofid_v2_crystalnets",
        "mofid2_crystalnets",
        "mofid_v2_crystalnets_group",
        "mofid2_crystalnets_group",
    ),
    "structure_matcher_strict": (
        "structure_matcher_strict",
        "structure_matcher",
        "sm",
        "structure_matcher_strict_group",
        "structure_matcher_group",
        "sm_group",
    ),
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
    "rac5_crystalnets": "rac_crystalnets",
    "mofid_v2_crystalnets": "mofid2_crystalnets",
    "structure_matcher_strict": "sm",
    "zeo": "zeo",
    "source_id": "source",
    "common_name": "name",
    "identity_union": "identity",
}


_MISSING_TEXT = {"", "NA", "N/A", "NONE", "NULL", "NOT_AVAILABLE", "UNAVAILABLE"}


PARENT_METHOD_CONTRACT_VERSION = "coremof-parent-method/1.5"
LEAKAGE_GUARD_CHOICES = ("auto", "main_union", "parent_only")

_PARENT_METHOD_DISPLAY_NAMES = {
    "priority_main": "conflict-aware RAC5-first parent components",
    "rac5": "exact RAC5 groups",
    "mofid_v2": "exact MOFid-v2 groups",
    "mofid_v1": "exact MOFid-v1 groups",
    "rac5_zeo": "exact RAC5-plus-Zeo++ groups",
    "rac5_crystalnets": "exact RAC5-plus-CrystalNets groups",
    "mofid_v2_crystalnets": "exact MOFid-v2-plus-CrystalNets groups",
    "structure_matcher_strict": "strict StructureMatcher component view",
    "zeo": "exact selected-Zeo++ groups",
    "source_id": "exact namespaced source-ID groups",
    "common_name": "exact canonicalized common-name groups",
    "identity_union": "provisional source-ID/MOFid transitive groups",
    "none": "independent-structure singleton control",
}

_LEAKAGE_GUARD_DISPLAY_NAMES = {
    "auto": "parent-method-aware leakage-guard selector",
    "main_union": "full-release main leakage components",
    "parent_only": "selected-parent-only leakage blocks",
}


def _freeze_contract(value):
    """Return an immutable, recursively frozen contract value."""

    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _freeze_contract(item) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_contract(item) for item in value)
    return value


CANONICALIZED_IDENTIFIER_TEXT_DEFINITION = _freeze_contract(
    {
        "term": "canonicalized_identifier_text",
        "display_name": "canonicalized identifier text",
        "project_defined_term": True,
        "purpose": (
            "Make release-authorized current identifier text comparable by exact "
            "equality; this is text cleanup only, not crystal-structure "
            "normalization."
        ),
        "current_release_text_steps_in_order": (
            "convert_the_release_authorized_current_value_to_text",
            "collapse_each_Unicode_whitespace_run_to_one_ASCII_space",
            "trim_leading_and_trailing_whitespace",
            "reject_a_whole_field_missing_or_execution_placeholder",
            "apply_Unicode_NFKC",
            "apply_Unicode_default_casefold",
        ),
        "current_release_case_insensitive_whole_field_placeholders": (
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
        ),
        "source_key": {
            "input_fields": ("metadata.source_database", "metadata.source_id"),
            "algorithm": (
                "canonicalize_each_field_separately_then_compare_the_ordered_"
                "source_database_source_id_pair_by_exact_equality"
            ),
            "database_namespace_prevents_cross_database_ID_matches": True,
            "scope": (
                "freshly_recomputed_over_every_current_row_of_each_named_release; "
                "no_prior_component_or_earlier-version_source_key_is_imported"
            ),
        },
        "mofid_key": {
            "input_fields": ("metadata.mofid_v2", "metadata.mofid_v1"),
            "algorithm": (
                "canonicalize_the_complete_release_authorized_current_MOFid_string_"
                "then_compare_"
                "by_exact_equality"
            ),
            "fuzzy_or_partial_string_matching": False,
        },
        "missing_behavior": (
            "a_null_or_rejected_whole-field_placeholder_supplies_no_equality_edge; "
            "two missing values never match"
        ),
        "excluded_operations": (
            "no_atom_bond_occupancy_coordinate_or_unit_cell_change",
            "no_CIF_coordinate_canonicalization",
            "no_general_chemical_punctuation_rewrite",
            "no_fuzzy_matching",
            "no_similarity_tolerance",
        ),
        "performed_by": (
            "the_release_builder; CoREMOF-tools consumes the release-authorized "
            "group/status/size columns and does not recalculate these text keys; "
            "the resolver's compatibility source fallback has its own narrower "
            "algorithm in the main_union definition"
        ),
        "user_facing_definition": (
            "Canonical identifier text is text processing only: convert a non-null "
            "value to text, collapse each Unicode whitespace run to one ASCII space, "
            "trim leading and trailing whitespace, reject empty or declared placeholder "
            "values, apply Unicode NFKC, and then casefold. It does not change a CIF, "
            "coordinates, atoms, occupancies, bonding, chemistry, topology, or unit cell."
        ),
    }
)


PARENT_GROUP_TRIAD_DEFINITION = _freeze_contract(
    {
        "term": "release_authorized_parent_group_triad",
        "display_name": "release-authorized parent status/group/size triad",
        "project_defined_term": True,
        "purpose": (
            "Expose one release-builder criterion without making CoREMOF-tools "
            "recalculate the underlying descriptor or identifier relation."
        ),
        "fields": ("<prefix>_status", "<prefix>_group", "<prefix>_size"),
        "status_semantics": {
            "MATCHED": "criterion_available_and_observed_group_size_at_least_2",
            "UNMATCHED": "criterion_available_and_observed_group_size_exactly_1",
            "NOT_AVAILABLE": (
                "criterion_unavailable; release_assigns_a_structure_specific_size_1_"
                "group_only_for_table_shape_and_it_supplies_no_scientific_edge"
            ),
        },
        "validation": (
            "status_must_be_MATCHED_UNMATCHED_or_NOT_AVAILABLE; group_must_be_a_nonblank_"
            "string_and_MATCHED_or_UNMATCHED_group_must_not_be_a_missing_placeholder; "
            "size_must_be_a_non-boolean_positive_integer_or_canonical_ASCII-decimal-"
            "text_without_sign_whitespace_or_leading_zero_"
            "and_equal_to_the_number_of_release_rows_with_"
            "that_group; MATCHED_requires_size_at_least_2; UNMATCHED_and_NOT_AVAILABLE_"
            "require_size_1"
        ),
        "availability_rule": (
            "MATCHED_and_UNMATCHED_are_available; NOT_AVAILABLE_is_missing_evidence"
        ),
        "missing_rule": (
            "NOT_AVAILABLE_never_matches_another_NOT_AVAILABLE_row_and_is_projected_"
            "to_a_unique_SINGLETON_structure_id_or_explicit_exclusion"
        ),
        "release_authorized_meaning": (
            "the triad is declared by parent_groups/parent_group_methods.json, appears "
            "in parent_groups/parent_groups.csv, and passes loader validation"
        ),
    }
)


RAC5_NUMERIC_FINGERPRINT_DEFINITION = _freeze_contract(
    {
        "term": "exact_RAC5_fingerprint",
        "display_name": "exact 264-value depth-5 RAC fingerprint",
        "purpose": "Compare the complete release-authorized RAC5 vector.",
        "input": (
            "all_264_ordered_descriptor_names_from_"
            "parent_group_methods.criteria.rac5.ordered_descriptors"
        ),
        "value_conversion": (
            "parse_each_trimmed_value_as_IEEE754_binary64; reject_nonfinite_or_"
            "unparseable_values; map_negative_zero_to_positive_zero; serialize_with_"
            "Python_float.hex"
        ),
        "comparison": (
            "exact_ordered_token_equality_with_rtol_0_and_atol_0; no_scaling_"
            "feature_deletion_imputation_rounding_or_similarity_tolerance"
        ),
        "missing_behavior": "one_missing_or_invalid_value_makes_RAC5_unavailable",
        "performed_by": (
            "the_release_builder; the_package_reads_the_validated_rac_status_group_size_"
            "triad"
        ),
    }
)


ZEO_NUMERIC_FINGERPRINT_DEFINITION = _freeze_contract(
    {
        "term": "exact_selected_Zeo_fingerprint",
        "display_name": "exact selected N2/He Zeo++ fingerprint",
        "purpose": "Compare the release-selected intensive pore geometry.",
        "probe_radii_A": {"N2": 1.655, "He": 1.32},
        "numeric_fields": (
            "density_g_cm3",
            "largest_cavity_diameter_A",
            "pore_limiting_diameter_A",
            "largest_free_path_diameter_A",
            "n2_accessible_surface_area_m2_cm3",
            "n2_accessible_surface_area_m2_g",
            "n2_nonaccessible_surface_area_m2_cm3",
            "n2_nonaccessible_surface_area_m2_g",
            "n2_accessible_volume_cm3_g",
            "n2_accessible_volume_fraction",
            "n2_nonaccessible_volume_cm3_g",
            "n2_nonaccessible_volume_fraction",
            "he_void_fraction",
        ),
        "hard_gates": (
            "n2_channel_dimension_exact_integer_equality",
            "periodicity_available_true_for_both_rows",
            "structure_periodic_dimension_exact_integer_equality",
        ),
        "value_conversion": (
            "parse_each_trimmed_numeric_value_as_IEEE754_binary64; reject_nonfinite_or_"
            "unparseable_values; map_negative_zero_to_positive_zero; serialize_with_"
            "Python_float.hex"
        ),
        "comparison": (
            "exact_ordered_token_and_hard_gate_equality_with_rtol_0_and_atol_0; "
            "no_scaling_imputation_rounding_or_similarity_tolerance"
        ),
        "missing_behavior": (
            "one_missing_or_invalid_numeric_value_or_hard_gate_makes_Zeo_unavailable"
        ),
        "excluded_fields": (
            "unit_cell_extensive_surface_areas_A2",
            "unit_cell_extensive_volumes_A3",
            "open_metal_sites",
            "topology_labels",
            "framework_component_counts",
            "zero_probe_properties",
        ),
        "performed_by": (
            "the_release_builder; the_package_reads_the_validated_zeo_status_group_size_"
            "triad"
        ),
    }
)


CURRENT_CRYSTALNETS_FINGERPRINT_DEFINITION = _freeze_contract(
    {
        "term": "current_CrystalNets_scientific_fingerprint",
        "display_name": "current CrystalNets scientific fingerprint",
        "project_defined_term": True,
        "purpose": (
            "Compare complete successful current CrystalNets results without using "
            "runtime, file, or diagnostic metadata."
        ),
        "current_meaning": (
            "the_current-result_topology_evidence_authorized_by_the_loaded_release_"
            "parent-method_file; the_package_does_not_search_for_a_newer_runtime_result"
        ),
        "performed_by": (
            "the_release_builder; the_package_reads_only_the_validated_optional_"
            "status_group_size_triad"
        ),
        "availability_gate": (
            "execution_status_SUCCESS_topology_available_true_error_null_and_nonempty_"
            "complete_SingleNodes_and_AllNodes_subnets"
        ),
        "included_fields": (
            "network_dimension",
            "interpenetrated_subnet_count",
            "catenation_degree",
            "subnet_count",
            "single_node_net",
            "all_node_net",
            "single_all_agree",
            "for_each_subnet_single_all_agree",
            "for_each_subnet_SingleNodes_status_dimension_topology_key_topology_name_"
            "topological_genome",
            "for_each_subnet_AllNodes_status_dimension_topology_key_topology_name_"
            "topological_genome",
        ),
        "canonicalization": (
            "canonical_JSON_with_sorted_object_keys; sort_complete_subnet_projections_"
            "by_their_canonical_JSON; retain_duplicate_subnets; hash_with_SHA256"
        ),
        "field_validation": {
            "count_fields": (
                "interpenetrated_subnet_count_and_catenation_degree_are_integers_"
                "at_least_1_and_each_equals_the_observed_subnet_count"
            ),
            "network_dimension": "null_or_one_of_0_1_2_3",
            "top_level_net_and_agreement": (
                "single_node_net_all_node_net_and_single_all_agree_may_be_null_for_"
                "heterogeneous_subnets_and_that_null_is_retained_as_scientific_data"
            ),
            "subnet_nodes": (
                "each_SingleNodes_and_AllNodes_view_has_status_SUCCESS_dimension_0_to_3_"
                "and_nonblank_topology_key; topology_name_and_topological_genome_may_"
                "be_null"
            ),
            "subnet_indices": "unique_contiguous_1_through_subnet_count_before_sorting",
        },
        "missing_behavior": (
            "an_incomplete_or_unsuccessful_result_supplies_no_topology_group_edge"
        ),
        "excluded_fields": (
            "runtime_seconds",
            "CIF_paths_or_hashes",
            "errors_or_diagnostics",
            "software_or_method_boilerplate",
            "original_subnet_index_or_order",
        ),
        "not_a_topology_similarity_tolerance": True,
    }
)


STRICT_STRUCTURE_MATCHER_DEFINITION = _freeze_contract(
    {
        "term": "pymatgen_structure_matcher_strict_v2",
        "display_name": "strict symmetric pymatgen StructureMatcher evidence",
        "project_defined_term": True,
        "purpose": (
            "Provide optional coordinate-and-lattice duplicate evidence as an audited "
            "direct pair-edge ledger and a non-transitive component view."
        ),
        "software": {
            "python": "3.9",
            "pymatgen": "2024.2.8",
            "numpy": "1.26.4",
        },
        "parser": {
            "implementation": "direct_pinned_pymatgen_CifParser",
            "site_tolerance": 0.0001,
            "frac_tolerance": 0.0001,
            "operations": (
                "expand_declared_symmetry_operations",
                "merge_generated_sites_within_site_tolerance",
                "round_fractional_coordinates_near_one_third_or_two_thirds_within_"
                "frac_tolerance",
                "check_occupancy",
                "sort_the_periodic_Structure",
            ),
            "disorder_policy": "preserve_disorder_natively",
            "forbidden_repairs": (
                "manual_CIF_repair",
                "occupancy_selection",
                "atom_deletion",
                "chemistry_edit",
            ),
        },
        "candidate_pairs": {
            "blocking_key": (
                "ElementComparator_fractional_composition_hash_of_the_parsed_structure"
            ),
            "enumeration": "all_pairs_within_each_equal_block_are_attempted",
            "not_exclusionary_blocks": (
                "routine_formula",
                "DOI",
                "source_database_or_source_id",
                "RAC5",
                "topology",
                "site_count",
                "cell_volume",
            ),
        },
        "matcher_settings": {
            "comparator": "ElementComparator",
            "ltol": 0.001,
            "stol": 0.001,
            "angle_tol": 0.01,
            "primitive_cell": True,
            "scale": False,
            "attempt_supercell": True,
            "allow_subset": False,
            "supercell_size": "num_sites",
            "ignored_species": (),
            "fit_symmetric_argument": True,
        },
        "direct_edge_rule": (
            "record_forward_and_reverse_fit_calls_with_symmetric_true; add_one_"
            "undirected_edge_only_when_both_directional_fits_succeed"
        ),
        "displacement_diagnostics": {
            "reported_values": (
                "forward_normalized_rms",
                "forward_normalized_max_displacement",
                "reverse_normalized_rms",
                "reverse_normalized_max_displacement",
            ),
            "normalization": (
                "periodic_site_displacement_divided_by_(V/Nsites)^(1/3)_for_the_"
                "corresponding_direction"
            ),
            "units": "dimensionless",
            "is_angstrom_RMSD": False,
            "uses_charnley_rmsd_Kabsch_package": False,
        },
        "failure_behavior": {
            "directional_disagreement": "ASYMMETRIC_RESULT_no_edge_and_unavailable",
            "parser_timeout_OOM_or_matcher_error": "NOT_AVAILABLE_not_unmatched",
            "incomplete_component": (
                "project_every_touched_structure_to_a_structure_specific_"
                "NOT_AVAILABLE_singleton"
            ),
        },
        "component_semantics": (
            "SM-_groups_are_connected_components_of_direct_edges_for_convenience; "
            "tolerance_matching_is_not_transitive_so_direct_edges_are_authoritative_"
            "and_component_degree_edge_count_possible_edge_count_completeness_and_"
            "clique_status_must_be_retained"
        ),
        "scope": (
            "optional_reference_only; excluded_from_priority_main_and_main_union"
        ),
        "excluded_method": (
            "the_historical_ltol_0.2_stol_0.3_angle_tol_5_scale_true_profile_is_"
            "documentation_only_and_is_neither_executed_nor_exposed"
        ),
    }
)


PRIORITY_MAIN_DEFINITION = _freeze_contract(
    {
        "schema_version": PARENT_METHOD_CONTRACT_VERSION,
        "identifier": "priority_main",
        "display_name": _PARENT_METHOD_DISPLAY_NAMES["priority_main"],
        "project_defined_identifier": True,
        "role": "explanatory_parent_resolution",
        "purpose": (
            "Choose one conflict-aware explanatory parent component for each "
            "structure from the three approved parent criteria. The word priority "
            "describes parent-evidence precedence; it does not rank, schedule, or "
            "recalculate failed scientific features."
        ),
        "summary": (
            "Conflict-aware hierarchy over release-authorized RAC5, MOFid v2, "
            "and MOFid v1 parent groups; it is not a row-by-row first-nonmissing "
            "fallback and it is separate from the leakage guard."
        ),
        "user_facing_definition": (
            "priority_main is the conflict-aware explanatory hierarchy over the full "
            "release: available exact RAC5 groups seed components, then exact MOFid-v2 "
            "groups and exact MOFid-v1 groups are processed in that precedence. A lower "
            "group touching zero stronger components creates a new component; one touching "
            "exactly one attaches only its unresolved rows; one touching two or more records "
            "PARENT_METHOD_CONFLICT, never merges the stronger components, and leaves "
            "lower-only rows unresolved. Remaining missing evidence becomes a structure-"
            "specific singleton or is explicitly excluded. "
            "It does not use Zeo++, CrystalNets, source ID, common name, CIF hash, "
            "identity_union, or StructureMatcher evidence."
        ),
        "not_feature_recalculation_queue": True,
        "input_fields": {
            "rac5": (
                "parent_groups.rac_status",
                "parent_groups.rac_group",
                "parent_groups.rac_size",
            ),
            "mofid_v2": (
                "parent_groups.mofid2_status",
                "parent_groups.mofid2_group",
                "parent_groups.mofid2_size",
            ),
            "mofid_v1": (
                "parent_groups.mofid1_status",
                "parent_groups.mofid1_group",
                "parent_groups.mofid1_size",
            ),
        },
        "parent_group_triad": PARENT_GROUP_TRIAD_DEFINITION,
        "resolution_scope": "complete_release_before_optional_subset",
        "priority_order": ("rac5", "mofid_v2", "mofid_v1"),
        "algorithm": (
            "anchor every available RAC5 group; process MOFid-v2 groups and then "
            "MOFid-v1 groups against already stronger components; attach unresolved "
            "members only when the lower group touches at most one stronger component"
        ),
        "availability_rule": (
            "use_a_nonmissing_release_authorized_group_and_treat_explicit_"
            "NOT_AVAILABLE_as_missing"
        ),
        "rules": (
            {
                "step": 1,
                "criterion": "rac5",
                "action": "anchor_each_available_rac5_group_as_a_component",
            },
            {
                "step": 2,
                "criterion": "mofid_v2",
                "action": "process_each_available_group_against_stronger_components",
            },
            {
                "step": 3,
                "criterion": "mofid_v1",
                "action": "process_each_available_group_against_stronger_components",
            },
        ),
        "lower_group_resolution": {
            "zero_stronger_components": "create_one_new_component_for_the_group",
            "one_stronger_component": "attach_unresolved_members_to_that_component",
            "multiple_stronger_components": (
                "keep_stronger_components_separate_record_parent_conflict_and_leave_"
                "unresolved_members_unassigned"
            ),
            "already_anchored_members": "remain_in_their_stronger_components",
            "conflict_unresolved_members_are_reconsidered_by_lower_steps": False,
        },
        "conflict": {
            "diagnostic": "PARENT_METHOD_CONFLICT",
            "stronger_components_are_merged": False,
            "unresolved_with_parent_only": "excluded",
            "unresolved_with_main_union": (
                "retained_only_if_main_union_supplies_a_leakage_block_and_marked_"
                "with_the_diagnostic"
            ),
            "ledger_is_retained_even_without_unresolved_members": True,
        },
        "conflict_behavior": (
            "never_merge_multiple_stronger_components_record_PARENT_METHOD_CONFLICT_"
            "and_leave_lower_only_members_unresolved"
        ),
        "missing_evidence": {
            "singleton": "assign_unique_SINGLETON_structure_id_group",
            "exclude": "exclude_as_MISSING_PARENT_EVIDENCE",
            "common_nulls_are_grouped": False,
        },
        "missing_behavior": (
            "unique_per_structure_singleton_by_default_or_explicit_"
            "MISSING_PARENT_EVIDENCE_exclusion"
        ),
        "output_group_ids": {
            "rac5_anchor": "RAC5:<release_authorized_current_rac5_group>",
            "mofid_v2_component": "MOFID_V2:<release_authorized_current_mofid_v2_group>",
            "mofid_v1_component": "MOFID_V1:<release_authorized_current_mofid_v1_group>",
            "missing_singleton": "SINGLETON:<structure_id>",
            "attached_member_rule": "keep_the_stronger_component_group_id",
            "published_group_text_is_preserved": True,
            "published_group_text_is_preserved_compatibility_note": (
                "the_key_name_is_retained_for_API_compatibility; it_means_the_group_"
                "text_stored_in_the_loaded_release_table_is_used_verbatim_even_when_"
                "that_release_is_a_non-published_stage-only_candidate"
            ),
        },
        "available_evidence_only": True,
        "excluded_inputs": (
            "rac5_zeo",
            "rac5_crystalnets",
            "mofid_v2_crystalnets",
            "structure_matcher_strict",
            "zeo",
            "source_id",
            "common_name",
            "identity_union",
            "cif_sha256",
        ),
        "recommended_leakage_guard": "main_union",
        "relation_to_other_terms": (
            "explains_one_parent_component_per_resolved_row; main_union is the "
            "separate broader partition-leakage graph and auto selects that graph"
        ),
    }
)


MAIN_UNION_DEFINITION = _freeze_contract(
    {
        "schema_version": PARENT_METHOD_CONTRACT_VERSION,
        "identifier": "main_union",
        "display_name": _LEAKAGE_GUARD_DISPLAY_NAMES["main_union"],
        "project_defined_identifier": True,
        "role": "leakage_guard",
        "purpose": (
            "Build indivisible train-validation-test blocks that prevent a path "
            "through any approved exact relation from crossing partitions."
        ),
        "summary": (
            "Full-release transitive union used to keep related structures in one "
            "split block; it is deliberately broader than priority_main and is not "
            "the explanatory parent assignment."
        ),
        "user_facing_definition": (
            "main_union is the split-leakage guard, not a parent or explanatory method "
            "and not proof of identity or parentage. Before any label, source, variant, "
            "metal, ID, or target filter, it takes transitive connected components over "
            "the complete release from exact full CIF SHA-256, database-namespaced source "
            "siblings, and available release-authorized RAC5, MOFid-v2, and MOFid-v1 "
            "edges. A complete lowercase 64-hex full CIF SHA-256 is mandatory for every "
            "row, and missing or malformed input fails closed; missing optional evidence "
            "adds no edge."
        ),
        "input_fields": {
            "cif_identity": "manifests/cif_manifest.csv.sha256",
            "source_siblings": (
                "parent_groups.source_status",
                "parent_groups.source_group",
                "parent_groups.source_size",
                "fallback_metadata.source_database",
                "fallback_metadata.source_id",
            ),
            "rac5": (
                "parent_groups.rac_status",
                "parent_groups.rac_group",
                "parent_groups.rac_size",
            ),
            "mofid_v2": (
                "parent_groups.mofid2_status",
                "parent_groups.mofid2_group",
                "parent_groups.mofid2_size",
            ),
            "mofid_v1": (
                "parent_groups.mofid1_status",
                "parent_groups.mofid1_group",
                "parent_groups.mofid1_size",
            ),
        },
        "parent_group_triad": PARENT_GROUP_TRIAD_DEFINITION,
        "normal_release_source_rule": (
            "a_loaded_release_must_supply_the_source_status_group_size_triad; the_"
            "resolver_uses_its_validated_group_and_does_not_recanonicalize_metadata"
        ),
        "compatibility_source_fallback": {
            "scope": (
                "only_a_manually_constructed_dataset_object_with_no_source_criterion_"
                "fields_for_that_structure; a_loaded_release_never_uses_this_path"
            ),
            "source_database_steps": (
                "require_exact_nonblank_string",
                "apply_upper",
            ),
            "source_id_steps": (
                "require_exact_nonblank_string",
                "apply_casefold",
                "split_on_Unicode_whitespace_and_join_with_one_ASCII_space",
            ),
            "database_namespace_retained": True,
            "Unicode_NFKC_applied": False,
            "nonblank_placeholder_rejection_applied": False,
            "explicit_NOT_AVAILABLE_behavior": "do_not_fallback_and_add_no_source_edge",
        },
        "explanatory_parent_method": False,
        "universe": "complete_unfiltered_release",
        "construction_order": "construct_before_all_label_identity_and_target_filters",
        "edge_sources": (
            {
                "criterion": "cif_sha256",
                "relation": "exact_full_lowercase_sha256_equality",
                "required_for_every_structure": True,
                "missing_action": "fail_closed",
            },
            {
                "criterion": "source_id",
                "relation": "same_database_namespaced_source_sibling_group",
                "required_for_every_structure": False,
                "fallback": (
                    "the_exact_compatibility_source_fallback_above_only_when_the_"
                    "source_criterion_is_absent; never_when_it_is_NOT_AVAILABLE"
                ),
            },
            {
                "criterion": "rac5",
                "relation": "same_available_release_authorized_group",
                "required_for_every_structure": False,
            },
            {
                "criterion": "mofid_v2",
                "relation": "same_available_release_authorized_group",
                "required_for_every_structure": False,
            },
            {
                "criterion": "mofid_v1",
                "relation": "same_available_release_authorized_group",
                "required_for_every_structure": False,
            },
        ),
        "graph_rule": "connected_components_of_the_transitive_union_of_all_edges",
        "algorithm": "connected_components_of_the_transitive_union_of_all_five_exact_edge_types",
        "priority_or_conflict_resolution": "none_all_listed_edges_are_union_edges",
        "conflict_behavior": "none_all_available_listed_edges_are_unioned_with_equal_status",
        "missing_behavior": (
            "missing_full_CIF_SHA256_fails_closed; each_missing_optional_relation_adds_"
            "no_edge; common_nulls_never_match"
        ),
        "filtered_rows_can_bridge_selected_rows": True,
        "missing_and_conflict_behavior": {
            "missing_cif_sha256": "fail_closed_before_splitting",
            "missing_optional_relation": "add_no_edge_for_that_relation",
            "priority_conflicts": (
                "irrelevant_to_union_construction_all_available_listed_edges_are_added"
            ),
            "common_nulls_are_grouped": False,
        },
        "group_name": "MAIN_plus_deterministic_sha256_of_sorted_component_members",
        "output_group_id": {
            "grammar": "MAIN-<lowercase_sha256_prefix>",
            "member_order": "ascending_structure_id",
            "member_serialization": "U+0000_NUL_joined_without_trailing_separator",
            "text_encoding": "UTF-8",
            "digest": "SHA-256",
            "hex_case": "lowercase",
            "initial_prefix_length": 16,
            "collision_rule": (
                "extend_the_digest_prefix_one_hex_character_at_a_time_until_unique"
            ),
        },
        "excluded_inputs": (
            "rac5_zeo",
            "rac5_crystalnets",
            "mofid_v2_crystalnets",
            "structure_matcher_strict",
            "zeo",
            "common_name",
            "identity_union",
        ),
        "relation_to_priority_main": (
            "priority_main_explains_parentage_main_union_only_constrains_partitioning"
        ),
        "relation_to_other_terms": (
            "broader_partition_guard_selected_by_auto_for_priority_main; not_an_"
            "explanatory_parent_method"
        ),
        "release_source_key_canonicalization": CANONICALIZED_IDENTIFIER_TEXT_DEFINITION,
    }
)


PARENT_ONLY_DEFINITION = _freeze_contract(
    {
        "schema_version": PARENT_METHOD_CONTRACT_VERSION,
        "identifier": "parent_only",
        "display_name": _LEAKAGE_GUARD_DISPLAY_NAMES["parent_only"],
        "project_defined_identifier": True,
        "role": "leakage_guard",
        "purpose": "Use no leakage relation beyond the selected parent method.",
        "summary": (
            "Use the selected explanatory parent groups themselves as split blocks; "
            "do not add the broader main_union relations."
        ),
        "user_facing_definition": (
            "parent_only uses only each group from the selected explanatory parent "
            "relation as one split block; it adds no relation from another method."
        ),
        "input_fields": "resolved_selected_parent_method.group_by_id",
        "algorithm": "use_each_resolved_selected_parent_group_as_one_split_block",
        "block_source": "resolved_selected_parent_method",
        "full_release_transitive_union_applied": False,
        "priority_main_conflict_members_without_a_parent_group": "excluded",
        "conflict_behavior": (
            "an_unresolved_priority_main_conflict_has_no_parent_block_and_is_excluded"
        ),
        "missing_behavior": (
            "follow_the_selected_parent_methods_singleton_or_exclude_policy"
        ),
        "excluded_inputs": (
            "all_relations_not_already_present_in_the_selected_parent_method",
        ),
        "relation_to_other_terms": (
            "narrow_sensitivity_guard; unlike_main_union_it_adds_no_cross_method_edges"
        ),
    }
)


AUTO_LEAKAGE_GUARD_DEFINITION = _freeze_contract(
    {
        "schema_version": PARENT_METHOD_CONTRACT_VERSION,
        "identifier": "auto",
        "display_name": _LEAKAGE_GUARD_DISPLAY_NAMES["auto"],
        "project_defined_identifier": True,
        "role": "leakage_guard_selector",
        "purpose": (
            "Select the package policy guard from the explanatory parent method "
            "before split blocks are constructed."
        ),
        "summary": (
            "Resolve to main_union for priority_main; resolve to parent_only for "
            "every other selectable explanatory parent method."
        ),
        "resolution_rules": (
            {
                "when_parent_method": "priority_main",
                "resolved_guard": "main_union",
            },
            {
                "when_parent_method": "any_other_explanatory_parent_method",
                "resolved_guard": "parent_only",
            },
        ),
        "resolution_time": "splitter_construction_before_parent_resolution",
        "input_fields": ("requested_parent_method",),
        "algorithm": (
            "if_parent_method_is_priority_main_select_main_union_else_select_parent_only"
        ),
        "missing_or_unknown_parent_method_behavior": "fail_validation",
        "missing_behavior": "a_missing_or_unknown_parent_method_fails_validation",
        "conflict_behavior": "does_not_resolve_parent_evidence_conflicts",
        "excluded_inputs": "does_not_read_scientific_evidence_or_construct_edges",
        "relation_to_other_terms": (
            "selector_only; the resolved main_union or parent_only definition controls blocks"
        ),
        "receipt_policy": (
            "record_auto_as_requested_guard_and_record_the_concrete_guard_separately"
        ),
        "does_not_construct_groups_itself": True,
    }
)


_PUBLISHED_METHOD_DETAILS = {
    "rac5": {
        "purpose": "Use exact depth-5 RAC fingerprint groups as the explanation.",
        "package_input_fields": (
            "parent_groups.rac_status",
            "parent_groups.rac_group",
            "parent_groups.rac_size",
        ),
        "release_relation": (
            "exact equality of all 264 ordered finite depth-5 RAC descriptors; "
            "IEEE-754 binary64 values are compared through float.hex after mapping "
            "negative zero to positive zero; rtol=0 and atol=0"
        ),
        "rac5_fingerprint": RAC5_NUMERIC_FINGERPRINT_DEFINITION,
        "excluded_inputs": "all_non_RAC5_evidence",
        "relation_to_default_methods": (
            "strongest anchor inside priority_main and one edge source inside main_union"
        ),
    },
    "mofid_v2": {
        "purpose": (
            "Use exact release-authorized current MOFid-v2 identifier groups as "
            "the explanation."
        ),
        "package_input_fields": (
            "parent_groups.mofid2_status",
            "parent_groups.mofid2_group",
            "parent_groups.mofid2_size",
        ),
        "release_relation": "exact equality of complete canonicalized MOFid-v2 text",
        "canonicalized_identifier_text": CANONICALIZED_IDENTIFIER_TEXT_DEFINITION,
        "excluded_inputs": "all_non_MOFid-v2_evidence",
        "relation_to_default_methods": (
            "second evidence level inside priority_main and one edge source inside main_union"
        ),
    },
    "mofid_v1": {
        "purpose": (
            "Use exact release-authorized current MOFid-v1 identifier groups as "
            "the explanation."
        ),
        "package_input_fields": (
            "parent_groups.mofid1_status",
            "parent_groups.mofid1_group",
            "parent_groups.mofid1_size",
        ),
        "release_relation": "exact equality of complete canonicalized MOFid-v1 text",
        "canonicalized_identifier_text": CANONICALIZED_IDENTIFIER_TEXT_DEFINITION,
        "excluded_inputs": "all_non_MOFid-v1_evidence",
        "relation_to_default_methods": (
            "last evidence level inside priority_main and one edge source inside main_union"
        ),
    },
    "rac5_zeo": {
        "purpose": "Inspect exact agreement of both RAC5 and selected Zeo++ fingerprints.",
        "package_input_fields": (
            "parent_groups.rac_zeo_status",
            "parent_groups.rac_zeo_group",
            "parent_groups.rac_zeo_size",
        ),
        "release_relation": (
            "exact equality of all 264 RAC5 descriptors and the complete authorized "
            "13-field intensive N2/He Zeo++ fingerprint, with equal N2 channel "
            "dimension and bonded-framework periodic dimension; both numeric vectors "
            "use binary64 float.hex after mapping negative zero to positive zero"
        ),
        "rac5_fingerprint": RAC5_NUMERIC_FINGERPRINT_DEFINITION,
        "zeo_fingerprint": ZEO_NUMERIC_FINGERPRINT_DEFINITION,
        "excluded_inputs": (
            "unit_cell_extensive_areas_or_volumes",
            "open_metal_sites",
            "topology",
            "framework_component_counts",
        ),
        "relation_to_default_methods": (
            "reference_sensitivity_only; excluded_from_priority_main_and_main_union"
        ),
    },
    "zeo": {
        "purpose": "Inspect exact agreement of the selected pore-geometry fingerprint.",
        "package_input_fields": (
            "parent_groups.zeo_status",
            "parent_groups.zeo_group",
            "parent_groups.zeo_size",
        ),
        "release_relation": (
            "exact equality of 13 finite intensive N2/He Zeo++ fields calculated "
            "with N2 radius 1.655 A and He radius 1.32 A, plus equal N2 channel "
            "dimension and bonded-framework periodic dimension; binary64 values use "
            "float.hex after mapping negative zero to positive zero with rtol=atol=0"
        ),
        "zeo_fingerprint": ZEO_NUMERIC_FINGERPRINT_DEFINITION,
        "excluded_inputs": (
            "RAC5",
            "unit_cell_extensive_areas_or_volumes",
            "open_metal_sites",
            "topology",
        ),
        "relation_to_default_methods": (
            "reference_sensitivity_only; excluded_from_priority_main_and_main_union"
        ),
    },
    "source_id": {
        "purpose": "Inspect source-database sibling records with the same source identifier.",
        "package_input_fields": (
            "parent_groups.source_status",
            "parent_groups.source_group",
            "parent_groups.source_size",
        ),
        "release_relation": (
            "exact equality of the ordered canonicalized source_database and source_id pair"
        ),
        "canonicalized_identifier_text": CANONICALIZED_IDENTIFIER_TEXT_DEFINITION,
        "excluded_inputs": "all_structural_and_descriptor_evidence",
        "relation_to_default_methods": (
            "reference explanation excluded from priority_main; source siblings are a "
            "separate edge source inside main_union"
        ),
    },
    "common_name": {
        "purpose": (
            "Inspect records with the same canonicalized current common name in "
            "the loaded release."
        ),
        "package_input_fields": (
            "parent_groups.name_status",
            "parent_groups.name_group",
            "parent_groups.name_size",
        ),
        "release_relation": (
            "exact equality after whitespace collapse, Unicode NFKC, and casefold"
        ),
        "canonicalized_identifier_text": CANONICALIZED_IDENTIFIER_TEXT_DEFINITION,
        "excluded_inputs": "all_structural_and_descriptor_evidence",
        "relation_to_default_methods": (
            "sparse_nonunique_reference_only; excluded_from_priority_main_and_main_union"
        ),
    },
    "rac5_crystalnets": {
        "purpose": "Inspect exact RAC5 plus current CrystalNets topology agreement.",
        "package_input_fields": (
            "parent_groups.rac_crystalnets_status",
            "parent_groups.rac_crystalnets_group",
            "parent_groups.rac_crystalnets_size",
        ),
        "release_relation": (
            "exact equality of all 264 finite RAC5 descriptors plus a complete "
            "successful current CrystalNets scientific fingerprint"
        ),
        "rac5_fingerprint": RAC5_NUMERIC_FINGERPRINT_DEFINITION,
        "crystalnets_fingerprint": CURRENT_CRYSTALNETS_FINGERPRINT_DEFINITION,
        "group_prefix_meaning": (
            "RT- identifies this RAC5-plus-CrystalNets criterion; the following "
            "compact digest labels a group and is not a CrystalNets topology name"
        ),
        "missing_behavior": (
            "any missing_or_nonfinite_RAC5_value_or_incomplete_unsuccessful_topology_"
            "supplies_no_group_evidence"
        ),
        "excluded_inputs": (
            "MOFid",
            "Zeo++",
            "source_ID",
            "StructureMatcher",
        ),
        "relation_to_default_methods": (
            "optional_reference_only; excluded_from_priority_main_and_main_union"
        ),
    },
    "mofid_v2_crystalnets": {
        "purpose": "Inspect exact MOFid-v2 plus current CrystalNets topology agreement.",
        "package_input_fields": (
            "parent_groups.mofid2_crystalnets_status",
            "parent_groups.mofid2_crystalnets_group",
            "parent_groups.mofid2_crystalnets_size",
        ),
        "release_relation": (
            "exact equality of complete canonicalized MOFid-v2 text plus a complete "
            "successful current CrystalNets scientific fingerprint"
        ),
        "eligible_mofid_v2_statuses": (
            "SUCCESS",
            "SUCCESS_TOPOLOGY_UNKNOWN",
            "SUCCESS_TOPOLOGY_ERROR",
            "SUCCESS_TOPOLOGY_TIMEOUT",
        ),
        "ineligible_mofid_v2_status_behavior": (
            "every_other_MOFid-v2_status_adds_no_edge_including_unresolved_timeout_"
            "error_no-MOF_unmatched-node_decomposition-error_and_ambiguous-node"
        ),
        "mofid_v2_crystalnets_rebuild_trigger": (
            "REBUILD_M2T_IF_AUTHORIZED_MOFID_V2_VALUES_CHANGE"
        ),
        "canonicalized_identifier_text": CANONICALIZED_IDENTIFIER_TEXT_DEFINITION,
        "crystalnets_fingerprint": CURRENT_CRYSTALNETS_FINGERPRINT_DEFINITION,
        "group_prefix_meaning": (
            "M2T- identifies this MOFid-v2-plus-CrystalNets criterion; the following "
            "compact digest labels a group and is not a CrystalNets topology name"
        ),
        "missing_behavior": (
            "missing_MOFid-v2_or_incomplete_unsuccessful_topology_supplies_no_group_evidence"
        ),
        "excluded_inputs": ("RAC5", "Zeo++", "source_ID", "StructureMatcher"),
        "relation_to_default_methods": (
            "optional_reference_only; excluded_from_priority_main_and_main_union; "
            "provisional_whenever_the_release-authorized_MOFid-v2_input_is_provisional"
        ),
        "provisional_scope": (
            "the_current_v26_method_remains_provisional_until_the_pinned_MOFid-v2_"
            "bundle_is_promoted_and_the_parent_table_is_rebuilt"
        ),
    },
    "structure_matcher_strict": {
        "purpose": (
            "Inspect connected components formed from audited direct symmetric strict "
            "pymatgen StructureMatcher pair edges."
        ),
        "package_input_fields": (
            "parent_groups.sm_status",
            "parent_groups.sm_group",
            "parent_groups.sm_size",
            "parent_group_methods.criteria.structure_matcher_strict",
            "parent_groups/structure_matcher_strict_evidence_receipt.json",
        ),
        "release_relation": (
            "add an undirected edge only for an audited pair whose forward and reverse "
            "strict fits both succeed, then report graph connected components"
        ),
        "direct_edge_authority": (
            "the_pinned_strict_pymatgen_symmetric_direct_pair_edge_ledger_is_authoritative"
        ),
        "component_interpretation": (
            "SM-_connected_components_are_convenience_views_for_a_non-transitive_"
            "tolerance_relation_and_are_not_duplicate_proof"
        ),
        "failure_projection": (
            "parser_timeout_OOM_matcher_or_execution_error_is_NOT_AVAILABLE_rather_than_"
            "unmatched"
        ),
        "strict_method": STRICT_STRUCTURE_MATCHER_DEFINITION,
        "strict_matcher_settings": STRICT_STRUCTURE_MATCHER_DEFINITION[
            "matcher_settings"
        ],
        "group_prefix_meaning": (
            "SM- identifies a connected-component convenience view; it does not assert "
            "that every member pair directly matches"
        ),
        "conflict_behavior": (
            "directional_fit_disagreement_or_unresolved_pair_creates_no_edge; any "
            "incomplete_component_is_projected_to_structure_specific_NOT_AVAILABLE_singletons"
        ),
        "missing_behavior": "unavailable evidence follows the selected singleton_or_exclude policy",
        "excluded_inputs": "historical_relaxed_StructureMatcher_evidence",
        "relation_to_default_methods": (
            "optional_reference_only; direct_edge_ledger_is_authoritative_and_the_method_"
            "is_excluded_from_priority_main_and_main_union"
        ),
    },
    "identity_union": {
        "display_name": "provisional source-ID/MOFid transitive groups",
        "purpose": (
            "Provide a broad provisional screening relation for records connected by "
            "source-ID or MOFid identifier evidence."
        ),
        "package_input_fields": (
            "parent_groups.identity_status",
            "parent_groups.identity_group",
            "parent_groups.identity_size",
        ),
        "release_construction": {
            "v26.0.1": (
                "freshly recompute equality edges and connected components over all "
                "and only the 36,628 current rows of this named release"
            ),
            "v26.0.2": (
                "freshly recompute equality edges and connected components over all "
                "42,574 current superset rows; do not seed or import a v26.0.1 component"
            ),
            "direct_edges": (
                "equal canonical database-namespaced source key; equal complete canonical "
                "MOFid-v2; or equal complete canonical MOFid-v1"
            ),
            "eligible_mofid_statuses": (
                "SUCCESS",
                "SUCCESS_TOPOLOGY_UNKNOWN",
                "SUCCESS_TOPOLOGY_ERROR",
                "SUCCESS_TOPOLOGY_TIMEOUT",
            ),
            "earlier_component_or_edge_import": False,
        },
        "algorithm": (
            "take transitive connected components of every available listed equality "
            "edge; there is no precedence and one criterion may bridge two others"
        ),
        "counted_unit": (
            "identity_size_is_the_number_of_structures_in_one_transitive_connected_"
            "component_not_a_count_of_edges_or_identifiers"
        ),
        "canonicalized_identifier_text": CANONICALIZED_IDENTIFIER_TEXT_DEFINITION,
        "missing_behavior": (
            "a missing_or_rejected_identifier_or_a_non-success_MOFid_status_adds_no_edge; "
            "unresolved_reconciliation_ambiguous_node_timeout_error_no-MOF_unmatched-node_"
            "and_decomposition-error_values_never_match; common_nulls_never_match; "
            "a release NOT_AVAILABLE row "
            "follows the selected singleton_or_exclude policy"
        ),
        "conflict_behavior": "none_all_available_listed_edges_are_unioned_transitively",
        "excluded_inputs": (
            "RAC5",
            "Zeo++",
            "CrystalNets_topology",
            "CIF_hash",
            "common_name",
            "StructureMatcher",
        ),
        "relation_to_default_methods": (
            "reference_screening_relation_only; it is not proof of structural identity "
            "and is excluded as an edge source from both priority_main and main_union; "
            "it_remains_provisional_until_the_pinned_MOFid_bundle_is_promoted_and_the_"
            "parent_table_is_rebuilt"
        ),
        "user_facing_definition": (
            "identity_union is the compatibility key for provisional source-ID/MOFid "
            "transitive groups, freshly recomputed over all current rows of each named "
            "release without importing an earlier component or edge. Each group is one "
            "connected component of direct edges from equal canonical database-namespaced "
            "source_database/source_id pairs or equal complete canonical MOFid-v2 and "
            "MOFid-v1 strings whose status is exactly SUCCESS, SUCCESS_TOPOLOGY_UNKNOWN, "
            "SUCCESS_TOPOLOGY_ERROR, or SUCCESS_TOPOLOGY_TIMEOUT. Every other status and "
            "missing or rejected evidence adds no edge. It does not use CIF coordinates "
            "or hashes, RAC5, Zeo++, CrystalNets, common names, DOI, or StructureMatcher; "
            "it is not proof of identity or common parentage, is candidate-only when its "
            "MOFid input is STAGE_ONLY, and is not main_union."
        ),
        "provisional_scope": (
            "current_v26_identity_components_are_provisional_while_their_MOFid_"
            "dependent_edges_use_unpromoted_input"
        ),
    },
    "none": {
        "purpose": "Run a control with one explanatory singleton per structure.",
        "package_input_fields": (),
        "release_relation": "no_equality_or_similarity_relation_is_read",
        "missing_behavior": "not_applicable_every_structure_is_available_as_its_own_singleton",
        "conflict_behavior": "none",
        "excluded_inputs": "all_parent_and_duplicate_evidence",
        "relation_to_default_methods": (
            "control_only; auto pairs it with parent_only so no broader edge is added"
        ),
    },
}


def parent_method_definition(method: str) -> Mapping[str, object]:
    """Return the machine-readable semantics of an explanatory parent method."""

    if not isinstance(method, str):
        raise TypeError("parent method must be a string")
    if method not in PARENT_METHODS or method == "main_union":
        raise ValueError("%r is not an explanatory parent method" % method)
    if method == "priority_main":
        return PRIORITY_MAIN_DEFINITION
    details = dict(_PUBLISHED_METHOD_DETAILS[method])
    details["input_fields"] = details.pop("package_input_fields")
    details["relation_to_other_terms"] = details.pop(
        "relation_to_default_methods"
    )
    details["parent_group_triad"] = (
        PARENT_GROUP_TRIAD_DEFINITION
        if method != "none"
        else "not_applicable_no_release_parent_criterion_is_read"
    )
    details.update(
        {
            "schema_version": PARENT_METHOD_CONTRACT_VERSION,
            "identifier": method,
            "display_name": _PARENT_METHOD_DISPLAY_NAMES[method],
            "project_defined_identifier": True,
            "role": "explanatory_parent_resolution",
            "summary": details.get("release_relation", details["purpose"]),
            "criterion": None if method == "none" else method,
            "resolution_scope": "complete_release_before_optional_subset",
            "algorithm": details.get(
                "algorithm",
                (
                    "use each available release-authorized group exactly as stored in "
                    "the loaded release table; do not combine it with another criterion"
                ),
            ),
            "conflict_behavior": details.get(
                "conflict_behavior",
                "none_single_criterion_groups_have_no_cross_criterion_precedence",
            ),
            "missing_behavior": details.get(
                "missing_behavior",
                (
                    "use_a_unique_per_structure_singleton_by_default_or_explicit_"
                    "MISSING_PARENT_EVIDENCE_exclusion"
                ),
            ),
            "missing_evidence": {
                "singleton": "assign_unique_SINGLETON_structure_id_group",
                "exclude": "exclude_as_MISSING_PARENT_EVIDENCE",
                "common_nulls_are_grouped": False,
            },
        }
    )
    return _freeze_contract(details)


def leakage_guard_definition(guard: str) -> Mapping[str, object]:
    """Return exact semantics for ``auto`` or a concrete leakage guard."""

    if not isinstance(guard, str):
        raise TypeError("leakage guard must be a string")
    if guard == "auto":
        return AUTO_LEAKAGE_GUARD_DEFINITION
    if guard == "main_union":
        return MAIN_UNION_DEFINITION
    if guard == "parent_only":
        return PARENT_ONLY_DEFINITION
    raise ValueError(
        "leakage guard must be 'auto', 'parent_only', or 'main_union'"
    )


def generated_output_terminology_definitions() -> Mapping[str, object]:
    """Return compact, locally complete definitions for generated receipts.

    Receipt writers place this closure before method-valued fields when keys
    are serialized in sorted order.  That makes the first machine-key use
    self-defining without duplicating the full scientific contracts.
    """

    return _freeze_contract(
        {
            "auto": (
                "auto is the parent-method-aware leakage-guard selector: it resolves "
                "priority_main to main_union and resolves every explicitly selected "
                "direct or reference explanatory parent method to parent_only before "
                "split blocks are constructed."
            ),
            "canonical_identifier_text": CANONICALIZED_IDENTIFIER_TEXT_DEFINITION[
                "user_facing_definition"
            ],
            "identity_union": parent_method_definition("identity_union")[
                "user_facing_definition"
            ],
            "main_union": MAIN_UNION_DEFINITION["user_facing_definition"],
            "mofid_v2_crystalnets": (
                "mofid_v2_crystalnets uses M2T- for an optional non-decisive reference "
                "group requiring exact equality of the complete canonical current "
                "MOFid-v2 string plus the same complete current-success CrystalNets "
                "fingerprint. Eligible MOFid-v2 statuses are exactly SUCCESS, "
                "SUCCESS_TOPOLOGY_UNKNOWN, SUCCESS_TOPOLOGY_ERROR, and "
                "SUCCESS_TOPOLOGY_TIMEOUT; the latter two are successful calculated "
                "identifiers whose embedded topology qualifier is ERROR or TIMEOUT, not "
                "MOFid execution failures. Every other status or incomplete input adds "
                "no edge. It is provisional/reference-only, must be rebuilt if authorized "
                "MOFid-v2 values change, and enters neither priority_main nor main_union."
            ),
            "parent_only": PARENT_ONLY_DEFINITION["user_facing_definition"],
            "priority_main": PRIORITY_MAIN_DEFINITION["user_facing_definition"],
            "rac5_crystalnets": (
                "rac5_crystalnets uses RT- for an optional non-decisive reference group "
                "requiring exact equality of all 264 available finite depth-5 RAC5 "
                "binary64 values plus the complete current-success CrystalNets fingerprint. "
                "Missing, nonfinite, partial, timed-out, failed, or otherwise incomplete "
                "input adds no match; it enters neither priority_main nor main_union."
            ),
            "structure_matcher_strict": (
                "structure_matcher_strict uses SM- for an optional convenience connected-"
                "component view over symmetric direct matches from the pinned strict "
                "pymatgen protocol. Direct edges are authoritative; components are not "
                "duplicate proof. Parser failures, timeout, OOM, matcher, and execution "
                "errors are NOT_AVAILABLE rather than unmatched. It enters neither "
                "priority_main nor main_union."
            ),
        }
    )


def resolve_leakage_guard(guard: str, parent_method: str) -> str:
    """Resolve ``auto`` without duplicating its project-defined policy."""

    # Validate both public identifiers through their authoritative lookups.
    leakage_guard_definition(guard)
    parent_method_definition(parent_method)
    if guard == "auto":
        return "main_union" if parent_method == "priority_main" else "parent_only"
    return guard


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
        # A validated release generation must still be byte-for-byte and
        # structure-for-structure identical to the generation sealed by the
        # loader.  Generic dataset-like inputs remain supported but carry no
        # release authority.
        from CoREMOF.dataset import (
            _reject_retired_or_reserved_keys,
            _validate_dataset_generation_if_present,
        )

        _validate_dataset_generation_if_present(dataset)
        self.dataset = dataset
        self.missing_parent = self._validate_missing_parent(missing_parent)
        rows = tuple(getattr(dataset, "metadata_rows"))
        self._rows = rows
        self._row_by_id: Dict[str, Mapping[str, object]] = {}
        self._ids: List[str] = []
        for row_index, row in enumerate(rows):
            if not isinstance(row, Mapping):
                raise TypeError("Every metadata row must be a mapping")
            _reject_retired_or_reserved_keys(
                row, "metadata_rows[%d]" % row_index
            )
            raw_structure_id = row.get("structure_id")
            if (
                not isinstance(raw_structure_id, str)
                or not raw_structure_id
                or raw_structure_id != raw_structure_id.strip()
            ):
                raise ValueError(
                    "Every metadata row must have an exact nonblank string structure_id"
                )
            structure_id = raw_structure_id
            if structure_id in self._row_by_id:
                raise ValueError("Duplicate structure_id in metadata_rows: %s" % structure_id)
            self._ids.append(structure_id)
            self._row_by_id[structure_id] = row
        parent_by_id = getattr(dataset, "parent_by_id", None)
        if parent_by_id is None:
            parent_by_id = {}
        if not isinstance(parent_by_id, Mapping):
            raise TypeError("parent_by_id must be a mapping when present")
        _reject_retired_or_reserved_keys(parent_by_id, "parent_by_id")
        for declaration_name in ("parent_group_methods", "dataset_info"):
            declaration = getattr(dataset, declaration_name, None)
            if declaration is not None:
                if not isinstance(declaration, Mapping):
                    raise TypeError("%s must be a mapping when present" % declaration_name)
                _reject_retired_or_reserved_keys(declaration, declaration_name)
        if any(not isinstance(key, str) or not key or key != key.strip() for key in parent_by_id):
            raise ValueError(
                "parent_by_id top-level keys must be exact nonblank strings"
            )
        allowed_entry_keys = {"structure_id"}
        for method, keys in _METHOD_KEYS.items():
            allowed_entry_keys.update(keys)
            base = _RELEASE_BASE[method]
            for stem in (base, method):
                allowed_entry_keys.update(
                    {
                        "%s_status" % stem,
                        "%s_group" % stem,
                        "%s_size" % stem,
                    }
                )
        for structure_id in self._ids:
            if structure_id not in parent_by_id:
                continue
            entry = parent_by_id[structure_id]
            if entry is not None and not isinstance(entry, (Mapping, str)):
                raise TypeError(
                    "parent_by_id entry for %s must be a mapping, a legacy "
                    "identity_union string, or None" % structure_id
                )
            if isinstance(entry, str) and (not entry or entry != entry.strip()):
                raise ValueError(
                    "legacy identity_union parent group must be an exact nonblank string"
                )
            if isinstance(entry, Mapping):
                if any(
                    not isinstance(key, str)
                    or not key
                    or key != key.strip()
                    for key in entry
                ):
                    raise ValueError(
                        "parent entry keys must be exact nonblank strings"
                    )
                unknown_entry_keys = set(entry).difference(allowed_entry_keys)
                if unknown_entry_keys:
                    raise ValueError(
                        "parent entry for %s contains unsupported, retired, or "
                        "reserved keys: %s"
                        % (
                            structure_id,
                            ", ".join(sorted(unknown_entry_keys)),
                        )
                    )
                if (
                    "structure_id" in entry
                    and entry["structure_id"] != structure_id
                ):
                    raise ValueError(
                        "parent entry structure_id differs from its table key"
                    )
        method_major_keys = {
            key for keys in _METHOD_KEYS.values() for key in keys
        }
        universe_ids = set(self._ids)
        ambiguous_top_level = universe_ids.intersection(method_major_keys).intersection(
            parent_by_id
        )
        if ambiguous_top_level:
            raise ValueError(
                "parent_by_id top-level keys are ambiguous between structure IDs "
                "and method-major aliases: %s"
                % ", ".join(sorted(ambiguous_top_level))
            )
        unknown_top_level = set(parent_by_id).difference(
            universe_ids.union(method_major_keys)
        )
        if unknown_top_level:
            raise ValueError(
                "parent_by_id contains unknown top-level ID/method keys: %s"
                % ", ".join(sorted(map(str, unknown_top_level)))
            )
        for key in method_major_keys.intersection(parent_by_id):
            method_map = parent_by_id[key]
            if method_map is not None and not isinstance(method_map, Mapping):
                raise TypeError(
                    "method-major parent table %s must be a mapping or None" % key
                )
            if isinstance(method_map, Mapping):
                unknown_ids = set(method_map).difference(universe_ids)
                if unknown_ids:
                    raise ValueError(
                        "method-major parent table %s contains unknown structure IDs: %s"
                        % (key, ", ".join(sorted(map(str, unknown_ids))))
                    )
        self._parent_by_id = parent_by_id

        cif_hashes = getattr(dataset, "cif_hashes", None)
        if cif_hashes is None:
            cif_hashes = {}
        if not isinstance(cif_hashes, Mapping):
            raise TypeError("cif_hashes must be a mapping when present")
        unknown_hash_ids = set(cif_hashes).difference(universe_ids)
        if unknown_hash_ids:
            raise ValueError(
                "cif_hashes contains unknown structure IDs: %s"
                % ", ".join(sorted(map(str, unknown_hash_ids)))
            )
        self._cif_hashes = cif_hashes
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
        from CoREMOF.dataset import _canonical_complete_mofid_v2

        for method, keys in _METHOD_KEYS.items():
            base = _RELEASE_BASE[method]
            group_keys = tuple(
                dict.fromkeys(("%s_group" % base, "%s_group" % method))
            )
            status_keys = tuple(
                dict.fromkeys(("%s_status" % base, "%s_status" % method))
            )
            size_keys = tuple(
                dict.fromkeys(("%s_size" % base, "%s_size" % method))
            )
            release_keys = set(group_keys + status_keys + size_keys)
            declared: List[Tuple[str, Mapping[str, object], str, str]] = []
            m2t_mofid_by_group: Dict[str, str] = {}
            for structure_id in self._ids:
                entry = self._parent_by_id.get(structure_id)
                if not isinstance(entry, Mapping):
                    continue
                present_release_keys = release_keys.intersection(entry)
                if not present_release_keys:
                    continue
                if not (
                    any(key in entry for key in group_keys)
                    and any(key in entry for key in status_keys)
                    and any(key in entry for key in size_keys)
                ):
                    raise ValueError(
                        "%s release parent entry for %s must contain a complete "
                        "status/group/size triad" % (method, structure_id)
                    )
                status = self._release_status(entry, method)
                group = self._release_group(entry, method, status)
                self._release_size(entry, method, status)
                if (
                    method == "mofid_v2_crystalnets"
                    and status in {"MATCHED", "UNMATCHED"}
                ):
                    mofid_status = self._row_by_id[structure_id].get(
                        "mofid_v2_status"
                    )
                    eligible = {
                        "SUCCESS",
                        "SUCCESS_TOPOLOGY_UNKNOWN",
                        "SUCCESS_TOPOLOGY_ERROR",
                        "SUCCESS_TOPOLOGY_TIMEOUT",
                    }
                    if mofid_status not in eligible:
                        raise ValueError(
                            "mofid_v2_crystalnets %s row for %s requires an exact "
                            "eligible mofid_v2_status" % (status, structure_id)
                        )
                    canonical_mofid = _canonical_complete_mofid_v2(
                        self._row_by_id[structure_id].get("mofid_v2"),
                        structure_id,
                    )
                    previous_mofid = m2t_mofid_by_group.setdefault(
                        group, canonical_mofid
                    )
                    if previous_mofid != canonical_mofid:
                        raise ValueError(
                            "mofid_v2_crystalnets group %s spans conflicting "
                            "canonical mofid_v2 values" % group
                        )
                for key in keys:
                    if key in entry and key not in release_keys:
                        raise ValueError(
                            "%s release parent entry for %s mixes a validated triad "
                            "with legacy alias %s" % (method, structure_id, key)
                        )
                declared.append((structure_id, entry, group, status))
            if not declared:
                continue
            if method in OPTIONAL_CRYSTALNETS_PARENT_METHODS and len(declared) != len(
                self._ids
            ):
                raise ValueError(
                    "%s release parent criterion must contain a complete validated "
                    "status/group/size triad for every structure" % method
                )
            counts: Dict[str, int] = {}
            for _, _, group, _ in declared:
                counts[group] = counts.get(group, 0) + 1
            for structure_id, entry, group, status in declared:
                self._validate_release_size(entry, method, status)
                advertised = self._release_size(entry, method, status)
                if advertised != counts[group]:
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
        Call :func:`parent_method_definition` for the complete machine-readable
        meaning of a selectable method.  ``main_union`` is a leakage graph, not
        an explanatory parent hierarchy; its contract is returned by
        :func:`leakage_guard_definition`.
        """

        if not isinstance(method, str):
            raise TypeError("parent method must be a string")
        if method not in PARENT_METHODS:
            raise ValueError(
                "Unknown parent method %r; choose one of %s"
                % (method, ", ".join(PARENT_METHODS))
            )
        missing = (
            self.missing_parent
            if missing_parent is None
            else self._validate_missing_parent(missing_parent)
        )
        if method in OPTIONAL_REFERENCE_PARENT_METHODS and not self._method_is_declared(
            method
        ):
            raise ValueError(
                "Parent method %r is not present in this release" % method
            )
        if method in OPTIONAL_CRYSTALNETS_PARENT_METHODS:
            # Import lazily so the lightweight parent module retains a simple
            # dependency boundary while sharing the loader's exact trusted
            # terminal RT/M2T contract.
            from CoREMOF.dataset import (
                _has_authenticated_crystalnets_reference_contract,
            )

            if not _has_authenticated_crystalnets_reference_contract(
                self.dataset, method
            ):
                raise ValueError(
                    "Parent method %r requires the exact validated release "
                    "criterion, definition, and integration contracts" % method
                )
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
        if isinstance(structure_ids, str):
            requested_values = (structure_ids,)
        elif isinstance(structure_ids, (list, tuple)):
            requested_values = tuple(structure_ids)
        else:
            raise TypeError(
                "structure_ids must be an exact string or ordered list/tuple "
                "of strings"
            )
        if any(
            not isinstance(value, str)
            or not value
            or value != value.strip()
            for value in requested_values
        ):
            raise ValueError(
                "structure_ids must contain exact nonblank strings"
            )
        requested = set(requested_values)
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

    def _method_is_declared(self, method: str) -> bool:
        """Return whether any release row declares the optional criterion."""

        if not isinstance(self._parent_by_id, Mapping):
            return False
        base = _RELEASE_BASE[method]
        triad_keys = {
            "%s_status" % base,
            "%s_group" % base,
            "%s_size" % base,
            "%s_status" % method,
            "%s_group" % method,
            "%s_size" % method,
        }
        for structure_id in self._ids:
            entry = self._parent_by_id.get(structure_id)
            if isinstance(entry, Mapping) and any(key in entry for key in triad_keys):
                return True
        return False

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
                base = _RELEASE_BASE[method]
                release_group_keys = tuple(
                    dict.fromkeys(("%s_group" % base, "%s_group" % method))
                )
                for key in release_group_keys:
                    if key in entry:
                        candidate = entry[key]
                        candidate_key = key
                        break
                if candidate_key is None:
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
            self._release_size(release_entry, method, status)
            candidate = self._release_group(release_entry, method, status)
            if status == "NOT_AVAILABLE":
                return None
        if isinstance(candidate, Mapping):
            for key in ("group", "group_id", "parent_group", "value"):
                if key in candidate:
                    candidate = candidate[key]
                    break
            else:
                raise ValueError("%s parent group must be a nonblank string" % method)
        if candidate is None:
            return None
        if not isinstance(candidate, str):
            raise ValueError("%s parent group must be a nonblank string" % method)
        if not candidate or candidate != candidate.strip():
            raise ValueError("%s parent group must be a nonblank string" % method)
        if candidate.upper() in _MISSING_TEXT:
            return None
        return candidate

    @staticmethod
    def _release_status(entry: Mapping[str, object], method: str) -> str:
        base = _RELEASE_BASE[method]
        values = [
            entry[key]
            for key in tuple(
                dict.fromkeys(("%s_status" % base, "%s_status" % method))
            )
            if key in entry
        ]
        if not values or any(not isinstance(value, str) for value in values):
            status = None
        else:
            if len(set(values)) != 1:
                raise ValueError("%s release parent status aliases conflict" % method)
            status = values[0]
        if status not in {"MATCHED", "UNMATCHED", "NOT_AVAILABLE"}:
            raise ValueError(
                "%s release parent status must be MATCHED, UNMATCHED, or "
                "NOT_AVAILABLE (got %r)" % (method, status)
            )
        return status

    @staticmethod
    def _validate_release_group(value: object, method: str, status: str) -> str:
        if not isinstance(value, str):
            raise ValueError("%s release parent group must be a nonblank string" % method)
        if not value or not value.strip():
            raise ValueError("%s release parent group must be a nonblank string" % method)
        if value != value.strip():
            raise ValueError(
                "%s release parent group must be an exact canonical string" % method
            )
        group = value
        if status in {"MATCHED", "UNMATCHED"} and group.upper() in _MISSING_TEXT:
            raise ValueError(
                "%s %s release parent group must not be a missing placeholder"
                % (method, status)
            )
        optional_prefixes = {
            "rac5_crystalnets": "RT-",
            "mofid_v2_crystalnets": "M2T-",
        }
        prefix = optional_prefixes.get(method)
        if prefix is not None and (
            not group.startswith(prefix)
            or len(group) < len(prefix) + 8
            or len(group) > len(prefix) + 64
            or any(character not in "0123456789ABCDEF" for character in group[len(prefix) :])
        ):
            raise ValueError(
                "%s release parent group must use %s plus 8 to 64 uppercase "
                "hexadecimal characters" % (method, prefix)
            )
        return group

    @classmethod
    def _release_group(
        cls, entry: Mapping[str, object], method: str, status: str
    ) -> str:
        base = _RELEASE_BASE[method]
        keys = tuple(dict.fromkeys(("%s_group" % base, "%s_group" % method)))
        groups = [
            cls._validate_release_group(entry[key], method, status)
            for key in keys
            if key in entry
        ]
        if not groups:
            raise ValueError("%s release parent group must be a nonblank string" % method)
        if len(set(groups)) != 1:
            raise ValueError("%s release parent group aliases conflict" % method)
        return groups[0]

    @classmethod
    def _release_size(
        cls, entry: Mapping[str, object], method: str, status: str
    ) -> int:
        base = _RELEASE_BASE[method]
        keys = tuple(dict.fromkeys(("%s_size" % base, "%s_size" % method)))
        values = []
        for key in keys:
            if key in entry:
                cls._validate_release_size_value(entry[key], method, status)
                value = entry[key]
                values.append(int(value.strip()) if isinstance(value, str) else value)
        if not values:
            raise ValueError("%s release parent size must be a positive integer" % method)
        if len(set(values)) != 1:
            raise ValueError("%s release parent size aliases conflict" % method)
        return values[0]

    @staticmethod
    def _validate_release_size(
        entry: Mapping[str, object], method: str, status: str
    ) -> None:
        base = _RELEASE_BASE[method]
        values = [
            entry[key]
            for key in tuple(dict.fromkeys(("%s_size" % base, "%s_size" % method)))
            if key in entry
        ]
        if not values:
            raise ValueError("%s release parent size must be a positive integer" % method)
        parsed = [
            ParentResolver._validate_release_size_value(value, method, status)
            for value in values
        ]
        if len(set(parsed)) != 1:
            raise ValueError("%s release parent size aliases conflict" % method)

    @staticmethod
    def _validate_release_size_value(value: object, method: str, status: str) -> int:
        if isinstance(value, bool):
            raise ValueError("%s release parent size must be a positive integer" % method)
        if isinstance(value, int):
            size = value
        elif isinstance(value, str):
            text = value
            if (
                not text
                or not text.isascii()
                or not text.isdigit()
                or (len(text) > 1 and text.startswith("0"))
            ):
                raise ValueError(
                    "%s release parent size must be a positive integer" % method
                )
            size = int(text)
        else:
            raise ValueError("%s release parent size must be a positive integer" % method)
        if size < 1:
            raise ValueError("%s release parent size must be a positive integer" % method)
        if status == "MATCHED" and size < 2:
            raise ValueError("%s MATCHED parent group must have size >= 2" % method)
        if status in {"UNMATCHED", "NOT_AVAILABLE"} and size != 1:
            raise ValueError("%s %s parent group must have size 1" % (method, status))
        return size

    def _source_evidence(self, structure_id: str) -> Optional[str]:
        """Return a database-namespaced source group."""

        row = self._row_by_id[structure_id]
        raw_database = row.get("source_database")
        if raw_database is not None and not isinstance(raw_database, str):
            raise ValueError("source_database must be a string when present")
        if raw_database in (None, ""):
            return None
        if raw_database != raw_database.strip():
            raise ValueError(
                "source_database must be an exact nonblank string when present"
            )
        database = raw_database.upper()
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
            raw_source_id = row.get("source_id")
            if raw_source_id is not None and not isinstance(raw_source_id, str):
                raise ValueError("source_id must be a string when present")
            if raw_source_id in (None, ""):
                return None
            if raw_source_id != raw_source_id.strip():
                raise ValueError(
                    "source_id must be an exact nonblank string when present"
                )
            source_id = raw_source_id
            source_group = " ".join(source_id.casefold().split())
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
        if not isinstance(candidate, str):
            raise ValueError(
                "CIF integrity hashes used for main_union must be full lowercase SHA-256"
            )
        text = candidate
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
    "AUTO_LEAKAGE_GUARD_DEFINITION",
    "CANONICALIZED_IDENTIFIER_TEXT_DEFINITION",
    "COMPUTED_PARENT_METHODS",
    "CURRENT_CRYSTALNETS_FINGERPRINT_DEFINITION",
    "DIRECT_PARENT_METHODS",
    "LEAKAGE_GUARD_CHOICES",
    "MAIN_UNION_DEFINITION",
    "OPTIONAL_REFERENCE_PARENT_METHODS",
    "OPTIONAL_CRYSTALNETS_PARENT_METHODS",
    "PARENT_METHODS",
    "PARENT_METHOD_CONTRACT_VERSION",
    "PARENT_GROUP_TRIAD_DEFINITION",
    "PARENT_ONLY_DEFINITION",
    "PRIORITY_MAIN_DEFINITION",
    "RAC5_NUMERIC_FINGERPRINT_DEFINITION",
    "REFERENCE_PARENT_METHODS",
    "SELECTABLE_PARENT_METHODS",
    "STRICT_STRUCTURE_MATCHER_DEFINITION",
    "ZEO_NUMERIC_FINGERPRINT_DEFINITION",
    "ParentConflict",
    "ParentResolution",
    "ParentResolver",
    "generated_output_terminology_definitions",
    "leakage_guard_definition",
    "parent_method_definition",
    "resolve_leakage_guard",
]
