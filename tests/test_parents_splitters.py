import json
import hashlib
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF.parents import (
    AUTO_LEAKAGE_GUARD_DEFINITION,
    CANONICALIZED_IDENTIFIER_TEXT_DEFINITION,
    CURRENT_CRYSTALNETS_FINGERPRINT_DEFINITION,
    LEAKAGE_GUARD_CHOICES,
    MAIN_UNION_DEFINITION,
    PARENT_GROUP_TRIAD_DEFINITION,
    PARENT_METHOD_CONTRACT_VERSION,
    PRIORITY_MAIN_DEFINITION,
    RAC5_NUMERIC_FINGERPRINT_DEFINITION,
    SELECTABLE_PARENT_METHODS,
    STRICT_STRUCTURE_MATCHER_DEFINITION,
    ZEO_NUMERIC_FINGERPRINT_DEFINITION,
    ParentResolver,
    leakage_guard_definition,
    parent_method_definition,
    resolve_leakage_guard,
)
import CoREMOF.splitters as splitters_module
from CoREMOF.splitters import (
    OfficialSplitUnavailableError,
    ParentGroupSplitter,
    split_release,
)


class _Dataset:
    def __init__(self, rows, parents=None, hashes=None):
        self.metadata_rows = tuple(rows)
        self.parent_by_id = parents or {}
        if hashes is None:
            self.cif_hashes = {
                row["structure_id"]: hashlib.sha256(
                    row["structure_id"].encode("utf-8")
                ).hexdigest()
                for row in self.metadata_rows
            }
        else:
            self.cif_hashes = hashes
        self.dataset_version = "v-test"
        self.input_hashes = {"metadata.csv": "abc123"}
        self.parent_group_methods = {
            "release_status": "PROVISIONAL_LATEST_AUDITED_SNAPSHOT"
        }


class _Classified:
    def __init__(self, dataset, labels):
        self.dataset = dataset
        self.label_by_id = labels
        self.checker_view = "3checker"


def _rows(ids):
    return [
        {
            "structure_id": structure_id,
            "source_database": "COD",
            "source_id": "SRC-%s" % structure_id,
            "structure_variant": "ASR",
            "metal_elements": "Cu;Zn" if index % 2 else "Cu",
        }
        for index, structure_id in enumerate(ids)
    ]


def _unique_hashes(ids):
    return {
        structure_id: hashlib.sha256(structure_id.encode("utf-8")).hexdigest()
        for structure_id in ids
    }


class ParentResolverTests(unittest.TestCase):
    def test_every_selectable_parent_and_guard_has_a_definition(self):
        self.assertEqual(len(SELECTABLE_PARENT_METHODS), len(set(SELECTABLE_PARENT_METHODS)))
        for method in SELECTABLE_PARENT_METHODS:
            with self.subTest(parent_method=method):
                definition = parent_method_definition(method)
                self.assertEqual(
                    definition["schema_version"], PARENT_METHOD_CONTRACT_VERSION
                )
                self.assertEqual(definition["identifier"], method)
                self.assertTrue(definition["display_name"])
                self.assertTrue(definition["project_defined_identifier"])
                self.assertTrue(definition["purpose"])
                self.assertTrue(definition["summary"])
                for field in (
                    "input_fields",
                    "algorithm",
                    "conflict_behavior",
                    "missing_behavior",
                    "excluded_inputs",
                    "relation_to_other_terms",
                ):
                    self.assertIn(field, definition)
                    self.assertTrue(definition[field] or field == "input_fields")
        for guard in LEAKAGE_GUARD_CHOICES:
            with self.subTest(leakage_guard=guard):
                definition = leakage_guard_definition(guard)
                self.assertEqual(
                    definition["schema_version"], PARENT_METHOD_CONTRACT_VERSION
                )
                self.assertEqual(definition["identifier"], guard)
                self.assertTrue(definition["display_name"])
                self.assertTrue(definition["project_defined_identifier"])
                self.assertTrue(definition["purpose"])
                self.assertTrue(definition["summary"])
                for field in (
                    "input_fields",
                    "algorithm",
                    "conflict_behavior",
                    "missing_behavior",
                    "excluded_inputs",
                    "relation_to_other_terms",
                ):
                    self.assertIn(field, definition)
                    self.assertTrue(definition[field])

    def test_identity_union_and_identifier_canonicalization_are_explicit(self):
        definition = parent_method_definition("identity_union")
        self.assertEqual(
            definition["input_fields"],
            (
                "parent_groups.identity_status",
                "parent_groups.identity_group",
                "parent_groups.identity_size",
            ),
        )
        self.assertEqual(
            definition["display_name"],
            "provisional source-ID/MOFid transitive groups",
        )
        self.assertIn("v26.0.1", definition["release_construction"])
        self.assertIn("v26.0.2", definition["release_construction"])
        self.assertIn("no precedence", definition["algorithm"])
        self.assertIn("common_nulls_never_match", definition["missing_behavior"])
        self.assertIn("RAC5", definition["excluded_inputs"])
        self.assertIn("excluded", definition["relation_to_other_terms"])
        canonical = definition["canonicalized_identifier_text"]
        self.assertEqual(canonical, CANONICALIZED_IDENTIFIER_TEXT_DEFINITION)
        self.assertEqual(
            canonical["current_release_text_steps_in_order"],
            (
                "convert_the_published_value_to_text",
                "collapse_each_Unicode_whitespace_run_to_one_ASCII_space",
                "trim_leading_and_trailing_whitespace",
                "reject_a_whole_field_missing_or_execution_placeholder",
                "apply_Unicode_NFKC",
                "apply_Unicode_default_casefold",
            ),
        )
        self.assertEqual(
            canonical["current_release_case_insensitive_whole_field_placeholders"],
            (
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
        )
        self.assertTrue(
            canonical["source_key"][
                "database_namespace_prevents_cross_database_ID_matches"
            ]
        )
        self.assertFalse(canonical["mofid_key"]["fuzzy_or_partial_string_matching"])
        inherited = canonical["inherited_v2601_identity_component_cleanup"]
        self.assertEqual(
            inherited["record_or_source_ID_steps"],
            (
                "convert_the_v11_value_to_text",
                "remove_everything_from_the_first_semicolon_record_suffix_onward",
                "replace_backslashes_with_forward_slashes",
                "remove_any_directory_path",
                "remove_one_terminal_.cif_or_.cif.gz_case_insensitively",
                "trim_leading_and_trailing_whitespace",
                "delete_every_Unicode_whitespace_character",
                "apply_Unicode_default_casefold_without_NFKC",
                "reject_a_whole_field_literal_placeholder",
                "remove_only_one_recognized_terminal_processing_bundle_"
                "ASR_pacman_FSR_pacman_ION_pacman_ION_ASR_pacman_or_ION_FSR_pacman",
                "reject_an_empty_or_whole_field_literal_placeholder_again",
            ),
        )
        self.assertEqual(
            inherited["mofid_steps"],
            (
                "convert_the_v11_value_to_text",
                "remove_the_semicolon_record_ID_suffix",
                "trim_leading_and_trailing_whitespace",
                "rewrite_only_the_literal_MOFidv2._marker_as_MOFid-v2.",
                "collapse_each_Unicode_whitespace_run_to_one_ASCII_space",
                "trim_leading_and_trailing_whitespace_again",
                "apply_Unicode_default_casefold_without_NFKC",
                "remove_terminal_.no_ref_only_from_a_MOFid-v1_value",
                "reject_literal_execution_and_leading_MOFid_NA_placeholders",
            ),
        )
        self.assertFalse(inherited["v11_refcode_edges_are_database_namespaced"])
        self.assertIn(
            "no_CIF_coordinate_canonicalization", canonical["excluded_operations"]
        )
        self.assertIn(
            "whole_field_asterisk",
            inherited["additional_rejected_placeholders"],
        )

    def test_optional_prefix_contracts_define_meaning_and_failure_behavior(self):
        expected = {
            "rac5_topology": "RT-",
            "mofid_v2_topology": "M2T-",
            "structure_matcher_strict": "SM-",
        }
        for method, prefix in expected.items():
            with self.subTest(method=method):
                definition = parent_method_definition(method)
                self.assertIn(prefix, definition["group_prefix_meaning"])
                self.assertTrue(definition["missing_behavior"])
                self.assertIn("excluded", definition["relation_to_other_terms"])
        topology = parent_method_definition("rac5_topology")["topology_fingerprint"]
        self.assertEqual(topology, CURRENT_CRYSTALNETS_FINGERPRINT_DEFINITION)
        for field in (
            "network_dimension",
            "catenation_degree",
            "single_node_net",
            "all_node_net",
            "for_each_subnet_SingleNodes_status_dimension_topology_key_topology_name_"
            "topological_genome",
        ):
            self.assertIn(field, topology["included_fields"])
        self.assertIn("runtime_seconds", topology["excluded_fields"])
        self.assertTrue(topology["not_a_topology_similarity_tolerance"])
        matcher = parent_method_definition("structure_matcher_strict")
        self.assertIn("forward and reverse", matcher["release_relation"])
        self.assertIn("not assert", matcher["group_prefix_meaning"])
        strict = matcher["strict_method"]
        self.assertEqual(strict, STRICT_STRUCTURE_MATCHER_DEFINITION)
        self.assertEqual(strict["software"]["pymatgen"], "2024.2.8")
        self.assertEqual(strict["software"]["numpy"], "1.26.4")
        self.assertEqual(strict["parser"]["site_tolerance"], 0.0001)
        self.assertEqual(strict["parser"]["frac_tolerance"], 0.0001)
        self.assertEqual(
            strict["candidate_pairs"]["blocking_key"],
            "ElementComparator_fractional_composition_hash_of_the_parsed_structure",
        )
        settings = strict["matcher_settings"]
        self.assertEqual(settings["ltol"], 0.001)
        self.assertEqual(settings["stol"], 0.001)
        self.assertEqual(settings["angle_tol"], 0.01)
        self.assertTrue(settings["fit_symmetric_argument"])
        self.assertEqual(settings["supercell_size"], "num_sites")
        self.assertIn("(V/Nsites)^(1/3)", strict["displacement_diagnostics"]["normalization"])
        self.assertEqual(strict["displacement_diagnostics"]["units"], "dimensionless")
        self.assertFalse(strict["displacement_diagnostics"]["is_angstrom_RMSD"])

    def test_parent_triad_numeric_and_current_topology_contracts_are_source_exact(self):
        triad = PARENT_GROUP_TRIAD_DEFINITION
        self.assertIn("size_at_least_2", triad["status_semantics"]["MATCHED"])
        self.assertIn("size_exactly_1", triad["status_semantics"]["UNMATCHED"])
        self.assertIn("no_scientific_edge", triad["status_semantics"]["NOT_AVAILABLE"])
        self.assertIn("passes loader validation", triad["release_authorized_meaning"])

        rac = parent_method_definition("rac5")["rac5_fingerprint"]
        self.assertEqual(rac, RAC5_NUMERIC_FINGERPRINT_DEFINITION)
        self.assertIn("IEEE754_binary64", rac["value_conversion"])
        self.assertIn("negative_zero_to_positive_zero", rac["value_conversion"])
        self.assertIn("Python_float.hex", rac["value_conversion"])
        self.assertIn("rtol_0_and_atol_0", rac["comparison"])

        zeo = parent_method_definition("zeo")["zeo_fingerprint"]
        self.assertEqual(zeo, ZEO_NUMERIC_FINGERPRINT_DEFINITION)
        self.assertEqual(len(zeo["numeric_fields"]), 13)
        self.assertEqual(dict(zeo["probe_radii_A"]), {"N2": 1.655, "He": 1.32})
        self.assertIn(
            "periodicity_available_true_for_both_rows", zeo["hard_gates"]
        )

        topology = CURRENT_CRYSTALNETS_FINGERPRINT_DEFINITION
        self.assertIn("loaded_release", topology["current_meaning"])
        self.assertIn("equals_the_observed_subnet_count", topology["field_validation"]["count_fields"])
        self.assertIn("may_be_null", topology["field_validation"]["top_level_net_and_agreement"])
        self.assertIn("topology_name_and_topological_genome_may_be_null", topology["field_validation"]["subnet_nodes"])

        priority_inputs = PRIORITY_MAIN_DEFINITION["input_fields"]
        self.assertEqual(priority_inputs["rac5"][-1], "parent_groups.rac_size")
        self.assertEqual(priority_inputs["mofid_v2"][-1], "parent_groups.mofid2_size")
        self.assertEqual(priority_inputs["mofid_v1"][-1], "parent_groups.mofid1_size")

    def test_auto_guard_lookup_and_resolution_are_explicit(self):
        definition = leakage_guard_definition("auto")
        self.assertIs(definition, AUTO_LEAKAGE_GUARD_DEFINITION)
        self.assertEqual(
            tuple(
                (rule["when_parent_method"], rule["resolved_guard"])
                for rule in definition["resolution_rules"]
            ),
            (
                ("priority_main", "main_union"),
                ("any_other_explanatory_parent_method", "parent_only"),
            ),
        )
        for method in SELECTABLE_PARENT_METHODS:
            with self.subTest(parent_method=method):
                expected = "main_union" if method == "priority_main" else "parent_only"
                self.assertEqual(resolve_leakage_guard("auto", method), expected)
                self.assertEqual(resolve_leakage_guard(expected, method), expected)

    def test_priority_main_has_an_explicit_immutable_contract(self):
        definition = parent_method_definition("priority_main")
        self.assertIs(definition, PRIORITY_MAIN_DEFINITION)
        self.assertTrue(definition["project_defined_identifier"])
        self.assertEqual(
            definition["priority_order"], ("rac5", "mofid_v2", "mofid_v1")
        )
        self.assertEqual(
            definition["lower_group_resolution"]["multiple_stronger_components"],
            "keep_stronger_components_separate_record_parent_conflict_and_leave_"
            "unresolved_members_unassigned",
        )
        self.assertFalse(definition["conflict"]["stronger_components_are_merged"])
        self.assertFalse(
            definition["missing_evidence"]["common_nulls_are_grouped"]
        )
        self.assertNotIn("structure_matcher_strict", definition["priority_order"])
        self.assertEqual(
            definition["output_group_ids"],
            {
                "rac5_anchor": "RAC5:<published_rac5_group>",
                "mofid_v2_component": "MOFID_V2:<published_mofid_v2_group>",
                "mofid_v1_component": "MOFID_V1:<published_mofid_v1_group>",
                "missing_singleton": "SINGLETON:<structure_id>",
                "attached_member_rule": "keep_the_stronger_component_group_id",
                "published_group_text_is_preserved": True,
            },
        )
        with self.assertRaises(TypeError):
            definition["summary"] = "changed"

    def test_main_union_contract_is_distinct_from_priority_explanation(self):
        definition = leakage_guard_definition("main_union")
        self.assertIs(definition, MAIN_UNION_DEFINITION)
        self.assertTrue(definition["project_defined_identifier"])
        self.assertEqual(definition["role"], "leakage_guard")
        self.assertFalse(definition["explanatory_parent_method"])
        self.assertEqual(definition["universe"], "complete_unfiltered_release")
        self.assertEqual(
            tuple(edge["criterion"] for edge in definition["edge_sources"]),
            ("cif_sha256", "source_id", "rac5", "mofid_v2", "mofid_v1"),
        )
        self.assertTrue(definition["edge_sources"][0]["required_for_every_structure"])
        self.assertEqual(
            definition["input_fields"]["source_siblings"][:3],
            (
                "parent_groups.source_status",
                "parent_groups.source_group",
                "parent_groups.source_size",
            ),
        )
        fallback = definition["compatibility_source_fallback"]
        self.assertIn("manually_constructed_dataset", fallback["scope"])
        self.assertFalse(fallback["Unicode_NFKC_applied"])
        self.assertFalse(fallback["nonblank_placeholder_rejection_applied"])
        self.assertEqual(
            fallback["explicit_NOT_AVAILABLE_behavior"],
            "do_not_fallback_and_add_no_source_edge",
        )
        self.assertNotIn(
            "structure_matcher_strict",
            tuple(edge["criterion"] for edge in definition["edge_sources"]),
        )
        self.assertEqual(
            definition["output_group_id"],
            {
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
        )
        with self.assertRaisesRegex(ValueError, "not an explanatory parent method"):
            parent_method_definition("main_union")

    def test_optional_topology_criterion_must_be_declared_by_release(self):
        resolver = ParentResolver(_Dataset(_rows(("a", "b"))))
        with self.assertRaisesRegex(ValueError, "not present in this release"):
            resolver.resolve("rac5_topology")

    def test_structure_matcher_is_optional_reference_not_default_hierarchy(self):
        ids = ("a", "b", "c")
        parents = {
            "a": {
                "sm_group": "SM-AAAA0001",
                "sm_status": "MATCHED",
                "sm_size": "2",
            },
            "b": {
                "sm_group": "SM-AAAA0001",
                "sm_status": "MATCHED",
                "sm_size": "2",
            },
            "c": {
                "sm_group": "SM-CCCC0001",
                "sm_status": "UNMATCHED",
                "sm_size": "1",
            },
        }
        resolver = ParentResolver(_Dataset(_rows(ids), parents))
        matched = resolver.resolve("structure_matcher_strict")
        self.assertEqual(matched.groups["a"], matched.groups["b"])
        self.assertNotEqual(matched.groups["a"], matched.groups["c"])
        priority = resolver.resolve("priority_main")
        self.assertEqual(priority.groups["a"], "SINGLETON:a")
        self.assertEqual(priority.groups["b"], "SINGLETON:b")
        union = resolver.resolve("main_union")
        self.assertNotEqual(union.groups["a"], union.groups["b"])

    def test_topology_combined_criteria_are_selectable_but_not_in_priority_main(self):
        ids = ("a", "b", "c")
        parents = {
            "a": {
                "rac_topology_group": "RT-AAAA0001",
                "rac_topology_status": "MATCHED",
                "rac_topology_size": "2",
                "mofid2_topology_group": "M2T-BBBB0001",
                "mofid2_topology_status": "MATCHED",
                "mofid2_topology_size": "2",
            },
            "b": {
                "rac_topology_group": "RT-AAAA0001",
                "rac_topology_status": "MATCHED",
                "rac_topology_size": "2",
                "mofid2_topology_group": "M2T-BBBB0001",
                "mofid2_topology_status": "MATCHED",
                "mofid2_topology_size": "2",
            },
            "c": {
                "rac_topology_group": "RT-CCCC0001",
                "rac_topology_status": "NOT_AVAILABLE",
                "rac_topology_size": "1",
                "mofid2_topology_group": "M2T-DDDD0001",
                "mofid2_topology_status": "NOT_AVAILABLE",
                "mofid2_topology_size": "1",
            },
        }
        resolver = ParentResolver(_Dataset(_rows(ids), parents))
        rac_topology = resolver.resolve("rac5_topology")
        mofid_topology = resolver.resolve("mofid_v2_topology")
        self.assertEqual(rac_topology.groups["a"], rac_topology.groups["b"])
        self.assertEqual(mofid_topology.groups["a"], mofid_topology.groups["b"])
        self.assertEqual(rac_topology.groups["c"], "SINGLETON:c")
        # The established RAC5 > MOFid-v2 > MOFid-v1 hierarchy is unchanged.
        priority = resolver.resolve("priority_main")
        self.assertEqual(priority.groups["a"], "SINGLETON:a")
        self.assertEqual(priority.groups["b"], "SINGLETON:b")

    def test_missing_direct_value_is_unique_singleton_by_default(self):
        dataset = _Dataset(
            _rows(("a", "b", "c")),
            {"a": {"rac5": "R1"}, "b": {"rac5": "R1"}, "c": {}},
        )
        resolution = ParentResolver(dataset).resolve("rac5")
        self.assertEqual(resolution.groups["a"], resolution.groups["b"])
        self.assertEqual(resolution.groups["c"], "SINGLETON:c")
        self.assertFalse(resolution.exclusions)
        with self.assertRaises(TypeError):
            resolution.groups["a"] = "changed"

        strict = ParentResolver(dataset, missing_parent="exclude").resolve("rac5")
        self.assertEqual(strict.exclusions["c"], "MISSING_PARENT_EVIDENCE")

    def test_release_not_available_group_is_not_treated_as_evidence(self):
        dataset = _Dataset(
            _rows(("missing", "x", "y")),
            {
                "missing": {
                    "rac_group": "R-SYNTHETIC",
                    "rac_status": "NOT_AVAILABLE",
                    "rac_size": "1",
                },
                "x": {"rac_group": "R1", "rac_status": "MATCHED", "rac_size": "2"},
                "y": {"rac_group": "R1", "rac_status": "MATCHED", "rac_size": "2"},
            },
        )
        resolution = ParentResolver(dataset).resolve("rac5")
        self.assertEqual(resolution.groups["missing"], "SINGLETON:missing")
        self.assertEqual(resolution.evidence_by_id["missing"], ())
        self.assertEqual(resolution.groups["x"], "R1")

    def test_release_status_and_size_are_closed_and_validated(self):
        bad_status = _Dataset(
            _rows(("a",)),
            {"a": {"rac_group": "R1", "rac_status": "AVAILABLE", "rac_size": "1"}},
        )
        with self.assertRaisesRegex(ValueError, "MATCHED, UNMATCHED, or NOT_AVAILABLE"):
            ParentResolver(bad_status).resolve("rac5")
        bad_size = _Dataset(
            _rows(("a",)),
            {"a": {"rac_group": "R1", "rac_status": "MATCHED", "rac_size": "1"}},
        )
        with self.assertRaisesRegex(ValueError, "size >= 2"):
            ParentResolver(bad_size).resolve("rac5")

    def test_priority_conflict_preserves_anchors_and_excludes_lower_only_row(self):
        ids = ("a", "b", "c", "d", "e", "f", "g")
        parents = {
            "a": {"rac5": "R1", "mofid_v2": "X"},
            "b": {"rac5": "R2", "mofid_v2": "X"},
            "c": {"mofid_v2": "X", "mofid_v1": "Z"},
            "d": {"mofid_v2": "Y", "mofid_v1": "Q"},
            "e": {"mofid_v2": "Y"},
            "f": {"mofid_v1": "Q"},
            "g": {},
        }
        resolution = ParentResolver(_Dataset(_rows(ids), parents)).resolve("priority_main")
        self.assertNotEqual(resolution.groups["a"], resolution.groups["b"])
        self.assertEqual(resolution.exclusions["c"], "PARENT_METHOD_CONFLICT")
        self.assertNotIn("c", resolution.groups)
        self.assertEqual(resolution.groups["d"], resolution.groups["e"])
        self.assertEqual(resolution.groups["d"], resolution.groups["f"])
        self.assertEqual(resolution.groups["g"], "SINGLETON:g")
        self.assertEqual(len(resolution.conflicts), 1)
        conflict = resolution.conflicts[0]
        self.assertEqual(conflict.lower_method, "mofid_v2")
        self.assertEqual(conflict.lower_group, "X")
        self.assertEqual(conflict.stronger_components, ("RAC5:R1", "RAC5:R2"))
        self.assertEqual(conflict.member_ids, ("a", "b", "c"))
        self.assertEqual(conflict.unresolved_ids, ("c",))

    def test_priority_records_anchor_only_conflicts_without_excluding_rows(self):
        parents = {
            "a": {"rac5": "R1", "mofid_v2": "M-X"},
            "b": {"rac5": "R2", "mofid_v2": "M-X"},
        }
        resolution = ParentResolver(_Dataset(_rows(("a", "b")), parents)).resolve(
            "priority_main"
        )
        self.assertFalse(resolution.exclusions)
        self.assertEqual(len(resolution.conflicts), 1)
        conflict = resolution.conflicts[0]
        self.assertEqual(conflict.unresolved_ids, ())
        self.assertEqual(
            conflict.component_members,
            (("RAC5:R1", ("a",)), ("RAC5:R2", ("b",))),
        )

    def test_main_union_keeps_transitive_full_universe_bridges(self):
        ids = ("a", "b", "c", "d", "e", "f")
        rows = _rows(ids)
        # b-c is a source bridge; c-d is a RAC bridge; d-e is MOFid v2.
        for row in rows:
            if row["structure_id"] in ("b", "c"):
                row["source_id"] = "SHARED-SOURCE"
        parents = {
            "c": {"rac5": "R-CD"},
            "d": {"rac5": "R-CD", "mofid_v2": "M-DE"},
            "e": {"mofid_v2": "M-DE"},
        }
        hashes = _unique_hashes(ids)
        hashes["a"] = hashes["b"] = "a" * 64
        resolution = ParentResolver(_Dataset(rows, parents, hashes)).resolve("main_union")
        self.assertEqual(len({resolution.groups[value] for value in ("a", "b", "c", "d", "e")}), 1)
        self.assertNotEqual(resolution.groups["f"], resolution.groups["a"])
        expected = hashlib.sha256(b"a\0b\0c\0d\0e").hexdigest()[:16]
        self.assertEqual(resolution.groups["a"], "MAIN-{}".format(expected))

    def test_main_union_rejects_truncated_or_placeholder_cif_hashes(self):
        dataset = _Dataset(_rows(("a",)), hashes={"a": "abc123"})
        with self.assertRaisesRegex(ValueError, "full lowercase SHA-256"):
            ParentResolver(dataset).resolve("main_union")

    def test_main_union_requires_a_hash_for_every_structure(self):
        dataset = _Dataset(_rows(("a", "b")), hashes={})
        with self.assertRaisesRegex(ValueError, "requires a full CIF SHA-256"):
            ParentResolver(dataset).resolve("main_union")

    def test_not_available_source_group_is_not_resurrected_from_metadata(self):
        rows = _rows(("a", "b"))
        for row in rows:
            row["source_id"] = "SAME-DISPLAY-ID"
        parents = {
            value: {
                "source_group": "S-SYNTHETIC-%s" % value,
                "source_status": "NOT_AVAILABLE",
                "source_size": "1",
            }
            for value in ("a", "b")
        }
        resolution = ParentResolver(_Dataset(rows, parents)).resolve("main_union")
        self.assertNotEqual(resolution.groups["a"], resolution.groups["b"])

    def test_mofid_v1_conflict_cannot_merge_two_mofid_v2_components(self):
        parents = {
            "a": {"mofid_v2": "M2-A", "mofid_v1": "M1-X"},
            "b": {"mofid_v2": "M2-B", "mofid_v1": "M1-X"},
            "c": {"mofid_v1": "M1-X"},
        }
        resolution = ParentResolver(_Dataset(_rows(("a", "b", "c")), parents)).resolve(
            "priority_main"
        )
        self.assertNotEqual(resolution.groups["a"], resolution.groups["b"])
        self.assertEqual(resolution.exclusions["c"], "PARENT_METHOD_CONFLICT")


class ParentGroupSplitterTests(unittest.TestCase):
    def _singleton_classified(self, count=100, order=None):
        ids = ["s%03d" % index for index in range(count)]
        if order is not None:
            ids = list(order)
        dataset = _Dataset(_rows(ids))
        labels = {
            structure_id: ("CR" if int(structure_id[1:]) % 2 else "NCR")
            for structure_id in ids
        }
        return _Classified(dataset, labels)

    def test_singleton_split_hits_requested_fractions_and_is_unpackable(self):
        result = ParentGroupSplitter(
            self._singleton_classified(), parent_method="none", random_state="stable-seed"
        ).train_valid_test_split((0.6, 0.2, 0.2))
        self.assertEqual(result.counts["train"], 60)
        self.assertEqual(result.counts["validation"], 20)
        self.assertEqual(result.counts["test"], 20)
        train, valid, test = result
        self.assertEqual(train, result.train_indices)
        self.assertEqual(valid, result.valid_indices)
        self.assertEqual(test, result.test_indices)
        self.assertEqual(result.leakage_audit["cross_split_block_count"], 0)
        self.assertEqual(
            result.achieved_fractions,
            {"train": 0.6, "validation": 0.2, "test": 0.2},
        )
        self.assertEqual(result.label_counts_by_split["train"], {"CR": 30, "NCR": 30})

    def test_fraction_contract_accepts_zero_and_rejects_invalid_values(self):
        splitter = ParentGroupSplitter(
            self._singleton_classified(count=10), parent_method="none"
        )
        all_train = splitter.train_valid_test_split((1.0, 0.0, 0.0))
        self.assertEqual(all_train.counts["train"], 10)
        self.assertEqual(all_train.counts["validation"], 0)
        for fractions in (
            (0.8, 0.2),
            (0.8, 0.3, -0.1),
            (0.8, 0.3, 0.0),
            (float("nan"), 0.5, 0.5),
            (0.0, 0.0, 0.0),
        ):
            with self.subTest(fractions=fractions):
                with self.assertRaises(ValueError):
                    splitter.train_valid_test_split(fractions)

    def test_assignment_is_independent_of_metadata_row_order(self):
        ids = ["s%03d" % index for index in range(80)]
        shuffled = list(ids)
        random.Random(991).shuffle(shuffled)
        first = ParentGroupSplitter(
            self._singleton_classified(order=ids),
            parent_method="none",
            random_state="order-test",
        ).train_valid_test_split()
        second = ParentGroupSplitter(
            self._singleton_classified(order=shuffled),
            parent_method="none",
            random_state="order-test",
        ).train_valid_test_split()
        self.assertEqual(first.assignments, second.assignments)
        self.assertEqual(
            first.receipt()["assignment_sha256"],
            second.receipt()["assignment_sha256"],
        )

    def test_unknown_stratification_field_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "Unknown stratify_by"):
            ParentGroupSplitter(
                self._singleton_classified(),
                parent_method="none",
                stratify_by=("label", "sorce_database"),
            )

    def test_external_label_maps_use_the_closed_canonical_vocabulary(self):
        classified = self._singleton_classified(count=4)
        classified.label_by_id["s000"] = "MAYBE"
        with self.assertRaisesRegex(ValueError, "unknown classification label"):
            ParentGroupSplitter(classified, parent_method="none")

    def test_main_union_is_a_guard_not_an_explanatory_parent_method(self):
        classified = self._singleton_classified(count=4)
        with self.assertRaisesRegex(TypeError, "parent method must be a string"):
            ParentResolver(classified.dataset).resolve(None)
        with self.assertRaisesRegex(TypeError, "parent method must be a string"):
            ParentGroupSplitter(classified, parent_method=None)
        with self.assertRaisesRegex(ValueError, "Unknown parent method"):
            ParentGroupSplitter(
                classified, parent_method="main_union"
            )

    def test_split_receipt_expands_parent_and_leakage_terms(self):
        result = ParentGroupSplitter(
            self._singleton_classified(count=8),
            parent_method="priority_main",
            leakage_guard="auto",
            random_state="contract-receipt",
        ).train_valid_test_split((0.5, 0.25, 0.25))
        receipt = result.receipt()
        self.assertEqual(receipt["parent_method"], "priority_main")
        self.assertEqual(
            receipt["parent_method_definition"]["priority_order"],
            ["rac5", "mofid_v2", "mofid_v1"],
        )
        self.assertIn("not a row-by-row", receipt["parent_method_definition"]["summary"])
        self.assertEqual(receipt["leakage_guard"], "main_union")
        self.assertEqual(receipt["requested_leakage_guard"], "auto")
        self.assertEqual(
            receipt["requested_leakage_guard_definition"]["identifier"], "auto"
        )
        self.assertEqual(
            receipt["leakage_guard_definition"]["graph_rule"],
            "connected_components_of_the_transitive_union_of_all_edges",
        )
        self.assertFalse(
            receipt["leakage_guard_definition"]["explanatory_parent_method"]
        )
        # The expanded explanation is ordinary JSON data, not a Python-only proxy.
        json.dumps(receipt, sort_keys=True, allow_nan=False)

        other = ParentGroupSplitter(
            self._singleton_classified(count=8),
            parent_method="none",
            leakage_guard="auto",
        ).train_valid_test_split((0.5, 0.25, 0.25))
        other_receipt = other.receipt()
        self.assertEqual(other_receipt["requested_leakage_guard"], "auto")
        self.assertEqual(other_receipt["leakage_guard"], "parent_only")

    def test_explicit_structure_id_filter_is_reported(self):
        classified = self._singleton_classified(count=10)
        result = ParentGroupSplitter(
            classified,
            parent_method="none",
            structure_ids=("s001", "s003", "s005"),
        ).train_valid_test_split((1.0, 0.0, 0.0))
        self.assertEqual(set(result.train_ids), {"s001", "s003", "s005"})
        self.assertEqual(result.exclusions["s000"], "STRUCTURE_ID_FILTER")
        with self.assertRaisesRegex(KeyError, "Unknown structure_id"):
            ParentGroupSplitter(
                classified,
                parent_method="none",
                structure_ids=("absent",),
            )

    def test_assignment_is_independent_of_python_hash_seed(self):
        script = r'''
import json
from types import SimpleNamespace
from CoREMOF.splitters import ParentGroupSplitter
ids = ["s%03d" % i for i in range(40)]
rows = [{"structure_id": x, "source_database": "COD", "source_id": x,
         "structure_variant": "ASR", "metal_elements": "Cu"} for x in ids]
dataset = SimpleNamespace(metadata_rows=rows, parent_by_id={}, cif_hashes={},
                          dataset_version="v", input_hashes={}, parent_group_methods={})
classified = SimpleNamespace(dataset=dataset,
    label_by_id={x: ("CR" if i % 2 else "NCR") for i, x in enumerate(ids)},
    checker_view="3checker")
result = ParentGroupSplitter(classified, parent_method="none",
                             random_state="hash-seed").train_valid_test_split()
print(json.dumps(dict(result.assignments), sort_keys=True))
'''
        outputs = []
        for seed in ("1", "918273"):
            environment = dict(os.environ)
            environment["PYTHONHASHSEED"] = seed
            outputs.append(
                subprocess.check_output(
                    [sys.executable, "-c", script],
                    cwd=str(Path(__file__).resolve().parents[1]),
                    env=environment,
                    text=True,
                )
            )
        self.assertEqual(outputs[0], outputs[1])

    def test_main_union_guard_blocks_mofid_bridge_across_rac_anchors(self):
        ids = ("a", "b", "c", "d", "e", "f")
        parents = {
            "a": {"rac5": "R1", "mofid_v2": "M-SHARED"},
            "b": {"rac5": "R2", "mofid_v2": "M-SHARED"},
            "c": {"rac5": "R3"},
            "d": {"rac5": "R4"},
            "e": {"rac5": "R5"},
            "f": {"rac5": "R6"},
        }
        classified = _Classified(_Dataset(_rows(ids), parents), {value: "CR" for value in ids})
        result = ParentGroupSplitter(
            classified,
            parent_method="rac5",
            leakage_guard="main_union",
            random_state=7,
        ).train_valid_test_split((0.5, 0.25, 0.25))
        self.assertEqual(result.assignments["a"], result.assignments["b"])
        self.assertEqual(result.leakage_groups["a"], result.leakage_groups["b"])
        self.assertTrue(result.leakage_audit["passed"])

    def test_priority_conflict_row_is_retained_when_main_union_guards_it(self):
        ids = ("a", "b", "bridge", "d", "e")
        parents = {
            "a": {"rac5": "R1", "mofid_v2": "M-SHARED"},
            "b": {"rac5": "R2", "mofid_v2": "M-SHARED"},
            "bridge": {"mofid_v2": "M-SHARED"},
            "d": {"rac5": "R3"},
            "e": {"rac5": "R4"},
        }
        classified = _Classified(
            _Dataset(_rows(ids), parents), {value: "CR" for value in ids}
        )
        result = ParentGroupSplitter(
            classified,
            parent_method="priority_main",
            leakage_guard="main_union",
        ).train_valid_test_split((0.6, 0.2, 0.2))
        self.assertIn("bridge", result.assignments)
        self.assertEqual(
            result.parent_diagnostics["bridge"], "PARENT_METHOD_CONFLICT"
        )
        self.assertNotIn("bridge", result.exclusions)
        self.assertEqual(result.assignments["a"], result.assignments["bridge"])
        self.assertIn("PARENT_METHOD_CONFLICTS_GUARDED", result.warnings)
        self.assertIn("PARENT_METHOD_CONFLICT_LEDGER_PRESENT", result.warnings)
        self.assertEqual(result.receipt()["parent_conflict_count"], 1)
        ledger = result.receipt()["parent_conflicts"][0]
        self.assertEqual(ledger["lower_method"], "mofid_v2")
        self.assertEqual(ledger["unresolved_ids"], ["bridge"])

    def test_filtered_middle_row_still_bridges_selected_rows(self):
        ids = ("a", "bridge", "c", "d", "e", "f")
        parents = {
            "bridge": {"mofid_v2": "M-BC"},
            "c": {"mofid_v2": "M-BC"},
        }
        hashes = _unique_hashes(ids)
        hashes["a"] = hashes["bridge"] = "b" * 64
        labels = {value: "CR" for value in ids}
        labels["bridge"] = "AMBIGUOUS"
        classified = _Classified(_Dataset(_rows(ids), parents, hashes), labels)
        result = ParentGroupSplitter(
            classified,
            parent_method="priority_main",
            leakage_guard="main_union",
            labels=("CR",),
        ).train_valid_test_split((0.5, 0.25, 0.25))
        self.assertEqual(result.assignments["a"], result.assignments["c"])
        self.assertEqual(result.exclusions["bridge"], "LABEL_FILTER")

    def test_prefiltered_classification_absence_is_an_explicit_exclusion(self):
        dataset = _Dataset(_rows(("a", "b", "c")))
        classified = _Classified(dataset, {"a": "CR", "b": "NCR"})
        result = split_release(classified, parent_method="none")
        self.assertEqual(result.exclusions["c"], "PRESELECTION_FILTER")
        self.assertTrue(result.filters["preselection"]["active"])
        self.assertTrue(
            all(
                index < 2
                for index in (
                    *result.train_indices,
                    *result.valid_indices,
                    *result.test_indices,
                )
            )
        )
        self.assertEqual(result.view_index_by_id, {"a": 0, "b": 1})
        self.assertEqual(result.index_by_id["c"], 2)
        self.assertEqual(result.filters["preselection"]["selected_count"], 2)

    def test_split_release_honours_or_rejects_explicit_preclassified_view(self):
        classified = self._singleton_classified(count=10)
        same = split_release(classified, checkers="3checker", parent_method="none")
        self.assertEqual(same.checker_view, "3checker")
        with self.assertRaisesRegex(ValueError, "preclassified dataset uses"):
            split_release(classified, checkers="5checker", parent_method="none")

    def test_split_release_loads_classifies_and_splits_a_path(self):
        from tests.test_dataset_labels import _make_release

        with tempfile.TemporaryDirectory() as directory:
            release_root = Path(directory, "release")
            release_root.mkdir()
            _make_release(release_root)
            result = split_release(
                release_root,
                checkers=3,
                parent_method="none",
                fractions=(0.5, 0.25, 0.25),
            )
        self.assertEqual(result.dataset_version, "vtest")
        self.assertEqual(result.checker_view, "3checker")
        self.assertEqual(sum(result.counts[name] for name in ("train", "validation", "test")), 2)

    def test_receipt_is_atomic_provisional_and_refuses_overwrite(self):
        result = ParentGroupSplitter(
            self._singleton_classified(count=20), parent_method="none"
        ).train_valid_test_split()
        with tempfile.TemporaryDirectory() as directory:
            csv_path, json_path = result.write(directory, stem="example")
            self.assertTrue(csv_path.exists())
            receipt = json.loads(Path(json_path).read_text(encoding="utf-8"))
            self.assertTrue(receipt["provisional_input"])
            self.assertFalse(receipt["official_split"])
            self.assertTrue(receipt["leakage_audit"]["passed"])
            self.assertEqual(receipt["algorithm"]["version"], "1.0")
            self.assertEqual(receipt["implementation"]["split_api_version"], "0.1.0")
            self.assertFalse(receipt["integrity"]["cif_files_verified"])
            self.assertEqual(set(receipt["implementation"]["source_sha256"]), {
                "dataset.py", "labels.py", "parents.py", "splitters.py"
            })
            self.assertIn("PROVISIONAL_PARENT_INPUT", receipt["warnings"])
            self.assertIn("achieved_fractions", receipt)
            self.assertIn("label_counts_by_split", receipt)
            self.assertEqual(receipt["parent_conflict_count"], 0)
            with self.assertRaises(FileExistsError):
                result.write(directory, stem="example")
            result.write(directory, stem="example", overwrite=True)
            with self.assertRaisesRegex(ValueError, "filename stem"):
                result.write(directory, stem="../escape")
        with self.assertRaises(TypeError):
            result.assignments["s000"] = "test"
        with self.assertRaises(TypeError):
            result.filters["labels"] = ["AMBIGUOUS"]

    def test_receipt_freezes_implementation_hashes_at_split_construction(self):
        first_hashes = {name: "a" * 64 for name in (
            "dataset.py", "labels.py", "parents.py", "splitters.py"
        )}
        second_hashes = {name: "b" * 64 for name in first_hashes}
        with patch("CoREMOF.splitters._implementation_hashes", return_value=first_hashes):
            result = ParentGroupSplitter(
                self._singleton_classified(count=10), parent_method="none"
            ).train_valid_test_split()
        with patch("CoREMOF.splitters._implementation_hashes", return_value=second_hashes):
            receipt = result.receipt()
        self.assertEqual(receipt["implementation"]["source_sha256"], first_hashes)

    def test_split_fails_if_imported_implementation_sources_drift(self):
        changed = dict(splitters_module._IMPORTED_BASE_IMPLEMENTATION_HASHES)
        changed["splitters.py"] = "0" * 64
        with patch.object(
            splitters_module,
            "_current_base_implementation_hashes",
            return_value=changed,
        ):
            with self.assertRaisesRegex(
                ValueError, "split implementation source changed after module import"
            ):
                ParentGroupSplitter(
                    self._singleton_classified(count=10), parent_method="none"
                ).train_valid_test_split()

    def test_pair_write_rolls_back_if_second_render_fails(self):
        result = ParentGroupSplitter(
            self._singleton_classified(count=20), parent_method="none"
        ).train_valid_test_split()
        with tempfile.TemporaryDirectory() as directory:
            original_to_json = result.to_json
            try:
                object.__setattr__(
                    result,
                    "to_json",
                    lambda *args, **kwargs: (_ for _ in ()).throw(
                        RuntimeError("synthetic render failure")
                    ),
                )
                with self.assertRaisesRegex(RuntimeError, "synthetic render failure"):
                    result.write(directory, stem="pair")
            finally:
                object.__setattr__(result, "to_json", original_to_json)
            self.assertFalse(Path(directory, "pair.csv").exists())
            self.assertFalse(Path(directory, "pair.json").exists())

    def test_pair_nonoverwrite_rollback_preserves_post_identity_replacement(self):
        result = ParentGroupSplitter(
            self._singleton_classified(count=20), parent_method="none"
        ).train_valid_test_split()
        with tempfile.TemporaryDirectory() as directory_text:
            directory = Path(directory_text)
            final_csv = directory / "pair.csv"
            real_link = splitters_module.os.link
            real_replace = splitters_module.os.replace
            real_samefile = splitters_module.os.path.samefile
            link_count = 0
            identity_checked = False

            def fail_second_pair_publication(source, destination):
                nonlocal link_count
                link_count += 1
                # Two links render the staged CSV/JSON; the third publishes
                # pair.csv and the fourth is the injected pair.json failure.
                if link_count == 4:
                    raise OSError("deterministic second publication failure")
                return real_link(source, destination)

            def replace_after_identity_check(left, right):
                nonlocal identity_checked
                answer = real_samefile(left, right)
                if answer and not identity_checked:
                    identity_checked = True
                    replacement = directory / ".concurrent-replacement"
                    replacement.write_text("foreign generation\n", encoding="utf-8")
                    real_replace(str(replacement), str(final_csv))
                return answer

            with patch.object(
                splitters_module.os,
                "link",
                side_effect=fail_second_pair_publication,
            ), patch.object(
                splitters_module.os.path,
                "samefile",
                side_effect=replace_after_identity_check,
            ):
                with self.assertRaisesRegex(
                    OSError, "deterministic second publication failure"
                ):
                    result.write(directory, stem="pair")

            self.assertTrue(identity_checked)
            self.assertEqual(
                final_csv.read_text(encoding="utf-8"), "foreign generation\n"
            )
            self.assertFalse((directory / "pair.json").exists())
            self.assertFalse(
                any(path.name.startswith(".pair.") for path in directory.iterdir())
            )

    def test_pair_overwrite_failure_removes_new_and_restores_existing_output(self):
        result = ParentGroupSplitter(
            self._singleton_classified(count=20), parent_method="none"
        ).train_valid_test_split()
        with tempfile.TemporaryDirectory() as directory_text:
            directory = Path(directory_text)
            final_csv = directory / "pair.csv"
            final_json = directory / "pair.json"
            old_json = b'{"old": true}\n'
            final_json.write_bytes(old_json)
            real_replace = splitters_module.os.replace
            replace_count = 0

            def fail_second_pair_replacement(source, destination):
                nonlocal replace_count
                replace_count += 1
                if replace_count == 2:
                    raise OSError("deterministic second replacement failure")
                return real_replace(source, destination)

            with patch.object(
                splitters_module.os,
                "replace",
                side_effect=fail_second_pair_replacement,
            ):
                with self.assertRaisesRegex(
                    OSError, "deterministic second replacement failure"
                ):
                    result.write(directory, stem="pair", overwrite=True)

            self.assertFalse(final_csv.exists())
            self.assertEqual(final_json.read_bytes(), old_json)
            self.assertFalse(
                any(path.name.startswith(".pair.") for path in directory.iterdir())
            )

    def test_official_split_request_fails_closed(self):
        with self.assertRaisesRegex(
            OfficialSplitUnavailableError, "No audited official split manifest"
        ):
            ParentGroupSplitter(
                self._singleton_classified(count=10),
                parent_method="none",
                official=True,
            )

    def test_custom_checker_view_is_marked_user_defined_in_receipt(self):
        classified = self._singleton_classified(count=10)
        classified.checker_view = "custom:MOFClassifier+MOSAEC"
        classified.checker_view_official = False
        result = ParentGroupSplitter(
            classified, parent_method="none"
        ).train_valid_test_split()
        self.assertEqual(
            result.receipt()["checker_view_kind"], "USER_DEFINED"
        )

    def test_final_input_state_uses_a_closed_two_artifact_contract(self):
        classified = self._singleton_classified(count=6)
        classified.dataset.parent_group_methods = {"release_status": "FINAL_CANDIDATE"}
        classified.dataset.dataset_info = {"release_status": "FINAL"}
        candidate = ParentGroupSplitter(
            classified, parent_method="none"
        ).train_valid_test_split()
        self.assertTrue(candidate.provisional_input)

        classified.dataset.parent_group_methods = {"release_status": "FINAL"}
        final = ParentGroupSplitter(
            classified, parent_method="none"
        ).train_valid_test_split()
        self.assertFalse(final.provisional_input)


if __name__ == "__main__":
    unittest.main()
