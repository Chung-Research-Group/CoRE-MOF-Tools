import copy
from dataclasses import replace
import json
import hashlib
import os
from pathlib import Path
import random
import pickle
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import CoREMOF.dataset as dataset_module

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
    generated_output_terminology_definitions,
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
        self.parent_group_methods = {}


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
            "mofid_v2": "MOF-SHARED",
            "mofid_v2_status": "SUCCESS",
        }
        for index, structure_id in enumerate(ids)
    ]


def _unique_hashes(ids):
    return {
        structure_id: hashlib.sha256(structure_id.encode("utf-8")).hexdigest()
        for structure_id in ids
    }


def _authorize_crystalnets(dataset, methods=(
    "rac5_crystalnets",
    "mofid_v2_crystalnets",
)):
    """Attach exact terminal contract copies to one synthetic dataset."""

    contracts = {
        method: copy.deepcopy(
            dataset_module._CRYSTALNETS_REFERENCE_CONTRACTS[method]
        )
        for method in methods
    }
    prefixes = {
        method: dataset_module._CRYSTALNETS_REFERENCE_METHODS[method][0]
        for method in methods
    }
    integration = dict(dataset_module._CRYSTALNETS_REFERENCE_INTEGRATION)
    dataset.parent_group_methods = {
        "csv_column_prefixes": prefixes,
        "criteria": copy.deepcopy(contracts),
        "definitions": copy.deepcopy(contracts),
        "crystalnets_reference_integration": copy.deepcopy(integration),
    }
    dataset.dataset_info = {
        "definitions": copy.deepcopy(contracts),
        "parent_grouping": {
            "criterion_contracts": copy.deepcopy(contracts),
            "crystalnets_reference_integration": copy.deepcopy(integration),
        },
    }
    return dataset


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

    def test_generated_output_terminology_closure_is_compact_and_immutable(self):
        definitions = generated_output_terminology_definitions()
        self.assertEqual(
            tuple(definitions),
            (
                "auto",
                "canonical_identifier_text",
                "identity_union",
                "main_union",
                "mofid_v2_crystalnets",
                "parent_only",
                "priority_main",
                "rac5_crystalnets",
                "structure_matcher_strict",
            ),
        )
        self.assertIn("priority_main to main_union", definitions["auto"])
        self.assertIn("parent_only", definitions["auto"])
        self.assertIn("text processing only", definitions["canonical_identifier_text"])
        self.assertIn("provisional source-ID/MOFid transitive groups", definitions["identity_union"])
        self.assertIn("freshly recomputed", definitions["identity_union"])
        self.assertIn("SUCCESS_TOPOLOGY_TIMEOUT", definitions["identity_union"])
        self.assertIn("not a parent or explanatory method", definitions["main_union"])
        self.assertIn("missing or malformed input fails closed", definitions["main_union"])
        self.assertIn("latter two are successful calculated", definitions["mofid_v2_crystalnets"])
        self.assertIn("must be rebuilt", definitions["mofid_v2_crystalnets"])
        self.assertIn("264", definitions["rac5_crystalnets"])
        self.assertIn("Direct edges are authoritative", definitions["structure_matcher_strict"])
        self.assertIn("selected explanatory parent relation", definitions["parent_only"])
        self.assertIn("MOFid-v2", definitions["priority_main"])
        self.assertIn("zero stronger components", definitions["priority_main"])
        self.assertIn("exactly one", definitions["priority_main"])
        self.assertIn("two or more", definitions["priority_main"])
        with self.assertRaises(TypeError):
            definitions["identity_union"] = "changed"

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
        self.assertFalse(
            definition["release_construction"]["earlier_component_or_edge_import"]
        )
        self.assertEqual(
            definition["release_construction"]["eligible_mofid_statuses"],
            (
                "SUCCESS",
                "SUCCESS_TOPOLOGY_UNKNOWN",
                "SUCCESS_TOPOLOGY_ERROR",
                "SUCCESS_TOPOLOGY_TIMEOUT",
            ),
        )
        self.assertIn(
            "freshly recompute", definition["release_construction"]["v26.0.1"]
        )
        self.assertIn(
            "do not seed", definition["release_construction"]["v26.0.2"]
        )
        self.assertIn("no precedence", definition["algorithm"])
        self.assertIn("common_nulls_never_match", definition["missing_behavior"])
        self.assertIn("RAC5", definition["excluded_inputs"])
        self.assertIn("excluded", definition["relation_to_other_terms"])
        canonical = definition["canonicalized_identifier_text"]
        self.assertEqual(canonical, CANONICALIZED_IDENTIFIER_TEXT_DEFINITION)
        self.assertEqual(
            canonical["current_release_text_steps_in_order"],
            (
                "convert_the_release_authorized_current_value_to_text",
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
        self.assertNotIn("inherited_v2601_identity_component_cleanup", canonical)
        self.assertIn("freshly_recomputed", canonical["source_key"]["scope"])
        self.assertIn("no_prior_component", canonical["source_key"]["scope"])
        self.assertIn(
            "no_CIF_coordinate_canonicalization", canonical["excluded_operations"]
        )

    def test_optional_prefix_contracts_define_meaning_and_failure_behavior(self):
        expected = {
            "rac5_crystalnets": ("RT-", "rac_crystalnets"),
            "mofid_v2_crystalnets": ("M2T-", "mofid2_crystalnets"),
            "structure_matcher_strict": ("SM-", "sm"),
        }
        for method, (prefix, column_prefix) in expected.items():
            with self.subTest(method=method):
                definition = parent_method_definition(method)
                self.assertIn(prefix, definition["group_prefix_meaning"])
                self.assertTrue(definition["missing_behavior"])
                self.assertIn("excluded", definition["relation_to_other_terms"])
                self.assertEqual(
                    definition["input_fields"][:3],
                    (
                        "parent_groups.{}_status".format(column_prefix),
                        "parent_groups.{}_group".format(column_prefix),
                        "parent_groups.{}_size".format(column_prefix),
                    ),
                )
                if method != "structure_matcher_strict":
                    self.assertEqual(len(definition["input_fields"]), 3)
        for method in ("rac5_crystalnets", "mofid_v2_crystalnets"):
            self.assertIn(method, PRIORITY_MAIN_DEFINITION["excluded_inputs"])
            self.assertIn(method, MAIN_UNION_DEFINITION["excluded_inputs"])
        self.assertEqual(
            PRIORITY_MAIN_DEFINITION["priority_order"],
            ("rac5", "mofid_v2", "mofid_v1"),
        )
        crystalnets = parent_method_definition("rac5_crystalnets")[
            "crystalnets_fingerprint"
        ]
        self.assertEqual(crystalnets, CURRENT_CRYSTALNETS_FINGERPRINT_DEFINITION)
        for field in (
            "network_dimension",
            "catenation_degree",
            "single_node_net",
            "all_node_net",
            "for_each_subnet_SingleNodes_status_dimension_topology_key_topology_name_"
            "topological_genome",
        ):
            self.assertIn(field, crystalnets["included_fields"])
        self.assertIn("runtime_seconds", crystalnets["excluded_fields"])
        self.assertTrue(crystalnets["not_a_topology_similarity_tolerance"])
        mofid_crystalnets = parent_method_definition("mofid_v2_crystalnets")
        self.assertEqual(
            mofid_crystalnets["eligible_mofid_v2_statuses"],
            (
                "SUCCESS",
                "SUCCESS_TOPOLOGY_UNKNOWN",
                "SUCCESS_TOPOLOGY_ERROR",
                "SUCCESS_TOPOLOGY_TIMEOUT",
            ),
        )
        self.assertIn(
            "every_other_MOFid-v2_status_adds_no_edge",
            mofid_crystalnets["ineligible_mofid_v2_status_behavior"],
        )
        self.assertEqual(
            mofid_crystalnets["mofid_v2_crystalnets_rebuild_trigger"],
            "REBUILD_M2T_IF_AUTHORIZED_MOFID_V2_VALUES_CHANGE",
        )
        self.assertEqual(
            mofid_crystalnets["crystalnets_fingerprint"],
            CURRENT_CRYSTALNETS_FINGERPRINT_DEFINITION,
        )
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
        self.assertIn("authoritative", matcher["direct_edge_authority"])
        self.assertIn("convenience", matcher["component_interpretation"])
        self.assertIn("not_duplicate_proof", matcher["component_interpretation"])
        self.assertIn(
            "NOT_AVAILABLE_rather_than_unmatched", matcher["failure_projection"]
        )

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
        self.assertTrue(definition["not_feature_recalculation_queue"])
        self.assertIn("does not rank, schedule, or recalculate", definition["purpose"])
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

    def test_optional_crystalnets_criterion_must_be_declared_by_release(self):
        resolver = ParentResolver(_Dataset(_rows(("a", "b"))))
        with self.assertRaisesRegex(ValueError, "not present in this release"):
            resolver.resolve("rac5_crystalnets")

    def test_crystalnets_reference_requires_full_release_triads_and_prefixes(self):
        cases = (
            (
                {"a": {"mofid_v2_crystalnets": "M2T-FAKE0001"}},
                "not present in this release",
            ),
            (
                {"a": {"mofid2_crystalnets_group": "M2T-ABCD0001"}},
                "complete status/group/size triad",
            ),
            (
                {
                    "a": {
                        "mofid2_crystalnets_status": "UNMATCHED",
                        "mofid2_crystalnets_group": "R-WRONG0001",
                        "mofid2_crystalnets_size": "1",
                    }
                },
                "must use M2T-",
            ),
            (
                {
                    "a": {
                        "mofid2_crystalnets_status": "ERROR",
                        "mofid2_crystalnets_group": "M2T-ABCD0001",
                        "mofid2_crystalnets_size": "1",
                    }
                },
                "MATCHED, UNMATCHED, or NOT_AVAILABLE",
            ),
        )
        for parents, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    ParentResolver(_Dataset(_rows(("a",)), parents)).resolve(
                        "mofid_v2_crystalnets"
                    )

        partial_rows = {
            "a": {
                "rac_crystalnets_status": "UNMATCHED",
                "rac_crystalnets_group": "RT-ABCD0001",
                "rac_crystalnets_size": "1",
            },
            "b": {},
        }
        with self.assertRaisesRegex(ValueError, "triad for every structure"):
            ParentResolver(_Dataset(_rows(("a", "b")), partial_rows))

        for method, field, prefix in (
            ("rac5_crystalnets", "rac_crystalnets", "RT-"),
            ("mofid_v2_crystalnets", "mofid2_crystalnets", "M2T-"),
        ):
            for suffix in (
                "A" * 7,
                "A" * 65,
                "abcdef12",
                "ABCDEFGZ",
            ):
                with self.subTest(method=method, suffix=suffix):
                    parents = {
                        "a": {
                            field + "_status": "UNMATCHED",
                            field + "_group": prefix + suffix,
                            field + "_size": "1",
                        }
                    }
                    with self.assertRaisesRegex(
                        ValueError, "8 to 64 uppercase hexadecimal"
                    ):
                        ParentResolver(_Dataset(_rows(("a",)), parents))

    def test_crystalnets_reference_requires_exact_release_contract_copies(self):
        parents = {
            "a": {
                "mofid2_crystalnets_status": "UNMATCHED",
                "mofid2_crystalnets_group": "M2T-DEADBEEF",
                "mofid2_crystalnets_size": "1",
            }
        }
        dataset = _Dataset(_rows(("a",)), parents)
        with self.assertRaisesRegex(ValueError, "exact validated release"):
            ParentResolver(dataset).resolve("mofid_v2_crystalnets")

        _authorize_crystalnets(dataset, ("mofid_v2_crystalnets",))
        resolution = ParentResolver(dataset).resolve("mofid_v2_crystalnets")
        self.assertEqual(resolution.groups["a"], "M2T-DEADBEEF")

        surfaces = (
            dataset.parent_group_methods["criteria"],
            dataset.parent_group_methods["definitions"],
            dataset.dataset_info["definitions"],
            dataset.dataset_info["parent_grouping"]["criterion_contracts"],
        )
        for surface in surfaces:
            with self.subTest(surface=id(surface)):
                decision = surface["mofid_v2_crystalnets"]["decision_contract"]
                original = decision["publication_authorized"]
                decision["publication_authorized"] = True
                with self.assertRaisesRegex(ValueError, "exact validated release"):
                    ParentResolver(dataset).resolve("mofid_v2_crystalnets")
                decision["publication_authorized"] = original

        for integration in (
            dataset.parent_group_methods["crystalnets_reference_integration"],
            dataset.dataset_info["parent_grouping"][
                "crystalnets_reference_integration"
            ],
        ):
            with self.subTest(integration=id(integration)):
                integration["publication_authorized"] = True
                with self.assertRaisesRegex(ValueError, "exact validated release"):
                    ParentResolver(dataset).resolve("mofid_v2_crystalnets")
                integration["publication_authorized"] = False

    def test_m2t_release_triads_require_eligible_explicit_mofid_status(self):
        parents = {
            "a": {
                "mofid2_crystalnets_status": "UNMATCHED",
                "mofid2_crystalnets_group": "M2T-DEADBEEF",
                "mofid2_crystalnets_size": "1",
            }
        }
        for value in (
            None,
            "ERROR",
            "TIMEOUT",
            "NOT_AVAILABLE",
            "success",
            " SUCCESS",
        ):
            with self.subTest(ineligible=value):
                rows = _rows(("a",))
                if value is None:
                    rows[0].pop("mofid_v2_status")
                else:
                    rows[0]["mofid_v2_status"] = value
                dataset = _authorize_crystalnets(
                    _Dataset(rows, copy.deepcopy(parents)),
                    ("mofid_v2_crystalnets",),
                )
                with self.assertRaisesRegex(ValueError, "eligible mofid_v2_status"):
                    ParentResolver(dataset).resolve("mofid_v2_crystalnets")

        for value in sorted(dataset_module.M2T_ELIGIBLE_MOFID_V2_STATUSES):
            with self.subTest(eligible=value):
                rows = _rows(("a",))
                rows[0]["mofid_v2_status"] = value
                dataset = _authorize_crystalnets(
                    _Dataset(rows, copy.deepcopy(parents)),
                    ("mofid_v2_crystalnets",),
                )
                self.assertEqual(
                    ParentResolver(dataset)
                    .resolve("mofid_v2_crystalnets")
                    .groups["a"],
                    "M2T-DEADBEEF",
                )

    def test_m2t_release_triads_bind_complete_canonical_mofid_values(self):
        singleton = {
            "a": {
                "mofid2_crystalnets_status": "UNMATCHED",
                "mofid2_crystalnets_group": "M2T-DEADBEEF",
                "mofid2_crystalnets_size": "1",
            }
        }
        for value in (
            None,
            "",
            "   ",
            "NOT_AVAILABLE",
            "error",
            False,
            1,
            [],
            {},
            b"MOF-A",
        ):
            with self.subTest(invalid=value):
                rows = _rows(("a",))
                if value is None:
                    rows[0].pop("mofid_v2")
                else:
                    rows[0]["mofid_v2"] = value
                dataset = _authorize_crystalnets(
                    _Dataset(rows, copy.deepcopy(singleton)),
                    ("mofid_v2_crystalnets",),
                )
                with self.assertRaisesRegex(
                    ValueError, "exact string|complete nonplaceholder mofid_v2"
                ):
                    ParentResolver(dataset).resolve("mofid_v2_crystalnets")

        matched = {
            structure_id: {
                "mofid2_crystalnets_status": "MATCHED",
                "mofid2_crystalnets_group": "M2T-ABCD0001",
                "mofid2_crystalnets_size": "2",
            }
            for structure_id in ("a", "b")
        }
        rows = _rows(("a", "b"))
        rows[0]["mofid_v2"] = "MOF-A"
        rows[1]["mofid_v2"] = "MOF-B"
        dataset = _authorize_crystalnets(
            _Dataset(rows, copy.deepcopy(matched)),
            ("mofid_v2_crystalnets",),
        )
        with self.assertRaisesRegex(ValueError, "conflicting canonical mofid_v2"):
            ParentResolver(dataset).resolve("mofid_v2_crystalnets")

        rows[0]["mofid_v2"] = "  ＭＯＦ-A  "
        rows[1]["mofid_v2"] = "mof-a"
        dataset = _authorize_crystalnets(
            _Dataset(rows, copy.deepcopy(matched)),
            ("mofid_v2_crystalnets",),
        )
        resolution = ParentResolver(dataset).resolve("mofid_v2_crystalnets")
        self.assertEqual(resolution.groups["a"], resolution.groups["b"])

    def test_retired_and_reserved_keys_are_closed_recursively_for_generic_inputs(self):
        for key in (
            "rac5_topology",
            "ＲＡＣ５＿Ｔｏｐｏｌｏｇｙ",
            "mofid_v2_topology_group",
            "_authority_generation_marker",
            "_official_checker_view_token",
        ):
            with self.subTest(parent_key=key):
                dataset = _Dataset(
                    _rows(("a",)),
                    {"a": {"rac5": {"nested": {key: True}}}},
                )
                with self.assertRaisesRegex(
                    ValueError, "retired or reserved declaration key"
                ):
                    ParentResolver(dataset)

        for key in ("mofid_v2_topology", "_validated_release_token"):
            with self.subTest(metadata_key=key):
                rows = _rows(("a",))
                rows[0]["nested"] = {key: "forged"}
                with self.assertRaisesRegex(
                    ValueError, "retired or reserved declaration key"
                ):
                    ParentResolver(_Dataset(rows))

        authority_keys = (
            "publication_authorized",
            "PUBLICATION.AUTHORIZED",
            "Publication-Authorized",
            "Ｐｕｂｌｉｃａｔｉｏｎ＿Ａｕｔｈｏｒｉｚｅｄ",
            "publication_status",
            "Publication-Status",
            "publication_authorization",
            "Publication Authorization",
            "publication_ready",
            "Publication-Ready",
            "publishable",
            "authoritative",
            "official_split",
            "Official-Split",
            "release_authority",
            "Release Authority",
            "release_status",
            "Release-Status",
            "release_state",
            "Release-State",
            "candidate_status",
            "candidate_state",
            "evidence_state",
            "priority_main_changed",
            "main_union_changed",
            "checker_view_official",
            "Checker.View-Official",
        )
        authority_claims = tuple(
            (key, value)
            for key in authority_keys
            for value in (
                True,
                False,
                "FINAL",
                "PROVISIONAL_LATEST_AUDITED_SNAPSHOT",
            )
        )
        for surface in ("metadata", "parent_group_methods", "dataset_info"):
            for key, value in authority_claims:
                with self.subTest(surface=surface, key=key, value=value):
                    rows = _rows(("a",))
                    dataset = _Dataset(rows)
                    if surface == "metadata":
                        rows[0]["nested"] = {key: value}
                    else:
                        setattr(dataset, surface, {"nested": {key: value}})
                    with self.assertRaisesRegex(
                        ValueError, "retired or reserved declaration key"
                    ):
                        ParentResolver(dataset)

        for declaration_name in ("parent_group_methods", "dataset_info"):
            for key in ("rac_topology_size", "_authority_locked"):
                with self.subTest(declaration=declaration_name, key=key):
                    dataset = _Dataset(_rows(("a",)))
                    setattr(dataset, declaration_name, {"x": {"y": {key: True}}})
                    with self.assertRaisesRegex(
                        ValueError, "retired or reserved declaration key"
                    ):
                        ParentResolver(dataset)

        staged = _Dataset(_rows(("a",)))
        staged.dataset_info = {
            "release_status": "PROVISIONAL_LATEST_AUDITED_SNAPSHOT"
        }
        self.assertEqual(
            ParentResolver(staged).resolve("none").groups["a"],
            "SINGLETON:a",
        )

        wrapped = _Dataset(_rows(("a",)))
        wrapped.parent_group_methods = {
            "criteria": {
                "rac5_crystalnets": [
                    {"release_state": dataset_module._RT_RELEASE_STATE}
                ]
            }
        }
        with self.assertRaisesRegex(
            ValueError, "retired or reserved declaration key"
        ):
            ParentResolver(wrapped)

    def test_release_triads_reject_legacy_and_conflicting_aliases(self):
        base = {
            "rac_status": "UNMATCHED",
            "rac_group": "R-GOOD0001",
            "rac_size": "1",
        }
        cases = (
            ({"rac5": "EVIL"}, "mixes a validated triad"),
            ({"rac5_group": "R-EVIL0001"}, "group aliases conflict"),
            ({"rac5_status": "MATCHED"}, "status aliases conflict"),
            ({"rac5_size": "2"}, "size aliases conflict|must have size 1"),
        )
        for addition, message in cases:
            with self.subTest(addition=addition):
                entry = dict(base)
                entry.update(addition)
                with self.assertRaisesRegex(ValueError, message):
                    ParentResolver(_Dataset(_rows(("a",)), {"a": entry}))

    def test_explicit_missing_parent_override_does_not_use_truthiness(self):
        resolver = ParentResolver(_Dataset(_rows(("a",))))
        for value in (False, 0, ""):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "missing_parent"):
                    resolver.resolve("rac5", missing_parent=value)

    def test_old_topology_method_names_are_not_selectable_aliases(self):
        resolver = ParentResolver(_Dataset(_rows(("a", "b"))))
        for method in (
            "rac5_topology",
            "mofid_v2_topology",
            "rac_topology",
            "mofid2_topology",
        ):
            with self.subTest(method=method):
                with self.assertRaisesRegex(ValueError, "[Uu]nknown parent method"):
                    resolver.resolve(method)
                with self.assertRaisesRegex(
                    ValueError, "not an explanatory parent method"
                ):
                    parent_method_definition(method)

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

    def test_crystalnets_combined_criteria_are_selectable_but_not_in_priority_main(self):
        ids = ("a", "b", "c")
        parents = {
            "a": {
                "rac_crystalnets_group": "RT-AAAA0001",
                "rac_crystalnets_status": "MATCHED",
                "rac_crystalnets_size": "2",
                "mofid2_crystalnets_group": "M2T-BBBB0001",
                "mofid2_crystalnets_status": "MATCHED",
                "mofid2_crystalnets_size": "2",
            },
            "b": {
                "rac_crystalnets_group": "RT-AAAA0001",
                "rac_crystalnets_status": "MATCHED",
                "rac_crystalnets_size": "2",
                "mofid2_crystalnets_group": "M2T-BBBB0001",
                "mofid2_crystalnets_status": "MATCHED",
                "mofid2_crystalnets_size": "2",
            },
            "c": {
                "rac_crystalnets_group": "RT-CCCC0001",
                "rac_crystalnets_status": "NOT_AVAILABLE",
                "rac_crystalnets_size": "1",
                "mofid2_crystalnets_group": "M2T-DDDD0001",
                "mofid2_crystalnets_status": "NOT_AVAILABLE",
                "mofid2_crystalnets_size": "1",
            },
        }
        resolver = ParentResolver(_authorize_crystalnets(_Dataset(_rows(ids), parents)))
        rac_crystalnets = resolver.resolve("rac5_crystalnets")
        mofid_crystalnets = resolver.resolve("mofid_v2_crystalnets")
        self.assertEqual(
            rac_crystalnets.groups["a"], rac_crystalnets.groups["b"]
        )
        self.assertEqual(
            mofid_crystalnets.groups["a"], mofid_crystalnets.groups["b"]
        )
        self.assertEqual(rac_crystalnets.groups["c"], "SINGLETON:c")
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

    def test_parent_resolver_subset_requires_an_ordered_string_container(self):
        resolver = ParentResolver(_Dataset(_rows(("a", "b"))))
        self.assertEqual(resolver.resolve("none", structure_ids="a").universe_ids, ("a",))
        self.assertEqual(
            resolver.resolve("none", structure_ids=["b", "a"]).universe_ids,
            ("a", "b"),
        )
        for value in (
            {"a", "b"},
            frozenset(("a", "b")),
            {"a": True},
            (item for item in ("a", "b")),
            b"ab",
        ):
            with self.subTest(container=type(value).__name__):
                with self.assertRaisesRegex(TypeError, "ordered list/tuple"):
                    resolver.resolve("none", structure_ids=value)

    def test_legacy_direct_parent_groups_require_nonblank_strings(self):
        for value in (True, False, 1, 1.5, [], {}, b"R1", ""):
            for wrapped in (value, {"group": value}):
                with self.subTest(value=wrapped):
                    with self.assertRaisesRegex(ValueError, "nonblank string"):
                        ParentResolver(
                            _Dataset(_rows(("a",)), {"a": {"rac5": wrapped}})
                        ).resolve("rac5")

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

    def test_release_parent_size_rejects_boolean_and_nonintegral_numeric_types(self):
        for value in (
            True,
            False,
            1.0,
            1.5,
            "1.0",
            "+1",
            "1e0",
            "١",
            "01",
        ):
            with self.subTest(value=value):
                dataset = _Dataset(
                    _rows(("a",)),
                    {
                        "a": {
                            "rac_group": "R1",
                            "rac_status": "UNMATCHED",
                            "rac_size": value,
                        }
                    },
                )
                with self.assertRaisesRegex(ValueError, "positive integer"):
                    ParentResolver(dataset)

    def test_release_parent_size_accepts_integer_and_exact_decimal_text(self):
        for value in (1, "1"):
            with self.subTest(value=value):
                dataset = _Dataset(
                    _rows(("a",)),
                    {
                        "a": {
                            "rac_group": "R1",
                            "rac_status": "UNMATCHED",
                            "rac_size": value,
                        }
                    },
                )
                resolution = ParentResolver(dataset).resolve("rac5")
                self.assertEqual(resolution.groups["a"], "R1")
        with self.assertRaisesRegex(ValueError, "positive integer"):
            ParentResolver(
                _Dataset(
                    _rows(("a",)),
                    {"a": {"rac_group": "R1", "rac_status": "UNMATCHED", "rac_size": " 1 "}},
                )
            )

    def test_release_parent_group_requires_a_nonmissing_string_when_available(self):
        for value in (True, False, 1, 1.5, [], {}, None, "", "  ", "NA", "null"):
            with self.subTest(value=value):
                dataset = _Dataset(
                    _rows(("a",)),
                    {
                        "a": {
                            "rac_group": value,
                            "rac_status": "UNMATCHED",
                            "rac_size": "1",
                        }
                    },
                )
                with self.assertRaisesRegex(
                    ValueError, "nonblank string|missing placeholder"
                ):
                    ParentResolver(dataset)

    def test_release_parent_status_and_group_are_exact_canonical_text(self):
        cases = (
            ({"rac_status": "matched", "rac_group": "R1"}, "status must be"),
            ({"rac_status": " MATCHED", "rac_group": "R1"}, "status must be"),
            ({"rac_status": "UNMATCHED", "rac_group": " R1"}, "exact canonical"),
            ({"rac_status": "UNMATCHED", "rac_group": "R1 "}, "exact canonical"),
        )
        for change, message in cases:
            with self.subTest(change=change):
                entry = {
                    "rac_status": "UNMATCHED",
                    "rac_group": "R1",
                    "rac_size": "1",
                }
                entry.update(change)
                with self.assertRaisesRegex(ValueError, message):
                    ParentResolver(_Dataset(_rows(("a",)), {"a": entry}))

    def test_identical_duplicate_release_aliases_are_accepted(self):
        entry = {
            "rac_status": "UNMATCHED",
            "rac5_status": "UNMATCHED",
            "rac_group": "R1",
            "rac5_group": "R1",
            "rac_size": "1",
            "rac5_size": 1,
        }
        resolution = ParentResolver(
            _Dataset(_rows(("a",)), {"a": entry})
        ).resolve("rac5")
        self.assertEqual(resolution.groups["a"], "R1")

    def test_metadata_structure_ids_must_be_nonblank_strings(self):
        for value in (True, False, 1, 1.5, [], {}, b"a", None, "", " a", "a "):
            with self.subTest(value=value):
                row = {
                    "structure_id": value,
                    "source_database": "COD",
                    "source_id": "SRC-a",
                }
                dataset = _Dataset((row,), hashes={})
                with self.assertRaisesRegex(
                    ValueError, "exact nonblank string structure_id"
                ):
                    ParentResolver(dataset)

    def test_parent_and_cif_tables_require_top_level_mappings(self):
        for value in (True, 1, 1.5, ["bad"], ("bad",), "bad", b"bad"):
            with self.subTest(table="parent_by_id", value=value):
                dataset = _Dataset(_rows(("a",)))
                dataset.parent_by_id = value
                with self.assertRaisesRegex(TypeError, "parent_by_id must be a mapping"):
                    ParentResolver(dataset)
            with self.subTest(table="cif_hashes", value=value):
                dataset = _Dataset(_rows(("a",)))
                dataset.cif_hashes = value
                with self.assertRaisesRegex(TypeError, "cif_hashes must be a mapping"):
                    ParentResolver(dataset)

    def test_explicit_per_id_parent_entries_fail_closed_except_legacy_identity_text(self):
        for value in (True, 1, 1.5, [], (), b"bad"):
            with self.subTest(value=value):
                dataset = _Dataset(_rows(("a",)), {"a": value})
                with self.assertRaisesRegex(TypeError, "parent_by_id entry for a"):
                    ParentResolver(dataset)
        for value in ("", " ", " ID", "ID "):
            with self.subTest(value=value):
                dataset = _Dataset(_rows(("a",)), {"a": value})
                with self.assertRaisesRegex(ValueError, "exact nonblank string"):
                    ParentResolver(dataset)

        resolver = ParentResolver(_Dataset(_rows(("a",)), {"a": "IDENTITY-A"}))
        self.assertEqual(resolver.resolve("identity_union").groups["a"], "IDENTITY-A")
        self.assertEqual(resolver.resolve("rac5").groups["a"], "SINGLETON:a")

    def test_absent_none_and_empty_parent_tables_mean_missing_evidence(self):
        for value in (None, {}):
            with self.subTest(value=value):
                dataset = _Dataset(_rows(("a",)))
                dataset.parent_by_id = value
                self.assertEqual(
                    ParentResolver(dataset).resolve("rac5").groups["a"],
                    "SINGLETON:a",
                )
        dataset = _Dataset(_rows(("a",)))
        del dataset.parent_by_id
        self.assertEqual(
            ParentResolver(dataset).resolve("rac5").groups["a"], "SINGLETON:a"
        )

    def test_method_major_parent_tables_fail_closed_on_malformed_containers(self):
        for value in (True, False, 1, 1.5, ["bad"], ("bad",), "bad", b"bad"):
            with self.subTest(value=value):
                dataset = _Dataset(_rows(("a",)), {"rac5": value})
                with self.assertRaisesRegex(
                    TypeError, "method-major parent table rac5 must be a mapping"
                ):
                    ParentResolver(dataset)
        for value in (None, {}):
            with self.subTest(legitimate_absence=value):
                resolution = ParentResolver(
                    _Dataset(_rows(("a",)), {"rac5": value})
                ).resolve("rac5")
                self.assertEqual(resolution.groups["a"], "SINGLETON:a")

    def test_parent_and_hash_tables_reject_unknown_or_ambiguous_keys(self):
        cases = (
            ({"typo": {"rac5": "R"}}, "unknown top-level"),
            ({"rac5": {"typo": "R"}}, "unknown structure IDs"),
            ({"unknown_method": {"a": "R"}}, "unknown top-level"),
            ({1: {"a": "R"}}, "non-string declaration key"),
        )
        for table, message in cases:
            with self.subTest(parent_by_id=table):
                with self.assertRaisesRegex(ValueError, message):
                    ParentResolver(_Dataset(_rows(("a",)), table))

        ambiguous = _Dataset(_rows(("rac5",)), {"rac5": {}})
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            ParentResolver(ambiguous)

        for hashes in ({"typo": "a" * 64}, {1: "a" * 64}):
            with self.subTest(cif_hashes=hashes):
                dataset = _Dataset(_rows(("a",)), hashes=hashes)
                with self.assertRaisesRegex(ValueError, "unknown structure IDs"):
                    ParentResolver(dataset)

    def test_main_union_source_namespace_fields_reject_nonstring_values(self):
        for field in ("source_database", "source_id"):
            for value in (True, 1, 1.5, [], {}, b"COD", " COD", "COD "):
                with self.subTest(field=field, value=value):
                    row = _rows(("a",))[0]
                    row[field] = value
                    dataset = _Dataset((row,))
                    with self.assertRaisesRegex(
                        ValueError, "%s must be (?:a string|an exact)" % field
                    ):
                        ParentResolver(dataset).resolve("main_union")

    def test_main_union_missing_source_identity_adds_no_source_edge(self):
        rows = _rows(("a", "b"))
        rows[0]["source_database"] = ""
        rows[0]["source_id"] = "SHARED"
        rows[1]["source_database"] = "COD"
        rows[1]["source_id"] = "SHARED"
        resolution = ParentResolver(
            _Dataset(rows, hashes=_unique_hashes(("a", "b")))
        ).resolve("main_union")
        self.assertNotEqual(resolution.groups["a"], resolution.groups["b"])

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

    def test_main_union_rejects_noncanonical_or_nonstring_cif_hashes(self):
        valid = "a" * 64
        for value in (
            "abc123",
            valid.upper(),
            " " + valid,
            valid + " ",
            True,
            1,
            b"a" * 64,
        ):
            with self.subTest(value=value):
                dataset = _Dataset(_rows(("a",)), hashes={"a": value})
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
        with self.assertRaisesRegex(ValueError, "not booleans"):
            splitter.train_valid_test_split((True, False, False))
        for fractions in (
            ("0.8", "0.1", "0.1"),
            {0.8, 0.1},
            frozenset((0.8, 0.1)),
            {"train": 0.8, "validation": 0.1, "test": 0.1},
            (value for value in (0.8, 0.1, 0.1)),
            b"fractions",
        ):
            with self.subTest(noncanonical=type(fractions).__name__):
                with self.assertRaisesRegex(
                    TypeError, "ordered list/tuple|non-boolean numeric"
                ):
                    splitter.train_valid_test_split(fractions)

    def test_random_state_accepts_only_stable_integer_or_string_forms(self):
        classified = self._singleton_classified(count=4)
        self.assertEqual(
            ParentGroupSplitter(
                classified, parent_method="none", random_state=-7
            ).random_state,
            "-7",
        )
        self.assertEqual(
            ParentGroupSplitter(
                classified, parent_method="none", random_state="stable"
            ).random_state,
            "stable",
        )
        for value in (True, False, 1.5, [], {}, None, object()):
            with self.subTest(value=type(value).__name__):
                with self.assertRaisesRegex(TypeError, "non-boolean integer or string"):
                    ParentGroupSplitter(
                        classified, parent_method="none", random_state=value
                    )
        with self.assertRaisesRegex(ValueError, "non-empty"):
            ParentGroupSplitter(classified, parent_method="none", random_state="")

    def test_split_filter_elements_must_be_exact_nonblank_strings(self):
        classified = self._singleton_classified(count=4)
        for name in (
            "labels",
            "sources",
            "variants",
            "metals",
            "structure_ids",
            "required_targets",
        ):
            for value in ((True,), (1,), (1.5,), ([],), ({},), (" x",), ("x ",)):
                with self.subTest(name=name, value=value):
                    with self.assertRaisesRegex(
                        (TypeError, ValueError), "strings|exact nonblank"
                    ):
                        ParentGroupSplitter(
                            classified,
                            parent_method="none",
                            **{name: value},
                        )
            for value in (
                {"x", "y"},
                frozenset(("x", "y")),
                {"x": True},
                (item for item in ("x", "y")),
                b"xy",
            ):
                with self.subTest(name=name, container=type(value).__name__):
                    with self.assertRaisesRegex(TypeError, "ordered list/tuple"):
                        ParentGroupSplitter(
                            classified,
                            parent_method="none",
                            **{name: value},
                        )

    def test_split_receipt_boolean_and_mapping_claims_require_exact_types(self):
        cases = (
            ("cif_files_verified", 1, "cif_files_verified must be a boolean"),
            ("cif_files_verified", "false", "cif_files_verified must be a boolean"),
            ("input_hashes", False, "input_hashes must be a mapping"),
            ("input_hashes", [], "input_hashes must be a mapping"),
            ("dataset_info", False, "dataset_info must be a mapping"),
            ("parent_group_methods", [], "parent_group_methods must be a mapping"),
        )
        for attribute, value, message in cases:
            with self.subTest(attribute=attribute, value=value):
                classified = self._singleton_classified(count=4)
                setattr(classified.dataset, attribute, value)
                with self.assertRaisesRegex((TypeError, ValueError), message):
                    ParentGroupSplitter(
                        classified, parent_method="none"
                    ).train_valid_test_split((1.0, 0.0, 0.0))

        for value in (False, 0, "", []):
            with self.subTest(selection_filters=value):
                classified = self._singleton_classified(count=4)
                classified.selection_filters = value
                with self.assertRaisesRegex(
                    TypeError, "selection_filters must be a mapping"
                ):
                    ParentGroupSplitter(classified, parent_method="none")

        for value in (False, 0, "", {}, {"target"}, [True], [" target"]):
            with self.subTest(target_columns=value):
                classified = self._singleton_classified(count=4)
                classified.dataset.target_columns = value
                with self.assertRaisesRegex(
                    (TypeError, ValueError), "target_columns"
                ):
                    ParentGroupSplitter(classified, parent_method="none")

        classified = self._singleton_classified(count=4)
        classified.selection_filters = {
            "steps": ({"sources": ("COD",)},),
            "selected_count": 4,
        }
        splitter = ParentGroupSplitter(classified, parent_method="none")
        self.assertEqual(
            splitter._selection_filters,
            {"steps": [{"sources": ["COD"]}], "selected_count": 4},
        )
        for value in ({1: "x"}, {"x": object()}, {"x": float("nan")}):
            with self.subTest(selection_filters_nested=value):
                classified = self._singleton_classified(count=4)
                classified.selection_filters = value
                with self.assertRaisesRegex(
                    TypeError, "selection_filters"
                ):
                    ParentGroupSplitter(classified, parent_method="none")

        for value in (1, "false", None):
            with self.subTest(checker_view_official=value):
                classified = self._singleton_classified(count=4)
                classified.checker_view_official = value
                with self.assertRaisesRegex(
                    TypeError, "checker_view_official must be a boolean"
                ):
                    ParentGroupSplitter(
                        classified, parent_method="none"
                    ).train_valid_test_split((1.0, 0.0, 0.0))

        for value in (True, 1, 1.5, [], {}, "", " v-test", "v-test "):
            with self.subTest(dataset_version=value):
                classified = self._singleton_classified(count=4)
                classified.dataset.dataset_version = value
                with self.assertRaisesRegex(
                    TypeError, "dataset_version must be an exact nonblank string"
                ):
                    ParentGroupSplitter(
                        classified, parent_method="none"
                    ).train_valid_test_split((1.0, 0.0, 0.0))

    def test_split_input_hash_keys_and_values_are_not_string_coerced(self):
        for hashes in ({1: "a"}, {"metadata.csv": True}, {"": "a"}, {"a": " "}):
            with self.subTest(hashes=hashes):
                classified = self._singleton_classified(count=4)
                classified.dataset.input_hashes = hashes
                with self.assertRaisesRegex(TypeError, "exact nonblank strings"):
                    ParentGroupSplitter(
                        classified, parent_method="none"
                    ).train_valid_test_split((1.0, 0.0, 0.0))

    def test_splitter_rejects_nonstring_metadata_structure_ids(self):
        for value in (True, 1, 1.5, [], {}, b"a", " a", "a "):
            with self.subTest(value=value):
                row = {
                    "structure_id": value,
                    "source_database": "COD",
                    "source_id": "SRC-a",
                    "structure_variant": "ASR",
                    "metal_elements": "Cu",
                }
                dataset = _Dataset((row,), hashes={})
                classified = _Classified(dataset, {})
                with self.assertRaisesRegex(
                    ValueError, "exact nonblank string structure_id"
                ):
                    ParentGroupSplitter(classified, parent_method="none")

    def test_splitter_rejects_noncanonical_selected_ids_labels_and_strata(self):
        for value in (True, 1, 1.5, [], {}, b"s000", " s000", "s000 "):
            with self.subTest(surface="selected_ids", value=value):
                classified = self._singleton_classified(count=4)
                classified.structure_ids = (value,)
                with self.assertRaisesRegex(
                    ValueError, "structure IDs must be exact nonblank strings"
                ):
                    ParentGroupSplitter(classified, parent_method="none")

            if not isinstance(value, (list, dict)):
                with self.subTest(surface="label_keys", value=value):
                    classified = self._singleton_classified(count=4)
                    classified.label_by_id = {value: "CR"}
                    with self.assertRaisesRegex(
                        ValueError, "label keys must be exact nonblank structure IDs"
                    ):
                        ParentGroupSplitter(classified, parent_method="none")

            with self.subTest(surface="stratify_by", value=value):
                classified = self._singleton_classified(count=4)
                with self.assertRaisesRegex(
                    ValueError, "stratify_by values must be exact nonblank strings"
                ):
                    ParentGroupSplitter(
                        classified, parent_method="none", stratify_by=(value,)
                    )

        classified = self._singleton_classified(count=4)
        classified.label_by_id["s000"] = True
        with self.assertRaisesRegex(ValueError, "labels must be strings or null"):
            ParentGroupSplitter(classified, parent_method="none")

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

        for value in (False, 0, 0.0, None):
            with self.subTest(stratify_by=value):
                with self.assertRaisesRegex(
                    TypeError, "string or iterable of strings"
                ):
                    ParentGroupSplitter(
                        self._singleton_classified(),
                        parent_method="none",
                        stratify_by=value,
                    )

        for value in (
            {"label", "source_database"},
            frozenset(("label", "source_database")),
            {"label": True},
            (item for item in ("label", "source_database")),
            b"label",
        ):
            with self.subTest(unordered_or_one_shot=type(value).__name__):
                with self.assertRaisesRegex(TypeError, "ordered list/tuple"):
                    ParentGroupSplitter(
                        self._singleton_classified(),
                        parent_method="none",
                        stratify_by=value,
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
        emitted_text = json.dumps(receipt, sort_keys=True, allow_nan=False)
        self.assertEqual(
            receipt["contract_definitions"],
            json.loads(
                json.dumps(
                    generated_output_terminology_definitions(),
                    sort_keys=True,
                    default=dict,
                )
            ),
        )
        for required in (
            "identity_union is the compatibility key for provisional source-ID/MOFid transitive groups",
            "priority_main is the conflict-aware explanatory hierarchy",
            "main_union is the split-leakage guard, not a parent or explanatory method",
            "Canonical identifier text is text processing only",
            "parent_only uses only each group from the selected explanatory parent relation",
        ):
            self.assertIn(required, emitted_text)
        # The expanded explanation is ordinary JSON data, not a Python-only proxy.
        json.loads(emitted_text)

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

    def test_assignment_and_ordered_receipt_are_independent_of_python_hash_seed(self):
        script = r'''
import hashlib
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
                             sources=["COD"],
                             stratify_by=["label", "source_database"],
                             random_state="hash-seed").train_valid_test_split()
payload = json.dumps(result.receipt(), sort_keys=True, separators=(",", ":"))
print(result.receipt()["assignment_sha256"],
      hashlib.sha256(payload.encode("utf-8")).hexdigest())
'''
        outputs = []
        for seed in ("1", "2", "3", "918273"):
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
        self.assertEqual(len(set(outputs)), 1)

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

        for value in (True, 1, 1.5, [], {}, " 3checker"):
            with self.subTest(existing_view=value):
                classified = self._singleton_classified(count=4)
                classified.checker_view = value
                with self.assertRaisesRegex(TypeError, "exact nonblank string"):
                    split_release(
                        classified, checkers="3checker", parent_method="none"
                    )

    def test_split_release_boolean_controls_are_type_exact(self):
        classified = self._singleton_classified(count=4)
        for value in (None, 0, 1, "false", [], {}):
            with self.subTest(verify_cif_files=value):
                with self.assertRaisesRegex(TypeError, "verify_cif_files"):
                    split_release(
                        classified,
                        parent_method="none",
                        verify_cif_files=value,
                    )
            with self.subTest(official=value):
                with self.assertRaisesRegex(TypeError, "official must be a boolean"):
                    split_release(
                        classified,
                        parent_method="none",
                        official=value,
                    )

        classified.dataset.cif_files_verified = "true"
        with self.assertRaisesRegex(TypeError, "cif_files_verified"):
            split_release(
                classified,
                parent_method="none",
                verify_cif_files=True,
            )

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
            self.assertEqual(receipt["implementation"]["split_api_version"], "0.2.0")
            self.assertEqual(receipt["schema_version"], "coremof-split-receipt/1.1")
            self.assertFalse(receipt["integrity"]["cif_files_verified"])
            self.assertEqual(set(receipt["implementation"]["source_sha256"]), {
                "_authority.py", "dataset.py", "labels.py", "parents.py", "splitters.py"
            })
            self.assertIn("PROVISIONAL_PARENT_INPUT", receipt["warnings"])
            self.assertIn("achieved_fractions", receipt)
            self.assertIn("label_counts_by_split", receipt)
            self.assertEqual(receipt["parent_conflict_count"], 0)
            with self.assertRaises(FileExistsError):
                result.write(directory, stem="example")
            result.write(directory, stem="example", overwrite=True)
            original_csv = csv_path.read_bytes()
            original_json = json_path.read_bytes()
            for invalid in ("false", 1, None):
                with self.subTest(overwrite=invalid):
                    with self.assertRaisesRegex(TypeError, "overwrite must be a boolean"):
                        result.write(directory, stem="example", overwrite=invalid)
                    with self.assertRaisesRegex(TypeError, "overwrite must be a boolean"):
                        result.to_csv(csv_path, overwrite=invalid)
                    with self.assertRaisesRegex(TypeError, "overwrite must be a boolean"):
                        result.to_json(json_path, overwrite=invalid)
                    self.assertEqual(csv_path.read_bytes(), original_csv)
                    self.assertEqual(json_path.read_bytes(), original_json)
            with self.assertRaisesRegex(ValueError, "filename stem"):
                result.write(directory, stem="../escape")
        with self.assertRaises(TypeError):
            result.assignments["s000"] = "test"
        with self.assertRaises(TypeError):
            result.filters["labels"] = ["AMBIGUOUS"]

    def test_receipt_revalidates_implementation_hashes_at_use(self):
        first_hashes = {name: "a" * 64 for name in (
            "_authority.py", "dataset.py", "labels.py", "parents.py", "splitters.py"
        )}
        second_hashes = {name: "b" * 64 for name in first_hashes}
        with patch("CoREMOF.splitters._implementation_hashes", return_value=first_hashes):
            result = ParentGroupSplitter(
                self._singleton_classified(count=10), parent_method="none"
            ).train_valid_test_split()
        with patch("CoREMOF.splitters._implementation_hashes", return_value=second_hashes):
            with self.assertRaisesRegex(ValueError, "implementation hash closure"):
                result.receipt()

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

    def test_split_result_authority_does_not_transfer_or_survive_mutation(self):
        result = ParentGroupSplitter(
            self._singleton_classified(count=20), parent_method="none"
        ).train_valid_test_split()
        for operation in (copy.copy, copy.deepcopy, pickle.dumps):
            with self.subTest(operation=operation.__name__):
                with self.assertRaisesRegex(TypeError, "cannot be"):
                    operation(result)

        for forged in (
            replace(result, checker_view="5checker", checker_view_official=True),
            replace(result, official_split=True),
        ):
            with self.subTest(forged=forged.checker_view_official):
                with self.assertRaisesRegex(
                    ValueError, "does not transfer|not produced by"
                ):
                    forged.receipt()

        object.__setattr__(result, "official_split", True)
        with self.assertRaisesRegex(ValueError, "changed after construction"):
            result.receipt()

    def test_official_split_request_fails_closed(self):
        with self.assertRaisesRegex(
            OfficialSplitUnavailableError, "No audited official split manifest"
        ):
            ParentGroupSplitter(
                self._singleton_classified(count=10),
                parent_method="none",
                official=True,
            )
        for value in (1, "false", None):
            with self.subTest(official=value):
                with self.assertRaisesRegex(TypeError, "official must be a boolean"):
                    ParentGroupSplitter(
                        self._singleton_classified(count=10),
                        parent_method="none",
                        official=value,
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

    def test_custom_checker_view_cannot_claim_official_status(self):
        classified = self._singleton_classified(count=10)
        classified.checker_view = "custom:evil"
        classified.checker_view_official = True
        with self.assertRaisesRegex(
            ValueError,
            "checker_view_official=True requires a canonical official checker preset",
        ):
            ParentGroupSplitter(
                classified, parent_method="none"
            ).train_valid_test_split()

        classified.checker_view = "3checker"
        with self.assertRaisesRegex(
            ValueError, "internally authenticated recomputed release view"
        ):
            ParentGroupSplitter(
                classified, parent_method="none"
            ).train_valid_test_split()

    def test_missing_official_marker_never_infers_an_official_view(self):
        classified = self._singleton_classified(count=10)
        classified.checker_view = "5checker"
        if hasattr(classified, "checker_view_official"):
            del classified.checker_view_official
        result = ParentGroupSplitter(
            classified, parent_method="none"
        ).train_valid_test_split()
        self.assertFalse(result.checker_view_official)
        self.assertEqual(result.receipt()["checker_view_kind"], "USER_DEFINED")

    def test_public_status_text_cannot_mint_nonprovisional_input(self):
        classified = self._singleton_classified(count=6)
        classified.dataset.dataset_info = {
            "release_status": "PROVISIONAL_LATEST_AUDITED_SNAPSHOT"
        }
        provisional = ParentGroupSplitter(
            classified, parent_method="none"
        ).train_valid_test_split()
        self.assertTrue(provisional.provisional_input)

        for parent_status, dataset_status in (
            ("FINAL_CANDIDATE", "FINAL"),
            ("FINAL", "FINAL"),
            ("PROVISIONAL_LATEST_AUDITED_SNAPSHOT", "PROVISIONAL_LATEST_AUDITED_SNAPSHOT"),
        ):
            with self.subTest(
                parent_status=parent_status,
                dataset_status=dataset_status,
            ):
                classified.dataset.parent_group_methods = {
                    "release_status": parent_status
                }
                classified.dataset.dataset_info = {
                    "release_status": dataset_status
                }
                with self.assertRaisesRegex(
                    ValueError, "retired or reserved declaration key"
                ):
                    ParentGroupSplitter(
                        classified, parent_method="none"
                    ).train_valid_test_split()


if __name__ == "__main__":
    unittest.main()
