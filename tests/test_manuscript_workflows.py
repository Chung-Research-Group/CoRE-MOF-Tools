"""Executable documentation checks, not a scientific production-data audit."""

import ast
import csv
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path
import re
import shlex
import tempfile
import unittest

from CoREMOF.benchmarks import normalize_group_criteria
from CoREMOF.cli import build_parser
from CoREMOF.targets import TargetSource
from test_benchmarks import _authenticated_classified


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT / "examples"
RECIPES = WORKSPACE / "grouped_workflow_recipes.md"


def _workflow_blocks():
    text = RECIPES.read_text(encoding="utf-8")
    blocks = re.findall(r"```python\n(.*?)\n```", text, flags=re.DOTALL)
    result = {}
    for block in blocks:
        name = re.match(r"# workflow: ([a-z_]+)\n", block)
        if name is None or name.group(1) in result:
            raise AssertionError("workflow blocks need unique executable names")
        result[name.group(1)] = block
    return result


class ManuscriptWorkflowTests(unittest.TestCase):
    def test_python_examples_compile_and_expected_functions_exist(self):
        blocks = _workflow_blocks()
        self.assertEqual(set(blocks), {
            "checker_views", "general_split", "frozen_benchmark", "attach_frozen"
        })
        for name, code in blocks.items():
            ast.parse(code, filename=RECIPES.name + ":" + name)
        self.assertEqual(
            normalize_group_criteria(("RT", "M2T")),
            ("rac5_crystalnets", "mofid_v2_crystalnets"),
        )

    def test_documented_workflow_executes_with_explicit_nonnumerical_profile(self):
        namespace = {}
        for name, code in _workflow_blocks().items():
            exec(compile(code, RECIPES.name + ":" + name, "exec"), namespace)
        with tempfile.TemporaryDirectory(prefix="coremof-manuscript-test-") as tmp:
            root = Path(tmp)
            dataset, _ = _authenticated_classified(
                root, cr=20, ncr=8, ambiguous=2, unchecked=1
            )
            views = list(namespace["checker_views"](dataset))
            self.assertEqual(len(views), 16)
            self.assertEqual([len(c) for c, _ in views].count(3), 10)
            self.assertEqual([len(c) for c, _ in views].count(4), 5)
            self.assertEqual([len(c) for c, _ in views].count(5), 1)
            for _, view in views:
                self.assertEqual(len(view), 31)
                self.assertTrue(view.checker_view.startswith("custom:"))

            # These synthetic rows have no optional reference evidence. This
            # tests the copyable API path, not representative clustering or
            # a production combined-reference release contract.
            availability, split, paths = namespace["write_general_split"](
                dataset, root / "general", group_criteria="priority_main",
                diversity="none",
            )
            self.assertIn("priority_main", availability)
            self.assertEqual(len(split.assignments), 28)
            self.assertTrue(all(path.is_file() for path in paths))
            with self.assertRaises(FileExistsError):
                split.write(root / "general")

            cohorts, suite = namespace["freeze_benchmark"](
                dataset, root / "suite", group_criteria="priority_main",
                diversity="none",
            )
            self.assertEqual(dict(cohorts.pool_counts), {"CR": 20, "NCR": 8})
            self.assertEqual(len(suite.runs), 12)
            self.assertFalse(suite.official_split)
            fixed = set(suite.fixed_test_ids)
            for run in suite.runs:
                expected = int((Decimal(run.requested_ncr_pool_fraction) * 8)
                               .quantize(Decimal(1), rounding=ROUND_HALF_UP))
                self.assertEqual(len(run.ncr_ids), expected)
                self.assertEqual(len(run.assignments), 20)
                self.assertEqual(
                    {sid for sid, part in run.assignments.items() if part == "test"},
                    fixed,
                )
                block_parts = {}
                for sid, part in run.assignments.items():
                    block_parts.setdefault(suite.effective_leakage_blocks[sid], set()).add(part)
                self.assertTrue(all(len(parts) == 1 for parts in block_parts.values()))

            target_path = root / "targets.csv"
            with target_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(("structure_id", "uptake"))
                for index, sid in enumerate(dataset.structure_ids):
                    writer.writerow((sid, "" if index == 0 else "1.5"))
            source = TargetSource(
                target_path, target_columns=("uptake",),
                value_types={"uptake": "float"}, units={"uptake": "mol/kg"},
                conditions={"uptake": {"temperature_K": 298, "pressure_bar": 65}},
            )
            before = suite.assignment_digest
            attached = namespace["attach_frozen"](
                suite, (source,), root / "attached"
            )
            self.assertEqual(attached.original_assignment_digest, before)
            self.assertEqual(len(attached.run_views), 12)
            for run in suite.runs:
                view = attached.run_views[run.run_key]
                self.assertEqual(dict(view.assignments), dict(run.assignments))
            self.assertTrue((root / "attached" / "coremof_cr_ncr_targets"
                             / "SHA256SUMS").is_file())

    def test_cli_examples_parse_with_current_parser(self):
        text = RECIPES.read_text(encoding="utf-8")
        commands = []
        for block in re.findall(r"```bash\n(.*?)\n```", text, flags=re.DOTALL):
            command = block.replace("\\\n", " ").strip()
            if command.startswith("coremof "):
                commands.append(build_parser().parse_args(shlex.split(command)[1:]))
        self.assertEqual([args.command for args in commands],
                         ["benchmark-cr-ncr", "attach-targets"])
        self.assertEqual(commands[0].group_criteria, ["RT", "M2T"])
        self.assertEqual(commands[1].missing, "keep")

    def test_local_links_and_sphinx_downloads_resolve(self):
        for path in (RECIPES, ROOT / "COMBINED_TARGET_DATASET.md",
                     ROOT / "ML_BENCHMARK_HANDOFF.md"):
            text = path.read_text(encoding="utf-8")
            for target in re.findall(r"\]\(([^\s)]+)\)", text):
                if "://" in target or target.startswith("#"):
                    continue
                candidate = path.parent / target.split("#", 1)[0]
                self.assertTrue(candidate.exists(), msg=str(candidate))
        landing = ROOT / "docs" / "source" / "research_workflows.rst"
        for target in re.findall(r":download:`[^`]*<([^>]+)>`",
                                 landing.read_text(encoding="utf-8")):
            self.assertTrue((landing.parent / target).is_file(), msg=target)
        self.assertIn("   research_workflows\n",
                      (landing.parent / "index.rst").read_text(encoding="utf-8"))

    def test_portable_skill_relative_links_resolve(self):
        skill = ROOT / ".agents" / "skills" / "coremof-dataset-use"
        self.assertTrue((skill / "SKILL.md").is_file())
        for path in sorted(skill.rglob("*.md")):
            text = path.read_text(encoding="utf-8")
            for target in re.findall(r"\]\(([^\s)]+)\)", text):
                if "://" in target or target.startswith("#"):
                    continue
                candidate = path.parent / target.split("#", 1)[0]
                self.assertTrue(candidate.is_file(), msg=str(candidate))


if __name__ == "__main__":
    unittest.main()
