"""Read-back checks for built distributions, enabled by an explicit artifact path."""

import hashlib
import os
from pathlib import Path
import tarfile
import unittest
import zipfile


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = os.environ.get("COREMOF_DISTRIBUTION_DIR")


@unittest.skipUnless(ARTIFACTS, "set COREMOF_DISTRIBUTION_DIR to inspect built artifacts")
class DistributionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = Path(ARTIFACTS)
        wheels = list(directory.glob("*.whl"))
        sources = list(directory.glob("*.tar.gz"))
        if len(wheels) != 1 or len(sources) != 1:
            raise AssertionError("provide exactly one wheel and one source distribution")
        with zipfile.ZipFile(wheels[0]) as archive:
            cls.wheel = {name: archive.read(name) for name in archive.namelist()}
        with tarfile.open(sources[0]) as archive:
            cls.source = {
                member.name.split("/", 1)[1]: archive.extractfile(member).read()
                for member in archive.getmembers() if member.isfile()
            }

    def test_python_code_matches_the_actual_source(self):
        for path in sorted((ROOT / "CoREMOF").rglob("*.py")):
            relative = path.relative_to(ROOT).as_posix()
            expected = hashlib.sha256(path.read_bytes()).hexdigest()
            for name, archive in (("wheel", self.wheel), ("sdist", self.source)):
                with self.subTest(artifact=name, file=relative):
                    self.assertIn(relative, archive)
                    self.assertEqual(hashlib.sha256(archive[relative]).hexdigest(), expected)

    def test_source_distribution_carries_executable_workflow_examples(self):
        examples = sorted((ROOT / "examples").rglob("*.py"))
        self.assertTrue(examples)
        for path in examples:
            relative = path.relative_to(ROOT).as_posix()
            with self.subTest(file=relative):
                self.assertIn(relative, self.source)
                self.assertEqual(self.source[relative], path.read_bytes())

    def test_source_distribution_carries_handbook_and_portable_skill(self):
        for relative in (
            "README_DATASET_SPLITTING.md",
            "ML_BENCHMARK_HANDOFF.md",
            "examples/README.md",
            "examples/CoREMOF_dataset_splitting_quickstart.ipynb",
            "COMBINED_TARGET_DATASET.md",
            "CoREMOF-COD_COMBINED_TARGET_COVERAGE_20260904.json",
            ".github/workflows/tests.yml",
            "docs/source/target_first_benchmark.rst",
            "manuscript/workflows.md",
            ".agents/skills/coremof-release-curation/SKILL.md",
            ".agents/skills/coremof-release-curation/references/dataset-splitting-ml-benchmark.md",
        ):
            with self.subTest(file=relative):
                self.assertIn(relative, self.source)
                self.assertEqual(self.source[relative], (ROOT / relative).read_bytes())

    def test_source_distribution_carries_runnable_docs(self):
        required = (
            "docs/Makefile",
            "docs/make.bat",
            "docs/requirements.txt",
            "docs/source/conf.py",
            "docs/source/_static/custom.css",
            "docs/source/_static/coremof-logo.png",
            "manuscript/README.md",
            "manuscript/manuscript.md",
            "manuscript/workflows.md",
            "manuscript/evidence.md",
            "manuscript/figures_and_tables.md",
            "manuscript/references.md",
            "manuscript/references.bib",
            "manuscript/verification.json",
            "manuscript/figures/workflow.svg",
            "manuscript/figures/workflow.pdf",
            "manuscript/figures/workflow.png",
            "verify_upload.py",
        )
        for relative in required:
            with self.subTest(file=relative):
                self.assertIn(relative, self.source)
                self.assertEqual(self.source[relative], (ROOT / relative).read_bytes())
        for path in sorted((ROOT / "docs" / "source").rglob("*.rst")):
            relative = path.relative_to(ROOT).as_posix()
            with self.subTest(file=relative):
                self.assertIn(relative, self.source)
                self.assertEqual(self.source[relative], path.read_bytes())

    def test_no_restricted_archives_or_build_caches(self):
        for name, archive in (("wheel", self.wheel), ("sdist", self.source)):
            for relative in archive:
                with self.subTest(artifact=name, file=relative):
                    self.assertNotIn("__pycache__", Path(relative).parts)
                    self.assertNotIn(".git", Path(relative).parts)
                    self.assertFalse(relative.startswith("CoREMOF/data/SI/"))
                    self.assertFalse(relative.startswith("CoREMOF/data/mosaec/"))
                    self.assertNotIn(Path(relative).name, {
                        "_release_checkers_protocol.py", "_release_checkers_worker.py",
                        "_release_mosaec_worker.py", "_release_setc_protocol.py", "_release_setc_worker.py",
                    })

    def test_checker_modules_contain_only_results_only_migration_notices(self):
        for archive in (self.wheel, self.source):
            for name in ('mosaec', 'release_mosaec', 'release_checkers', 'release_setc'):
                code = archive['CoREMOF/' + name + '.py'].decode()
                self.assertIn('results_only', code)
                for forbidden in ('import ccdc', 'import mofchecker', 'import subprocess', 'import torch'):
                    self.assertNotIn(forbidden, code)


if __name__ == "__main__":
    unittest.main()
