"""Static regression checks for the dataset-splitting README handbook."""

import re
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
HANDBOOK = ROOT / "README_DATASET_SPLITTING.md"


def _heading_slug(value):
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9 _-]", "", value)
    return re.sub(r"[ _]+", "-", value)


class HandbookTests(unittest.TestCase):
    def test_main_readme_links_the_handbook(self):
        self.assertTrue(HANDBOOK.is_file())
        main_readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn(
            "https://github.com/Chung-Research-Group/CoRE-MOF-Tools/"
            "blob/main/README_DATASET_SPLITTING.md",
            main_readme,
        )

    def test_python_examples_are_syntactically_valid(self):
        text = HANDBOOK.read_text(encoding="utf-8")
        self.assertEqual(text.count("```") % 2, 0)
        blocks = re.findall(r"```python\n(.*?)\n```", text, flags=re.DOTALL)
        self.assertGreaterEqual(len(blocks), 40)
        for index, block in enumerate(blocks, start=1):
            with self.subTest(block=index):
                compile(block, "handbook-python-block-{}".format(index), "exec")

    def test_table_of_contents_anchors_resolve(self):
        text = HANDBOOK.read_text(encoding="utf-8")
        headings = {
            _heading_slug(match.group(1))
            for match in re.finditer(r"^#{1,6}\s+(.+?)\s*$", text, re.MULTILINE)
        }
        anchors = re.findall(r"\]\(#([^)]+)\)", text)
        self.assertGreaterEqual(len(anchors), 20)
        self.assertEqual(sorted(set(anchors).difference(headings)), [])


if __name__ == "__main__":
    unittest.main()
