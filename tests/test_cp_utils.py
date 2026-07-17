from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from CoREMOF.models.cp_app.utils import (
    cv_from_frequencies,
    read_frequencies_from_mesh,
    select_structures,
)


class CpUtilsTests(unittest.TestCase):
    def test_mesh_reader_uses_current_pyyaml_api(self):
        with tempfile.TemporaryDirectory() as directory:
            mesh = Path(directory, "mesh.yaml")
            mesh.write_text(
                "phonon:\n  - band:\n      - frequency: 1.0\n      - frequency: 2.0\n",
                encoding="utf-8",
            )
            frequencies = read_frequencies_from_mesh(mesh)
            self.assertEqual(frequencies.tolist(), [33.35641, 66.71282])

    def test_single_member_structure_type_does_not_raise_index_error(self):
        frame = pd.DataFrame(
            {
                "structure_type": ["only", "another"],
                "atom_types": ["C", "N"],
            }
        )
        selected = select_structures(2, frame)
        self.assertEqual(selected, {0, 1})

    def test_select_structures_rejects_impossible_sample_size(self):
        frame = pd.DataFrame({"structure_type": ["a"], "atom_types": ["C"]})
        with self.assertRaisesRegex(ValueError, "Cannot select 2"):
            select_structures(2, frame)

    def test_heat_capacity_remains_finite_at_high_frequency(self):
        value = cv_from_frequencies(1.0, np.array([1.0e8]))
        self.assertTrue(np.isfinite(value))
        self.assertGreaterEqual(value, 0.0)


if __name__ == "__main__":
    unittest.main()
