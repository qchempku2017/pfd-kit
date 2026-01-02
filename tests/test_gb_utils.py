import os
import unittest
from pymatgen.core import Structure
from pymatgen.core.interface import GrainBoundary
from pfd.utils import gb_utils


class TestGbUtils(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(test_dir, "data")
        cls.si = Structure.from_file(os.path.join(data_dir, "Si.cif"))
        cls.tio2 = Structure.from_file(os.path.join(data_dir, "TiO2.cif"))

    def test_get_shifted_grain_boundaries(self):
        # Test for Si, (110) twist GB, sigma 3.
        gbs = gb_utils.get_shifted_grain_boundaries(
            self.si, (1, 1, 0), 70.52877936550931,
            grain_plane=(1, 1, 1),  # Tilt + twist.
            n_shifts_per_ab_direction=2, min_ab_size=8
        )
        self.assertEqual(len(gbs), 4)
        for gb in gbs:
            self.assertIsInstance(gb, GrainBoundary)
            self.assertGreaterEqual(gb.lattice.a, 8)
            self.assertGreaterEqual(gb.lattice.b, 8)

    def test_remove_random_gb_sites(self):
        # Use a generated GB for Si, sigma 3 (110) pure twist.
        gbs = gb_utils.get_shifted_grain_boundaries(
            self.si, (1, 1, 0), 70.52877936550931, n_shifts_per_ab_direction=2, min_ab_size=8
        )
        gb = gbs[0]
        names, vac_gbs = gb_utils.remove_random_gb_sites(
            gb, remove_ratio=0.5, n_sample_max=2
        )
        self.assertEqual(len(names), len(vac_gbs))
        self.assertEqual(len(names), 2)
        for vac in vac_gbs:
            self.assertLess(len(vac), len(gb))
            self.assertGreater(len(vac), 0)

    def test_generate_gbs_with_random_vacancies(self):
        # Use Si, test vacancy generation
        names, vac_gbs, gbs = gb_utils.generate_gbs_with_random_vacancies(
            self.si, (1, 1, 0), 70.52877936550931,
            n_shifts_per_ab_direction=2, min_ab_size=8,
            min_vacancy_ratio=0, max_vacancy_ratio=0.5,
            num_vacancy_ratios=2, n_sample_per_ratio=1
        )
        self.assertTrue(len(names) >= 2)
        self.assertEqual(len(names), len(vac_gbs))
        self.assertTrue(all(isinstance(gb, GrainBoundary) for gb in vac_gbs))
        self.assertEqual(len(gbs), 4)  # 2x2 shifts
        self.assertEqual(len(names), 8)  # 4 GBs * 2 vacancy ratios, 1 sample each.
        self.assertTrue(all("gb_" in name for name in names))
        self.assertFalse(any("/" in name for name in names))

        # Use TiO2, test vacancy generation, (001) sigma 5, should capture warning.
        with self.assertLogs(level="WARNING") as cm:
            names, vac_gbs, gbs = gb_utils.generate_gbs_with_random_vacancies(
                self.tio2, (0, 0, 1), 36.86989764584402,
                grain_plane=(1, 0, 0),  # Pure tilt.
                n_shifts_per_ab_direction=3, min_ab_size=8,
                min_vacancy_ratio=0, max_vacancy_ratio=0.5,
                num_vacancy_ratios=3, n_sample_per_ratio=2
            )
        self.assertTrue(any('warn' in m.lower() for m in cm.output))
        self.assertEqual(len(names), len(vac_gbs))
        self.assertTrue(all(isinstance(gb, GrainBoundary) for gb in vac_gbs))
        self.assertEqual(len(gbs), 9)  # 3x3 shifts
        self.assertEqual(len(names), 45)  # 9 GBs * (1 + 2 * 2)
        self.assertTrue(all("gb_" in name for name in names))
        self.assertFalse(any("/" in name for name in names))


if __name__ == "__main__":
    unittest.main()
