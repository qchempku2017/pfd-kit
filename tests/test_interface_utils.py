import os
import unittest
from pymatgen.core import Structure
from pymatgen.core.interface import Interface
from pfd.utils import interface_utils


def check_isolated_site(
        structure: Structure,
        r: float=3.0,
) -> bool:
    """Check whether a structure contains isolated sites."""
    for site in structure.sites:
        nns = structure.get_neighbors(site, r)
        if len(nns) == 0:
            return True
    return False


class TestInterfaceUtils(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(test_dir, "data")
        cls.si = Structure.from_file(os.path.join(data_dir, "Si.cif"))
        cls.tio2 = Structure.from_file(os.path.join(data_dir, "TiO2.cif"))

    def test_get_interfaces_si_tio2(self):
        # Test basic interface generation between Si and TiO2
        names, interfaces = interface_utils.get_interfaces(
            film_struct=self.si,
            substrate_struct=self.tio2,
            film_miller=(1, 0, 0),
            substrate_miller=(0, 0, 1),
            film_thickness=8.0,
            substrate_thickness=10.0,
            min_ab=8.0,
            gap=2.0,
            vacuum_over_film=10.0,
            filter_out_sym_slabs=True,  # Prevent too many generated.
        )
        self.assertTrue(len(names) > 0)
        self.assertEqual(len(names), len(interfaces))
        for iface in interfaces:
            self.assertIsInstance(iface, Interface)
            self.assertGreaterEqual(iface.lattice.a, 8)
            self.assertGreaterEqual(iface.lattice.b, 8)
            self.assertGreaterEqual(iface.lattice.c, 18)

    def test_remove_random_interface_sites(self):
        # Generate a simple interface and remove random sites
        _, interfaces = interface_utils.get_interfaces(
            film_struct=self.si,
            substrate_struct=self.tio2,
            film_miller=(1, 0, 0),
            substrate_miller=(0, 0, 1),
            film_thickness=8.0,
            substrate_thickness=10.0,
            min_ab=8.0,
            gap=2.0,
            vacuum_over_film=10.0,
            filter_out_sym_slabs=True,  # Prevent too many generated.
        )
        iface = interfaces[0]
        names, vac_ifaces = interface_utils.remove_random_interface_sites(
            iface, remove_ratio=0.2, vacancy_depth=2.0, n_sample_max=2
        )
        self.assertEqual(len(names), len(vac_ifaces))
        self.assertEqual(len(names), 2)
        for vac in vac_ifaces:
            self.assertIsInstance(vac, Interface)
            self.assertLess(len(vac), len(iface))
            self.assertGreater(len(vac), 0)
            self.assertFalse(check_isolated_site(vac))

        # Test remove only O2-.
        names, vac_ifaces = interface_utils.remove_random_interface_sites(
            iface, remove_ratio=0.3, vacancy_depth=2.0, n_sample_max=3, remove_atom_types=["O2-"],
        )
        for vac in vac_ifaces:
            self.assertGreater(iface.composition["Si0+"], 0)
            self.assertGreater(iface.composition["O2-"], 0)
            self.assertGreater(iface.composition["Ti4+"], 0)
            self.assertEqual(vac.composition["Si0+"], iface.composition["Si0+"])
            self.assertLess(vac.composition["O2-"], iface.composition["O2-"])
            self.assertEqual(vac.composition["Ti4+"], iface.composition["Ti4+"])


    def test_generate_interfaces_with_random_vacancies(self):
        # Generate interfaces with random vacancies for Si/TiO2
        # Should get warning about removing ionic species.
        with self.assertLogs(level="WARNING") as cm:
            names, vac_ifaces, orig_ifaces = interface_utils.generate_interfaces_with_random_vacancies(
                film_struct=self.tio2,
                substrate_struct=self.si,
                film_miller=(0, 0, 1),
                substrate_miller=(1, 0, 0),
                film_thickness=8.0,
                max_area=300,  # Prevent too many generated.
                substrate_thickness=10.0,
                min_ab=8.0,
                gap=2.0,
                vacuum_over_film=10.0,
                min_vacancy_ratio=0.0,
                max_vacancy_ratio=0.1,
                num_vacancy_ratios=2,
                n_sample_per_ratio=1,
                max_return_interfaces=5,
                seed=42,
                filter_out_sym_slabs=True,  # Prevent too many generated.
            )
        self.assertTrue(any('warn' in m.lower() for m in cm.output))
        self.assertTrue(len(names) > 0)
        self.assertEqual(len(names), len(vac_ifaces))
        self.assertTrue(all(isinstance(v, Interface) for v in vac_ifaces))
        self.assertTrue(all(isinstance(o, Interface) for o in orig_ifaces))
        self.assertTrue(all("remove_" in name for name in names))
        self.assertFalse(any("/" in name for name in names))
        self.assertLessEqual(len(names), 5)

        # Now return all interfaces without filtering.
        names, vac_ifaces, orig_ifaces = interface_utils.generate_interfaces_with_random_vacancies(
            film_struct=self.tio2,
            substrate_struct=self.si,
            film_miller=(0, 0, 1),
            substrate_miller=(1, 0, 0),
            film_thickness=8.0,
            max_area=300,  # Prevent too many generated.
            substrate_thickness=10.0,
            min_ab=8.0,
            gap=2.0,
            vacuum_over_film=10.0,
            min_vacancy_ratio=0.0,
            max_vacancy_ratio=0.1,
            num_vacancy_ratios=2,
            n_sample_per_ratio=1,
            max_return_interfaces=1000000,
            seed=42,
            filter_out_sym_slabs=True,  # Prevent too many generated.
        )
        self.assertGreater(len(names), 0)
        self.assertEqual(len(names), len(vac_ifaces))
        self.assertEqual(len(names), 2*len(orig_ifaces))  # 2 vacancy ratios per original interface, 1 sample each.


if __name__ == "__main__":
    unittest.main()
