from typing import Optional, Dict, List

from pymatgen.core import Structure
from pymatgen.core.interface import Interface
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder
from pymatgen.analysis.interfaces.zsl import ZSLGenerator
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
from numpy.random import Generator
import math
from tqdm import tqdm
import logging
from logging import Logger

from pfd.utils.structure_utils import (
    deduplicate_structures,
    get_site_charge,
    remove_isolated_atoms,
    get_z_range_indices,
)


def get_interfaces(
        film_struct: Structure,
        substrate_struct: Structure,
        film_miller: tuple[int, int, int],
        substrate_miller: tuple[int, int, int],
        symprec: float=0.1,
        angle_tol: float=8,
        max_area: float=400,
        max_area_ratio_tol: float=0.09,
        max_length_tol: float=0.03,
        max_angle_tol: float=0.05,
        zsl_bidirectional: bool=False,
        termination_ftol: float=0.01,
        filter_out_sym_slabs: bool=False,
        film_thickness: float=10.0,
        substrate_thickness: float=15.0,
        min_ab: float=10.0,
        gap: float=2.0,
        vacuum_over_film: float=20.0,
        max_abs_von_mises_strain: Optional[float]=None,
        max_abs_volume_strain: Optional[float]=None,
        deduplicate_interfaces: bool=False,
        structure_matcher_kwargs: Optional[Dict]=None,
        logger: Optional[Logger]=None,
) -> tuple[list[str], list[Interface]]:
    """Generate coherent interfaces between film and substrate structures.

    Generates all possible coherent interfaces given film and substrate structures, and their respective
    miller indices. The resulting interfaces are represented as two stacked slabs with coherent in-plane lattices,
    separated by a gap in the z-direction and wrapped in vacuum.

    Args:
        film_struct (Structure): Film structure. Better be converted to conventional cell beforehand.
        substrate_struct (Structure): Substrate structure. Better be converted to conventional cell beforehand.
        film_miller (tuple[int, int, int]): Miller indices for film surface.
        substrate_miller (tuple[int, int, int]): Miller indices for substrate surface.
        symprec (float, optional): Symmetry precision for SpacegroupAnalyzer. Defaults to 0.1.
        angle_tol (float, optional): Angle tolerance for SpacegroupAnalyzer. Defaults to 8.
        max_area (float, optional): Max interface area for ZSLGenerator. Defaults to 400.
        max_area_ratio_tol (float, optional): Max relative area ratio tolerance for ZSLGenerator. Defaults to 0.09.
        max_length_tol (float, optional): Max relative length tolerance for ZSLGenerator. Defaults to 0.03.
        max_angle_tol (float, optional): Max relative angle tolerance for ZSLGenerator. Defaults to 0.05 (5%).
        zsl_bidirectional (bool, optional): Whether to use bidirectional match of surface lattice in ZSLGenerator.
            Defaults to False. Refer to ZSLGenerator docs for details.
        termination_ftol (float, optional): Fractional tolerance for termination matching. Defaults to 0.01.
        filter_out_sym_slabs (bool, optional): Whether to filter out symmetric slabs. Defaults to False.
            Set False to ensure all inequivalent termination combinations are included. See docs of
            CoherentInterfaceBuilder for details.
        film_thickness (float, optional): Thickness of film in Angstroms. Defaults to 10.0.
        substrate_thickness (float, optional): Thickness of substrate in Angstroms. Defaults to 15.0.
        min_ab (float, optional): Minimum in-plane lattice constants a and b in Angstroms.
            Defaults to 10.0.
        gap (float, optional): Gap between film and substrate in Angstroms. Defaults to 2.0, usually the optimal.
        vacuum_over_film (float, optional): Vacuum thickness over film in Angstroms. Defaults to 20.0.
        max_abs_von_mises_strain (Optional[float], optional): Max absolute von Mises strain allowed for interfaces.
            Interfaces with higher von Mises strain will be filtered out. Defaults to None, meaning no filtering.
        max_abs_volume_strain (Optional[float], optional): Max absolute volume strain allowed for interfaces.
            Interfaces with higher volume strain will be filtered out. Defaults to None, meaning no filtering.
        deduplicate_interfaces (bool, optional): Whether to deduplicate interfaces using StructureMatcher.
            Defaults to False, because it can be time-consuming.
        structure_matcher_kwargs (Optional[Dict], optional): Additional kwargs for StructureMatcher if
            deduplication is enabled. Defaults to None.
        logger (Optional[Logger], optional): Logger for logging messages. Defaults to None.
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    logger.info(
        "** Notice: Coherence of given miller indices will not be checked here!"
        " We recommend using pymatgen.analysis.interfaces.substrate_analyzer.SubstrateAnalyzer"
        " and its calculate() method to perform coherence check before interface construction, and select"
        " appropriate miller indices, volume strain cutoff and von-mises strain cutoff accordingly."
        " To enumerate all valid miller indices, use"
        " pymatgen.core.surface.get_symmetrically_distinct_miller_indices."
    )
    zsl = ZSLGenerator(
        max_area=max_area,
        max_area_ratio_tol=max_area_ratio_tol,
        max_length_tol=max_length_tol,
        max_angle_tol=max_angle_tol,
        bidirectional=zsl_bidirectional,
    )
    # Convert to conventional cell.
    film_sa = SpacegroupAnalyzer(film_struct, symprec=symprec, angle_tolerance=angle_tol)
    film_conv = film_sa.get_conventional_standard_structure()
    subs_sa = SpacegroupAnalyzer(substrate_struct, symprec=symprec, angle_tolerance=angle_tol)
    subs_conv = subs_sa.get_conventional_standard_structure()

    # Build interfaces.
    builder = CoherentInterfaceBuilder(
        film_structure=film_conv,
        substrate_structure=subs_conv,
        film_miller=film_miller,
        substrate_miller=substrate_miller,
        termination_ftol=termination_ftol,
        zslgen=zsl,
        filter_out_sym_slabs=filter_out_sym_slabs,
        label_index=True,  # Must label_index to ensure all inequivalent termination combos are included.
    )
    all_interfaces = []
    all_names = []
    max_von_mises = max_abs_von_mises_strain if max_abs_von_mises_strain is not None else np.inf
    max_vol = max_abs_volume_strain if max_abs_volume_strain is not None else np.inf

    structure_matcher_kwargs = structure_matcher_kwargs or {}
    sm = StructureMatcher(**structure_matcher_kwargs)
    for termination in builder.terminations:
        logger.info(f"Generating interfaces for termination: {termination}.")
        interfaces = list(
            builder.get_interfaces(
                termination=termination,
                film_thickness=film_thickness,
                substrate_thickness=substrate_thickness,
                in_layers=False,  # Use thickness in Angstroms instead of number of layers.
                vacuum_over_film=vacuum_over_film,
                gap=gap
            )
        )
        if deduplicate_interfaces:
            interfaces = deduplicate_structures(interfaces, matcher=sm, progress_bar=True)
        logger.info(f"Number of generated interfaces: {len(interfaces)}.")
        for ii, interface in enumerate(interfaces):
            # Filter by von_mises strain and relative volume expansion/compression.
            von_mises = interface.interface_properties["von_mises_strain"]
            vol = np.trace(interface.interface_properties["strain"])
            if (von_mises <= max_von_mises) and (vol <= max_vol):
                interface_cp = interface.copy()
                # Make supercell to ensure min a,b lattice constants.
                mult_a = int(np.ceil(min_ab / interface.lattice.a))
                mult_b = int(np.ceil(min_ab / interface.lattice.b))
                interface_cp.make_supercell([mult_a, mult_b, 1])  # When optimizing with VASP, ISIF=4 recommended.
                all_interfaces.append(interface_cp)
                all_names.append("-".join(termination) + f"_interface_{ii}")
    return all_names, all_interfaces


def remove_random_interface_sites(
        interface: Interface,
        remove_ratio: float,
        vacancy_depth: float=3.0,
        rng: Optional[Generator]=None,
        remove_atom_types: List[str]=None,
        n_sample_max: int=3,
        detect_isolated_atom_range: float=3.0,
        remove_isolated_atom: bool=True
) -> tuple[list[str], list[Interface]]:
    """Remove random atoms from the interface within a specified z-range.

    Args:
        interface (Interface): Input interface structure.
        remove_ratio (float): Ratio of atoms to remove within the vacancy region.
        vacancy_depth (float, optional): Depth of the vacancy region in Angstroms. Defaults to 3.0.
        rng (Optional[Generator], optional): Numpy random generator. If None, a default generator
            is created. Defaults to None.
        remove_atom_types (List[str], optional): List of atom types to consider for removal.
            If None, all atom types are considered. Defaults to None.
        n_sample_max (int, optional): Maximum number of samples to generate. Defaults to 3.
        detect_isolated_atom_range (float, optional): Radius to detect isolated atoms after removal.
            Defaults to 3.0.
        remove_isolated_atom (bool, optional): Whether to remove isolated atoms after vacancy creation.
            Defaults to True.
    Returns:
        tuple[list[str], list[Interface]]: Tuple of vacancy interface names and corresponding interface
            structures.
    """
    if rng is None:
        rng = np.random.default_rng()

    subs_bot = interface.frac_coords[interface.substrate_indices, 2].min()
    subs_top = interface.frac_coords[interface.substrate_indices, 2].max()
    film_bot = interface.frac_coords[interface.film_indices, 2].min()
    film_top = interface.frac_coords[interface.film_indices, 2].max()

    # Get fractional z range for vacancy region.
    if subs_bot >= film_top:  # substrate above film. Rare case.
        zmid = (film_top + subs_bot) / 2
    else:  # normal case: film above substrate
        zmid = (subs_top + film_bot) / 2

    dz = vacancy_depth / interface.lattice.c
    zmin = zmid - dz
    zmax = zmid + dz

    if remove_atom_types is None:
        remove_atom_types = [str(site.specie) for site in interface]
    remove_atom_types = set(remove_atom_types)

    interfaced_site_inds = get_z_range_indices(interface, zmin, zmax)
    interfaced_site_inds = [
        ii for ii in interfaced_site_inds if str(interface[ii].specie) in remove_atom_types
    ]
    n_sites = len(interfaced_site_inds)
    n_remove = min(int(np.ceil(remove_ratio * n_sites)), n_sites)
    n_sample = min(math.comb(n_sites, n_remove), n_sample_max)

    vac_names = []
    vac_interfaces = []
    for samp_id in range(n_sample):
        remove_inds = rng.choice(interfaced_site_inds, n_remove, replace=False).tolist()
        interface_cp = interface.copy()
        interface_cp.remove_sites(remove_inds)
        # Detect and delete isolated atoms.
        if remove_isolated_atom:
            interface_cp = remove_isolated_atoms(interface_cp, radius=detect_isolated_atom_range)
        real_remove_ratio = float(len(interface) - len(interface_cp)) / n_sites
        vac_interfaces.append(interface_cp)
        vac_names.append(f"remove_{remove_ratio:.4f}_sample_{samp_id}_actual_{real_remove_ratio:.4f}")

    return vac_names, vac_interfaces


def generate_interfaces_with_random_vacancies(
        film_struct: Structure,
        substrate_struct: Structure,
        film_miller: tuple[int, int, int],
        substrate_miller: tuple[int, int, int],
        symprec: float=0.1,
        angle_tol: float=8,
        max_area: float=400,
        max_area_ratio_tol: float=0.09,
        max_length_tol: float=0.03,
        max_angle_tol: float=0.05,
        zsl_bidirectional: bool=False,
        termination_ftol: float=0.01,
        filter_out_sym_slabs: bool=False,
        film_thickness: float=10.0,
        substrate_thickness: float=15.0,
        min_ab: float=10.0,
        gap: float=2.0,
        vacuum_over_film: float=20.0,
        max_abs_volume_strain: Optional[float]=None,
        max_abs_von_mises_strain: Optional[float]=None,
        deduplicate_interfaces: bool=False,
        structure_matcher_kwargs: Optional[Dict]=None,
        remove_atom_types: Optional[List[str]]=None,
        vacancy_depth: float=3.0,  # Vacancy region computed as z = (film_bottom_z + substrate_top_z) / 2 +- vacancy_depth / c
        min_vacancy_ratio: float=0.0,
        max_vacancy_ratio: float=0.2,
        num_vacancy_ratios: int=1,
        n_sample_per_ratio: int=3,
        detect_isolated_atom_range: float=3.0,
        remove_isolated_atom: bool=True,
        seed: Optional[int]=None,
        max_return_interfaces: int=500,
        # If more than this, randomly sample this number. Prevents overly many computations.
        logger: Optional[Logger]=None,
) -> tuple[list[str], list[Interface], list[Interface]]:
    """Generate coherent interfaces with random vacancies at the contact area.

    Args:
        film_struct (Structure): Film structure. Better be converted to conventional cell beforehand.
        substrate_struct (Structure): Substrate structure. Better be converted to conventional cell beforehand.
        film_miller (tuple[int, int, int]): Miller indices for film surface.
        substrate_miller (tuple[int, int, int]): Miller indices for substrate surface.
        symprec (float, optional): Symmetry precision for SpacegroupAnalyzer. Defaults to 0.1.
        angle_tol (float, optional): Angle tolerance for SpacegroupAnalyzer. Defaults to 8.
        max_area (float, optional): Max interface area for ZSLGenerator. Defaults to 400.
        max_area_ratio_tol (float, optional): Max relative area ratio tolerance for ZSLGenerator. Defaults to 0.09.
        max_length_tol (float, optional): Max relative length tolerance for ZSLGenerator. Defaults to 0.03.
        max_angle_tol (float, optional): Max relative angle tolerance for ZSLGenerator. Defaults to 0.05 (5%).
        zsl_bidirectional (bool, optional): Whether to use bidirectional match of surface lattice in ZSLGenerator.
            Defaults to False. Refer to ZSLGenerator docs for details.
        termination_ftol (float, optional): Fractional tolerance for termination matching. Defaults to 0.01.
        filter_out_sym_slabs (bool, optional): Whether to filter out symmetric slabs. Defaults to False.
            Set False to ensure all inequivalent termination combinations are included. See docs of
            CoherentInterfaceBuilder for details. We recommend trying CoherentInterfaceBuilder separately first
            to see how many interfaces are generated for the given miller indices before setting this.
        film_thickness (float, optional): Thickness of film in Angstroms. Defaults to 10.0.
        substrate_thickness (float, optional): Thickness of substrate in Angstroms. Defaults to 15.0.
        min_ab (float, optional): Minimum in-plane lattice constants a and b in Angstroms.
            Defaults to 10.0.
        gap (float, optional): Gap between film and substrate in Angstroms. Defaults to 2.0, usually the optimal.
        vacuum_over_film (float, optional): Vacuum thickness over film in Angstroms. Defaults to 20.0.
        max_abs_volume_strain (Optional[float], optional): Max absolute volume strain allowed for interfaces.
            Interfaces with higher volume strain will be filtered out. Defaults to None, meaning no filtering.
        max_abs_von_mises_strain (Optional[float], optional): Max absolute von Mises strain allowed for interfaces.
            Interfaces with higher von Mises strain will be filtered out. Defaults to None, meaning no filtering.
        deduplicate_interfaces (bool, optional): Whether to deduplicate interfaces using StructureMatcher.
            Defaults to False, because it can be time-consuming.
        structure_matcher_kwargs (Optional[Dict], optional): Additional kwargs for StructureMatcher if
            deduplication is enabled. Defaults to None.
        remove_atom_types (Optional[List[str]], optional): List of atom types to consider for removal.
            If None, all atom types are considered. Defaults to None.
        vacancy_depth (float, optional): Depth of the vacancy region into film and substrate in Angstroms.
            Will remove from region (middle - depth, middle + depth).Defaults to 3.0.
        min_vacancy_ratio (float, optional):
            Minimum ratio of atoms to remove within the vacancy region. Defaults to 0.0.
        max_vacancy_ratio (float, optional):
            Maximum ratio of atoms to remove within the vacancy region. Defaults to 0.2.
        num_vacancy_ratios (int, optional):
            Number of vacancy ratios to sample between min and max. Defaults to 1, means only to use
            min_vacancy_ratio.
        n_sample_per_ratio (int, optional): Number of samples to generate per vacancy ratio. Defaults to 3.
        detect_isolated_atom_range (float, optional): Radius to detect isolated atoms after removal.
            Defaults to 3.0.
        remove_isolated_atom (bool, optional): Whether to remove isolated atoms after vacancy creation.
            Defaults to True.
        seed (Optional[int], optional): Random seed for reproducibility. Defaults to None.
        max_return_interfaces (int, optional): Maximum number of interfaces to return. If more are generated,
            random sampling is performed to limit the number. Defaults to 500.
        logger (Optional[Logger], optional): Logger for logging messages. Defaults to None.
    Returns:
        tuple[list[str], list[Interface], list[Interface]]: Tuple of vacancy interface names,
            corresponding interface structures with vacancies, and original interfaces without vacancies.
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.info("**** Generating valid interfaces.")
    interface_names, interfaces = get_interfaces(
        film_struct=film_struct, substrate_struct=substrate_struct,
        film_miller=film_miller, substrate_miller=substrate_miller,
        symprec=symprec, angle_tol=angle_tol,
        max_area=max_area, max_area_ratio_tol=max_area_ratio_tol,
        max_length_tol=max_length_tol, max_angle_tol=max_angle_tol,
        zsl_bidirectional=zsl_bidirectional,
        termination_ftol=termination_ftol,
        filter_out_sym_slabs=filter_out_sym_slabs,
        film_thickness=film_thickness,
        substrate_thickness=substrate_thickness,
        min_ab=min_ab,
        gap=gap, vacuum_over_film=vacuum_over_film,
        max_abs_von_mises_strain=max_abs_von_mises_strain,
        max_abs_volume_strain=max_abs_volume_strain,
        deduplicate_interfaces=deduplicate_interfaces,
        structure_matcher_kwargs=structure_matcher_kwargs,
        logger=logger,
    )
    # Create random vacancies.
    if max_vacancy_ratio > 0 and (
            any(get_site_charge(site) != 0 for site in film_struct) or
            any(get_site_charge(site) != 0 for site in substrate_struct)
    ):
        logger.warning(
            "Removing atoms from ionic interface. May result in charge imbalance."
            " Make sure you know what you are doing!"
        )

    rng = np.random.default_rng(seed)
    remove_ratios = np.linspace(min_vacancy_ratio, max_vacancy_ratio, num_vacancy_ratios)
    all_names = []
    all_vacancy_interfaces = []
    iterator = tqdm(
        zip(interface_names, interfaces), total=len(interface_names), desc="Enumerating random vacancies",
    )
    for interface_name, interface in iterator:
        for remove_ratio in remove_ratios:
            for vac_name, vac_interface in zip(*remove_random_interface_sites(
                    interface, remove_ratio,
                    vacancy_depth=vacancy_depth,
                    remove_atom_types=remove_atom_types,
                    rng=rng,
                    n_sample_max=n_sample_per_ratio,
                    detect_isolated_atom_range=detect_isolated_atom_range,
                    remove_isolated_atom=remove_isolated_atom
            )):
                all_names.append(f"{interface_name}_{vac_name}")
                all_vacancy_interfaces.append(vac_interface)
    n_total = len(all_names)
    if n_total > max_return_interfaces:  # Require sub-sampling.
        logger.warning(
            f"More structures ({n_total}) generated than required ({max_return_interfaces})."
            f" Down sampling at random."
        )
        sample_inds = rng.choice(n_total, size=max_return_interfaces, replace=False).astype(int).tolist()
        all_names = [all_names[i] for i in range(n_total) if i in sample_inds]
        all_vacancy_interfaces = [all_vacancy_interfaces[i] for i in range(n_total) if i in sample_inds]

    # Replace slash with underscore to prevent filename saving issues.
    all_names = [name.replace("/", "_") for name in all_names]
    logger.info(f"Total vacancy interfaces generated: {len(all_names)}")
    return all_names, all_vacancy_interfaces, interfaces
