"""Utility functions for building grain boundary models in PFD."""
from itertools import product
from typing import Optional, List

from pymatgen.core import Structure
from pymatgen.core.interface import GrainBoundaryGenerator, GrainBoundary
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
from numpy.random import Generator
import math

from pfd.utils.structure_utils import (
    get_site_charge,
    get_z_range_indices,
    remove_isolated_atoms
)

import logging
from logging import Logger


def get_shifted_grain_boundaries(
        prim: Structure,
        rotation_axis: tuple[int, int, int],
        rotation_angle: float,
        grain_plane: Optional[tuple[int, int, int]]=None,
        n_shifts_per_ab_direction: int=5,
        expand_times: int=4,
        min_ab_size: float=10.0,
        vacuum_thickness: float=0.0,
        c_normal: bool=False,
        ratio: Optional[List[int]]=None,
        symprec: float=0.1,
        angle_tol: float=5,
        max_search: int=20,
        coincidence_tol: float=1e-6,
        rm_ratio: float=0.7,
) -> List[GrainBoundary]:
    """Generate grain boundaries with different shifts on a and b directions of the plane.

    In PFD-kit, we follow the grain plane - rotation axis - rotation angle convention for grain
    boundary definition.

    Notice: the users are responsible for selecting the grain boundary building parameters.
    We recommend playing with GranBoundaryGenerator and its utility functions, particularly those for
    searching coincidence site lattices rotation parameters from Sigma values.

    Also, be aware that some grain boundary structure may be generated to be too large. You may want to
    confirm this before submitting to PFD.

    Args:
        prim (Structure): The input primitive structure. Recommended to be conventional standard structure.
        rotation_axis (tuple[int, int, int]): The rotation axis of the other grain with respect to the reference grain.
        rotation_angle (float): The rotation angle of the other in degrees.
        grain_plane (tuple[int, int, int], optional):
            The grain boundary plane normal vector. If None, will choose to be perpendicular to rotation_axis,
            i.e, twist GB.
        n_shifts_per_ab_direction (int, optional): Number of shifts to sample along a and b directions. Default is 5,
            i.e, will sample [0.0, 0.2, 0.4, 0.6, 0.8].
        expand_times (int, optional): Number of times to expand the grain boundary structure along c direction.
            Default is 4. Used to prevent too thin grain boundary structures that may cause interactions
            between periodic images.
        min_ab_size (float, optional): Minimum size of the grain boundary structure along a and b directions.
            Default is 10.0 Angstrom. The grain boundary structure will be expanded accordingly.
        vacuum_thickness (float, optional): Thickness of vacuum layer to add along c direction. Default is 0.0.
        c_normal (bool, optional): Whether to search proper supercell matrix such that the grain boundary plane will be
            normal along c direction. Default is False, as setting True may lead to overly large structures.
        ratio (list[int], optional): lattice axial ratio. Typically not needed.
            See GrainBoundaryGenerator.gb_from_parameters documentation for details.
        symprec (float, optional): Symmetry precision for space group analysis in SymmetryAnalyzer. Default is 0.1.
        angle_tol (float, optional): Angle tolerance for space group analysis in SymmetryAnalyzer. Default is 5 degrees.
        max_search (int, optional): Maximum number of attempts to find coincidence site lattice. Default is 20.
        coincidence_tol (float, optional): Tolerance for coincidence site lattice search. Default is 1e-6.
        rm_ratio (float): Ratio of thresholding minimum interatomic distance relative to the shortest interatomic
            distance in bulk structure. Atoms closer than this threshold will be removed. Default is 0.7.
    Returns:
        List[GrainBoundary]: List of generated grain boundary structures with different ab shifts.
    """
    sa = SpacegroupAnalyzer(prim, symprec, angle_tolerance=angle_tol)
    conv = sa.get_conventional_standard_structure()
    gbg = GrainBoundaryGenerator(
        initial_structure=conv, symprec=symprec, angle_tolerance=angle_tol,
    )

    dxs = np.linspace(0.0, 1.0 - 1.0 / n_shifts_per_ab_direction, n_shifts_per_ab_direction)
    gbs = []
    # generate with grain boundary parameters. the users are responsible for selecting these building parameters.
    for shift in list(product(dxs, dxs)):
        gb = (
            gbg.gb_from_parameters(
                rotation_axis=rotation_axis,
                rotation_angle=rotation_angle,
                plane=grain_plane,
                expand_times=expand_times,
                vacuum_thickness=vacuum_thickness,
                ab_shift=shift,
                normal=c_normal,
                ratio=ratio,
                max_search=max_search,
                tol_coi=coincidence_tol,
                rm_ratio=rm_ratio,
                quick_gen=False,
            )
        )
        mult_a = int(np.ceil(min_ab_size / gb.lattice.a))
        mult_b = int(np.ceil(min_ab_size / gb.lattice.b))
        gb.make_supercell([mult_a, mult_b, 1])
        gbs.append(gb)

    return gbs


def remove_random_gb_sites(
        gb: GrainBoundary,
        remove_ratio: float,
        vacancy_depth: float=3.0,
        rng: Optional[Generator]=None,
        remove_atom_types: Optional[List[str]]=None,
        n_sample_max: int=5,
        detect_isolated_atom_range: float=3.0,
        remove_isolated_atom: bool=True
) -> tuple[List[str], List[GrainBoundary]]:
    """Remove random grain boundary sites to form vacancies.

    Args:
        gb (GrainBoundary): The input grain boundary structure.
        remove_ratio (float): Ratio of grain boundary sites to remove.
        vacancy_depth (float, optional): Depth from the grain boundary plane to consider for vacancy formation.
            Default is 3.0 Angstrom, below and above the grain boundary plane.
        rng (Generator, optional): Numpy random number generator. If None, will create a new one.
        remove_atom_types (List[str], optional): List of atom types to consider for removal.
            If None, will consider all atom types.
        n_sample_max (int, optional): Maximum number of samples to generate. Default is 5.
        detect_isolated_atom_range (float, optional): Distance threshold to detect isolated atoms after removal.
            Default is 3.0 Angstrom.
        remove_isolated_atom (bool, optional): Whether to remove isolated atoms after vacancy formation.
            Default is True.
    Returns:
        tuple[List[str], List[GrainBoundary]]:
            A tuple of list of names and list of grain boundary structures with vacancies.
    """
    if rng is None:
        rng = np.random.default_rng()

    dz = vacancy_depth / gb.lattice.c

    if remove_atom_types is None:
        remove_atom_types = [str(site.specie) for site in gb]
    remove_atom_types = set(remove_atom_types)

    # Assume grain boundary plane is at z = 0.5 and z = 0.0 in fractional coordinates.
    gb_site_inds = (
            get_z_range_indices(gb, 0.5 - dz, 0.5 + dz)
            + get_z_range_indices(gb, 0.0, dz)
            + get_z_range_indices(gb, 1.0 - dz, 1.0)
    )

    gb_site_inds = [ii for ii in gb_site_inds if str(gb[ii].specie) in remove_atom_types]
    n_sites = len(gb_site_inds)
    n_remove = min(int(np.ceil(remove_ratio * n_sites)), n_sites)
    n_sample = min(math.comb(n_sites, n_remove), n_sample_max)

    vac_names = []
    vac_gbs = []
    for samp_id in range(n_sample):
        remove_inds = rng.choice(gb_site_inds, n_remove, replace=False).tolist()
        gb_cp = gb.copy()
        gb_cp.remove_sites(remove_inds)
        # Detect and delete isolated atoms.
        if remove_isolated_atom:
            gb_cp = remove_isolated_atoms(gb_cp, radius=detect_isolated_atom_range)
        real_remove_ratio = float(len(gb) - len(gb_cp)) / n_sites
        vac_gbs.append(gb_cp)
        vac_names.append(f"remove_{remove_ratio:.4f}_sample_{samp_id}_actual_{real_remove_ratio:.4f}")

    return vac_names, vac_gbs


def generate_gbs_with_random_vacancies(
        prim: Structure,
        rotation_axis: tuple[int, int, int],
        rotation_angle: float,
        grain_plane: Optional[tuple[int, int, int]]=None,
        n_shifts_per_ab_direction: int=5,
        expand_times: int=4,
        min_ab_size: float=10.0,
        vacuum_thickness: float=0.0,
        c_normal: bool=False,
        ratio: Optional[List[int]]=None,
        symprec: float=0.1,
        angle_tol: float=5,
        max_search: int=20,
        coincidence_tol: float=1e-6,
        rm_ratio: float=0.7,
        remove_atom_types: Optional[List[str]]=None,
        min_vacancy_ratio: float=0.0,
        max_vacancy_ratio: float=0.3,
        num_vacancy_ratios: int=1,
        n_sample_per_ratio: int=5,
        vacancy_depth: float=3.0,
        seed: Optional[int]=None,
        detect_isolated_atom_range: float=3.0,
        remove_isolated_atom: bool=True,
        max_return_gbs: int=500,
        logger: Optional[Logger]=None,
) -> tuple[List[str], List[GrainBoundary], List[GrainBoundary]]:
    """Generate grain boundary structures with random vacancies.

    In PFD-kit, we follow the grain plane - rotation axis - rotation angle convention for grain
    boundary definition.

    Notice: the users are responsible for selecting the grain boundary building parameters.
    We recommend playing with GranBoundaryGenerator and its utility functions, particularly those for
    searching coincidence site lattices rotation parameters from Sigma values.

    Also, be aware that some grain boundary structure may be generated to be too large. You may want to
    confirm this before submitting to PFD.

    Args:
        prim (Structure): The input primitive structure. Recommended to be conventional standard structure.
        rotation_axis (tuple[int, int, int]): The rotation axis of the other grain with respect to the reference grain.
        rotation_angle (float): The rotation angle of the other in degrees.
        grain_plane (tuple[int, int, int], optional):
            The grain boundary plane normal vector. If None, will choose to be perpendicular to rotation_axis,
            i.e, twist GB.
        n_shifts_per_ab_direction (int, optional): Number of shifts to sample along a and b directions. Default is 5,
            i.e, will sample [0.0, 0.2, 0.4, 0.6, 0.8].
        expand_times (int, optional): Number of times to expand the grain boundary structure along c direction.
            Default is 4, used to prevent too thin grain boundary structures that may cause interactions
            between periodic images.
        min_ab_size (float, optional): Minimum size of the grain boundary structure along a and b directions.
            Default is 10.0 Angstrom. The grain boundary structure will be expanded accordingly.
        vacuum_thickness (float, optional): Thickness of vacuum layer to add along c direction. Default is 0.0.
        c_normal (bool, optional): Whether to search proper supercell matrix such that the grain boundary plane will be
            normal along c direction. Default is False, as setting True may lead to overly large structures.
        ratio (list[int], optional): lattice axial ratio. Typically not needed.
            See GrainBoundaryGenerator.gb_from_parameters documentation for details.
        symprec (float, optional): Symmetry precision for space group analysis in SymmetryAnalyzer. Default is 0.1.
        angle_tol (float, optional): Angle tolerance for space group analysis in SymmetryAnalyzer. Default is 5 degrees.
        max_search (int, optional): Maximum number of attempts to find coincidence site lattice. Default is 20.
        coincidence_tol (float, optional): Tolerance for coincidence site lattice search. Default is 1e-6.
        rm_ratio (float): Ratio of thresholding minimum interatomic distance relative to the shortest interatomic
            distance in bulk structure. Atoms closer than this threshold will be removed. Default is 0.7.
        remove_atom_types (List[str], optional): List of atom types to consider for removal.
            If None, will consider all atom types.
        min_vacancy_ratio (float, optional): Minimum ratio of grain boundary sites to remove. Default is 0.0.
        max_vacancy_ratio (float, optional): Maximum ratio of grain boundary sites to remove. Default is 0.3.
        num_vacancy_ratios (int, optional): Number of vacancy ratios to sample between min and max. Default is 1.
        n_sample_per_ratio (int, optional): Number of samples to generate per vacancy ratio. Default is 5.
        vacancy_depth (float, optional): Depth from the grain boundary plane to consider for vacancy formation.
            Default is 3.0 Angstrom, below and above the grain boundary plane.
        seed (int, optional): Random seed for reproducibility. If None, will use random seed.
        detect_isolated_atom_range (float, optional): Distance threshold to detect isolated atoms after removal.
            Default is 3.0 Angstrom.
        remove_isolated_atom (bool, optional): Whether to remove isolated atoms after vacancy formation.
            Default is True.
        max_return_gbs (int, optional): Maximum number of grain boundary structures to return.
            If more structures are generated, will randomly sample this amount to return. Default is 500.
        logger (Logger, optional): Logger for logging messages. If None, will create a new logger.
    Returns:
        tuple[List[str], List[GrainBoundary], List[GrainBoundary]]:
            A tuple of list of names, list of grain boundary structures with vacancies,
            and list of original grain boundary structures before vacancy formation.
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    gbs = get_shifted_grain_boundaries(
        prim,
        rotation_axis,
        rotation_angle,
        grain_plane=grain_plane,
        n_shifts_per_ab_direction=n_shifts_per_ab_direction,
        expand_times=expand_times,
        min_ab_size=min_ab_size,
        vacuum_thickness=vacuum_thickness,
        c_normal=c_normal,
        ratio=ratio,
        symprec=symprec,
        angle_tol=angle_tol,
        max_search=max_search,
        coincidence_tol=coincidence_tol,
        rm_ratio=rm_ratio,
    )
    rng = np.random.default_rng(seed)
    remove_ratios = np.linspace(min_vacancy_ratio, max_vacancy_ratio, num_vacancy_ratios)
    vac_gbs = []
    vac_names = []

    if max_vacancy_ratio > 0 and any(get_site_charge(site) != 0 for site in prim):
        logger.warning(
            "Removing atoms from ionic grain boundary. May result in charge imbalance."
            " Make sure you know what you are doing!"
        )

    for gb_id, gb in enumerate(gbs):
        logger.info(f"Generating random vacancies, GB: {gb_id + 1}/{len(gbs)}.")
        prefix = f"gb_{gb_id}_"
        for remove_ratio in remove_ratios:
            sub_names, sub_gbs = remove_random_gb_sites(
                gb, remove_ratio, vacancy_depth=vacancy_depth,
                rng=rng,
                remove_atom_types=remove_atom_types,
                n_sample_max=n_sample_per_ratio,
                detect_isolated_atom_range=detect_isolated_atom_range,
                remove_isolated_atom=remove_isolated_atom
            )
            sub_names = [prefix + name for name in sub_names]
            vac_gbs.extend(sub_gbs)
            vac_names.extend(sub_names)

    n_total = len(vac_names)
    if n_total > max_return_gbs:  # Require sub-sampling.
        logger.warning(f"More structures ({n_total}) generated than required ({max_return_gbs}). Down sampling.")
        sample_inds = rng.choice(n_total, size=max_return_gbs, replace=False).astype(int).tolist()
        vac_names = [vac_names[i] for i in range(n_total) if i in sample_inds]
        vac_gbs = [vac_gbs[i] for i in range(n_total) if i in sample_inds]
    # Replace slash with underscore to prevent filename saving issues.
    vac_names = [name.replace("/", "_") for name in vac_names]
    return vac_names, vac_gbs, gbs  # gbs are to be saved into json.
