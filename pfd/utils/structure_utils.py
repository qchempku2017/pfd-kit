"""Utility functions for handling pymatgen structures in PFD."""
from typing import List, Optional
import numpy as np
from pymatgen.core import Structure, Site, Element
from pymatgen.analysis.structure_matcher import StructureMatcher
from tqdm import tqdm
import sys


def remove_isolated_atoms(
        struct: Structure,
        radius: float=3.0
) -> Structure:
    """Remove isolated atoms from a structure.

    Args:
        struct (Structure): Input structure.
        radius (float, optional): Radius to consider neighbors. Defaults to 3.0.
    Returns:
        Structure: New structure with isolated atoms removed.
    """
    isolated_inds = []
    for i in range(len(struct)):
        if len(struct.get_neighbors(struct[i], radius)) == 0:
            isolated_inds.append(i)
    struct_cp = struct.copy()
    struct_cp.remove_sites(isolated_inds)
    return struct_cp


def get_site_charge(site: Site):
    """Get the charge of a site. Returns 0 for elemental species."""
    if isinstance(site.specie, Element):
        return 0
    else:
        return site.specie.oxi_state


def get_z_range_indices(
        structure: Structure,
        zmin: float,
        zmax: float
):
    """Get indices of atoms within a specified z-coordinate range in fractional coordinates."""
    frac_coords = structure.frac_coords
    frac_coords = frac_coords - np.floor(frac_coords)
    return np.where((frac_coords[:, 2] >= zmin) & (frac_coords[:, 2] < zmax))[0].tolist()


def deduplicate_structures(
        structs: List[Structure],
        matcher: Optional[StructureMatcher]=None,
        progress_bar: bool=False
) -> List[Structure]:
    """Remove duplicate structures from a list using a StructureMatcher.

    Args:
        structs (List[Structure]): List of structures to deduplicate.
        matcher (Optional[StructureMatcher], optional): StructureMatcher to use.
            If None, a default StructureMatcher is used. Defaults to None.
        progress_bar (bool, optional): Whether to show a progress bar. Defaults to False.
    Returns:
        List[Structure]: List of deduplicated structures.
    """
    if matcher is None:
        matcher = StructureMatcher()
    non_dupe_inds = []
    if not progress_bar:
        iterator = range(len(structs))
    else:
        iterator = tqdm(
            range(len(structs)), total=len(structs), desc="Deduplicating structures",
            file=sys.stdout,  # Ensures tqdm works in Jupyter notebooks.
        )
    for i in iterator:
        dupe = False
        for j in non_dupe_inds:
            if matcher.fit(structs[i], structs[j]):
                dupe = True
                break
        if not dupe:
            non_dupe_inds.append(i)
    return [structs[i] for i in non_dupe_inds]
