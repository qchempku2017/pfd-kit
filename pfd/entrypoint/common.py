
from pathlib import (
    Path,
)
from typing import (
    Dict,
    List,
    Optional,
    Union,
    Tuple
)
import os
import dflow

from monty.serialization import dumpfn
from pymatgen.io.ase import AseAtomsAdaptor

from pfd.utils import (
    bohrium_config_from_dict,
    workflow_config_from_dict,
    perturb
)
from pfd.utils.slab_utils import generate_slabs_with_random_vacancies
from pfd.utils.interface_utils import generate_interfaces_with_random_vacancies

from ase import Atoms
from ase.io import read,write


def global_config_workflow(
    wf_config,
):
    # dflow_config, dflow_s3_config
    workflow_config_from_dict(wf_config)

    if os.getenv("DFLOW_DEBUG"):
        dflow.config["mode"] = "debug"
        return None

    # bohrium configuration
    if wf_config.get("bohrium_config") is not None:
        bohrium_config_from_dict(wf_config["bohrium_config"])


def expand_sys_str(root_dir: Union[str, Path]) -> List[str]:
    root_dir = Path(root_dir)
    matches = [str(d) for d in root_dir.rglob("*") if (d / "type.raw").is_file()]
    if (root_dir / "type.raw").is_file():
        matches.append(str(root_dir))
    return matches


def expand_idx(in_list) -> List[int]:
    ret = []
    for ii in in_list:
        if isinstance(ii, int):
            ret.append(ii)
        elif isinstance(ii, str):
            # e.g., 0-41:1
            step_str = ii.split(":")
            if len(step_str) > 1:
                step = int(step_str[1])
            else:
                step = 1
            range_str = step_str[0].split("-")
            if len(range_str) == 2:
                ret += range(int(range_str[0]), int(range_str[1]), step)
            elif len(range_str) == 1:
                ret += [int(range_str[0])]
            else:
                raise RuntimeError("not expected range string", step_str[0])
    ret = sorted(list(set(ret)))
    return ret


def perturb_cli(
    atoms_path_ls: List[Union[str,Path]], 
    pert_num: int, 
    cell_pert_fraction: float, 
    atom_pert_distance: float, 
    atom_pert_style: str, 
    atom_pert_prob: float, 
    supercell: Optional[Union[int, Tuple[int,int,int]]] = None,
    ):
    """A CLI function to perturb structures from file paths.
    """
    #pert_atoms_ls = []
    for atoms_path in atoms_path_ls:
        atoms_ls = read(atoms_path,index=':')
        pert_atom_ls = perturb(
            atoms_ls,
            pert_num,
            cell_pert_fraction,
            atom_pert_distance,
            atom_pert_style,
            atom_pert_prob=atom_pert_prob,
            supercell=supercell
        )
        write("pert_"+Path(atoms_path).stem+'.extxyz',pert_atom_ls,format='extxyz')


def slab_cli(
        atoms_path_ls: List[Union[str,Path]],
        miller_indices: List[Tuple[int,int,int]],
        **kwargs,
    ):
    """A CLI function to create slabs from file paths.

    Args:
        atoms_path_ls: List of file paths containing structures.
        miller_indices: List of Miller indices for slab generation.
        **kwargs: Additional arguments for slab generation. See `generate_slabs_with_random_vacancies`
        in `pfd.utils.slab_utils` for details.
    """
    slab_data = {}
    for atoms_path in atoms_path_ls:
        atoms_ls = read(atoms_path,index=':')
        name = Path(atoms_path).stem
        for atoms_id, atoms in enumerate(atoms_ls):
            for miller_index in miller_indices:
                vac_names, vac_slabs, slabs = generate_slabs_with_random_vacancies(
                    AseAtomsAdaptor.get_structure(atoms),
                    miller_index=miller_index,
                    **kwargs,
                )
                keyname = f"{name}_{atoms_id}_miller_{miller_index[0]}_{miller_index[1]}_{miller_index[2]}"
                slab_data[keyname] = {}
                slab_data[keyname]["slabs"] = slabs
                slab_data[keyname]["vac_names"] = vac_names
                vac_slabs_atoms = [AseAtomsAdaptor.get_atoms(slab) for slab in vac_slabs]
                write(keyname+".extxyz", vac_slabs_atoms, format='extxyz')
    dumpfn(slab_data, "slab_data.json")  # save slab data for audit.


def interface_cli(
        film_atoms_path: Union[str, Path],
        substrate_atoms_path: Union[str, Path],
        film_miller: Tuple[int, int, int],
        substrate_miller: Tuple[int, int, int],
        **kwargs,
    ):
    """A CLI function to create interfaces from file paths.

    Currently, only support single film and single substrate structure, and single miller
    index for each.

    Args:
        film_atoms_path: File path containing film structure.
        substrate_atoms_path: File path containing substrate structure.
        film_miller: Miller index for film slab generation.
        substrate_miller: Miller index for substrate slab generation.
        **kwargs: Additional arguments for interface generation. See `generate_interfaces_with_random_vacancies`
        in `pfd.utils.interface_utils` for details.
    """
    interface_data = {}
    film_atoms = read(film_atoms_path, index=0)
    substrate_atoms = read(substrate_atoms_path, index=0)
    film_name = Path(film_atoms_path).stem
    substrate_name = Path(substrate_atoms_path).stem
    vac_names, vac_interfaces, interfaces = generate_interfaces_with_random_vacancies(
        AseAtomsAdaptor.get_structure(film_atoms),
        AseAtomsAdaptor.get_structure(substrate_atoms),
        film_miller=film_miller,
        substrate_miller=substrate_miller,
        **kwargs,
    )
    keyname = (
        f"film_{film_name}_miller"
        f"_{film_miller[0]}_{film_miller[1]}_{film_miller[2]}"
        f"_substrate_{substrate_name}_miller"
        f"_{substrate_miller[0]}_{substrate_miller[1]}_{substrate_miller[2]}"
    )
    interface_data[keyname] = {}
    interface_data[keyname]["interfaces"] = interfaces
    interface_data[keyname]["vac_names"] = vac_names
    vac_interface_atoms = [AseAtomsAdaptor.get_atoms(interface) for interface in vac_interfaces]
    write(keyname + ".extxyz", vac_interface_atoms, format='extxyz')
    dumpfn(interface_data, "interface_data.json")  # save interface data for audit.
