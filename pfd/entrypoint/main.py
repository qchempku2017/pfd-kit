import json
import logging

from pfd.utils.download_pfd_artifacts import (
    print_op_download_setting,
)
from .download import download, download_by_def, download_end_result
from .status import status
from .submit import FlowGen, resubmit_workflow
from .common import (
    expand_idx,
    perturb_cli,
    slab_cli,
    interface_cli,
    gb_cli,
)

from pfd.entrypoint.parsers import parse_args


def main():
    # logging
    logging.basicConfig(level=logging.INFO)
    logo = r"""
    ____  _____ ____    _  __ _ _   
   |  _ \|  ___|  _ \  | |/ /(_) |_ 
   | |_) | |_  | | | | | / / | | __|
   |  __/|  _| | |_| | |  /\ | | |_ 
   |_|   |_|   |____/  |_|\_\|_|\__|
    """
    print(logo)
    args = parse_args()

    if args.command == "submit":
        print("Submitting workflow")
        with open(args.CONFIG) as fp:
            config = json.load(fp)
        FlowGen(config).submit(only_submit=args.monitering)

    elif args.command == "resubmit":
        with open(args.CONFIG) as fp:
            config = json.load(fp)
        wfid = args.ID
        resubmit_workflow(
            wf_config=config,
            wfid=wfid,
            list_steps=args.list,
            reuse=args.reuse,
            fold=args.fold,
            only_submit=args.monitering,
        )

    elif args.command == "download":
        with open(args.CONFIG) as fp:
            config = json.load(fp)
        wfid = args.ID
        if args.list_supported is not None and args.list_supported:
            print(print_op_download_setting())
        elif args.keys is not None:
            download(
                wfid,
                config,
                wf_keys=args.keys,
                prefix=args.prefix,
                chk_pnt=args.no_check_point,
            )
        elif args.step_definitions:
            download_by_def(
                wfid,
                config,
                iterations=(
                    expand_idx(args.iterations) if args.iterations is not None else None
                ),
                step_defs=args.step_definitions,
                prefix=args.prefix,
                chk_pnt=args.no_check_point,
            )
        else:
            download_end_result(wfid, config, prefix=args.prefix)
    elif args.command == "status":
        with open(args.CONFIG) as fp:
            config = json.load(fp)
        wfid = args.ID
        status(wfid, config)
    elif args.command == "perturb":
        perturb_cli(
            atoms_path_ls=args.ATOMS,
            pert_num=args.pert_num,
            cell_pert_fraction=args.cell_pert_fraction,
            atom_pert_distance=args.atom_pert_distance,
            atom_pert_style=args.atom_pert_style,
            atom_pert_prob=args.atom_pert_prob,
            supercell=args.supercell
        )
    elif args.command == "slab":
        slab_cli(
            atoms_path_ls=args.ATOMS,
            miller_indices=args.miller_indices,
            symprec=args.symprec,
            angle_tol=args.angle_tol,
            min_slab_ab=args.min_slab_ab,
            min_slab=args.min_slab,
            min_vac=args.min_vac,
            max_normal_search=args.max_normal_search,
            symmetrize_slab=args.symmetrize_slab,
            tasker2_modify_polar=args.tasker2_modify_polar,
            drop_polar=args.drop_polar,
            remove_atom_types=args.remove_atom_types,
            min_vacancy_ratio=args.min_vacancy_ratio,
            max_vacancy_ratio=args.max_vacancy_ratio,
            num_vacancy_ratios=args.num_vacancy_ratios,
            n_sample_per_ratio=args.n_sample_per_ratio,
            surface_mapping_fractol=args.surface_mapping_fractol,
            seed=args.seed,
            detect_isolated_atom_range=args.detect_isolated_atom_range,
            remove_isolated_atom=args.remove_isolated_atom,
            max_return_slabs=args.max_return_slabs,
        )
    elif args.command == "interface":
        interface_cli(
            film_atoms_path=args.FILM,
            substrate_atoms_path=args.SUBSTRATE,
            film_miller=args.film_miller,
            substrate_miller=args.substrate_miller,
            symprec=args.symprec,
            angle_tol=args.angle_tol,
            max_area=args.max_area,
            max_area_ratio_tol=args.max_area_ratio_tol,
            max_length_tol=args.max_length_tol,
            max_angle_tol=args.max_angle_tol,
            zsl_bidirectional=args.zsl_bidirectional,
            termination_ftol=args.termination_ftol,
            filter_out_sym_slabs=args.filter_out_sym_slabs,
            film_thickness=args.film_thickness,
            substrate_thickness=args.substrate_thickness,
            min_ab=args.min_ab,
            gap=args.gap,
            vacuum_over_film=args.vacuum_over_film,
            max_abs_volume_strain=args.max_abs_volume_strain,
            max_abs_von_mises_strain=args.max_abs_von_mises_strain,
            deduplicate_interfaces=args.deduplicate_interfaces,
            structure_matcher_kwargs=args.structure_matcher_kwargs,
            remove_atom_types=args.remove_atom_types,
            vacancy_depth=args.vacancy_depth,
            min_vacancy_ratio=args.min_vacancy_ratio,
            max_vacancy_ratio=args.max_vacancy_ratio,
            num_vacancy_ratios=args.num_vacancy_ratios,
            n_sample_per_ratio=args.n_sample_per_ratio,
            detect_isolated_atom_range=args.detect_isolated_atom_range,
            remove_isolated_atom=args.remove_isolated_atom,
            seed=args.seed,
            max_return_interfaces=args.max_return_interfaces,
        )
    elif args.command == "gb":
        gb_cli(
            prim_path=args.ATOM,
            rotation_axis=args.rotation_axis,
            rotation_angle=args.rotation_angle,
            grain_plane=args.grain_plane,
            n_shifts_per_ab_direction=args.n_shifts_per_ab_direction,
            expand_times=args.expand_times,
            min_ab_size=args.min_ab_size,
            vacuum_thickness=args.vacuum_thickness,
            c_normal=args.c_normal,
            ratio=args.ratio,
            symprec=args.symprec,
            angle_tol=args.angle_tol,
            max_search=args.max_search,
            coincidence_tol=args.coincidence_tol,
            rm_ratio=args.rm_ratio,
            remove_atom_types=args.remove_atom_types,
            min_vacancy_ratio=args.min_vacancy_ratio,
            max_vacancy_ratio=args.max_vacancy_ratio,
            num_vacancy_ratios=args.num_vacancy_ratios,
            n_sample_per_ratio=args.n_sample_per_ratio,
            vacancy_depth=args.vacancy_depth,
            seed=args.seed,
            detect_isolated_atom_range=args.detect_isolated_atom_range,
            remove_isolated_atom=args.remove_isolated_atom,
            max_return_gbs=args.max_return_gbs,
        )
    else:
        raise RuntimeError(f"unknown command {args.command}")
