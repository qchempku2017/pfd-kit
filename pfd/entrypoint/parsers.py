import argparse
import textwrap
from typing import (
    List,
    Optional,
)


from pfd import __version__


def add_parser_run(subparsers: argparse._SubParsersAction):
    ##########################################
    # submit
    parser_run = subparsers.add_parser(
        "submit",
        help="Submit workflows",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_run.add_argument(
        "CONFIG", help="the config file in json format defining the workflow."
    )
    parser_run.add_argument(
        "-m",
        "--monitering",
        action="store_false",
        help="Keep monitering the progress",
    )
    return parser_run


def add_parser_resubmit(subparsers: argparse._SubParsersAction):
    ##########################################
    # resubmit
    parser_resubmit = subparsers.add_parser(
        "resubmit",
        help="Submit workflows",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_resubmit.add_argument(
        "CONFIG", help="the config file in json format defining the workflow."
    )
    parser_resubmit.add_argument("ID", help="the ID of existing workflow")
    parser_resubmit.add_argument(
        "-l",
        "--list",
        action="store_true",
        help="list the Steps of the existing workflow.",
    )
    parser_resubmit.add_argument(
        "-u",
        "--reuse",
        type=str,
        nargs="+",
        default=None,
        help="specify which Steps to reuse. e.g., 0-41,\
            the first to the 41st steps.",
    )
    parser_resubmit.add_argument(
        "-f",
        "--fold",
        action="store_true",
        help="if set then super OPs are folded to be reused in the new workflow",
    )
    parser_resubmit.add_argument(
        "-m",
        "--monitering",
        action="store_false",
        help="Keep monitering the progress",
    )
    return parser_resubmit


def add_parser_download(subparsers: argparse._SubParsersAction):
    ##########################################
    # download
    parser_download = subparsers.add_parser(
        "download",
        help=(
            "Download the artifacts of PFD workflow steps. User needs to provide the input json file as well as the workflow ID. The command would then download the end model if workflow is successfully completed.\n"
        ),
        description=(
            textwrap.dedent(
                """
            """
            )
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser_download.add_argument("CONFIG", help="the config file in json format.")
    parser_download.add_argument("ID", help="the ID of the existing workflow.")
    parser_download.add_argument(
        "-l",
        "--list-supported",
        action="store_true",
        help="list all supported steps artifacts",
    )
    parser_download.add_argument(
        "-k",
        "--keys",
        type=str,
        nargs="+",
        help="the keys of the downloaded steps. If not provided download all artifacts",
    )
    parser_download.add_argument(
        "-i",
        "--iterations",
        type=str,
        nargs="+",
        help="the iterations to be downloaded, support ranging expression as 0-10.",
    )
    parser_download.add_argument(
        "-d",
        "--step-definitions",
        type=str,
        nargs="+",
        help="the definition for downloading step artifacts",
    )
    parser_download.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="the prefix of the path storing the download artifacts",
    )
    parser_download.add_argument(
        "-n",
        "--no-check-point",
        action="store_false",
        help="if specified, download regardless whether check points exist.",
    )
    return parser_download


def add_parser_perturb(subparsers: argparse._SubParsersAction):
    #########################################
    # perturb
    parser_perturb = subparsers.add_parser(
        "perturb",
        help="Perturb structures from files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_perturb.add_argument(
        "ATOMS",
        type=str,
        nargs="+",
        help="the structure files to be perturbed, support multiple files.",
    )
    parser_perturb.add_argument(
        "-n",
        "--pert-num",
        type=int,
        default=1,
        help="the number of perturbed structures to be generated for each input structure.",
    )
    parser_perturb.add_argument(
        "-c",
        "--cell-pert-fraction",
        type=float,
        default=0.05,
        help="the fraction of cell perturbation.",
    )
    parser_perturb.add_argument(
        "-d",
        "--atom-pert-distance",
        type=float,
        default=0.2,
        help="the distance to perturb the atom.",
    )
    parser_perturb.add_argument(
        "-s",
        "--atom-pert-style",
        type=str,
        default="normal",
        help="the style of perturbation.",
    )
    parser_perturb.add_argument(
        "-a",
        "--atom-pert-prob",
        type=float,
        default=1.0,
        help="the probability of perturbing each atom.",
    )
    parser_perturb.add_argument(
        "-r",
        "--supercell",
        type=int,
        nargs="+",
        default=None,
        help="the supercell replication, support int or 3 ints.",
    )
    return parser_perturb


def add_parser_slab(subparsers: argparse._SubParsersAction):
    #########################################
    # slab
    parser_slab = subparsers.add_parser(
        "slab",
        help="Generate slabs from structures with optional random surface vacancies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_slab.add_argument(
        "ATOMS",
        type=str,
        nargs="+",
        help="the structure files to generate slabs, support multiple files.",
    )
    parser_slab.add_argument(
        "-m",
        "--miller-indices",
        type=str,
        nargs="+",
        required=True,
        help="Miller indices for slab generation."
             " Each set as a quoted string, e.g. --miller-indices '1 2 3' '1 0 0'",
    )
    parser_slab.add_argument(
        "--symprec",
        type=float,
        default=0.1,
        help="Symmetry precision for SpacegroupAnalyzer.",
    )
    parser_slab.add_argument(
        "--angle-tol",
        type=float,
        default=8,
        help="Angle tolerance for SpacegroupAnalyzer.",
    )
    parser_slab.add_argument(
        "--min-slab-ab",
        type=float,
        default=12.0,
        help="Minimum slab size in a and b directions after supercell construction.",
    )
    parser_slab.add_argument(
        "--min-slab",
        type=float,
        default=12.0,
        help="Minimum slab thickness in c direction.",
    )
    parser_slab.add_argument(
        "--min-vac",
        type=float,
        default=20.0,
        help="Minimum vacuum thickness in c direction.",
    )
    parser_slab.add_argument(
        "--max-normal-search",
        type=int,
        default=20,
        help="Maximum integer supercell factor to search for a normal c direction.",
    )
    parser_slab.add_argument(
        "--symmetrize-slab",
        action="store_true",
        help="Whether to symmetrize the slab.",
    )
    parser_slab.add_argument(
        "--no-tasker2-modify-polar",
        action="store_false",
        dest="tasker2_modify_polar",
        help="Whether not to apply Tasker 2 modification to polar slabs.",
    )
    parser_slab.add_argument(
        "--no-drop-polar",
        action="store_false",
        dest="drop_polar",
        help="Whether not to drop polar slabs after Tasker modification.",
    )
    parser_slab.add_argument(
        "--remove-atom-types",
        type=str,
        nargs="+",
        default=None,
        help="List of atom types (as strings) that can be removed. If None, all atom types can be removed.",
    )
    parser_slab.add_argument(
        "--min-vacancy-ratio",
        type=float,
        default=0.0,
        help="Minimum vacancy ratio on surface sites.",
    )
    parser_slab.add_argument(
        "--max-vacancy-ratio",
        type=float,
        default=0.3,
        help="Maximum vacancy ratio on surface sites.",
    )
    parser_slab.add_argument(
        "--num-vacancy-ratios",
        type=int,
        default=1,
        help="Number of vacancy ratios to sample between min and max.",
    )
    parser_slab.add_argument(
        "--n-sample-per-ratio",
        type=int,
        default=5,
        help="Number of random samples to generate per vacancy ratio.",
    )
    parser_slab.add_argument(
        "--surface-mapping-fractol",
        type=float,
        default=1e-5,
        help="Fractional coordinate tolerance for symmetry mapping.",
    )
    parser_slab.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility.",
    )
    parser_slab.add_argument(
        "--detect-isolated-atom-range",
        type=float,
        default=3.0,
        help="Distance range to detect isolated atoms after removal.",
    )
    parser_slab.add_argument(
        "--no-remove-isolated-atom",
        action="store_false",
        dest="remove_isolated_atom",
        help="Whether not to remove isolated atoms after site removal.",
    )
    parser_slab.add_argument(
        "--max-return-slabs",
        type=int,
        default=500,
        help="Maximum number of slabs to return. If more structures are generated, random sampling is applied.",
    )
    return parser_slab


def add_parser_interface(subparsers: argparse._SubParsersAction):
    #########################################
    # interface
    parser_interface = subparsers.add_parser(
        "interface",
        help="Generate interfaces from structures with optional random surface vacancies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_interface.add_argument(
        "FILM",
        type=str,
        help="Path to the film structure file (single file only), recommend conventional cell.",
    )
    parser_interface.add_argument(
        "SUBSTRATE",
        type=str,
        help="Path to the substrate structure file (single file only), recommend conventional cell.",
    )
    parser_interface.add_argument(
        "--film-miller",
        type=str,
        required=True,
        help="Miller index for film (string of three integers, e.g. '1 0 0')",
    )
    parser_interface.add_argument(
        "--substrate-miller",
        type=str,
        required=True,
        help="Miller index for substrate (string of three integers, e.g. '1 0 0')",
    )
    parser_interface.add_argument(
        "--symprec",
        type=float,
        default=0.1,
        help="Symmetry precision for SpacegroupAnalyzer.",
    )
    parser_interface.add_argument(
        "--angle-tol",
        type=float,
        default=8,
        help="Angle tolerance for SpacegroupAnalyzer.",
    )
    parser_interface.add_argument(
        "--max-area",
        type=float,
        default=400,
        help="Max interface area for ZSLGenerator.",
    )
    parser_interface.add_argument(
        "--max-area-ratio-tol",
        type=float,
        default=0.09,
        help="Max relative area ratio tolerance for ZSLGenerator.",
    )
    parser_interface.add_argument(
        "--max-length-tol",
        type=float,
        default=0.03,
        help="Max relative length tolerance for ZSLGenerator.",
    )
    parser_interface.add_argument(
        "--max-angle-tol",
        type=float,
        default=0.05,
        help="Max relative angle tolerance for ZSLGenerator.",
    )
    parser_interface.add_argument(
        "--zsl-bidirectional",
        action="store_true",
        help="Use bidirectional match of surface lattice in ZSLGenerator.",
    )
    parser_interface.add_argument(
        "--termination-ftol",
        type=float,
        default=0.01,
        help="Fractional tolerance for termination matching.",
    )
    parser_interface.add_argument(
        "--filter-out-sym-slabs",
        action="store_true",
        help="Filter out symmetric slabs.",
    )
    parser_interface.add_argument(
        "--film-thickness",
        type=float,
        default=10.0,
        help="Thickness of film in Angstroms.",
    )
    parser_interface.add_argument(
        "--substrate-thickness",
        type=float,
        default=15.0,
        help="Thickness of substrate in Angstroms.",
    )
    parser_interface.add_argument(
        "--min-ab",
        type=float,
        default=10.0,
        help="Minimum in-plane lattice constants a and b in Angstroms.",
    )
    parser_interface.add_argument(
        "--gap",
        type=float,
        default=2.0,
        help="Gap between film and substrate in Angstroms.",
    )
    parser_interface.add_argument(
        "--vacuum-over-film",
        type=float,
        default=20.0,
        help="Vacuum thickness over film in Angstroms.",
    )
    parser_interface.add_argument(
        "--max-abs-volume-strain",
        type=float,
        default=None,
        help="Max absolute volume strain allowed for interfaces.",
    )
    parser_interface.add_argument(
        "--max-abs-von-mises-strain",
        type=float,
        default=None,
        help="Max absolute von Mises strain allowed for interfaces.",
    )
    parser_interface.add_argument(
        "--deduplicate-interfaces",
        action="store_true",
        help="Deduplicate interfaces using StructureMatcher.",
    )
    parser_interface.add_argument(
        "--structure-matcher-kwargs",
        type=str,
        default=None,
        help="JSON string of kwargs for StructureMatcher.",
    )
    parser_interface.add_argument(
        "--remove-atom-types",
        type=str,
        nargs='+',
        default=None,
        help="List of atom types to consider for removal. If not set, all atom types are considered.",
    )
    parser_interface.add_argument(
        "--vacancy-depth",
        type=float,
        default=3.0,
        help="Depth of the vacancy region into film and substrate in Angstroms.",
    )
    parser_interface.add_argument(
        "--min-vacancy-ratio",
        type=float,
        default=0.0,
        help="Minimum ratio of atoms to remove within the vacancy region.",
    )
    parser_interface.add_argument(
        "--max-vacancy-ratio",
        type=float,
        default=0.2,
        help="Maximum ratio of atoms to remove within the vacancy region.",
    )
    parser_interface.add_argument(
        "--num-vacancy-ratios",
        type=int,
        default=1,
        help="Number of vacancy ratios to sample between min and max.",
    )
    parser_interface.add_argument(
        "--n-sample-per-ratio",
        type=int,
        default=3,
        help="Number of samples to generate per vacancy ratio.",
    )
    parser_interface.add_argument(
        "--detect-isolated-atom-range",
        type=float,
        default=3.0,
        help="Radius to detect isolated atoms after removal.",
    )
    parser_interface.add_argument(
        "--no-remove-isolated-atom",
        action="store_false",
        dest="remove_isolated_atom",
        help="Do not remove isolated atoms after vacancy creation.",
    )
    parser_interface.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility.",
    )
    parser_interface.add_argument(
        "--max-return-interfaces",
        type=int,
        default=500,
        help="Maximum number of interfaces to return. If more are generated, random sampling is applied.",
    )
    return parser_interface


def add_parser_status(subparsers: argparse._SubParsersAction):
    #########################################
    # status
    parser_status = subparsers.add_parser(
        "status",
        help="Check exploration status",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_status.add_argument(
        "CONFIG", help="the config file in json format defining the workflow."
    )
    parser_status.add_argument("ID", help="the ID of existing workflow")
    return parser_status


def main_parser() -> argparse.ArgumentParser:
    """PFD-kit commandline options argument parser.

    Notes
    -----
    This function is used by documentation.

    Returns
    -------
    argparse.ArgumentParser
        the argument parser
    """
    parser = argparse.ArgumentParser(
        description="PFD-kit: fine-tune and distillation from pre-trained atomic models"
                    "machine learning potential energy models.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(title="Valid subcommands", dest="command")

    _ = add_parser_run(subparsers)
    _ = add_parser_resubmit(subparsers)
    _ = add_parser_download(subparsers)
    _ = add_parser_perturb(subparsers)
    _ = add_parser_slab(subparsers)
    _ = add_parser_status(subparsers)
    _ = add_parser_interface(subparsers)

    # Add the version argument
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="show the version number and exit",
    )
    return parser


def parse_args(args: Optional[List[str]] = None):
    """PFD-kit commandline options argument parsing.

    Parameters
    ----------
    args : List[str]
        list of command line arguments, main purpose is testing default option None
        takes arguments from sys.argv
    """
    parser = main_parser()

    parsed_args = parser.parse_args(args=args)
    if parsed_args.command is None:
        parser.print_help()

    # Process miller_indices for slab command
    if getattr(args, 'command', None) == 'slab':
        try:
            processed = []
            for s in args.miller_indices:
                if isinstance(s, str):
                    processed.append(tuple(map(int, s.split())))
                elif isinstance(s, (tuple, list)) and len(s) == 3:
                    processed.append(tuple(map(int, s)))
                else:
                    raise ValueError(f"Invalid miller index: {s}")
            args.miller_indices = processed
        except Exception as e:
            raise ValueError(
                f"Failed to parse --miller-indices: {args.miller_indices}."
                f" Each must be a string of three integers, e.g. '1 2 3'."
            ) from e

    if getattr(args, 'command', None) == 'interface':
        # Process film_miller and substrate_miller for interface command
        try:
            args.film_miller = tuple(map(int, args.film_miller.split()))
            args.substrate_miller = tuple(map(int, args.substrate_miller.split()))
            if len(args.film_miller) != 3 or len(args.substrate_miller) != 3:
                raise ValueError("Miller indices must have exactly three integers.")
        except Exception as e:
            raise ValueError(
                f"Failed to parse --film-miller: {args.film_miller} or --substrate-miller: {args.substrate_miller}."
                f" Each must be a string of three integers, e.g. '1 0 0'."
            ) from e
        # Process structure_matcher_kwargs if provided
        if args.structure_matcher_kwargs is not None:
            import json
            try:
                args.structure_matcher_kwargs = json.loads(args.structure_matcher_kwargs)
            except json.JSONDecodeError as e:
                raise ValueError(
                    f"Failed to parse --structure-matcher-kwargs: {args.structure_matcher_kwargs}."
                    " It must be a valid JSON string."
                ) from e

    return parsed_args
