"""
This script projects sparse data from a ROOT file based on the provided configuration.

Run:
    python project_data_from_sparse.py config.yaml
"""

from pathlib import Path
import argparse
from itertools import product
from concurrent.futures import ProcessPoolExecutor
import sys
sys.path.append(str(Path(__file__).resolve().parents[1] / 'utils'))
from utils import pt_folders_exist, get_root_files_in_directory
import numpy as np
import yaml
import ROOT
# pylint: disable=no-member


def get_file_names_from_string(
        input_files_names_str: str,
        config_input: dict,
        pt_cuts: tuple):
    """
    Get input ROOT file names from a directory string.

    Args:
    - input_files_names_str (str): Directory path containing ROOT files.

    Returns:
    - input_files_names (list): A list of input ROOT file names.
    """
    if Path(input_files_names_str).is_dir():
        # if a directory is provided, get all ROOT files in it
        if config_input["is_preprocessed"]:
            # check if pt folders exist
            if not pt_folders_exist(input_files_names_str, pt_cuts):
                raise ValueError("Not all pt folders exist in the provided directory.")

            # get all root files in the pt folders
            folders = [
                Path(input_files_names_str) / f'pt_{int(pt_min*10)}_{int(pt_max*10)}'
                for pt_min, pt_max in zip(pt_cuts[0], pt_cuts[1])
            ]
            files = []
            for f in folders:
                files += get_root_files_in_directory(f)
            return files

        return get_root_files_in_directory(input_files_names_str)

    return [input_files_names_str]

def get_file_names(config_input, pt_cuts):
    """
    Get input ROOT file names based on the provided configuration.

    Args:
    - config_input (dict): Configuration dictionary containing input information.
    - pt_cuts (tuple): Tuple containing pt minimums and maximums.

    Returns:
    - input_files_names (list): A list of input ROOT file names.
    """
    input_files_names = config_input["names"]

    if isinstance(input_files_names, str):
        return get_file_names_from_string(input_files_names, config_input, pt_cuts)

    if isinstance(input_files_names, list):
        all_files = []
        for item in input_files_names:
            all_files += get_file_names_from_string(item, config_input, pt_cuts)
        return all_files

    raise ValueError("Input names must be a string or a list of strings.")

def get_inputs(config_input, pt_cuts):
    """
    Get input ROOT files based on the provided configuration.

    Args:
    - config_input (dict): Configuration dictionary containing input information.

    Returns:
    - input_files (list): A list of TFile objects opened from the ROOT files.
    """
    input_files_names = get_file_names(config_input, pt_cuts)
    sparse_names = config_input["sparses"]
    if not isinstance(sparse_names, list):
        sparse_names = [sparse_names]
    if len(sparse_names) != len(input_files_names):
        if len(sparse_names) == 1:
            # single sparse, multiple files
            # -> get the same sparse from all files
            sparse_names = sparse_names * len(input_files_names)
        elif len(input_files_names) == 1:
            # single input file, multiple sparses
            # -> get all the sparses from the same file
            input_files_names = input_files_names * len(sparse_names)
        else:
            raise ValueError("Number of input files and sparse names do not match.")

    return input_files_names, sparse_names

def load_sparse_from_task(config_input, pt_cuts):
    """
    Load sparse data from a ROOT file based on the provided configuration.

    Args:
    - config_input (dict): Configuration dictionary containing input information.
    - pt_cuts (tuple): Tuple containing pt minimums and maximums.
    
    Returns:
    - sparses (list): A list of THnSparse objects loaded from the ROOT file.
    """
    input_files_names, sparse_names = get_inputs(config_input, pt_cuts)

    sparses = []
    for input_file_name, sparse_name in zip(input_files_names, sparse_names):
        with ROOT.TFile.Open(input_file_name, "READ") as input_file:
            sparses.append(input_file.Get(sparse_name))
    print(sparses)
    return sparses


def load_event_collision_histograms_from_task(config_input, cutset):
    """
    Load event histogram from a ROOT file based on the provided configuration.

    Args:
    - config_input (dict): Configuration dictionary containing input information.

    Returns:
    - h_ev (TH1F): Event histogram loaded from the ROOT file.
    """
    input_files_names, sparse_names = get_inputs(config_input, cutset)
    h_evs = []
    h_colls = []
    for file_name in input_files_names:
        with ROOT.TFile.Open(file_name, "READ") as input_file:
            h_evs.append(input_file.Get(config_input["event_histogram"]))
            h_evs[-1].SetDirectory(0)
            h_colls.append(input_file.Get(config_input["collision_histogram"]))
            h_colls[-1].SetDirectory(0)
    if len(input_files_names) == 1 and len(sparse_names) > 1:
        h_evs = h_evs * len(sparse_names)
        h_colls = h_colls * len(sparse_names)
    return h_evs, h_colls


def get_cuts(cuts_config):
    """
    Generate a list of cuts based on the provided configuration.

    Parameters:
        - cuts_config (dict): A dictionary containing the cutset.

    Returns:
        - cuts (list): A list of cutsets, where each cutset is represented as a
            tuple of dictionaries.
            Each dictionary contains 'min', 'max', and 'axis' keys.
            If centrality selections are provided, they are included in the cuts.
    """
    if "cent" in cuts_config:
        cent_cfg = cuts_config.pop("cent")
        centrality_sels = [
            {'varname': 'cent', 'axis': cent_cfg['axisnum'], 'min': mn, 'max': mx}
            for mn, mx in zip(cent_cfg['min'], cent_cfg['max'])
        ]
    else:
        centrality_sels = None

    var_names = [v for v in cuts_config if v != 'mass']
    cuts_list  = [cuts_config[v] for v in var_names]
    axes       = [c['axisnum'] for c in cuts_list]
    mins_list  = [c['min'] for c in cuts_list]
    maxs_list  = [c['max'] for c in cuts_list]

    # Check for list-of-lists (centrality-dependent cuts per pt bin)
    has_cent_dependent = any(isinstance(mins[0], list) for mins in mins_list)

    if has_cent_dependent:
        if not centrality_sels:
            raise ValueError(
                "Centrality selections must be provided when any cut is a list of lists."
            )
        n_cent = len(centrality_sels)

        per_var_per_cent = []
        for var_name, axis, mins, maxs in zip(var_names, axes, mins_list, maxs_list):
            if isinstance(mins[0], list):
                assert len(mins) == n_cent, (
                    f"'{var_name}' has {len(mins)} cent bins but {n_cent} were provided."
                )
                per_var_per_cent.append([
                    [{'varname': var_name, 'min': mn, 'max': mx, 'axis': axis}
                     for mn, mx in zip(mins_cent, maxs_cent)]
                    for mins_cent, maxs_cent in zip(mins, maxs)
                ])
            else:
                cut_list = [
                    {'varname': var_name, 'min': mn, 'max': mx, 'axis': axis}
                    for mn, mx in zip(mins, maxs)
                ]
                per_var_per_cent.append([cut_list] * n_cent)

        cuts = []
        for i_cent, cent_sel in enumerate(centrality_sels):
            pt_cuts = list(zip(*[per_var_per_cent[v][i_cent] for v in range(len(var_names))]))
            for cut_tuple in pt_cuts:
                cuts.append((*cut_tuple, cent_sel))

    else:
        per_var = [
            [{'varname': var_name, 'min': mn, 'max': mx, 'axis': axis}
             for mn, mx in zip(mins, maxs)]
            for var_name, axis, mins, maxs in zip(var_names, axes, mins_list, maxs_list)
        ]

        pt_cuts = list(zip(*per_var))

        if centrality_sels:
            cuts = [
                (*pt_cut, cent_sel)
                for pt_cut, cent_sel in product(pt_cuts, centrality_sels)
            ]
        else:
            cuts = pt_cuts

    return cuts


def project_sparse_worker(sparse, cut):
    """
    Projects sparse data based on given cuts and returns histograms for mass and pt.

    Parameters:
        sparse (ROOT.THnSparse): The sparse data to be projected.
        cut (list of dict): List of dictionaries specifying the cuts. Each dictionary should have:
            - "varname" (str): The variable name (e.g., "Pt", "Cent").
            - "axis" (int): The axis index to apply the cut on.
            - "min" (float): The minimum value for the cut.
            - "max" (float): The maximum value for the cut.

    Returns:
        tuple: A tuple containing:
            - h_mass (ROOT.TH1): The projected mass histogram.
            - h_pt (ROOT.TH1): The projected pt histogram.
    """
    pt_min, pt_max = None, None
    cent_min, cent_max = None, None

    for cut_var in cut:
        sparse.GetAxis(cut_var["axis"]).SetRangeUser(cut_var["min"], cut_var["max"])
        if cut_var["varname"] == "pt":
            pt_min, pt_max = cut_var["min"], cut_var["max"]
        if cut_var["varname"] == "cent":
            cent_min, cent_max = cut_var["min"], cut_var["max"]

    # Project the sparse data
    h_mass = sparse.Projection(0)
    h_mass.SetDirectory(0)
    h_pt = sparse.Projection(1)
    h_pt.SetDirectory(0)

    if cent_min is not None and cent_max is not None:
        suffix = f'_{pt_min * 10:.0f}_{pt_max * 10:.0f}_cent_{cent_min:.0f}_{cent_max:.0f}'
    else:
        suffix = f'_{pt_min * 10:.0f}_{pt_max * 10:.0f}'
    h_mass.SetName(f'h_mass{suffix}')
    h_pt.SetName(f'h_pt{suffix}')

    return h_mass, h_pt

def to_be_skipped(cut, pt_mins, pt_maxs, isparse):
    """
    Determine if the current cut should be skipped based on pt cuts and sparse index.

    Parameters:
        cut (list of dict): List of dictionaries specifying the cuts.
        pt_mins (list): List of pt minimums from the cuts configuration.
        pt_maxs (list): List of pt maximums from the cuts configuration.
        isparse (int): Index of the current sparse being processed.

    Returns:
        bool: True if the cut should be skipped, False otherwise.
    """
    pt_min, pt_max = None, None
    for cut_var in cut:
        if cut_var["varname"] == "pt":
            pt_min, pt_max = cut_var["min"], cut_var["max"]
    if (pt_min, pt_max) != (pt_mins[isparse], pt_maxs[isparse]):
        return True
    return False

def project_sparse(config_file_name):  # pylint: disable=too-many-locals, too-many-statements, too-many-branches
    """
    Project THnSparse data based on the configuration provided in the YAML file.

    Parameters:
        config_file_name (str): Path to the configuration YAML file.
    """
    with open(config_file_name, 'r', encoding="utf8") as config_file:
        config = yaml.load(config_file, yaml.FullLoader)

    with open(config["inputs"]["cutset"], 'r', encoding="utf8") as cuts_config_file:
        cuts_config = yaml.load(cuts_config_file, yaml.FullLoader)

    cent_mins, cent_maxs = None, None
    if "cent" in cuts_config:
        cent_mins = cuts_config["cent"]["min"].copy()
        cent_maxs = cuts_config["cent"]["max"].copy()

        # remove 0 - 100 centrality selection from the cuts config (for collision histogram)
        for idx, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
            if cent_min == 0 and cent_max == 100:
                cent_mins.pop(idx)
                cent_maxs.pop(idx)
                break

    pt_mins, pt_maxs = cuts_config["pt"]["min"], cuts_config["pt"]["max"]
    sparses = load_sparse_from_task(config["inputs"], (pt_mins, pt_maxs))
    h_evs, h_colls = load_event_collision_histograms_from_task(config["inputs"], (pt_mins, pt_maxs))

    cuts = get_cuts(cuts_config)
    if isinstance(config["inputs"]["names"], str) and config["inputs"]["is_preprocessed"]:
        # automatically set output labels based on pt cuts
        output_labels = []
        for f in sorted(Path(config["inputs"]["names"]).iterdir()):
            if not f.name.startswith("pt_"):
                continue
            suffix = f.name[3:]  # remove 'pt_' prefix
            pt_min_str, pt_max_str = suffix.split('_')[0:2]
            output_labels.append(f'data_{int(pt_min_str):.0f}_{int(pt_max_str):.0f}.root')
    else:
        output_labels = config["output"]["file_names"]
        if not isinstance(output_labels, list):
            output_labels = [output_labels]

    for isparse, (sparse, h_ev, h_coll, output_label) in enumerate(
        zip(sparses, h_evs, h_colls, output_labels)
    ):
        results = []
        with ProcessPoolExecutor(max_workers=min(len(cuts), 64)) as executor:
            for cut in cuts:
                if config["inputs"]["is_preprocessed"]:
                    if to_be_skipped(cut, pt_mins, pt_maxs, isparse):
                        continue
                results.append(executor.submit(project_sparse_worker, sparse, cut))
        out_file_name = Path(config["output"]["directory"]) / output_label
        out_file_name.parent.mkdir(parents=True, exist_ok=True)
        with ROOT.TFile(str(out_file_name), "RECREATE") as _:
            for result in results:
                h_mass, h_pt = result.result()
                h_mass.Write()
                h_pt.Write()
            if config["inputs"]["is_preprocessed"]:
                # rescale the histograms since they have already been merged
                h_ev.Scale(1.0 / len(sparses))
                h_coll.Scale(1.0 / len(sparses))
            h_ev.Write("h_ev")
            h_coll.Write("h_coll")

            h_coll_proj = h_coll.ProjectionX()
            h_coll_proj.SetName('h_coll_proj')
            h_coll_proj.Write()

            if cent_mins and cent_maxs:
                cent_bins = cent_mins + [cent_maxs[-1]]
                h_coll_rebinned = h_coll_proj.Rebin(
                    len(cent_bins) - 1,
                    'h_coll_rebinned',
                    np.asarray(cent_bins, "d")
                )
                h_coll_rebinned.SetName('h_coll_rebinned')
                h_coll_rebinned.Write()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('config_file_name', metavar='text', default='config.yaml',
                        help='configuration file name')
    args = parser.parse_args()

    project_sparse(args.config_file_name)
