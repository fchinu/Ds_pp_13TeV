"""
Script to produce multiplicity weights for the Ds/D+ ratio analysis.

Run:
    python produce_mult_weights.py config.yaml
"""

import itertools
import os
import argparse
import yaml
import ROOT
# pylint: disable=no-member

CENT_MINS = [0,     0,  1,  10, 30, 50, 70]
CENT_MAXS = [100,   1,  10, 30, 50, 70, 100]


def __find_matching_mc(info_data, proj_mc, particle, origin):
    """
    Find a matching Monte Carlo (MC) histogram based on provided data information.

    Parameters:
    - info_data (dict): Dictionary containing information about the data histogram.
        Expected keys are 'task_name' and 'h_name'.
    - proj_mc (list): List of tuples, where each tuple contains a dictionary with
        MC histogram information and the corresponding histogram object.
    - particle (str): The particle type to match.
    - origin (str): The origin to match.

    Returns:
    - tuple: A tuple containing the matching MC histogram information
        dictionary and the histogram object.

    Raises:
    FileNotFoundError: If no matching MC histogram is found.
    """
    for info_mc, h_mc in proj_mc:
        if info_mc['as_mc'] and info_data['task_name'] == info_mc['task_name'] and\
                info_data['h_name'] == info_mc['h_name'] and info_mc['particle'] == particle and\
                info_mc['origin'] == origin:

            return info_mc, h_mc
    # If not found, try to find a MC histogram with the same task_name and h_name,
    # but treated as data (e.g., for histos_names_as_data_only)
    for info_mc, h_mc in proj_mc:
        if not info_mc['as_mc'] and info_data['task_name'] == info_mc['task_name'] and\
                info_data['h_name'] == info_mc['h_name']:

            return info_mc, h_mc

    raise FileNotFoundError(
        f"Cannot find matching MC histogram for {info_data['task_name']} {info_data['h_name']}"
    )


def get_weights(proj_mc, proj_data, cfg):
    """
    Calculate weights for particle origins based on Monte Carlo (MC) and data projections.

    Parameters:
    - proj_mc (list): A list of tuples containing MC projection information and histograms.
    - proj_data (list): A list of tuples containing data projection information and histograms.
    - cfg (dict): Configuration dictionary containing 'particles' and 'origins' keys.

    Returns:
    - dict: A dictionary where keys are tuples of (particle, origin, task_name, h_name)
        and values are lists containing:
        - The weights histogram.
        - The MC projection tuple.
        - The data projection tuple.
    """
    particles = cfg['particles']
    origins = cfg['origins']
    particles_origins = itertools.product(particles, origins)
    h_weights = {}

    for particle, origin in particles_origins:
        for p_data in proj_data:
            info_data, h_data = p_data
            p_mc = __find_matching_mc(info_data, proj_mc, particle, origin)
            _, h_mc = p_mc
            key = (particle, origin, info_data['task_name'], info_data['h_name'])
            h_weights[key] = h_data.Clone(f"{info_data['h_name']}_weights")
            h_weights[key].Divide(h_mc)
            h_weights[key] = [
                h_weights[particle, origin, info_data['task_name'], info_data['h_name']],
                p_mc, p_data
            ]
    return h_weights


def __get_histogram(input_info):
    """
    Retrieve a histogram from a ROOT file based.

    Parameters:
    - input_info (dict): A dictionary containing the following keys:
        - 'task_name' (str): The task name.
        - 'as_mc' (bool): Flag indicating whether the data should be considered
            as from Monte Carlo simulation.
        - 'particle' (str): The particle type (required if 'as_mc' is True).
        - 'origin' (str): The origin of the particle (required if 'as_mc' is True).
        - 'h_name' (str): The name of the histogram.
        - 'file_name' (str): The path to the ROOT file.

    Returns:
        ROOT.TH1: The requested histogram.

    Raises:
        Exception: If the histogram cannot be found in the specified file.
    """
    path = input_info['task_name']
    if input_info['as_mc']:
        path += f"/MC/{input_info['particle']}/{input_info['origin']}"
    else:
        path += "/Data"
    path += f"/{input_info['h_name']}"

    with ROOT.TFile.Open(input_info['file_name']) as infile:
        h = infile.Get(path)
    if not h:
        raise FileNotFoundError(f"Cannot find histogram {path} in file {input_info['file_name']}")
    return h


def get_histos(cfg, is_mc=False):  # pylint: disable=too-many-locals
    """
    Retrieve histograms based on the provided configuration.

    Args:
    - cfg (dict): Configuration dictionary containing input file paths, task names,
        histogram names, particles, and origins.
    - is_mc (bool, optional): Flag indicating whether to process Monte Carlo (MC) data.
        Defaults to False.

    Returns:
    - list: A list of tuples with input information and the corresponding histogram.
    """
    task_names = cfg['inputs']['task_names']
    histos_names = cfg['inputs']['histos_names']
    histos_names_as_data_only = cfg['inputs']['histos_names_as_data_only']
    particles = cfg['particles']
    origins = cfg['origins']
    if is_mc:
        input_file = cfg['inputs']['mc']
    else:
        input_file = cfg['inputs']['data']

    histos = []

    input_combinations_data = itertools.product(
        task_names, histos_names_as_data_only + histos_names
    )
    for task_name, h_name in input_combinations_data:
        input_info_data = {
            'file_name': input_file,
            'task_name': task_name,
            'h_name': h_name,
            'as_mc': False
        }
        histos.append((input_info_data, __get_histogram(input_info_data)))
    if is_mc:
        input_combinations_mc = itertools.product(task_names, histos_names, particles, origins)
        for task_name, h_name, particle, origin in input_combinations_mc:
            input_info_mc = {
                'file_name': input_file,
                'task_name': task_name,
                'h_name': h_name,
                'as_mc': True,
                'particle': particle,
                'origin': origin
            }
            histos.append((input_info_mc, __get_histogram(input_info_mc)))
    return histos


def get_projections(histos, cent_min=None, cent_max=None):
    """
    Generate projections of histograms.

    Parameters:
    - histos (list of tuples): A list of tuples where each tuple contains:
        - info (dict): Dictionary containing metadata about the histogram.
        - h (ROOT.TH1): The histogram object to be projected.
    - cent_min (float, optional): Minimum centrality value.
    - cent_max (float, optional): Maximum centrality value.

    Returns:
    - list of tuples: A list of tuples where each tuple contains:
        - out_info (dict): Updated metadata dictionary with projection details.
        - h_proj (ROOT.TH1): The projected and scaled histogram.
    """
    h_projections = []
    for info, h in histos:
        out_info = info.copy()
        h_name = f"{info['h_name']}_projected"
        if cent_min is not None and cent_max is not None:
            h.GetYaxis().SetRangeUser(cent_min, cent_max)
            h_name += f"_{cent_min}_{cent_max}"
            out_info['cent_min'] = cent_min
            out_info['cent_max'] = cent_max
        h_proj = h.ProjectionX(h_name, 0, -1, "e")
        h_proj.Scale(1./h_proj.Integral())
        if cent_min is not None and cent_max is not None:
            h_proj.SetMarkerStyle(ROOT.kFullCircle)
            h_proj.SetMarkerColor(ROOT.kRed)
            h_proj.SetLineColor(ROOT.kRed)
        else:
            h_proj.SetMarkerStyle(ROOT.kFullDiamond)
            h_proj.SetMarkerColor(ROOT.kAzure + 3)
            h_proj.SetLineColor(ROOT.kAzure + 3)
        h_projections.append((out_info, h_proj))
    return h_projections


def dump_to_root(cfg, h_weights, cent_min, cent_max):  # pylint: disable=too-many-locals
    """
    Dump histograms and weights to a ROOT file organized by centrality and particle origin.

    Parameters:
    - cfg (dict): Configuration dictionary containing input task names and output directory/suffix.
    - h_weights (dict): Dictionary containing histogram weights and associated metadata.
    - cent_min (int): Minimum centrality value.
    - cent_max (int): Maximum centrality value.
    """
    for task_name, suffix in zip(cfg['inputs']['task_names'], cfg['output']['suffix']):
        out_file_name = os.path.join(cfg['output']['dir'], f"weights_{suffix}.root")
        with ROOT.TFile(out_file_name, "UPDATE") as out_file:
            # create a directory for the centrality range
            if not out_file.Get(f"centrality_{cent_min}_{cent_max}"):
                out_file.mkdir(f"centrality_{cent_min}_{cent_max}")
            canvases = {}
            for (particle, origin, _, _), weights in h_weights.items():
                out_file.cd(f"centrality_{cent_min}_{cent_max}")
                cent_dir = out_file.GetDirectory(f"centrality_{cent_min}_{cent_max}")
                # create a directory for the centrality range
                if not cent_dir.Get(f"{particle}_{origin}"):
                    cent_dir.mkdir(f"{particle}_{origin}")
                cent_dir.cd(f"{particle}_{origin}")
                h_weight, (_, h_mc), (info_data, h_data) = weights
                if info_data['task_name'] == task_name:
                    h_weight.Write(f"{info_data['h_name']}_weights_{particle}_{origin}")
                    h_mc.Write(f"{info_data['h_name']}_mc_{particle}_{origin}")
                    h_data.Write(f"{info_data['h_name']}_data_{particle}_{origin}")
                    if (particle, origin) not in canvases:
                        def get_title(h_name):
                            return {
                                "hNPvContribAll": "All events",
                                "hNPvContribCands": "Events with a candidate",
                                "hNPvContribCandsInSignalRegionDs": "Events with cand in mass",
                                "hNPvContribCandsInSignalRegionDplus": "Events with cand in mass"
                            }[h_name]

                        def get_col(h_name):
                            return {
                                "hNPvContribAll": 1,
                                "hNPvContribCands": 2,
                                "hNPvContribCandsInSignalRegionDs": 3,
                                "hNPvContribCandsInSignalRegionDplus": 3
                            }[h_name]

                        canvases[particle, origin] = ROOT.TCanvas(
                            f"canvas_{task_name}_{particle}_{origin}_{cent_min}_{cent_max}",
                            f"canvas_{task_name}_{particle}_{origin}_{cent_min}_{cent_max}",
                            800, 600
                        )
                        canvases[particle, origin].Divide(3, 2)
                    if "InSignalRegion" in info_data['h_name']:
                        if particle not in info_data['h_name']:
                            continue
                    canvases[particle, origin].cd(get_col(info_data['h_name']))
                    min_y = min(h_data.GetMinimum(), h_mc.GetMinimum())
                    max_y = max(h_data.GetMaximum(), h_mc.GetMaximum())
                    min_x = h_data.GetXaxis().GetBinLowEdge(1)
                    max_x = h_data.GetXaxis().GetBinLowEdge(h_data.GetNbinsX()+1)
                    ROOT.gPad.DrawFrame(
                        min_x, min_y*0.9, max_x, max_y*1.1,
                        f"{get_title(info_data['h_name'])};N_{{PV}};Counts"
                    )
                    h_data.Draw("same")
                    h_mc.Draw("same")
                    canvases[particle, origin].cd(get_col(info_data['h_name']) + 3)
                    h_weight.Draw()
            for key, canvas in canvases.items():
                out_file.cd(f"centrality_{cent_min}_{cent_max}/{key[0]}_{key[1]}")
                canvas.Write()


def main(cfg_file_name):
    """
    Produce multiplicity weights.

    Parameters:
    - cfg_file_name (str): Path to the configuration file in YAML format.
    """
    with open(cfg_file_name, 'r', encoding="utf-8") as cfg_file:
        cfg = yaml.safe_load(cfg_file)

    ROOT.gStyle.SetOptLogy()

    histos_data = get_histos(cfg, is_mc=False)
    histos_mc = get_histos(cfg, is_mc=True)
    histos_projections_mc = get_projections(histos_mc)

    # Create output file
    for suffix in cfg['output']['suffix']:
        out_file_name = os.path.join(cfg['output']['dir'], f"weights_{suffix}.root")
        out_file = ROOT.TFile(out_file_name, "RECREATE")
        out_file.Close()

    for cent_min, cent_max in zip(CENT_MINS, CENT_MAXS):
        histos_projection_data = get_projections(histos_data, cent_min, cent_max)
        h_weights = get_weights(histos_projections_mc, histos_projection_data, cfg)

        dump_to_root(cfg, h_weights, cent_min, cent_max)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce multiplicity weights')
    parser.add_argument('config_file', type=str, help='Config file name')
    args = parser.parse_args()

    ROOT.TH1.AddDirectory(False)
    ROOT.TH2.AddDirectory(False)

    main(args.config_file)
