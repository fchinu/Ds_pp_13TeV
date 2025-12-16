'''
This script is used to pre-process a/multi large AnRes.root for the BDT training:
    - split the input by pT
    - obtain the sigma from prompt enhance sample
python3 pre_process.py config.yml AnRes_1.root AnRes_2.root --pre --sigma  
'''
import os
import sys
import yaml
import array
from ROOT import TFile, TObject
import argparse
import gc
import concurrent.futures
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f"{script_dir}/")
sys.path.append(f"{script_dir}/../utils/")
from utils import make_dir_root_file, logger

def check_existing_outputs(file_path):
    """
    Check ih_ev.GetName()ut file already exists, and open it accordingly.

    Args:
        file_path (str): path to the output ROOT file.

    Returns:
        tuple: (out_file, write_opt) where out_file is the opened TFile and write_opt is the write option.
    """
    if os.path.exists(file_path):
        logger(f"    Updating file: {file_path}", level='WARNING')
        out_file = TFile(file_path, 'update')
        write_opt = TObject.kOverwrite
    else:
        logger(f"    Creating file: {file_path}", level='WARNING')
        out_file = TFile.Open(file_path, 'recreate')
        write_opt = TObject.kSingleKey # Standard

    return out_file, write_opt

def get_inputs(file, full_cfg, sparse_cfg, debug=False):
    """Load a single sparse and axes info
    
    Args:
        file (TFile): input ROOT file
        full_cfg (dict): full configuration dictionary
        sparse_cfg (dict): sparse configuration dictionary
        debug (bool, optional): print debug info. Defaults to False.
    
    Returns:
        tuple: (sparse, axes) where sparse is the loaded sparse histogram and axes is the dictionary of axes information
    """
    axes_reco = {
        'Mass': 0,
        'Pt': 1,
        'Cent': 2,
        'ScoreBkg': 3,
        'ScorePrompt': 4,
        'ScoreFD': 5,
        'NPvContributors': 6,
        'PtB': 7,
        'BHadronFlag': 8
    }
    axes_gen = {
        'Pt': 0,
        'Y': 1,
        'NPvContributors': 2,
        'Centrality': 3,
        'PtB': 4,
        'BHadronFlag': 5
    }
    sparse = file.Get(sparse_cfg['path'])
    if sparse is None:
        logger(f"Sparse {sparse_cfg['type']} not found in file {file.GetName()} at path {sparse_cfg['path']}", level='ERROR')
    else:
        logger(f"Sparse {sparse} loaded from {sparse_cfg['path']}", level='INFO')

    return sparse, axes_reco if sparse_cfg['type'] == 'reco' else axes_gen

def load_event_collision_histograms_from_task(input_cfg, infile):
    """
    Load event histogram from a ROOT file based on the provided configuration.

    Args:
    - input_cfg (dict): input configuration dictionary
    - infile (TFile): input ROOT file

    Returns:
    - h_ev (TH1F): Event histogram loaded from the ROOT file.
    """
    h_ev = infile.Get(input_cfg["event_histogram"])
    h_coll = infile.Get(input_cfg["collision_histogram"])
    return h_ev, h_coll

def process_sparse(i_file, infile, full_cfg, sparse_cfg, prep_out_dir, input_cfg):
    """
    Process a single sparse from an input file for all pt bins according to the configuration.
    
    Args:
        i_file (int): index of the input file
        infile (TFile): input ROOT file
        full_cfg (dict): full configuration dictionary
        sparse_cfg (dict): sparse configuration dictionary
        prep_out_dir (str): output directory for pre-processed files
        input_cfg (dict): input configuration dictionary
    """

    logger(f'Processing file {i_file}, {infile.GetName()}')
    sparse, axes = get_inputs(infile, full_cfg, sparse_cfg, True)

    h_ev, h_coll = load_event_collision_histograms_from_task(input_cfg, infile)

    # Only for flow with SP, not for correlations (applied in O2Physics)
    if axes.get('Cent') is not None:
        cent_min, cent_max = full_cfg['centrality']
        logger(f"Applying cent cut to sparse {sparse} with value {cent_min} -- {cent_max}", "INFO")
        sparse.GetAxis(axes['Cent']).SetRangeUser(cent_min, cent_max)

    pt_mins, pt_maxs = full_cfg['ptbins'][:-1], full_cfg['ptbins'][1:]
    bkg_maxs = full_cfg['preprocess']['bkg_cuts']
    axes_to_keep, rebin = [list(axes.values())[:sparse.GetNdimensions()]], [1]*sparse.GetNdimensions()
    sparse_type, sparse_path = sparse_cfg['type'], sparse_cfg['path']
    sparse_dir, sparse_name = sparse_path.rsplit('/', 1)[0], sparse_path.split('/')[-1]

    logger(f"Projecting sparse {sparse_cfg['type']} for file {i_file} into pT bins ({pt_mins} - {pt_maxs}) with bkg cuts {bkg_maxs}", level='INFO')
    for pt_min, pt_max, bkg_max in zip(pt_mins, pt_maxs, bkg_maxs):
        logger(f"Processing pT bin {pt_min} - {pt_max} with bkg max {bkg_max}", level='INFO')
        # Create output file
        out_file_dir = f"{prep_out_dir}/Preprocess/{input_cfg['outdir']}/pt_{int(pt_min*10)}_{int(pt_max*10)}/jobs"
        os.makedirs(out_file_dir, exist_ok=True)
        out_file_path = f'{out_file_dir}/AnalysisResults_{i_file}.root'
        out_file, write_opt = check_existing_outputs(out_file_path)

        sparse.GetAxis(axes.get('PtTrig', axes.get('Pt'))).SetRangeUser(pt_min, pt_max) # PtTrig for correlations, Pt for SP flow
        if axes.get('ScoreBkg') is not None: # Skip sparses for generated info
            sparse.GetAxis(axes['ScoreBkg']).SetRangeUser(0, bkg_max)
        make_dir_root_file(sparse_dir, out_file)
        out_file.cd(sparse_dir)
        sparse.Write(sparse_name, write_opt)
        h_ev.Write(h_ev.GetName().split('/')[-1], write_opt)
        h_coll.Write(h_coll.GetName(), write_opt)
        out_file.Delete(sparse_name + ";*")

        gc.collect()

        out_file.Close()
        logger(f'----> Finished processing pT bin {pt_min} - {pt_max} for {i_file}, sparse: {sparse_cfg["type"]}\n', "INFO")

    del sparse
    gc.collect()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument('config', metavar='text', default='config.yml', help='configuration file')
    parser.add_argument("--workers", "-w", type=int, default=1, help="number of workers")
    args = parser.parse_args()

    with open(args.config, 'r') as cfg_pre:
        full_cfg = yaml.safe_load(cfg_pre)

    output_dir = full_cfg['outdir']

    for input_cfg in full_cfg['preprocess']['inputs']:
        files = [TFile.Open(fp, 'r') for fp in input_cfg['files']]
        for sparse_cfg in input_cfg['sparses']:
            logger(f"##### Skimming {sparse_cfg['type']} #####", "WARNING")
            with concurrent.futures.ThreadPoolExecutor(args.workers) as executor:
                tasks_data = [executor.submit(process_sparse, i_file, file, full_cfg, sparse_cfg, output_dir, input_cfg) for i_file, file in enumerate(files)]
        [TFile.Close(file) for file in files]

        for task in concurrent.futures.as_completed(tasks_data):
            task.result()

        # use hadd to merge the results in the directories for the different pT bins
        for directory in os.listdir(f"{output_dir}/Preprocess/{input_cfg['outdir']}"):
            logger(f"Merging files in {output_dir}/Preprocess/{input_cfg['outdir']}/{directory}/jobs ...", "INFO")
            prep_dir = f"{output_dir}/Preprocess/{input_cfg['outdir']}/{directory}"
            files_to_merge = [f"./jobs/{file}" for file in os.listdir(f"{prep_dir}/jobs") if file.startswith("AnalysisResults_") and file.endswith(".root")]
            files_to_merge_str = ' '.join(files_to_merge)
            merged_file = f"{prep_dir}/AnalysisResults_{directory}.root"
            log_merge = f"{prep_dir}/log_merge.txt"
            hadd = f"hadd -f {merged_file} {files_to_merge_str} > {log_merge}"
            logger(f"\n\nRunning command: {hadd}", "INFO")
            try:
                os.system(f"cd {prep_dir} && {hadd}")
            except Exception as e:
                logger(f"Error while creating hadd command: {e}", "ERROR")
                os.system(f"rm {merged_file}") # Remove the partially merged file
                os.system(f"mv {log_merge} {prep_dir}/log_merge_error.txt")

        logger(f"Finished processing\n\n", "INFO")
