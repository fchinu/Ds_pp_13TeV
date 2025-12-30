"""Script to generate variations of FD fraction calculations based on different trial selections."""
import argparse
from pathlib import Path
import json
import sys
import ROOT

sys.path.append(str(Path(__file__).parent.parent.parent.parent / "FD_Fraction" / "data_driven"))
sys.path.append(str(Path(__file__).parent.parent.parent.parent / "utils"))
from compute_fraction_cutvar import main as compute_fraction_main  # pylint: disable=wrong-import-position,import-error
from plot_utils import get_discrete_matplotlib_palette  # pylint: disable=wrong-import-position,import-error

# pylint: disable=no-member # ROOT module members

ROOT.TH1.AddDirectory(False)
colors, _ = get_discrete_matplotlib_palette("tab10")

def get_trials(inlist, option):
    """
    Get the trials from a list based on the provided option.

    Args:
    - inlist (list): list of items to filter
    - option (str): option to filter by

    Returns:
    - trials (list): filtered list of trials
    """
    switcher = {
        'all': inlist,
        'odd': inlist[::2],
        'even': inlist[1::2],
        'last_trials': inlist[len(inlist)//3:],
        'central_trials': inlist[len(inlist)//6:5*len(inlist)//6],
        'first_trials': inlist[:2*len(inlist)//3]
    }
    return switcher.get(option, inlist)

def create_config(infile: str, trial: str):
    """
    Create a configuration file for a specific trial.
    Args:
        infile (str): path to the input configuration file
        trial (str): trial option to filter by
    Returns:
        outfile (str): path to the output configuration file
    """
    # Create output directory
    out_dir = Path(__file__).parent / trial
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(infile, 'r', encoding='utf-8') as f:
        config = json.load(f)

    raw_yields_cfg = config['rawyields']['inputfiles']
    efficiency_cfg = config['efficiencies']['inputfiles']
    raw_yield_trials = get_trials(raw_yields_cfg, trial)
    efficiency_trials = get_trials(efficiency_cfg, trial)

    config['rawyields']['inputfiles'] = raw_yield_trials
    config['efficiencies']['inputfiles'] = efficiency_trials

    particle = infile.parent.name
    outdir = out_dir / particle
    outdir.mkdir(parents=True, exist_ok=True)

    config['output']['directory'] = str(outdir)

    outfile = outdir / infile.name
    with open(outfile, 'w', encoding='utf-8') as f:
        json.dump(config, f, indent=4)

    return outfile

def run(indir):
    """
    Run the trial selection based on the input configuration file.

    Args:
        infile (str): path to the input configuration file
    """
    trials = ['all', 'odd', 'even', 'last_trials', 'central_trials', 'first_trials']

    infiles = Path(indir).rglob("*.json")
    for infile in infiles:
        for trial in trials:
            outfile = create_config(infile, trial)
            compute_fraction_main(str(outfile))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get trials based on option")
    parser.add_argument("trial_dir", type=str, help="Input configuration folder")
    args = parser.parse_args()

    run(args.trial_dir)
