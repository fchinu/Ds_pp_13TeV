import argparse
from pathlib import Path
import json
import sys
sys.path.append(str(Path(__file__).parent.parent.parent.parent / "FD_Fraction" / "data_driven"))
sys.path.append(str(Path(__file__).parent.parent.parent.parent / "utils"))
from compute_fraction_cutvar import main as compute_fraction_main
import ROOT
ROOT.TH1.AddDirectory(False)
from plot_utils import get_discrete_matplotlib_palette

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

    with open(infile, 'r') as f:
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
    with open(outfile, 'w') as f:
        json.dump(config, f, indent=4)

    return outfile

def dump_outputs(indir):
    """
    Produce a summary of outputs in the input directory.
    Args:
        indir (Path): path to the input directory
    """
    infiles = indir.rglob("*.root")
    fracs = {"Ds": {}, "Dplus": {}}
    for file in infiles:
        particle = file.parent.name
        trial = file.parent.parent.name
        if file.stem.split('_')[-1] == "MB":
            cent_min, cent_max = "0", "100"
        else:
            cent_min, cent_max = file.stem.split('_')[-2], file.stem.split('_')[-1]

        with ROOT.TFile.Open(str(file)) as f:
            h_frac = f.Get("hRawFracPrompt_cent_0_10")
        if (cent_min, cent_max) not in fracs[particle]:
            fracs[particle][(cent_min, cent_max)] = {}
        fracs[particle][(cent_min, cent_max)][trial] = h_frac

    # set all as last
    for particle in fracs:
        for cent_range in fracs[particle]:
            all_trial = fracs[particle][cent_range].pop('all')
            fracs[particle][cent_range] = {**fracs[particle][cent_range], 'all': all_trial}
        

    for particle in fracs:
        for cent_range in fracs[particle]:
            c = ROOT.TCanvas(f"c_{particle}_{cent_range[0]}_{cent_range[1]}", f"c_{particle}_{cent_range[0]}_{cent_range[1]}")
            c.DrawFrame(0, 0.5, 24, 1, f";#it{{p}}_{{T}} (GeV/#it{{c}});Prompt {'D_{s}^{+}' if particle == 'Ds' else 'D^{+}'} fraction")
            leg = ROOT.TLegend(0.6, 0.2, 0.89, 0.59)
            leg.SetTextSize(0.03)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            for i_trial, (trial, h_frac) in enumerate(fracs[particle][cent_range].items()):
                h_frac.SetLineWidth(2)
                h_frac.SetMarkerStyle(ROOT.kFullCircle)
                h_frac.SetLineColor(colors[i_trial])
                h_frac.SetMarkerColor(colors[i_trial])
                h_frac.Draw("same")
                leg.AddEntry(h_frac, trial, "lp")
            leg.Draw()
            outname = f"fraction_{particle}_cent_{cent_range[0]}_{cent_range[1]}"
            c.SaveAs(str(indir / f"{outname}.pdf"))




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


    dump_outputs(Path(__file__).parent)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get trials based on option")
    parser.add_argument("trial_dir", type=str, help="Input configuration folder")
    args = parser.parse_args()

    run(args.trial_dir)