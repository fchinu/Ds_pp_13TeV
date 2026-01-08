"""Compute FD fraction systematics vs multiplicity."""

from pathlib import Path
import argparse
import itertools
import sys
import yaml
import ROOT
import numpy as np

sys.path.append(str(Path(__file__).parent.parent.parent.parent / "FD_Fraction" / "data_driven"))
sys.path.append(str(Path(__file__).parent.parent.parent.parent / "utils"))
from plot_utils import get_discrete_from_continuous_matplotlib_palette  # pylint: disable=wrong-import-position,import-error

# pylint: disable=no-member # ROOT module members

colors, _ = get_discrete_from_continuous_matplotlib_palette("magma", 60)
colors = colors[5::]  # Take some sparse colors

ROOT.TH1.AddDirectory(False)

def draw_figures(fracs, indir: Path):
    """
    Draw figures for the fractions stored in the fracs dictionary.
    Args:
        fracs (dict): dictionary containing fraction histograms
        indir (Path): input directory path
    """
    for particle, cents in fracs.items():
        for cent_range, trials in cents.items():
            c = ROOT.TCanvas(
                f"c_{particle}_{cent_range[0]}_{cent_range[1]}",
                f"c_{particle}_{cent_range[0]}_{cent_range[1]}"
            )
            c.DrawFrame(
                0, 0.5, 24, 1,
                f";#it{{p}}_{{T}} (GeV/#it{{c}});"
                f"Prompt {'D_{s}^{+}' if particle == 'Ds' else 'D^{+}'} fraction"
            )
            leg = ROOT.TLegend(0.6, 0.2, 0.89, 0.59)
            leg.SetTextSize(0.03)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            for i_trial, (trial, h_frac) in enumerate(trials.items()):
                h_frac.SetLineWidth(2)
                h_frac.SetMarkerStyle(ROOT.kFullCircle)
                h_frac.SetLineColor(colors[i_trial])
                h_frac.SetMarkerColor(colors[i_trial])
                h_frac.Draw("same")
                leg.AddEntry(h_frac, trial, "lp")
            leg.Draw()
            outname = f"fraction_{particle}_cent_{cent_range[0]}_{cent_range[1]}"
            c.SaveAs(str(indir / f"{outname}.pdf"))

def load_inputs(cfg):
    """
    Load input histograms based on the configuration.

    Args:
        cfg (dict): configuration dictionary

    Returns:
        fracs (dict): dictionary containing loaded fraction histograms
        indir (Path): input directory path
    """
    indir = Path(cfg["trial_dir"])
    infiles = indir.rglob("*.root")
    fracs = {"Ds": {}, "Dplus": {}}
    for file in infiles:
        particle = file.parent.name
        trial = file.parent.parent.name
        if file.stem.split('_')[-1] == "MB":
            cent_min, cent_max = "0", "100"
        else:
            cent_min, cent_max = file.stem.split('_')[-2], file.stem.split('_')[-1]

        with ROOT.TFile.Open(str(file)) as infile:
            h_frac = infile.Get("hRawFracPrompt_cent_0_10")
        if (cent_min, cent_max) not in fracs[particle]:
            fracs[particle][(cent_min, cent_max)] = {}
        fracs[particle][(cent_min, cent_max)][trial] = h_frac
    return fracs

def dump_uncertainties(fracs, cfg): # pylint: disable=too-many-locals
    """
    Dump uncertainties and assigned systematic errors based on the fraction histograms.
    Args:
        fracs (dict): dictionary containing fraction histograms
        cfg (dict): configuration dictionary
    """
    outdir = Path(cfg["output_dir"])
    outdir.mkdir(parents=True, exist_ok=True)
    with ROOT.TFile.Open(str(outdir / "syst.root"), "RECREATE") as f_out:
        pass  # Just to create/clear the file

    for ds_items in fracs["Ds"].items():
        cent_range, ds_trials = ds_items
        dplus_trials = fracs["Dplus"][cent_range]

        if cent_range == ('0', '100'):
            # We don't look at the fraction
            continue

        ratios = {}
        for trial in ds_trials.keys():
            ds_hist = ds_trials[trial]
            dplus_hist = dplus_trials[trial]
            ratio_hist = ds_hist.Clone(f"ratio_{trial}_cent_{cent_range[0]}_{cent_range[1]}")
            ratio_hist.Divide(dplus_hist)
            ratios[trial] = ratio_hist

        rms_pt = []
        for i_pt in range(1, ratios['all'].GetNbinsX() + 1):
            rms_pt.append([])
            for _, trial in ratios.items():
                rms_pt[-1].append(trial.GetBinContent(i_pt) / ratios['central'].GetBinContent(i_pt))
        rms_pt = np.std(np.array(rms_pt), axis=1)

        with ROOT.TFile.Open(str(outdir / "syst.root"), "UPDATE") as f_out:
            h_rms = ratios['all'].Clone(
                f"hRMS_DsDplus_cent_{cent_range[0]}_{cent_range[1]}"
            )
            for i_pt in range(1, ratios['all'].GetNbinsX() + 1):
                h_rms.SetBinContent(i_pt, rms_pt[i_pt - 1])
                h_rms.SetBinError(i_pt, 0)
            h_rms.Write()

            h_assigned_syst = ratios['all'].Clone(
                f"assigned_syst_{cent_range[0]}_{cent_range[1]}"
            )

            for i_pt in range(1, ratios['all'].GetNbinsX() + 1):
                assigned_syst = cfg["assigned_syst"][f"{cent_range[0]}-{cent_range[1]}"][i_pt - 1]
                h_assigned_syst.SetBinContent(i_pt, assigned_syst)
                h_assigned_syst.SetBinError(i_pt, 0)
            h_assigned_syst.Write()

            c_trials = ROOT.TCanvas(
                f"c_ratios_cent_{cent_range[0]}_{cent_range[1]}",
                f"c_ratios_cent_{cent_range[0]}_{cent_range[1]}"
            )
            c_trials.DrawFrame(
                0, 0, 24, 2,
                ";#it{p}_{T} (GeV/#it{c});"
                "D_{s}^{+}/D^{+} prompt fraction ratio"
            )
            leg = ROOT.TLegend(0.6, 0.7, 0.89, 0.89)
            leg.SetTextSize(0.03)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            for i_trial, (trial, ratio_hist) in enumerate(ratios.items()):
                ratio_hist.SetLineWidth(2)
                ratio_hist.SetMarkerStyle(ROOT.kFullCircle)
                ratio_hist.SetLineColor(colors[i_trial])
                ratio_hist.SetMarkerColor(colors[i_trial])
                ratio_hist.Draw("same")
                leg.AddEntry(ratio_hist, trial, "lp")
            leg.Draw()
            c_trials.Write()

        ratios_all_comb = {}
        for trial_ds, trial_dplus in itertools.product(ds_trials.keys(), dplus_trials.keys()):
            ds_hist = ds_trials[trial_ds]
            dplus_hist = dplus_trials[trial_dplus]
            ratio_hist = ds_hist.Clone(f"ratio_{trial_ds}_{trial_dplus}_cent_{cent_range[0]}_{cent_range[1]}")
            ratio_hist.Divide(dplus_hist)
            ratios_all_comb[f"{trial_ds}_{trial_dplus}"] = ratio_hist

        rms_pt_all_comb = []
        for i_pt in range(1, ratios_all_comb['all_all'].GetNbinsX() + 1):
            rms_pt_all_comb.append([])
            for _, trial in ratios_all_comb.items():
                rms_pt_all_comb[-1].append(trial.GetBinContent(i_pt))
        central_array = np.array([ratios_all_comb['central_central'].GetBinContent(i_pt) for i_pt in range(1, ratios_all_comb['all_all'].GetNbinsX() + 1)])
        rms_pt_all_comb = np.sqrt(np.std(np.array(rms_pt_all_comb), axis=1)**2 + (np.mean(np.array(rms_pt_all_comb), axis=1)-central_array)**2) / central_array

        with ROOT.TFile.Open(str(outdir / "syst.root"), "UPDATE") as f_out:
            h_rms = ratios_all_comb['all_all'].Clone(
                f"hRMS_all_comb_DsDplus_cent_{cent_range[0]}_{cent_range[1]}"
            )
            for i_pt in range(1, ratios_all_comb['all_all'].GetNbinsX() + 1):
                h_rms.SetBinContent(i_pt, rms_pt_all_comb[i_pt - 1])
                h_rms.SetBinError(i_pt, 0)
            h_rms.Write()

            c_trials = ROOT.TCanvas(
                f"c_ratios_cent_{cent_range[0]}_{cent_range[1]}",
                f"c_ratios_cent_{cent_range[0]}_{cent_range[1]}"
            )
            c_trials.DrawFrame(
                0, 0, 24, 2,
                ";#it{p}_{T} (GeV/#it{c});"
                "D_{s}^{+}/D^{+} prompt fraction ratio"
            )
            leg = ROOT.TLegend(0.6, 0.7, 0.89, 0.89)
            leg.SetTextSize(0.03)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            for i_trial, (trial, ratio_hist) in enumerate(ratios_all_comb.items()):
                ratio_hist.SetLineWidth(2)
                ratio_hist.SetMarkerStyle(ROOT.kFullCircle)
                ratio_hist.SetLineColor(colors[i_trial])
                ratio_hist.SetMarkerColor(colors[i_trial])
                ratio_hist.Draw("same")
                leg.AddEntry(ratio_hist, trial, "lp")

            ratios_all_comb["all_all"].SetLineWidth(3)
            ratios_all_comb["all_all"].SetMarkerStyle(ROOT.kFullCircle)
            ratios_all_comb["all_all"].SetLineColor(ROOT.kBlack)
            ratios_all_comb["all_all"].SetMarkerColor(ROOT.kBlack)
            ratios_all_comb["all_all"].Draw("same")
            leg.Draw()
            c_trials.Write()

            c_trials.SaveAs(str(outdir / f"ratio_all_comb_cent_{cent_range[0]}_{cent_range[1]}.pdf"))

def dump_outputs(cfg):  # pylint: disable=too-many-locals
    """
    Produce a summary of outputs in the input directory.
    Args:
        cfg (dict): configuration dictionary
    """
    fracs = load_inputs(cfg)

    # set all as last
    for particle, cents in fracs.items():
        for cent_range, trials in cents.items():
            all_trial = trials.pop('all')
            cents[cent_range] = {**trials, 'all': all_trial}

    # We add the central case as well
    for particle, cents in fracs.items():
        for cent_range in cents.keys():  # pylint: disable=consider-using-dict-items
            if cent_range == ('0', '100'):
                continue
            filename = (
                Path(cfg["central_frac"]) / particle /
                f"CutVar{particle}_pp13TeV_{cent_range[0]}_{cent_range[1]}.root"
            )
            filename = str(filename)
            with ROOT.TFile.Open(filename) as infile:
                h_frac = infile.Get(f"hRawFracPrompt_cent_{cent_range[0]}_{cent_range[1]}")
                cents[cent_range]['central'] = h_frac

    draw_figures(fracs, Path(cfg["trial_dir"]))
    dump_uncertainties(fracs, cfg)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute FD fraction systematics vs multiplicity"
    )
    parser.add_argument("config", type=str, help="Path to the config_syst.yml file")
    args = parser.parse_args()

    with open(args.config, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)

    dump_outputs(config)
