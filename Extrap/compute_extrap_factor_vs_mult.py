import argparse
import sys
import numpy as np
import uproot
import ROOT

def set_graph_style(graph, color):
    """
    """
    graph.SetLineWidth(2)
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(ROOT.kFullCircle)


def compute_extrap_factor(infile_name, outfile_name, model):
    """
    Main function for extrapolation
    """

    if model not in ["epos4hq", "epos4hq_oo", "pythia"]:
        print(f"Model {model} not supported. Exit")
        sys.exit()

    hist_corr_mult = None
    if model == "epos4hq":
        df = uproot.open(infile_name)["tree_dmesons"].arrays(library="pd")
        df_ds_midy = df.query("abs(y) < 0.5 and abs(id) == 340")
        df_dp_midy = df.query("abs(y) < 0.5 and abs(id) == 240")
        multiplicity = [0, 75, 145, 265, 437, 681, 1037, 1369, 2047]
    elif model == "epos4hq_oo":
        df = uproot.open(infile_name)["tree_dmesons"].arrays(library="pd")
        df_ds_midy = df.query("abs(y) < 0.5 and abs(id) == 340")
        df_dp_midy = df.query("abs(y) < 0.5 and abs(id) == 240")
        multiplicity = [0, 40, 80, 110, 150, 200]
    else:
        df = uproot.open(infile_name)["treeD"].arrays(library="pd")
        df_ds_midy = df.query("abs(yD) < 0.5 and abs(pdgD) == 431")
        df_dp_midy = df.query("abs(yD) < 0.5 and abs(pdgD) == 411")
        multiplicity = [0, 10, 20, 30, 40, 50, 100]
        infile_root = ROOT.TFile.Open(infile_name)
        hist_corr_mult = infile_root.Get("histMultFT0MVsMultMid05")
        hist_corr_mult.SetDirectory(0)
        infile_root.Close()

    graph_extrapfactor_ds, graph_extrapfactor_dp, graph_extrapfactor_ds_over_dp = (ROOT.TGraphErrors(1) for _ in range(3))
    graph_extrapfactor_ds.SetNameTitle("graph_extrapfactor_ds", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus}) = D_{s}^{#plus}(#it{p}_{T} integrated)/D_{s}^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_dp.SetNameTitle("graph_extrapfactor_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D^{#plus}) = D^{#plus}(#it{p}_{T} integrated)/D^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_ds_over_dp.SetNameTitle("graph_extrapfactor_ds_over_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus})/#alpha(D^{#plus})")
    set_graph_style(graph_extrapfactor_ds, ROOT.kBlack)
    set_graph_style(graph_extrapfactor_dp, ROOT.kBlack)
    set_graph_style(graph_extrapfactor_ds_over_dp, ROOT.kBlack)
    for imult, (mult_min, mult_max) in enumerate(zip(multiplicity[:-1], multiplicity[1:])):
        mult_name = "nch_mid"
        pt_name = "pt"
        if model == "pythia":
            mult_name = "multiplicity_ft0m"
            pt_name = "ptD"
        df_ds_midy_mult = df_ds_midy.query(f"{mult_min} < {mult_name} < {mult_max}")
        df_dp_midy_mult = df_dp_midy.query(f"{mult_min} < {mult_name} < {mult_max}")
        df_ds_midy_mult_pt = df_ds_midy_mult.query(f"{pt_name} > 1")
        df_dp_midy_mult_pt = df_dp_midy_mult.query(f"{pt_name} > 1")

        n_ds = len(df_ds_midy_mult)
        n_ds_ptgtr1 = len(df_ds_midy_mult_pt)
        n_dp = len(df_dp_midy_mult)
        n_dp_ptgtr1 = len(df_dp_midy_mult_pt)

        extrap_factor_ds = float(n_ds)/n_ds_ptgtr1
        extrap_factor_dp = float(n_dp)/n_dp_ptgtr1
        rho_ds = np.sqrt(float(n_ds_ptgtr1)/n_ds)
        rho_dp = np.sqrt(float(n_dp_ptgtr1)/n_dp)

        unc_extrap_factor_ds = np.sqrt(n_ds/n_ds_ptgtr1**2 + n_ds_ptgtr1**3/n_ds**4 - 2*rho_ds * n_ds_ptgtr1/n_ds**3)
        unc_extrap_factor_dp = np.sqrt(n_dp/n_dp_ptgtr1**2 + n_dp_ptgtr1**3/n_dp**4 - 2*rho_dp * n_dp_ptgtr1/n_dp**3)
        ratio = extrap_factor_ds/extrap_factor_dp
        unc_ratio = np.sqrt(unc_extrap_factor_ds**2/extrap_factor_ds**2 + unc_extrap_factor_dp**2/extrap_factor_dp**2) * ratio

        if "epos4hq" in model:
            mult_cent = (mult_max + mult_min) / 2
            mult_unc = (mult_max - mult_min) / 2
        else:
            mult_bin_min = hist_corr_mult.GetYaxis().FindBin(mult_min)
            mult_bin_max = hist_corr_mult.GetYaxis().FindBin(mult_max)
            hist_mult = hist_corr_mult.ProjectionX(f"hist_mult_{imult}", mult_bin_min, mult_bin_max)
            mult_cent = hist_mult.GetMean()
            mult_unc  = hist_mult.GetRMS()

        graph_extrapfactor_ds.SetPoint(imult, mult_cent, extrap_factor_ds)
        graph_extrapfactor_ds.SetPointError(imult, mult_unc, unc_extrap_factor_ds)
        graph_extrapfactor_dp.SetPoint(imult, mult_cent, extrap_factor_dp)
        graph_extrapfactor_dp.SetPointError(imult, mult_unc, unc_extrap_factor_dp)
        graph_extrapfactor_ds_over_dp.SetPoint(imult, mult_cent, ratio)
        graph_extrapfactor_ds_over_dp.SetPointError(imult, mult_unc, unc_ratio)

    outfile = ROOT.TFile(outfile_name, "recreate")
    graph_extrapfactor_ds.Write()
    graph_extrapfactor_dp.Write()
    graph_extrapfactor_ds_over_dp.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute extrapolation factor for Ds and D+ vs imppar')
    parser.add_argument('--input_file', '-i', metavar='text', default="EPOS4HQ_HF_VsMult_PbPb_skimmed.root", help='Input file')
    parser.add_argument('--output_file', '-o', metavar='text', default="ds_dp_extrap_factors_epos4hq.root", help='Output file')
    parser.add_argument('--model', '-m', metavar='text', default="epos4hq", help='Model (options: pythia, epos4hq)')
    args = parser.parse_args()

    compute_extrap_factor(args.input_file, args.output_file, args.model)
