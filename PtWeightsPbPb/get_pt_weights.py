"""
Script for evaluation of pt weights (interpolation)
"""

import argparse
import sys
import pandas as pd
import ROOT
from utils import read_tamu, read_lido, read_fonll, get_centrality_interpolation, get_fonll_times_raa

def compute_weights(hists_dndpt, hist_mcgen, base_name):
    """
    Function to compute FONLL x RAA
    """
    hists = {}
    for cent_key, hist_dndpt in hists_dndpt.items():
        hists[cent_key] = hist_dndpt.Clone(f"{base_name}_{cent_key}")
        hists[cent_key].Divide(hist_mcgen)

    return hists


def get_pt_weights(model):
    """
    Main function
    """

    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetPadTopMargin(0.035)
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetTitleSize(0.045)
    ROOT.gStyle.SetLabelSize(0.045)
    ROOT.gStyle.SetOptStat(0)

    # read MC gen
    infile_mc = ROOT.TFile.Open("GenPtShapes_LHC25a5.root")
    hist_ds_y_vs_pt = infile_mc.Get("hf-task-mc-gen-pt-rap-shapes/CharmHadrons/hRapVsPtPrompt431")
    hist_d_y_vs_pt = infile_mc.Get("hf-task-mc-gen-pt-rap-shapes/CharmHadrons/hRapVsPtPrompt411")
    hist_bs_y_vs_pt = infile_mc.Get("hf-task-mc-gen-pt-rap-shapes/BeautyHadrons/hRapVsPt531")
    hist_b_y_vs_pt = infile_mc.Get("hf-task-mc-gen-pt-rap-shapes/BeautyHadrons/hRapVsPt511")
    hist_ds_y_vs_pt.SetDirectory(0)
    hist_d_y_vs_pt.SetDirectory(0)
    hist_bs_y_vs_pt.SetDirectory(0)
    hist_b_y_vs_pt.SetDirectory(0)
    y_min_d = hist_ds_y_vs_pt.GetYaxis().FindBin(-0.4999)
    y_max_d = hist_ds_y_vs_pt.GetYaxis().FindBin(0.4999)
    hist_ds_mcgen = hist_ds_y_vs_pt.ProjectionX("hist_ds_mcgen", y_min_d, y_max_d)
    hist_d_mcgen = hist_d_y_vs_pt.ProjectionX("hist_d_mcgen", y_min_d, y_max_d)
    hist_bs_mcgen = hist_bs_y_vs_pt.ProjectionX("hist_bs_mcgen")
    hist_b_mcgen = hist_b_y_vs_pt.ProjectionX("hist_b_mcgen")
    for hist in [hist_ds_mcgen, hist_d_mcgen, hist_bs_mcgen, hist_b_mcgen]:
        hist.Scale(1./hist.Integral())
        hist.SetLineColor(ROOT.kBlack)
        hist.SetMarkerColor(ROOT.kBlack)
        hist.SetLineWidth(2)

    # read fonll
    hist_fonll_ds = read_fonll("FONLL_Ds_5dot36TeV.txt", "hist_fonll_ds")
    hist_fonll_d = read_fonll("FONLL_Dp_5dot36TeV.txt", "hist_fonll_d")
    hist_fonll_b = read_fonll("FONLL_B_5dot36TeV.txt", "hist_fonll_b")

    # read RAA models
    gr_ds_raa, gr_d_raa, gr_b_raa, gr_bs_raa = {}, {}, {}, {}
    if model == "tamu":
        gr_d_raa["0_10"] = read_tamu("PromptD0_TAMU_RAA_5TeV_010.txt", "gr_d_raa_tamu_0_10")
        gr_d_raa["30_50"] = read_tamu("PromptD0_TAMU_RAA_5TeV_3050.txt", "gr_d_raa_tamu_30_50")
        gr_ds_raa["0_10"] = read_tamu("PromptDs_TAMU_RAA_5TeV_010.txt", "gr_ds_raa_tamu_0_10")
        gr_ds_raa["30_50"] = read_tamu("PromptDs_TAMU_RAA_5TeV_3050.txt", "gr_ds_raa_tamu_30_50")
        gr_b_raa["0_10"] = read_tamu("B_TAMU_RAA_5TeV_010.txt", "gr_b_raa_tamu_0_10", True)
        gr_b_raa["30_50"] = read_tamu("B_TAMU_RAA_5TeV_3050.txt", "gr_b_raa_tamu_30_50", True)
        gr_bs_raa["0_10"] = read_tamu("Bs_TAMU_RAA_5TeV_010.txt", "gr_bs_raa_tamu_0_10", True)
        gr_bs_raa["30_50"] = read_tamu("Bs_TAMU_RAA_5TeV_3050.txt", "gr_bs_raa_tamu_30_50", True)

        # assume most peripheral to be = 0.8
        gr_ds_raa["80_90"] = gr_ds_raa["0_10"].Clone("gr_ds_raa_tamu_80_90")
        for ipt in range(gr_ds_raa["80_90"].GetN()):
            gr_ds_raa["80_90"].SetPoint(ipt, gr_ds_raa["80_90"].GetPointX(ipt), 0.8)
        gr_d_raa["80_90"] = gr_d_raa["0_10"].Clone("gr_d_raa_tamu_80_90")
        for ipt in range(gr_d_raa["80_90"].GetN()):
            gr_d_raa["80_90"].SetPoint(ipt, gr_d_raa["80_90"].GetPointX(ipt), 0.8)
        gr_b_raa["80_90"] = gr_b_raa["0_10"].Clone("gr_b_raa_tamu_80_90")
        for ipt in range(gr_b_raa["80_90"].GetN()):
            gr_b_raa["80_90"].SetPoint(ipt, gr_b_raa["80_90"].GetPointX(ipt), 0.8)
        gr_bs_raa["80_90"] = gr_bs_raa["0_10"].Clone("gr_bs_raa_tamu_80_90")
        for ipt in range(gr_bs_raa["80_90"].GetN()):
            gr_bs_raa["80_90"].SetPoint(ipt, gr_bs_raa["80_90"].GetPointX(ipt), 0.8)

    if model == "lido":
        gr_d_raa["0_10"] = read_lido("lido/D-meson-Raa-010.dat", "gr_d_raa_lido_0_10")
        gr_d_raa["30_50"] = read_lido("lido/D-meson-Raa-3050.dat", "gr_d_raa_lido_30_50")
        gr_b_raa["0_10"] = read_lido("lido/B-meson-Raa-010.dat", "gr_b_raa_lido_0_10")
        gr_b_raa["30_50"] = read_lido("lido/B-meson-Raa-3050.dat", "gr_b_raa_lido_30_50")

        # assume most peripheral to be = 0.8
        gr_d_raa["80_90"] = gr_d_raa["0_10"].Clone("gr_d_raa_lido_80_90")
        for ipt in range(gr_d_raa["80_90"].GetN()):
            gr_d_raa["80_90"].SetPoint(ipt, gr_d_raa["80_90"].GetPointX(ipt), 0.8)

        gr_b_raa["80_90"] = gr_b_raa["0_10"].Clone("gr_b_raa_lido_80_90")
        for ipt in range(gr_b_raa["80_90"].GetN()):
            gr_b_raa["80_90"].SetPoint(ipt, gr_b_raa["80_90"].GetPointX(ipt), 0.8)

    # interpolate to obtain all centralities
    get_centrality_interpolation(gr_d_raa, f"gr_d_raa_{model}")
    get_centrality_interpolation(gr_b_raa, f"gr_b_raa_{model}")
    if model != "lido": # only D and B for lido
        get_centrality_interpolation(gr_ds_raa, f"gr_ds_raa_{model}")
        get_centrality_interpolation(gr_bs_raa, f"gr_bs_raa_{model}")

    hist_dndpt_d = get_fonll_times_raa(hist_fonll_d, gr_d_raa, f"hist_d_dndpt_{model}")
    hist_dndpt_b = get_fonll_times_raa(hist_fonll_b, gr_b_raa, f"hist_b_dndpt_{model}")
    if model != "lido": # only D and B for lido
        hist_dndpt_ds = get_fonll_times_raa(hist_fonll_ds, gr_ds_raa, f"hist_ds_dndpt_{model}")
        hist_dndpt_bs = get_fonll_times_raa(hist_fonll_b, gr_bs_raa, f"hist_bs_dndpt_{model}")

    hist_weights_d = compute_weights(hist_dndpt_d, hist_d_mcgen, f"hist_d_weights_{model}")
    hist_weights_b = compute_weights(hist_dndpt_b, hist_b_mcgen, f"hist_b_weights_{model}")
    if model != "lido": # only D and B for lido
        hist_weights_ds = compute_weights(hist_dndpt_ds, hist_ds_mcgen, f"hist_ds_weights_{model}")
        hist_weights_bs = compute_weights(hist_dndpt_bs, hist_bs_mcgen, f"hist_bs_weights_{model}")

    colors = {
        "0_10": ROOT.kRed+2,
        "10_20": ROOT.kRed+1,
        "20_30": ROOT.kOrange+9,
        "30_40": ROOT.kOrange+7,
        "40_50": ROOT.kSpring+2,
        "50_60": ROOT.kGreen+2, 
        "60_70": ROOT.kAzure+2,
        "70_80": ROOT.kAzure+4,
        "80_90": ROOT.kBlue+2
    }

    outfile = ROOT.TFile(f"ptweights_{model}.root", "recreate")
    for hist in [hist_ds_mcgen, hist_d_mcgen, hist_bs_mcgen, hist_b_mcgen]:
        hist.Write()

    outfile.mkdir("raa")
    canv_raa = ROOT.TCanvas("canv_raa", "", 1000, 1000)
    canv_raa.Divide(2, 2)
    canv_raa.cd(1).DrawFrame(0., 0., 12., 2., ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}(D)")
    leg = ROOT.TLegend(0.2, 0.6, 0.9, 0.9)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(2)

    for cent, graph in gr_d_raa.items():
        if cent == "30_50":
            continue
        graph.SetLineColor(colors[cent])
        graph.SetLineWidth(2)
        graph.SetMarkerColor(colors[cent])
        cent_min = int(cent.split(sep="_")[0])
        cent_max = int(cent.split(sep="_")[1])
        leg.AddEntry(graph, f"{cent_min}#minus{cent_max}%", "l")
        graph.Draw()
        outfile.cd("raa")
        graph.Write()
    leg.Draw()
    canv_raa.cd(2).DrawFrame(0., 0., 12., 2., ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}(D_{s}^{#plus})")
    if model != "lido": # only D and B for lido
        for cent, graph in gr_ds_raa.items():
            if cent == "30_50":
                continue
            graph.SetLineColor(colors[cent])
            graph.SetLineWidth(2)
            graph.SetMarkerColor(colors[cent])
            graph.Draw()
            outfile.cd("raa")
            graph.Write()
    canv_raa.cd(3).DrawFrame(0., 0., 24., 3., ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}(B)")
    for cent, graph in gr_b_raa.items():
        if cent == "30_50":
            continue
        graph.SetLineColor(colors[cent])
        graph.SetLineWidth(2)
        graph.SetMarkerColor(colors[cent])
        graph.Draw()
        outfile.cd("raa")
        graph.Write()
    canv_raa.cd(4).DrawFrame(0., 0., 24., 3., ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}(B_{s}^{0})")
    if model != "lido": # only D and B for lido
        for cent, graph in gr_bs_raa.items():
            if cent == "30_50":
                continue
            graph.SetLineColor(colors[cent])
            graph.SetLineWidth(2)
            graph.SetMarkerColor(colors[cent])
            graph.Draw()
            outfile.cd("raa")
            graph.Write()
    canv_raa.Write()

    outfile.cd()
    outfile.mkdir("dndpt")
    canv_dndpt = ROOT.TCanvas("canv_dndpt", "", 1000, 1000)
    canv_dndpt.Divide(2, 2)
    canv_dndpt.cd(1).DrawFrame(0., 1.e-7, 24., 1.e-1, ";#it{p}_{T} (GeV/#it{c});d#it{N}(D)/d#it{p}_{T} (a.u.)")
    canv_dndpt.cd(1).SetLogy()
    hist_d_mcgen.Draw("histsame")
    for cent, hist in hist_dndpt_d.items():
        hist.SetLineColor(colors[cent])
        hist.SetLineWidth(2)
        hist.SetMarkerColor(colors[cent])
        hist.Draw("histsame")
        outfile.cd("dndpt")
        hist.Write()
    canv_dndpt.cd(2).DrawFrame(0., 1.e-7, 24., 1.e-1, ";#it{p}_{T} (GeV/#it{c});d#it{N}(D_{s}^{#plus})/d#it{p}_{T} (a.u.)")
    canv_dndpt.cd(2).SetLogy()
    hist_ds_mcgen.Draw("histsame")
    if model != "lido": # only D and B for lido
        for cent, hist in hist_dndpt_ds.items():
            hist.SetLineColor(colors[cent])
            hist.SetLineWidth(2)
            hist.SetMarkerColor(colors[cent])
            hist.Draw("histsame")
            outfile.cd("dndpt")
            hist.Write()
    canv_dndpt.cd(3).DrawFrame(0., 1.e-7, 50., 1.e-1, ";#it{p}_{T} (GeV/#it{c});d#it{N}(B)/d#it{p}_{T} (a.u.)")
    canv_dndpt.cd(3).SetLogy()
    hist_b_mcgen.Draw("histsame")
    for cent, hist in hist_dndpt_b.items():
        hist.SetLineColor(colors[cent])
        hist.SetLineWidth(2)
        hist.SetMarkerColor(colors[cent])
        hist.Draw("histsame")
        outfile.cd("dndpt")
        hist.Write()
    leg_dndpt = leg.Clone()
    leg_dndpt.SetY1NDC(0.2)
    leg_dndpt.SetY2NDC(0.45)
    leg_dndpt.AddEntry(hist_ds_mcgen, "MC gen", "l")
    leg_dndpt.Draw()
    canv_dndpt.cd(4).DrawFrame(0., 1.e-7, 50., 1.e-1, ";#it{p}_{T} (GeV/#it{c});d#it{N}(B_{s}^{0})/d#it{p}_{T} (a.u.)")
    canv_dndpt.cd(4).SetLogy()
    hist_bs_mcgen.Draw("histsame")
    if model != "lido": # only D and B for lido
        for cent, hist in hist_dndpt_bs.items():
            hist.SetLineColor(colors[cent])
            hist.SetLineWidth(2)
            hist.SetMarkerColor(colors[cent])
            hist.Draw("histsame")
            outfile.cd("dndpt")
            hist.Write()
    canv_dndpt.Write()

    outfile.cd()
    outfile.mkdir("weights")
    canv_weights = ROOT.TCanvas("canv_weights", "", 1000, 1000)
    canv_weights.Divide(2, 2)
    canv_weights.cd(1).DrawFrame(0., 1.e-1, 24., 1.e1, ";#it{p}_{T} (GeV/#it{c});#it{p}_{T}-weights (D)")
    canv_weights.cd(1).SetLogy()
    for cent, hist in hist_weights_d.items():
        hist.SetLineColor(colors[cent])
        hist.SetLineWidth(2)
        hist.SetMarkerColor(colors[cent])
        hist.Draw("histsame")
        outfile.cd("weights")
        hist.Write()
    canv_weights.cd(2).DrawFrame(0., 1.e-1, 24., 1.e1, ";#it{p}_{T} (GeV/#it{c});#it{p}_{T}-weights (D_{s}^{#plus})")
    canv_weights.cd(2).SetLogy()
    if model != "lido": # only D and B for lido
        for cent, hist in hist_weights_ds.items():
            hist.SetLineColor(colors[cent])
            hist.SetLineWidth(2)
            hist.SetMarkerColor(colors[cent])
            hist.Draw("histsame")
        outfile.cd("weights")
        hist.Write()
    canv_weights.cd(3).DrawFrame(0., 1.e-1, 50., 1.e1, ";#it{p}_{T} (GeV/#it{c});#it{p}_{T}-weights (B)")
    canv_weights.cd(3).SetLogy()
    for cent, hist in hist_weights_b.items():
        hist.SetLineColor(colors[cent])
        hist.SetLineWidth(2)
        hist.SetMarkerColor(colors[cent])
        hist.Draw("histsame")
        outfile.cd("weights")
        hist.Write()
    leg_weights = leg.Clone()
    leg_weights.SetY1NDC(0.2)
    leg_weights.SetY2NDC(0.45)
    leg_weights.Draw()
    canv_weights.cd(4).DrawFrame(0., 1.e-1, 50., 1.e1, ";#it{p}_{T} (GeV/#it{c});#it{p}_{T}-weights (B_{s}^{0})")
    canv_weights.cd(4).SetLogy()
    if model != "lido": # only D and B for lido
        for cent, hist in hist_weights_bs.items():
            hist.SetLineColor(colors[cent])
            hist.SetLineWidth(2)
            hist.SetMarkerColor(colors[cent])
            hist.Draw("histsame")
            outfile.cd("weights")
            hist.Write()
    canv_weights.Write()
    outfile.Close()

    canv_raa.SaveAs(f"ptweights_{model}_raa.pdf")
    canv_dndpt.SaveAs(f"ptweights_{model}_dndpt.pdf")
    canv_weights.SaveAs(f"ptweights_{model}.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate pT weights")
    parser.add_argument("--model", "-m", metavar="text", help="Enabled model", default="tamu")
    args = parser.parse_args()

    models = ["tamu", "lido"]
    if args.model not in models:
        print(f"Model {args.model} not yet implemented. Available options: {models}") 
        sys.exit()

    get_pt_weights(args.model)
