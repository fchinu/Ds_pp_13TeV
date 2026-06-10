"""
Script for interpolation of RAA
"""

import argparse
import sys
import pandas as pd
import ROOT
sys.path.append(__file__.rsplit("/", 2)[0]) # to import utils
from utils import read_tamu, read_lido, get_centrality_interpolation

def get_interpolated_raa(model):
    """
    Main function
    """

    gr_ds_raa, gr_d_raa, gr_npds_raa, gr_npd_raa, gr_b_raa, gr_bs_raa = {}, {}, {}, {}, {}, {}
    if model == "tamu":
        gr_d_raa["0_10"] = read_tamu("PromptD0_TAMU_RAA_5TeV_010.txt", "gr_d_raa_tamu_0_10")
        gr_d_raa["30_50"] = read_tamu("PromptD0_TAMU_RAA_5TeV_3050.txt", "gr_d_raa_tamu_30_50")
        gr_ds_raa["0_10"] = read_tamu("PromptDs_TAMU_RAA_5TeV_010.txt", "gr_ds_raa_tamu_0_10")
        gr_ds_raa["30_50"] = read_tamu("PromptDs_TAMU_RAA_5TeV_3050.txt", "gr_ds_raa_tamu_30_50")
        gr_npd_raa["0_10"] = read_tamu("NonPromptD0_TAMU_RAA_5TeV_010.txt", "gr_npd_raa_tamu_0_10", True)
        gr_npd_raa["30_50"] = read_tamu("NonPromptD0_TAMU_RAA_5TeV_3050.txt", "gr_npd_raa_tamu_30_50", True)
        gr_npds_raa["0_10"] = read_tamu("NonPromptDs_TAMU_RAA_5TeV_010.txt", "gr_npds_raa_tamu_0_10", True)
        gr_npds_raa["30_50"] = read_tamu("NonPromptDs_TAMU_RAA_5TeV_3050.txt", "gr_npds_raa_tamu_30_50", True)
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

        gr_npds_raa["80_90"] = gr_npds_raa["0_10"].Clone("gr_npds_raa_tamu_80_90")
        for ipt in range(gr_npds_raa["80_90"].GetN()):
            gr_npds_raa["80_90"].SetPoint(ipt, gr_npds_raa["80_90"].GetPointX(ipt), 0.8)
        gr_npd_raa["80_90"] = gr_npd_raa["0_10"].Clone("gr_npd_raa_tamu_80_90")
        for ipt in range(gr_npd_raa["80_90"].GetN()):
            gr_npd_raa["80_90"].SetPoint(ipt, gr_npd_raa["80_90"].GetPointX(ipt), 0.8)

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
    if model != "lido": # non prompt and strange mesons not available for lido
        get_centrality_interpolation(gr_ds_raa, f"gr_ds_raa_{model}")
        get_centrality_interpolation(gr_npds_raa, f"gr_npds_raa_{model}")
        get_centrality_interpolation(gr_npd_raa, f"gr_npd_raa_{model}")
        get_centrality_interpolation(gr_bs_raa, f"gr_bs_raa_{model}")

    outfile = ROOT.TFile(f"raa{model}_interpolated.root", "recreate")
    for graph in gr_d_raa.values():
        graph.Write()
    for graph in gr_b_raa.values():
        graph.Write()
    if model != "lido":
        for graph in gr_ds_raa.values():
            graph.Write()
        for graph in gr_npds_raa.values():
            graph.Write()
        for graph in gr_npd_raa.values():
            graph.Write()
        for graph in gr_bs_raa.values():
            graph.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate pT weights")
    parser.add_argument("--model", "-m", metavar="text", help="Enabled model", default="tamu")
    args = parser.parse_args()

    models = ["tamu", "lido"]
    if args.model not in models:
        print(f"Model {args.model} not yet implemented. Available options: {models}") 
        sys.exit()

    get_interpolated_raa(args.model)
