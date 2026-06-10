"""
Script for interpolation of RAA
"""

import argparse
import sys
import pandas as pd
import ROOT
sys.path.append(__file__.rsplit("/", 2)[0]) # to import utils
from utils import read_tamu, read_lido, read_langevin, get_centrality_interpolation

def get_interpolated_raa(model):
    """
    Main function
    """

    gr_ds_raa, gr_d_raa, gr_npds_raa, gr_npd_raa, gr_b_raa, gr_bs_raa = {}, {}, {}, {}, {}, {}
    if model == "langevin":
        gr_d_raa["0_10"] = read_langevin("langevin/OO_RAA_D0_010.dat", "gr_d_raa_langevin_0_10")
        gr_d_raa["30_50"] = read_langevin("langevin/OO_RAA_D0_3050.dat", "gr_d_raa_langevin_30_50")
        gr_d_raa["50_80"] = read_langevin("langevin/OO_RAA_D0_5080.dat", "gr_d_raa_langevin_50_80")
        gr_b_raa["0_10"] = read_langevin("langevin/OO_RAA_D0_010.dat", "gr_b_raa_langevin_0_10")
        gr_b_raa["30_50"] = read_langevin("langevin/OO_RAA_D0_3050.dat", "gr_b_raa_langevin_30_50")
        gr_b_raa["50_80"] = read_langevin("langevin/OO_RAA_D0_5080.dat", "gr_b_raa_langevin_50_80")

        # # assume most peripheral to be = 0.8
        # gr_ds_raa["80_90"] = gr_ds_raa["0_10"].Clone("gr_ds_raa_langevin_80_90")
        # for ipt in range(gr_ds_raa["80_90"].GetN()):
        #     gr_ds_raa["80_90"].SetPoint(ipt, gr_ds_raa["80_90"].GetPointX(ipt), 0.8)
        # gr_d_raa["80_90"] = gr_d_raa["0_10"].Clone("gr_d_raa_langevin_80_90")
        # for ipt in range(gr_d_raa["80_90"].GetN()):
        #     gr_d_raa["80_90"].SetPoint(ipt, gr_d_raa["80_90"].GetPointX(ipt), 0.8)

        # gr_npds_raa["80_90"] = gr_npds_raa["0_10"].Clone("gr_npds_raa_langevin_80_90")
        # for ipt in range(gr_npds_raa["80_90"].GetN()):
        #     gr_npds_raa["80_90"].SetPoint(ipt, gr_npds_raa["80_90"].GetPointX(ipt), 0.8)
        # gr_npd_raa["80_90"] = gr_npd_raa["0_10"].Clone("gr_npd_raa_langevin_80_90")
        # for ipt in range(gr_npd_raa["80_90"].GetN()):
        #     gr_npd_raa["80_90"].SetPoint(ipt, gr_npd_raa["80_90"].GetPointX(ipt), 0.8)

        # gr_b_raa["80_90"] = gr_b_raa["0_10"].Clone("gr_b_raa_langevin_80_90")
        # for ipt in range(gr_b_raa["80_90"].GetN()):
        #     gr_b_raa["80_90"].SetPoint(ipt, gr_b_raa["80_90"].GetPointX(ipt), 0.8)
        # gr_bs_raa["80_90"] = gr_bs_raa["0_10"].Clone("gr_bs_raa_langevin_80_90")
        # for ipt in range(gr_bs_raa["80_90"].GetN()):
        #     gr_bs_raa["80_90"].SetPoint(ipt, gr_bs_raa["80_90"].GetPointX(ipt), 0.8)


    # interpolate to obtain all centralities
    get_centrality_interpolation(gr_d_raa, f"gr_d_raa_{model}")
    get_centrality_interpolation(gr_b_raa, f"gr_b_raa_{model}")
    # if model != "lido": # non prompt and strange mesons not available for lido
    #     get_centrality_interpolation(gr_ds_raa, f"gr_ds_raa_{model}")
    #     get_centrality_interpolation(gr_npds_raa, f"gr_npds_raa_{model}")
    #     get_centrality_interpolation(gr_npd_raa, f"gr_npd_raa_{model}")
    #     get_centrality_interpolation(gr_bs_raa, f"gr_bs_raa_{model}")

    outfile = ROOT.TFile(f"raa{model}_interpolated.root", "recreate")
    for graph in gr_d_raa.values():
        graph.Write()
    for graph in gr_b_raa.values():
        graph.Write()
    # if model != "lido":
    #     for graph in gr_ds_raa.values():
    #         graph.Write()
    #     for graph in gr_npds_raa.values():
    #         graph.Write()
    #     for graph in gr_npd_raa.values():
    #         graph.Write()
    #     for graph in gr_bs_raa.values():
    #         graph.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate pT weights")
    parser.add_argument("--model", "-m", metavar="text", help="Enabled model", default="langevin")
    args = parser.parse_args()

    models = ["langevin"]
    if args.model not in models:
        print(f"Model {args.model} not yet implemented. Available options: {models}") 
        sys.exit()

    get_interpolated_raa(args.model)
