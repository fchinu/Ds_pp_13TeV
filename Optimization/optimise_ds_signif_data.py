"""
Script to perform D0 signal extraction as a function of BDT
"""

import os
import sys
import argparse
import numpy as np
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
import ROOT


def set_graph_style(graph, color=ROOT.kBlack):
    """
    Method for style setting
    """

    graph.SetMarkerStyle(ROOT.kFullCircle)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetLineWidth(2)

def perform_fit(file_name, histo_name, fitter_name, sgn_funcs, bkg_funcs,
                mass_min, mass_max):
    """
    Method for fit
    """
    labels_signal_pdf = [r"$\mathrm{D_{s}}^{+}$", r"$\mathrm{D}^{+}$"]
    labels_bkg_pdf = ["Comb. bkg."]

    data_hdl = DataHandler(file_name, histoname=histo_name, limits=[mass_min, mass_max])

    fitter = F2MassFitter(data_hdl,
                          name_signal_pdf=sgn_funcs,
                          name_background_pdf=bkg_funcs,
                          name=fitter_name, tol=0.1,
                          label_signal_pdf=labels_signal_pdf,
                          label_bkg_pdf=labels_bkg_pdf)
    fitter.set_particle_mass(0, pdg_id=431, limits=[1.95, 1.98])
    fitter.set_particle_mass(1, pdg_id=411, limits=[1.85, 1.88])
    fitter.set_signal_initpar(0, "sigma", 0.005, limits=[0.0001, 0.06])
    fitter.set_signal_initpar(1, "sigma", 0.005, limits=[0.0001, 0.06])
    fitter.set_signal_initpar(0, "frac", 0.1, limits=[0., 0.5])
    fitter.set_signal_initpar(1, "frac", 0.1, limits=[0., 0.5])
    fitter.set_background_initpar(0, "c0", 0.4)
    fitter.set_background_initpar(0, "c1", -0.2)
    fitter.set_background_initpar(0, "c2", -0.01)
    fitter.set_background_initpar(0, "c3", 0.01)
    fitter.set_background_initpar(0, "lam", -1.)

    fitter.mass_zfit()

    return fitter

# main function to perform fits
def optimise_signif(infile_name, force_project, project_only, bdt_min, bdt_max, bdt_step,
                    pt_min, pt_max, cent_min, cent_max, outputdir):
    """
    Method for fitting
    """

    ROOT.gStyle.SetTitleSize(0.045, "xy")
    ROOT.gStyle.SetLabelSize(0.04, "xy")
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadTopMargin(0.035)
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    infile_name_histos = os.path.join(
        outputdir, f"hist_mass_pt{pt_min:.0f}_{pt_max:.0f}_cent{cent_min:.0f}_{cent_max:.0f}.root")
    if force_project or not os.path.isfile(infile_name_histos):
        project(infile_name, bdt_min, bdt_max, bdt_step, pt_min, pt_max, cent_min, cent_max, outputdir)
    else:
        print("WARNING: using already available projection")

    if project_only:
        sys.exit()

    cuts_to_test = np.arange(bdt_min, bdt_max + bdt_step, bdt_step)
    ncuts = len(cuts_to_test)

    # we first fit the high significance cases
    outfile_name = os.path.join(
        outputdir, f"significance_scan_d0_{pt_min:.0f}_{pt_max:.0f}_cent{cent_min:.0f}_{cent_max:.0f}.root")
    outfile = ROOT.TFile(outfile_name, "recreate")
    outfile.Close()

    graph_rawyield, graph_sigma, graph_mean, graph_significance, graph_soverb = (
        [] for _ in range(5))
    part_names = ["ds", "dplus"]
    part_colors = [ROOT.kRed+1, ROOT.kAzure+4]
    for isig in range(2):
        graph_rawyield.append(ROOT.TGraphErrors(ncuts))
        graph_sigma.append(ROOT.TGraphErrors(ncuts))
        graph_mean.append(ROOT.TGraphErrors(ncuts))
        graph_significance.append(ROOT.TGraphErrors(ncuts))
        graph_soverb.append(ROOT.TGraphErrors(ncuts))
        graph_rawyield[isig].SetNameTitle(f"graph_rawyield_{part_names[isig]}",
                                          ";BDT bkg < ; raw yield")
        graph_sigma[isig].SetNameTitle(f"graph_sigma_{part_names[isig]}",
                                       ";BDT bkg < ; #sigma (GeV/#it{c}^{2})")
        graph_mean[isig].SetNameTitle(f"graph_mean_{part_names[isig]}",
                                      ";BDT bkg < ; #mu (GeV/#it{c}^{2})")
        graph_significance[isig].SetNameTitle(f"graph_significance_{part_names[isig]}",
                                              ";BDT bkg < ; significance(3#sigma)")
        graph_soverb[isig].SetNameTitle(f"graph_soverb_{part_names[isig]}",
                                        ";BDT bkg < ; S/B(3#sigma)")
        set_graph_style(graph_rawyield[isig], part_colors[isig])
        set_graph_style(graph_sigma[isig], part_colors[isig])
        set_graph_style(graph_mean[isig], part_colors[isig])
        set_graph_style(graph_significance[isig], part_colors[isig])
        set_graph_style(graph_soverb[isig], part_colors[isig])

    max_rawy, max_signif, max_soverb = 0., 0., 0.
    for icut, cut in enumerate(cuts_to_test):
        fitter_data = perform_fit(infile_name_histos,
                                  f"hist_mass_bdt{cut:.3f}",
                                  f"ds_bdt{cut:.3f}",
                                  ["gaussian", "gaussian"],
                                  ["chebpol2"],
                                  1.75,
                                  2.15)

        if fitter_data.get_fit_result.converged:
            for isig in range(2):
                rawyield = fitter_data.get_raw_yield(isig)
                max_rawy = max(max_rawy, rawyield[0])
                sigma = fitter_data.get_signal_parameter(isig, "sigma")
                mean = fitter_data.get_mass(isig)
                significance = fitter_data.get_significance(isig, nsigma=3)
                max_signif = max(max_signif, significance[0])
                soverb = fitter_data.get_signal_over_background(isig, nsigma=3)
                max_soverb = max(max_soverb, soverb[0])
                graph_rawyield[isig].SetPoint(icut, cut, rawyield[0])
                graph_rawyield[isig].SetPointError(icut, 0, rawyield[1])
                graph_sigma[isig].SetPoint(icut, cut, sigma[0])
                graph_sigma[isig].SetPointError(icut, 0, sigma[1])
                graph_mean[isig].SetPoint(icut, cut, mean[0])
                graph_mean[isig].SetPointError(icut, 0, mean[1])
                graph_soverb[isig].SetPoint(icut, cut, soverb[0])
                graph_soverb[isig].SetPointError(icut, 0, soverb[1])
                graph_significance[isig].SetPoint(icut, cut, significance[0])
                graph_significance[isig].SetPointError(icut, 0, significance[1])

            fitter_data.dump_to_root(outfile_name,
                                     option="update",
                                     suffix=f"_bdt{cut:.3f}")

            xaxis = r"$M(\mathrm{KK\pi})$ (GeV/$c^{2}$)"
            fig, _ = fitter_data.plot_mass_fit(style="ATLAS", figsize=(8, 8), axis_title=xaxis)
            figres = fitter_data.plot_raw_residuals(figsize=(8, 8), style="ATLAS", axis_title=xaxis)

            fig.savefig(
                os.path.join(outputdir, f"massfit_pt{pt_min:.0f}_{pt_max:.0f}_bdt{cut:.3f}.pdf")
            )
            figres.savefig(
                os.path.join(outputdir, f"massfitres_pt{pt_min:.0f}_{pt_max:.0f}_bdt{cut:.3f}.pdf")
            )

    canv = ROOT.TCanvas("canv_signif", "", 1920, 1080)
    canv.Divide(3, 2)
    canv.cd(1).DrawFrame(
        bdt_min, 0., bdt_max, max_rawy*1.2, ";BDT bkg < ; raw yield")
    for isig in range(2):
        graph_rawyield[isig].Draw("p")
    canv.cd(2).DrawFrame(
        bdt_min, 0., bdt_max, 0.05, ";BDT bkg < ; #sigma (GeV/#it{c}^{2})")
    for isig in range(2):
        graph_sigma[isig].Draw("p")
    canv.cd(3).DrawFrame(
        bdt_min, 1.8, bdt_max, 2.0, ";BDT bkg < ; #mu (GeV/#it{c}^{2})")
    for isig in range(2):
        graph_mean[isig].Draw("p")
    canv.cd(4).DrawFrame(
        bdt_min, 0., bdt_max, max_signif*1.2, ";BDT bkg < ; significance(3#sigma)")
    for isig in range(2):
        graph_significance[isig].Draw("p")
    canv.cd(5).DrawFrame(
        bdt_min, 0., bdt_max, max_soverb*1.2, ";BDT bkg < ; S/B(3#sigma)")
    for isig in range(2):
        graph_soverb[isig].Draw("p")
    canv.cd(6)
    leg = ROOT.TLegend(0.2, 0.5, 0.8, 0.8)
    leg.SetTextSize(0.05)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetHeader(f"#splitline{{{pt_min:.0f} < #it{{p}}_{{T}} < {pt_max:.0f} GeV/#it{{c}}}}{{Centrality {cent_min:.0f}#minus{cent_max:.0f}%}}")
    leg.AddEntry(graph_significance[0], "D_{s}^{+}", "lp")
    leg.AddEntry(graph_significance[1], "D^{+}", "lp")
    leg.Draw()

    canv.SaveAs(outfile_name.replace(".root", ".pdf"))

    outfile_nocut = ROOT.TFile(outfile_name, "update")
    canv.Write()
    for isig in range(2):
        graph_rawyield[isig].Write()
        graph_sigma[isig].Write()
        graph_mean[isig].Write()
        graph_soverb[isig].Write()
        graph_significance[isig].Write()
    outfile_nocut.Close()


# function to project the sparse
def project(infile_name, bdt_min, bdt_max, bdt_step, pt_min, pt_max, cent_min, cent_max, outputdir):
    """
    Method for sparse projection
    """

    infile = ROOT.TFile.Open(infile_name)
    sparse = infile.Get("hf-task-ds/Data/hSparseMass")
    infile.Close()

    outfile = ROOT.TFile(
        os.path.join(outputdir, f"hist_mass_pt{pt_min:.0f}_{pt_max:.0f}_cent{cent_min:.0f}_{cent_max:.0f}.root"), "recreate")
    pt_bin_min = sparse.GetAxis(1).FindBin(pt_min*1.001)
    pt_bin_max = sparse.GetAxis(1).FindBin(pt_max*0.999)
    cent_bin_min = sparse.GetAxis(2).FindBin(cent_min*1.001)
    cent_bin_max = sparse.GetAxis(2).FindBin(cent_max*0.999)
    sparse.GetAxis(1).SetRange(pt_bin_min, pt_bin_max)
    sparse.GetAxis(2).SetRange(cent_bin_min, cent_bin_max)
    histos_pt = []
    for icut, cut in enumerate(np.arange(bdt_min, bdt_max+bdt_step, bdt_step)):
        bdt_bkg_bin_max = sparse.GetAxis(3).FindBin(cut*0.999)
        sparse.GetAxis(3).SetRange(1, bdt_bkg_bin_max)
        histos_pt.append(sparse.Projection(0))
        histos_pt[icut].SetName(f"hist_mass_bdt{cut:.3f}")
        outfile.cd()
        histos_pt[icut].Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infile", "-i", metavar="text",
                        default="AnalysisResults.root", help="input file")
    parser.add_argument("--force_project", "-f", action="store_true",
                        default=False, help="force projection of THnSparse")
    parser.add_argument("--project_only", "-p", action="store_true",
                        default=False, help="do projection only, no fit")
    parser.add_argument("--bdt_min", "-mi", type=float,
                        default=0.002, help="minimum BDT cut to test")
    parser.add_argument("--bdt_max", "-ma", type=float,
                        default=0.02, help="maximum BDT cut to test")
    parser.add_argument("--bdt_step", "-s", type=float,
                        default=0.002, help="step between one BDT cut and another")
    parser.add_argument("--pt_min", "-ptmi", type=float,
                        default=0., help="minimum pT bin")
    parser.add_argument("--pt_max", "-ptma", type=float,
                        default=1., help="maximum pT bin")
    parser.add_argument("--cent_min", "-cmi", type=float,
                        default=0., help="minimum centrality bin")
    parser.add_argument("--cent_max", "-cma", type=float,
                        default=1., help="maximum centrality bin")
    parser.add_argument("--output_dir", "-o", metavar="text",
                        default=".", help="output directory")
    args = parser.parse_args()

    optimise_signif(args.infile, args.force_project, args.project_only,
                    args.bdt_min, args.bdt_max, args.bdt_step,
                    args.pt_min, args.pt_max, args.cent_min, args.cent_max,
                    args.output_dir)
