"""

"""
import argparse
import ROOT
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from plot_utils import get_discrete_matplotlib_palette
ci = ROOT.TColor.GetFreeColorIndex()

def set_data_style(hist, isnorm=False):
    """

    """

    hist.SetDirectory(0)
    hist.Sumw2()
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetMarkerSize(0.8)
    hist.SetLineWidth(2)
    hist.SetLineColor(ROOT.kBlack)
    if isnorm:
        hist.GetYaxis().SetTitle(
            f"Normalised counts per {hist.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}")
        hist.GetYaxis().SetRangeUser(0., hist.GetMaximum() * 1.3)
        hist.SetMarkerSize(1.6)
    else:
        hist.GetYaxis().SetTitle(f"Counts per {hist.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}")
        hist.GetYaxis().SetRangeUser(0.1, hist.GetMaximum() * 1.3)
        # set x range
        hist.GetXaxis().SetRangeUser(hist.GetBinLowEdge(1), hist.GetBinLowEdge(hist.GetNbinsX()) + hist.GetBinWidth(hist.GetNbinsX())-0.001)
    hist.GetXaxis().SetTitle("#it{M}(D^{#mp}#pi^{#pm}) (GeV/#it{c}^{2})")

    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetNdivisions(508)
    hist.GetYaxis().SetNdivisions(508)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetDecimals()
    hist.GetYaxis().SetMaxDigits(3)


def get_function_fromhist(hist, func_type="totfunc", norm_fact = 1., colors=None):
    """

    """
    hist.Scale(norm_fact)
    hist.SetFillStyle(1001)
    hist.SetLineWidth(0)

    spline = ROOT.TSpline3(hist)
    spline.SetLineWidth(3)
    spline.SetNpx(500)
    if func_type == "totfunc":
        spline.SetLineColor(ROOT.kBlue + 2)
        hist.SetFillColorAlpha(ROOT.kBlue + 2, 0.5)
    elif func_type == "bkg":
        spline.SetLineColor(ROOT.kRed + 1)
        hist.SetLineColor(0)
        hist.SetFillColor(10)
        hist.SetLineStyle(9)
        spline.SetLineStyle(9)
    elif "bkg_corr" in func_type:
        i_bkg = int(func_type.split("_")[-1])
        spline.SetLineColor(colors[2])
        hist.SetFillColorAlpha(colors[1], 0.75)
        spline.SetFillColorAlpha(colors[1], 0.75)
        spline.SetFillStyle(1000)
        spline.SetLineStyle(i_bkg+2)
    elif func_type == "signal_0":
        spline.SetLineColor(ROOT.kAzure + 4)
        hist.SetFillColorAlpha(ROOT.kAzure + 4, 0.5)
    elif func_type == "signal_1":
        spline.SetLineColor(ROOT.kTeal - 7)
        hist.SetFillColorAlpha(ROOT.kTeal - 7, 0.5)

    return spline

def plot(infile_name, colors, pt_min, pt_max, cent_min, cent_max):
    """

    """
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetPadTopMargin(0.065)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.TGaxis.SetMaxDigits(2)

    suffix = f"{pt_min*10}_{pt_max*10}_cent_{cent_min}_{cent_max}_"

    infile = ROOT.TFile.Open(infile_name)
    hist_mass = infile.Get(f"hdata_{suffix}")
    set_data_style(hist_mass)

    norm = 1./hist_mass.Integral()
    hist_mass_norm = hist_mass.Clone()
    hist_mass_norm.Scale(norm)
    set_data_style(hist_mass_norm, True)
    hist_signal_ds = infile.Get(f"signal_0_{suffix}")
    hist_signal_ds_norm = hist_signal_ds.Clone("hist_signal_ds")
    hist_signal_dplus = infile.Get(f"signal_1_{suffix}")
    hist_signal_dplus_norm = hist_signal_dplus.Clone("hist_signal_dplus")
    hist_bkg = infile.Get(f"bkg_1_{suffix}")
    hist_bkg_corr_0 = infile.Get(f"bkg_0_{suffix}")
    hist_tot_func = infile.Get(f"total_func_{suffix}")
    func_bkg = get_function_fromhist(hist_bkg, "bkg")
    func_bkg_corr_0 = get_function_fromhist(hist_bkg_corr_0, "bkg_corr_0", 1., colors)
    corr_bkg_stack = ROOT.THStack("corr_bkg_stack", "")
    corr_bkg_stack.Add(hist_bkg)
    corr_bkg_stack.Add(hist_bkg_corr_0)

    func_signal_ds = get_function_fromhist(hist_signal_ds, "signal_0")
    func_signal_dplus = get_function_fromhist(hist_signal_dplus, "signal_1")
    func_totfunc = get_function_fromhist(infile.Get(f"total_func_{suffix}"), "totfunc")

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextColor(ROOT.kBlack)
    lat.SetTextSize(0.045)

    lat_small = ROOT.TLatex()
    lat_small.SetNDC()
    lat_small.SetTextFont(42)
    lat_small.SetTextColor(ROOT.kBlack)
    lat_small.SetTextSize(0.04)

    leg = ROOT.TLegend(0.55, 0.6, 0.85, 0.93)
    leg.SetTextSize(0.035)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(hist_mass, "Data", "p")
    leg.AddEntry(hist_signal_ds, "#splitline{#lower[0.2]{D_{s}^{+}#rightarrow#phi#pi^{#plus}#rightarrow K^{#plus}K^{#minus}#pi^{#plus}}}{#lower[-0.2]{and charge conj.}}", "f")
    leg.AddEntry(hist_signal_dplus, "#splitline{#lower[0.2]{D^{+}#rightarrow#phi#pi^{#plus}#rightarrow K^{#plus}K^{#minus}#pi^{#plus}}}{#lower[-0.2]{and charge conj.}}", "f")
    leg.AddEntry(func_bkg, "Comb. background", "l")
    leg.AddEntry(hist_bkg_corr_0, "#splitline{#lower[0.2]{D^{+}#rightarrow#pi^{+}K^{#minus}#pi^{+}}}{#lower[-0.2]{correlated background}}", "f")
    leg.AddEntry(func_totfunc, "Total fit function", "l")

    leg2 = ROOT.TLegend(0.55, 0.6, 0.85, 0.93)
    leg2.SetTextSize(0.035)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.AddEntry(hist_mass, "Data", "p")
    leg2.AddEntry(func_signal_ds, "#splitline{#lower[0.2]{D_{s}^{+}#rightarrow#phi#pi^{#plus}#rightarrow K^{#plus}K^{#minus}#pi^{#plus}}}{#lower[-0.2]{and charge conj.}}", "f")
    leg2.AddEntry(func_signal_dplus, "#splitline{#lower[0.2]{D^{+}#rightarrow#phi#pi^{#plus}#rightarrow K^{#plus}K^{#minus}#pi^{#plus}}}{#lower[-0.2]{and charge conj.}}", "f")
    leg2.AddEntry(func_bkg, "Comb. background", "l")
    leg2.AddEntry(hist_bkg_corr_0, "#splitline{#lower[0.2]{D^{+}#rightarrow#pi^{+}K^{#minus}#pi^{+}}}{#lower[-0.2]{correlated background}}", "f")
    leg2.AddEntry(func_totfunc, "Total fit function", "l")

    canv_masses = ROOT.TCanvas("canv_masses", "", 500, 500)
    h_frame = canv_masses.DrawFrame(
        hist_mass.GetBinLowEdge(1),
        0.,
        hist_mass.GetBinLowEdge(hist_mass.GetNbinsX()) + hist_mass.GetBinWidth(hist_mass.GetNbinsX())-0.0002,
        hist_mass.GetMaximum() * 1.3,
        ";#it{M}(K^{#pm}K^{#mp}#pi^{#pm}) (GeV/#it{c}^{2})"f";Counts per {hist_mass.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}"
    )
    h_frame.GetXaxis().SetTitleSize(0.045)
    h_frame.GetXaxis().SetNdivisions(508)
    h_frame.GetYaxis().SetNdivisions(508)
    h_frame.GetXaxis().SetLabelSize(0.045)
    h_frame.GetYaxis().SetTitleSize(0.045)
    h_frame.GetYaxis().SetTitleOffset(1.3)
    h_frame.GetYaxis().SetLabelSize(0.045)
    h_frame.GetYaxis().SetDecimals()
    h_frame.GetYaxis().SetMaxDigits(3)
    h_frame.GetYaxis().ChangeLabel(1, 1, 0)
    hist_mass.DrawCopy("same")
    corr_bkg_stack.Draw("NOCLEAR hist same")
    hist_bkg.Draw("histsame")
    func_bkg.Draw("same")
    hist_mass.DrawCopy("][ same")
    func_totfunc.Draw("lsame")
    hist_signal_ds.DrawCopy("][ histsame")
    hist_signal_dplus.DrawCopy("][ histsame")
    func_signal_ds.Draw("lsame")
    func_signal_dplus.Draw("lsame")
    lat.DrawLatex(0.15, 0.86, "ALICE Preliminary")
    lat.DrawLatex(0.15, 0.8, "pp,#sqrt{#it{s}} = 13.6 TeV")
    lat.DrawLatex(0.15, 0.74, "#font[132]{#it{L}}_{int} = 6.5 pb^{#minus1}")
    lat.DrawLatex(0.15, 0.68, f"{pt_min}#kern[0.25]{{<}}#kern[0.25]{{#it{{p}}_{{T}}}}#kern[0.25]{{<}}#kern[0.25]{{{pt_max}}} GeV/#it{{c}}")
    lat.DrawLatex(0.15, 0.62, f"{cent_min}#kern[0.25]{{#font[122]{{-}}}}{cent_max}% FT0M")

    leg.Draw()
    leg2.Draw()

    canv_masses.RedrawAxis()

    canv_masses.SaveAs(f"/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/VsMult/D_mass_{suffix}.pdf")
    canv_masses.SaveAs(f"/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/VsMult/D_mass_{suffix}.root")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_file", "-i", metavar="text",
                        default="Projections_RawYields/data/VsMult/FT0M/2023_h1_fix_template_params_sigma_for_figure/fits/fits_.root",
                        help="input root file", required=False)
    parser.add_argument("--palette", "-p", metavar="text",
                        default="tab10",
                        help="matplotlib palette name for corr. bkg", required=False)
    parser.add_argument("--ptmin", type=int,
                        default=None,
                        help="minimum pt", required=True)
    parser.add_argument("--ptmax", type=int,
                        default=None,
                        help="maximum pt", required=True)
    parser.add_argument("--centmin", type=int,
                        default=None,
                        help="minimum centrality", required=True)
    parser.add_argument("--centmax", type=int,
                        default=None,
                        help="maximum centrality", required=True)
    args = parser.parse_args()

    ROOT.TH1.AddDirectory(False)
    COLORS, _ = get_discrete_matplotlib_palette(args.palette)

    plot(args.input_file, COLORS, args.ptmin, args.ptmax, args.centmin, args.centmax)
