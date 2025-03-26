import ROOT
import numpy as np
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from plot_utils import get_discrete_matplotlib_palette


infile = "/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/VsMult/efficiency_LHC24h1.root"

cent_mins = [0, 1, 10, 30, 50, 70]
cent_maxs = [1, 10, 30, 50, 70, 100]
pt_mins = [1, 2, 4, 6, 8, 12]
pt_maxs = [2, 4, 6, 8, 12, 24]

colors, _ = get_discrete_matplotlib_palette("tab10")
markers_full = [
    ROOT.kFullCircle,
    ROOT.kFullDiamond,
    ROOT.kFullSquare,
    ROOT.kFullDoubleDiamond,
    ROOT.kFullCross,
    ROOT.kFullFourTrianglesPlus
]
markers_open = [
    ROOT.kOpenCircle,
    ROOT.kOpenDiamond,
    ROOT.kOpenSquare,
    ROOT.kOpenDoubleDiamond,
    ROOT.kOpenCross,
    ROOT.kOpenFourTrianglesPlus
]
markers_size = [
    1.5,
    2,
    1.5,
    2,
    1.5,
    2
]

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetLabelSize(0.05, "XYZ")
ROOT.gStyle.SetTitleSize(0.05, "XYZ")
ROOT.gStyle.SetTitleOffset(1.1, "X")
ROOT.gStyle.SetTitleOffset(1.4, "Y")
ROOT.gStyle.SetLabelOffset(0.01, "Y")
ROOT.gStyle.SetLabelOffset(-0.01, "X")

h_effs_ds_vs_pt = []
h_effs_dplus = []
with ROOT.TFile(infile) as f:
    for cent_min, cent_max in zip(cent_mins, cent_maxs):
        h_effs_ds_vs_pt.append(f.Get(f"eff_DsPrompt_cent_{cent_min}_{cent_max}"))
        h_effs_dplus.append(f.Get(f"eff_DplusPrompt_cent_{cent_min}_{cent_max}"))

h_effs_ds_vs_cent = []
h_effs_dplus_vs_cent = []

for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
    h_effs_ds_vs_cent.append(ROOT.TH1D(f"eff_DsPrompt_pt_{pt_min}_{pt_max}", f";FT0M percentile;D_{{s}}^{{+}} efficiency", len(cent_mins), np.asarray([1.e-2 if x==0 else x for x in cent_mins] + [cent_maxs[-1]], "d")))
    h_effs_dplus_vs_cent.append(ROOT.TH1D(f"eff_DplusPrompt_pt_{pt_min}_{pt_max}", f";FT0M percentile;D_{{s}}^{{+}} efficiency", len(cent_mins), np.asarray([1.e-2 if x==0 else x for x in cent_mins] + [cent_maxs[-1]], "d")))

    for i_cent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
        h_effs_ds_vs_cent[i_pt].SetBinContent(i_cent + 1, h_effs_ds_vs_pt[i_cent].GetBinContent(i_pt + 1))
        h_effs_ds_vs_cent[i_pt].SetBinError(i_cent + 1, h_effs_ds_vs_pt[i_cent].GetBinError(i_pt + 1))
        h_effs_dplus_vs_cent[i_pt].SetBinContent(i_cent + 1, h_effs_dplus[i_cent].GetBinContent(i_pt + 1))
        h_effs_dplus_vs_cent[i_pt].SetBinError(i_cent + 1, h_effs_dplus[i_cent].GetBinError(i_pt + 1))

for i_pt, (h_ds, h_dplus) in enumerate(zip(h_effs_ds_vs_cent, h_effs_dplus_vs_cent)):
    h_ds.SetMarkerColor(colors[i_pt])
    h_ds.SetLineColor(colors[i_pt])
    h_ds.SetMarkerStyle(markers_full[i_pt])
    h_ds.SetMarkerSize(markers_size[i_pt])
    h_ds.SetLineWidth(2)

    h_dplus.SetMarkerColor(colors[i_pt])
    h_dplus.SetLineColor(colors[i_pt])
    h_dplus.SetMarkerStyle(markers_full[i_pt])
    h_dplus.SetMarkerSize(markers_size[i_pt])
    h_dplus.SetLineWidth(2)
c = ROOT.TCanvas("c", "c", 1000, 600)

# Create sub pads
pad_ds = ROOT.TPad("pad_ds", "pad_ds", 0.0, 0.0, 0.527777, 1.0)
pad_dplus = ROOT.TPad("pad_dplus", "pad_dplus", 0.527777, 0.0, 1.0, 1.0)

# Set margins for sub pads
pad_ds.SetLeftMargin(0.15)
pad_ds.SetRightMargin(0.0)
pad_ds.SetTopMargin(0.05)
pad_ds.SetBottomMargin(0.15)

pad_dplus.SetLeftMargin(0.0)
pad_dplus.SetRightMargin(0.05)
pad_dplus.SetTopMargin(0.05)
pad_dplus.SetBottomMargin(0.15)

# Draw sub pads
pad_ds.Draw()
pad_dplus.Draw()

# Draw on pad_ds
pad_ds.cd()
ROOT.gPad.SetLogy()
ROOT.gPad.SetLogx()
h_frame_ds = ROOT.gPad.DrawFrame(2.e-1, 1.e-3, 99.9, 1., ";FT0M percentile;Acceptance#kern[0.5]{#times} Efficiency")
leg_ds = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
leg_ds.SetBorderSize(0)
leg_ds.SetFillStyle(0)
for i_pt, h_ds in enumerate(h_effs_ds_vs_cent):
    h_ds.Draw("same p")
    leg_ds.AddEntry(h_ds, f"{pt_mins[i_pt]} < #it{{p}}_{{T}} < {pt_maxs[i_pt]} GeV/#it{{c}}", "pl")
#leg_ds.Draw()

latex_alice = ROOT.TLatex()
latex_alice.SetTextFont(42)
latex_alice.SetTextSize(0.06)
latex_alice.SetNDC()
latex_alice.DrawLatex(0.18, 0.88, "ALICE Preliminary")

latex_ds = ROOT.TLatex()
latex_ds.SetTextFont(42)
latex_ds.SetTextSize(0.045)
latex_ds.DrawLatex(0.3, 2.1e-3, "D_{s}^{+}#rightarrow#phi#pi^{+}#rightarrowK^{+}K^{#minus}#pi^{+} and charge conj.")
latex_ds.DrawLatex(0.3, 1.4e-3, "Prompt D_{s}^{+}")

pad_ds.RedrawAxis()

# Draw on pad_dplus
pad_dplus.cd()
ROOT.gPad.SetLogy()
ROOT.gPad.SetLogx()
h_frame_dplus = ROOT.gPad.DrawFrame(2.e-1, 1.e-3, 100, 1., ";FT0M percentile;D_{s}^{+} efficiency")
h_frame_dplus.GetXaxis().SetTitleSize(0.05*1.117)
h_frame_dplus.GetXaxis().SetTitleOffset(1.1/1.117)
h_frame_dplus.GetXaxis().SetLabelSize(0.05*1.117)
h_frame_dplus.GetXaxis().SetLabelOffset(-0.015*1.117)
leg_dplus = ROOT.TLegend(0.03, 0.3, 0.91, 0.45)
leg_dplus.SetBorderSize(0)
leg_dplus.SetFillStyle(0)
leg_dplus.SetTextSize(0.045)
leg_dplus.SetNColumns(2)
leg_dplus.SetColumnSeparation(0.035)
leg_dplus.SetMargin(0.15)
for i_pt, h_dplus in enumerate(h_effs_dplus_vs_cent):
    h_dplus.Draw("same p")
    leg_dplus.AddEntry(h_dplus, f"{pt_mins[i_pt]}#kern[0.6]{{<}}#kern[0.4]{{#it{{p}}_{{T}}}}#kern[0.25]{{<}} {pt_maxs[i_pt]} GeV/#it{{c}}", "pl")
leg_dplus.Draw()

latex_dplus = ROOT.TLatex()
latex_dplus.SetTextFont(42)
latex_dplus.SetTextSize(0.05)
latex_dplus.DrawLatex(0.27, 2.1e-3, "D^{+}#rightarrow#phi#pi^{+}#rightarrowK^{+}K^{#minus}#pi^{+} and charge conj.")
latex_dplus.DrawLatex(0.27, 1.4e-3, "Prompt D^{+}")

latex_syst = ROOT.TLatex()
latex_syst.SetTextFont(42)
latex_syst.SetTextSize(0.05)
latex_syst.SetNDC()
latex_syst.DrawLatex(0.05, 0.88, "pp,#kern[0.25]{#sqrt{#it{s}}} = 13.6 TeV, |#it{y}| < 0.5")

pad_ds.RedrawAxis()
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/efficiency/efficiency.pdf")



