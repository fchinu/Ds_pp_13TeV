import ROOT
import numpy as np

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetAxisMaxDigits(3)
ROOT.gStyle.SetTitleSize(0.045, "XY")
ROOT.gStyle.SetLabelSize(0.045, "Y")

infile_ds = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/mult_differential/2023/ds/doublecb/CutVarDs_pp13TeV_10_30.root"

with ROOT.TFile.Open(infile_ds) as file_ds:
    h_data = file_ds.Get("hRawYieldVsCut_pt10_20")
    h_prompt = file_ds.Get("hRawYieldPromptVsCut_pt10_20")
    h_nonprompt = file_ds.Get("hRawYieldNonPromptVsCut_pt10_20")
    h_sum = file_ds.Get("hRawYieldSumVsCut_pt10_20")

cutsets = np.linspace(0.05, 0.75, 15)
cutsets = np.append(cutsets, 0.775)
h_data.GetXaxis().Set(len(cutsets), 0, len(cutsets))
h_prompt.GetXaxis().Set(len(cutsets), 0, len(cutsets))
h_nonprompt.GetXaxis().Set(len(cutsets), 0, len(cutsets))
h_sum.GetXaxis().Set(len(cutsets), 0, len(cutsets))


h_data.SetMarkerStyle(ROOT.kFullCircle)
h_data.SetMarkerSize(1.5)
h_data.SetMarkerColor(ROOT.kBlack)
h_data.SetLineColor(ROOT.kBlack)
h_data.SetLineWidth(2)
h_data.SetTitle(";Minimum BDT score for non-prompt D_{#lower[-0.3]{s}}^{+};Raw yield")
h_data.GetXaxis().SetLabelSize(0.06)
h_data.GetYaxis().SetRangeUser(0, 1.25*h_data.GetMaximum())
h_data.GetYaxis().SetTitleOffset(1.1)
h_data.GetXaxis().SetTitleOffset(1.1)

for i, cut in enumerate(cutsets):
    if i%2 == 0:
        h_data.GetXaxis().SetBinLabel(i+1, f"{cut:.2f}")
    else:
        h_data.GetXaxis().SetBinLabel(i+1, "")

c_ds = ROOT.TCanvas("c_ds", "c_ds", 600, 600)
c_ds.SetLeftMargin(0.12)
c_ds.SetRightMargin(0.05)
c_ds.SetTopMargin(0.05)
c_ds.SetBottomMargin(0.12)

h_data.Draw("same p")
h_prompt.Draw("same hist")
h_nonprompt.Draw("same hist")
h_sum.Draw("same hist")

latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextSize(0.045)
latex.SetTextAlign(12)
latex.DrawLatexNDC(0.15, 0.9, "ALICE Preliminary")
latex.DrawLatexNDC(0.15, 0.85, "pp,#kern[0.5]{#sqrt{#it{s}}} = 13.6 TeV")

latex.DrawLatexNDC(0.5, 0.75, "D_{s}^{+}#rightarrow#phi#pi^{+}#rightarrowK^{+}K^{#minus}#pi^{+}")
latex.DrawLatexNDC(0.5, 0.7, "and charge conjugate")
latex.DrawLatexNDC(0.5, 0.65, "1 < #it{p}_{T} < 2 GeV/#it{c}")
latex.DrawLatexNDC(0.5, 0.6, "10#minus30% FT0M")

legend = ROOT.TLegend(0.5, 0.35, 0.8, 0.55)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.AddEntry(h_data, "Data", "p")
legend.AddEntry(h_prompt, "Prompt", "f")
legend.AddEntry(h_nonprompt, "Non-prompt", "f")
legend.AddEntry(h_sum, "Total", "l")
legend.Draw()

c_ds.RedrawAxis()

c_ds.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/fprompt/VsMult/cut_var_1_2_10_30.pdf")