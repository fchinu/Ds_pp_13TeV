import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette

colors, _ = get_discrete_matplotlib_palette("tab10")

# Get inputs
infileDsData = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root")
hRawFracDs = infileDsData.Get("hRawFracPrompt")
hRawFracDs.SetDirectory(0)
infileDsData.Close()

infileDplusData = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root")
hRawFracDplus = infileDplusData.Get("hRawFracPrompt")
hRawFracDplus.SetDirectory(0)
infileDplusData.Close()

infileDsTheory = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/theory_driven/FD_fraction_Ds_LHC24d3a.root")
gRawFracDsTheory = infileDsTheory.Get("gfraction")
infileDsTheory.Close()

infileDplusTheory = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/theory_driven/FD_fraction_Dplus_LHC24d3a.root")
gRawFracDplusTheory = infileDplusTheory.Get("gfraction")
infileDplusTheory.Close()

# Set style
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

hRawFracDs.SetMarkerStyle(ROOT.kFullCircle)
hRawFracDs.SetMarkerColor(colors[0])
hRawFracDs.SetLineColor(colors[0])
hRawFracDs.SetMarkerSize(1.5)
hRawFracDs.SetLineWidth(2)

gRawFracDsTheory.SetMarkerStyle(ROOT.kOpenCircle)
gRawFracDsTheory.SetMarkerColor(colors[0])
gRawFracDsTheory.SetLineColor(colors[0])
gRawFracDsTheory.SetMarkerSize(1.5)
gRawFracDsTheory.SetLineWidth(2)

hRawFracDplus.SetMarkerStyle(ROOT.kFullDiamond)
hRawFracDplus.SetMarkerColor(colors[1])
hRawFracDplus.SetLineColor(colors[1])
hRawFracDplus.SetMarkerSize(2)
hRawFracDplus.SetLineWidth(2)

gRawFracDplusTheory.SetMarkerStyle(ROOT.kOpenDiamond)
gRawFracDplusTheory.SetMarkerColor(colors[1])
gRawFracDplusTheory.SetLineColor(colors[1])
gRawFracDplusTheory.SetMarkerSize(2)
gRawFracDplusTheory.SetLineWidth(2)

# Draw comparison
c = ROOT.TCanvas("c", "c", 800, 600)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.12)


hFrame = c.DrawFrame(0, 0.7, 24, 1.2, ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c}); Prompt fraction")
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(0.9)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetXaxis().SetLabelSize(0.04)
hFrame.GetYaxis().SetLabelSize(0.04)

ROOT.gStyle.SetLineStyleString(11,"50 15")
lineAtOne = ROOT.TLine(0, 1, 24, 1)
lineAtOne.SetLineStyle(11)
lineAtOne.SetLineColor(ROOT.kGray+1)
lineAtOne.SetLineWidth(2)
lineAtOne.Draw("same")

hRawFracDs.Draw("same")
hRawFracDplus.Draw("same")
gRawFracDsTheory.Draw("pez same")
gRawFracDplusTheory.Draw("pez same")

legend = ROOT.TLegend(0.6, 0.65, 0.8, 0.9)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.045)
legend.AddEntry(hRawFracDs, "D_{s}^{+} data#font[122]{-}driven", "lp")
legend.AddEntry(gRawFracDsTheory, "D_{s}^{+} theory#font[122]{-}driven", "lp")
legend.AddEntry(hRawFracDplus, "D^{+} data#font[122]{-}driven", "lp")
legend.AddEntry(gRawFracDplusTheory, "D^{+} theory#font[122]{-}driven", "lp")
legend.Draw("same")

thesisText = ROOT.TLatex(0.17, 0.85, "This Thesis")
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.07)
thesisText.Draw()

ppText = ROOT.TLatex(0.17, 0.8, "pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
ppText.SetNDC()
ppText.SetTextFont(42)
ppText.SetTextSize(0.05)
ppText.Draw()

DecayText = ROOT.TLatex(0.17, 0.7, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.17, 0.65, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ROOT.gPad.RedrawAxis()

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/fprompt/Compare_data_theory_frac.pdf")

