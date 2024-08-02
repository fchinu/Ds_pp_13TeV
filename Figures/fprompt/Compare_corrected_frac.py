import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette

colors, _ = get_discrete_matplotlib_palette("tab10")

# Get inputs
infileDs = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root")
hCorrFracDs = infileDs.Get("hCorrFracPrompt")
hCorrFracDs.SetDirectory(0)
infileDs.Close()

infileDplus = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root")
hCorrFracDplus = infileDplus.Get("hCorrFracPrompt")
hCorrFracDplus.SetDirectory(0)
infileDplus.Close()

# Set style
hCorrFracDs.SetMarkerStyle(ROOT.kFullCircle)
hCorrFracDs.SetMarkerColor(colors[0])
hCorrFracDs.SetLineColor(colors[0])
hCorrFracDs.SetMarkerSize(1.5)
hCorrFracDs.SetLineWidth(2)

hCorrFracDplus.SetMarkerStyle(ROOT.kFullDiamond)
hCorrFracDplus.SetMarkerColor(colors[1])
hCorrFracDplus.SetLineColor(colors[1])
hCorrFracDplus.SetMarkerSize(2)
hCorrFracDplus.SetLineWidth(2)

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

# Draw comparison
c = ROOT.TCanvas("c", "c", 800, 600)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.12)


hFrame = c.DrawFrame(0, 0.6, 24, 1.2666666, ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c}); Corrected prompt fraction")
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

hCorrFracDs.Draw("same")
hCorrFracDplus.Draw("same")

legend = ROOT.TLegend(0.7, 0.75, 0.9, 0.85)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.04)
legend.AddEntry(hCorrFracDs, "D_{s}^{+}", "lp")
legend.AddEntry(hCorrFracDplus, "D^{+}", "lp")
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

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/fprompt/Compare_prompt_frac.pdf")
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/fprompt/Compare_prompt_frac.png")

