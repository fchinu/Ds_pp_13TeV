import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import set_matplotlib_palette, get_discrete_matplotlib_palette

set_matplotlib_palette('coolwarm')
#ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)

inFile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Optimization/Significance_ptDiff_2_3.root")
hSignif = inFile.Get("hSignificanceScan_2_2.5")
hSignif.SetDirectory(0)
hSignif.SetTitle("Significance scan;Minimum background score;Minimum prompt score;Significance")
hSignif.GetZaxis().SetTitleSize(0.04)
hSignif.GetZaxis().SetTitleOffset(1)
inFile.Close()

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
ROOT.gPad.SetRightMargin(0.13)
ROOT.gPad.SetLeftMargin(0.1)
ROOT.gPad.SetBottomMargin(0.1)
ROOT.gPad.SetTopMargin(0.05)

hFrame = canvas.DrawFrame(0.1, 0.1, 0.65, 0.65, ";Maximum background score;Minimum prompt score;Significance")
hFrame.GetXaxis().SetTitleSize(0.04)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetYaxis().SetTitleSize(0.04)
hFrame.GetYaxis().SetTitleOffset(1)

hSignif.GetZaxis().SetRangeUser(22,24.5)
hSignif.Smooth()
hSignif.Draw("colz same")
ROOT.gPad.Modified()
ROOT.gPad.Update()
ROOT.gPad.SetGrid()
ROOT.gPad.RedrawAxis()

gWP = ROOT.TGraph(1)
gWP.SetPoint(0, 0.25, 0.2)
gWP.SetMarkerStyle(ROOT.kFullCrossX)
gWP.SetMarkerSize(3)
gWP.SetMarkerColor(ROOT.kTeal-7)
gWP.Draw("P same")

gWPBorder = gWP.Clone()
gWPBorder.SetMarkerStyle(67)
gWPBorder.SetMarkerColor(ROOT.kBlack)
gWPBorder.Draw("P same")

thesisText = ROOT.TLatex(0.4, 0.7, "This Thesis")
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.07)
thesisText.Draw()

ppText = ROOT.TLatex(0.4, 0.65, "pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
ppText.SetNDC()
ppText.SetTextFont(42)
ppText.SetTextSize(0.05)
ppText.Draw()

DecayText = ROOT.TLatex(0.4, 0.5, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.4, 0.45, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ptText = ROOT.TLatex(0.4, 0.37, '2.0 < #it{p}_{T} < 2.5 Ge#kern[-0.03]{V}/#it{c}')
ptText.SetNDC()
ptText.SetTextFont(42)
ptText.SetTextSize(0.05)
ptText.Draw("same")

legend = ROOT.TLegend(0.34, 0.37, 0.5, 0.42)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.05)
legend.AddEntry(gWP, "Working point", "p")
#legend.Draw()

legendBorder = ROOT.TLegend(0.34, 0.37, 0.5, 0.42)
legendBorder.SetBorderSize(0)
legendBorder.SetFillStyle(0)
legendBorder.SetTextFont(42)
legendBorder.SetTextSize(0.05)
legendBorder.AddEntry(gWPBorder, "Working point", "p")
#legendBorder.Draw()


canvas.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Significance_Scan/Significance_Scan_2_2p5.pdf")
canvas.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Significance_Scan/Significance_Scan_2_2p5.png")




