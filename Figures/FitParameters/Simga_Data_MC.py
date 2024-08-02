import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette


infile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/SigmaDsOverDplusMC.root")
hSigmaDsMC = infile.Get("hSigmaDsMC")
hSigmaDsMC.SetDirectory(0)
hSigmaDplusMC = infile.Get("hSigmaDplusMC")
hSigmaDplusMC.SetDirectory(0)
infile.Close()

infile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Fits_No_fixed_Dplus_sigma.root")
hSigmaDplusData = infile.Get("hRawYieldsSigmaSecondPeak")
hSigmaDplusData.SetDirectory(0)
hSigmaDplusData.Scale(1000)
hSigmaDsData = infile.Get("hRawYieldsSigma")
hSigmaDsData.SetDirectory(0)
hSigmaDsData.Scale(1000)
infile.Close()

# set style
colors, _ = get_discrete_matplotlib_palette("tab10")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

hSigmaDsMC.SetLineColor(colors[0])
hSigmaDsMC.SetMarkerColor(colors[0])
hSigmaDsMC.SetMarkerStyle(ROOT.kOpenCircle)
hSigmaDsMC.SetMarkerSize(1.5)
hSigmaDsMC.SetLineWidth(2)

hSigmaDplusMC.SetLineColor(colors[1])
hSigmaDplusMC.SetMarkerColor(colors[1])
hSigmaDplusMC.SetMarkerStyle(ROOT.kOpenDoubleDiamond)
hSigmaDplusMC.SetMarkerSize(2)
hSigmaDplusMC.SetLineWidth(2)

hSigmaDsData.SetLineColor(colors[0])
hSigmaDsData.SetMarkerColor(colors[0])
hSigmaDsData.SetMarkerStyle(ROOT.kFullCircle)
hSigmaDsData.SetMarkerSize(1.5)
hSigmaDsData.SetLineWidth(2)

hSigmaDplusData.SetLineColor(colors[1])
hSigmaDplusData.SetMarkerColor(colors[1])
hSigmaDplusData.SetMarkerStyle(ROOT.kFullDoubleDiamond)
hSigmaDplusData.SetMarkerSize(2)
hSigmaDplusData.SetLineWidth(2)

# Draw histograms

c = ROOT.TCanvas("c", "c", 800, 600)
c.SetRightMargin(0.05)
c.SetLeftMargin(0.1)
c.SetBottomMargin(0.12)
c.SetTopMargin(0.05)

hFrame = c.DrawFrame(0, 0, 24, 40, ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});Peak width (Me#kern[-0.05]{V}/#it{c}^{2})")
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(0.9)
hFrame.GetXaxis().SetLabelSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.05)

hSigmaDsMC.Draw("same")
hSigmaDplusMC.Draw("same")
hSigmaDsData.Draw("same")
hSigmaDplusData.Draw("same")

leg = ROOT.TLegend(0.6, 0.15, 0.9, 0.4)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.045)
leg.AddEntry(hSigmaDsData, "Data D_{s}^{+}", "lp")
leg.AddEntry(hSigmaDplusData, "Data D^{+}", "lp")
leg.AddEntry(hSigmaDsMC, "MC D_{s}^{+}", "lp")
leg.AddEntry(hSigmaDplusMC, "MC D^{+}", "lp")
leg.Draw()

thesisText = ROOT.TLatex(0.15, 0.85, 'This Thesis')
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.06)
thesisText.Draw("same")

collisionText = ROOT.TLatex(0.15, 0.79, 'pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.05]{V}')
collisionText.SetNDC()
collisionText.SetTextFont(42)
collisionText.SetTextSize(0.05)
collisionText.Draw("same")

DecayText = ROOT.TLatex(0.15, 0.7, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.15, 0.64, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ROOT.gPad.RedrawAxis()

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/SigmaDs_Dplus.pdf")

# ratios
hRatioData = hSigmaDsData.Clone("hRatioData")
hRatioData.Divide(hSigmaDplusData)

hRatioMC = hSigmaDsMC.Clone("hRatioMC")
hRatioMC.Divide(hSigmaDplusMC)

hRatioMC.SetLineColor(colors[2])
hRatioMC.SetMarkerColor(colors[2])
hRatioMC.SetMarkerStyle(ROOT.kOpenCircle)
hRatioMC.SetMarkerSize(1.5)
hRatioMC.SetLineWidth(2)

colors, a = get_discrete_matplotlib_palette("tab20b")


hRatioData.SetLineColor(colors[18])
hRatioData.SetMarkerColor(colors[18])
hRatioData.SetMarkerStyle(ROOT.kFullCircle)
hRatioData.SetMarkerSize(1.5)
hRatioData.SetLineWidth(2)


cRatio = ROOT.TCanvas("cRatio", "cRatio", 800, 600)
cRatio.SetRightMargin(0.05)
cRatio.SetLeftMargin(0.1)
cRatio.SetBottomMargin(0.12)
cRatio.SetTopMargin(0.05)

hFrameRatio = cRatio.DrawFrame(0, 0.7, 24, 2., ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});#sigma(D_{s}^{+})/#sigma(D^{+})")
hFrameRatio.GetXaxis().SetTitleSize(0.05)
hFrameRatio.GetXaxis().SetTitleOffset(1)
hFrameRatio.GetYaxis().SetTitleSize(0.05)
hFrameRatio.GetYaxis().SetTitleOffset(0.9)
hFrameRatio.GetXaxis().SetLabelSize(0.05)
hFrameRatio.GetYaxis().SetLabelSize(0.05)

hRatioData.Draw("same")
hRatioMC.Draw("same")

legRatio = ROOT.TLegend(0.3, 0.15, 0.6, 0.3)
legRatio.SetBorderSize(0)
legRatio.SetFillStyle(0)
legRatio.SetTextSize(0.045)
legRatio.AddEntry(hRatioData, "Data", "lp")
legRatio.AddEntry(hRatioMC, "MC", "lp")
legRatio.Draw()

thesisText = ROOT.TLatex(0.15, 0.85, 'This Thesis')
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.06)
thesisText.Draw("same")

collisionText = ROOT.TLatex(0.15, 0.79, 'pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.05]{V}')
collisionText.SetNDC()
collisionText.SetTextFont(42)
collisionText.SetTextSize(0.05)
collisionText.Draw("same")

DecayText = ROOT.TLatex(0.15, 0.7, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.15, 0.64, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ROOT.gPad.RedrawAxis()

cRatio.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/SigmaDsOverDplusMC.pdf")



