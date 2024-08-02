import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette

infile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/Data/RawYields_Data.root")
hRawYieldsMean = infile.Get("hRawYieldsMean")
hRawYieldsMean.SetDirectory(0)
hRawYieldsSigma = infile.Get("hRawYieldsSigma")
hRawYieldsSigma.SetDirectory(0)
hRawYieldsMeanSecondPeak = infile.Get("hRawYieldsMeanSecondPeak")
hRawYieldsMeanSecondPeak.SetDirectory(0)
hRawYieldsSigmaSecondPeak = infile.Get("hRawYieldsSigmaSecondPeak")
hRawYieldsSigmaSecondPeak.SetDirectory(0)
hRawYieldsSignificance = infile.Get("hRawYieldsSignificance")
hRawYieldsSignificance.SetDirectory(0)
hRawYieldsSignificanceSecondPeak = infile.Get("hRawYieldsSignificanceSecondPeak")
hRawYieldsSignificanceSecondPeak.SetDirectory(0)
infile.Close()

# Set style
colors, _ = get_discrete_matplotlib_palette("tab10")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

hRawYieldsMean.SetLineColor(colors[0])
hRawYieldsMean.SetMarkerColor(colors[0])
hRawYieldsMean.SetMarkerStyle(ROOT.kFullCircle)
hRawYieldsMean.SetMarkerSize(1.5)
hRawYieldsMean.SetLineWidth(2)

hRawYieldsSigma.SetLineColor(colors[0])
hRawYieldsSigma.SetMarkerColor(colors[0])
hRawYieldsSigma.SetMarkerStyle(ROOT.kFullCircle)
hRawYieldsSigma.SetMarkerSize(1.5)
hRawYieldsSigma.SetLineWidth(2)

hRawYieldsMeanSecondPeak.SetLineColor(colors[1])
hRawYieldsMeanSecondPeak.SetMarkerColor(colors[1])
hRawYieldsMeanSecondPeak.SetMarkerStyle(ROOT.kFullDoubleDiamond)
hRawYieldsMeanSecondPeak.SetMarkerSize(2)
hRawYieldsMeanSecondPeak.SetLineWidth(2)

hRawYieldsSigmaSecondPeak.SetLineColor(colors[1])
hRawYieldsSigmaSecondPeak.SetMarkerColor(colors[1])
hRawYieldsSigmaSecondPeak.SetMarkerStyle(ROOT.kFullDoubleDiamond)
hRawYieldsSigmaSecondPeak.SetMarkerSize(2)
hRawYieldsSigmaSecondPeak.SetLineWidth(2)

hRawYieldsSignificance.SetLineColor(colors[0])
hRawYieldsSignificance.SetMarkerColor(colors[0])
hRawYieldsSignificance.SetMarkerStyle(ROOT.kFullCircle)
hRawYieldsSignificance.SetMarkerSize(1.5)
hRawYieldsSignificance.SetLineWidth(2)

hRawYieldsSignificanceSecondPeak.SetLineColor(colors[1])
hRawYieldsSignificanceSecondPeak.SetMarkerColor(colors[1])
hRawYieldsSignificanceSecondPeak.SetMarkerStyle(ROOT.kFullDoubleDiamond)
hRawYieldsSignificanceSecondPeak.SetMarkerSize(2)
hRawYieldsSignificanceSecondPeak.SetLineWidth(2)

# Draw histograms
cMean = ROOT.TCanvas("cMean", "cMean", 800, 600)
cMean.SetRightMargin(0.05)
cMean.SetLeftMargin(0.13)
cMean.SetBottomMargin(0.12)
cMean.SetTopMargin(0.05)

hFrame = cMean.DrawFrame(0, 1.85, 24, 1.98, ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});Peak mean (Ge#kern[-0.05]{V}/#it{c}^{2})")
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetXaxis().SetLabelSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.05)

hRawYieldsMean.Draw("same")
hRawYieldsMeanSecondPeak.Draw("same")

PDGLineDplus = ROOT.TLine(0, 1.86966, 24, 1.86966)
PDGLineDplus.SetLineStyle(2)
PDGLineDplus.SetLineWidth(2)
PDGLineDplus.SetLineColor(colors[1])
PDGLineDplus.Draw("same")

ROOT.gStyle.SetLineStyleString(11, "30 10 10 10 ")
PDGLineDs = ROOT.TLine(0, 1.96835, 24, 1.96835)
PDGLineDs.SetLineStyle(11)
PDGLineDs.SetLineWidth(2)
PDGLineDs.SetLineColor(colors[0])
PDGLineDs.Draw("same")

legMean = ROOT.TLegend(0.7, 0.4, 0.9, 0.7)
legMean.SetBorderSize(0)
legMean.SetFillStyle(0)
legMean.SetTextSize(0.05)
legMean.AddEntry(hRawYieldsMean, "D_{s}^{+}", "lp")
legMean.AddEntry(hRawYieldsMeanSecondPeak, "D^{+}", "lp")
legMean.AddEntry(PDGLineDs, "PDG D_{s}^{+}", "l")
legMean.AddEntry(PDGLineDplus, "PDG D^{+}", "l")
legMean.Draw()


thesisText = ROOT.TLatex(0.2, 0.65, 'This Thesis')
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.06)
thesisText.Draw("same")

collisionText = ROOT.TLatex(0.2, 0.59, 'pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.05]{V}')
collisionText.SetNDC()
collisionText.SetTextFont(42)
collisionText.SetTextSize(0.05)
collisionText.Draw("same")

DecayText = ROOT.TLatex(0.2, 0.5, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.2, 0.44, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ROOT.gPad.RedrawAxis()

cMean.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Mean.png")
cMean.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Mean.pdf")

cSigma = ROOT.TCanvas("cSigma", "cSigma", 800, 600)
cSigma.SetRightMargin(0.05)
cSigma.SetLeftMargin(0.1)
cSigma.SetBottomMargin(0.12)
cSigma.SetTopMargin(0.05)

hFrame = cSigma.DrawFrame(0, 0, 24, 40, ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});Peak width (Me#kern[-0.05]{V}/#it{c}^{2})")
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(0.9)
hFrame.GetXaxis().SetLabelSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.05)

hRawYieldsSigma.Scale(1000)
hRawYieldsSigma.Draw("same")
hRawYieldsSigmaSecondPeak.Scale(1000)
hRawYieldsSigmaSecondPeak.Draw("same")

legSigma = ROOT.TLegend(0.7, 0.2, 0.9, 0.35)
legSigma.SetBorderSize(0)
legSigma.SetFillStyle(0)
legSigma.SetTextSize(0.05)
legSigma.AddEntry(hRawYieldsSigma, "D_{s}^{+}", "lp")
legSigma.AddEntry(hRawYieldsSigmaSecondPeak, "D^{+}", "lp")
legSigma.Draw()

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

cSigma.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Sigma.png")
cSigma.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Sigma.pdf")

cSignificance = ROOT.TCanvas("cSignificance", "cSignificance", 800, 600)
cSignificance.SetRightMargin(0.05)
cSignificance.SetLeftMargin(0.1)
cSignificance.SetBottomMargin(0.12)
cSignificance.SetTopMargin(0.05)

hFrame = cSignificance.DrawFrame(0, 0, 24, 125, ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});Significance")
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.)
hFrame.GetXaxis().SetLabelSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.05)

hRawYieldsSignificance.Draw("same")
hRawYieldsSignificanceSecondPeak.Draw("same")

legSignificance = ROOT.TLegend(0.65, 0.45, 0.9, 0.6)
legSignificance.SetBorderSize(0)
legSignificance.SetFillStyle(0)
legSignificance.SetTextSize(0.05)
legSignificance.AddEntry(hRawYieldsSignificance, "D_{s}^{+}", "lp")
legSignificance.AddEntry(hRawYieldsSignificanceSecondPeak, "D^{+}", "lp")
legSignificance.Draw()

thesisText = ROOT.TLatex(0.5, 0.85, 'This Thesis')
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.06)
thesisText.Draw("same")

collisionText = ROOT.TLatex(0.5, 0.79, 'pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.05]{V}')
collisionText.SetNDC()
collisionText.SetTextFont(42)
collisionText.SetTextSize(0.05)
collisionText.Draw("same")

DecayText = ROOT.TLatex(0.5, 0.7, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.5, 0.64, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ROOT.gPad.RedrawAxis()

cSignificance.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Significance.png")
cSignificance.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Significance.pdf")

infile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Fits_with_frac.root")
hFrac_Bkg = infile.Get("hFracs_background0")
hFrac_Bkg.SetDirectory(0)
infile.Close()

for i in range(1, hFrac_Bkg.GetNbinsX()+1):
    hFrac_Bkg.SetBinError(i, 0)

hFrac_Bkg.SetLineColor(colors[0])
hFrac_Bkg.SetMarkerColor(colors[0])
hFrac_Bkg.SetMarkerStyle(ROOT.kFullCircle)
hFrac_Bkg.SetMarkerSize(1.5)
hFrac_Bkg.SetLineWidth(2)
hFrac_Bkg.SetLineStyle(11)


cFrac = ROOT.TCanvas("cFrac", "cFrac", 800, 600)
cFrac.SetRightMargin(0.05)
cFrac.SetLeftMargin(0.12)
cFrac.SetBottomMargin(0.12)
cFrac.SetTopMargin(0.05)
#cFrac.SetLogy()

hFrame = cFrac.DrawFrame(0, 0, 24, 0.06, ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});Relative contribution of D^{+}#rightarrow #pi^{+}K^{#font[122]{-}}#pi^{+}")
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetXaxis().SetLabelSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.05)

hFrac_Bkg.Draw("hist,same")
hFrac_Bkg.DrawClone("p0,same")

thesisText = ROOT.TLatex(0.2, 0.85, 'This Thesis')
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.06)
thesisText.Draw("same")

collisionText = ROOT.TLatex(0.2, 0.79, 'pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.05]{V}')
collisionText.SetNDC()
collisionText.SetTextFont(42)
collisionText.SetTextSize(0.05)
collisionText.Draw("same")

DecayText = ROOT.TLatex(0.2, 0.7, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.2, 0.64, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ROOT.gPad.RedrawAxis()

cFrac.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Frac.png")
cFrac.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FitParameters/Frac.pdf")

