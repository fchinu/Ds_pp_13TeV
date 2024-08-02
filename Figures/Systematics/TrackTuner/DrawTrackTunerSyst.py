import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette

# Get inputs
infile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/CompareResults.root")
hCentral = infile.Get("hRatioCentral")
hCentral.SetDirectory(0)
hPtSmear = infile.Get("hRatioPtSmearing")
hPtSmear.SetDirectory(0)
hWorseReso = infile.Get("hRatioWorseResoPulls")
hWorseReso.SetDirectory(0)

# Set style
colors, _ = get_discrete_matplotlib_palette("tab10")

hCentral.SetLineColor(colors[3])
hCentral.SetMarkerColor(colors[3])
hCentral.SetMarkerStyle(ROOT.kFullCircle)
hCentral.SetMarkerSize(1.5)
hCentral.SetLineWidth(2)

hPtSmear.SetLineColor(colors[0])
hPtSmear.SetMarkerColor(colors[0])
hPtSmear.SetMarkerStyle(ROOT.kFullDiamond)
hPtSmear.SetMarkerSize(2)
hPtSmear.SetLineWidth(2)

hWorseReso.SetLineColor(colors[1])
hWorseReso.SetMarkerColor(colors[1])
hWorseReso.SetMarkerStyle(ROOT.kFullDoubleDiamond)
hWorseReso.SetMarkerSize(2)
hWorseReso.SetLineWidth(2)

# Draw
c = ROOT.TCanvas("c", "c", 1000, 600)
c.Divide(2,1,0.0001,0.0001)

c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetBottomMargin(0.12)
ROOT.gPad.SetTopMargin(0.05)

hFrame = ROOT.gPad.DrawFrame(0., 0.15, 24, 0.9,";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});D_{s}^{+}/D^{+}")
hFrame.GetXaxis().SetTitleOffset(1.1)
hFrame.GetYaxis().SetTitleOffset(1.3)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetLabelSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.05)
ROOT.gPad.SetTickx(1)
ROOT.gPad.SetTicky(1)

hCentral.Draw("same")
hPtSmear.Draw("same")
hWorseReso.Draw("same")

thesisText = ROOT.TLatex(0.2, 0.85, "This Thesis")
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.07)
thesisText.Draw()

ppText = ROOT.TLatex(0.2, 0.79, "pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
ppText.SetNDC()
ppText.SetTextFont(42)
ppText.SetTextSize(0.05)
ppText.Draw()

DecayText = ROOT.TLatex(0.2, 0.71, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.2, 0.66, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ROOT.gPad.RedrawAxis()

c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetBottomMargin(0.12)
ROOT.gPad.SetTopMargin(0.05)

hFrame = ROOT.gPad.DrawFrame(0., 0.6, 24, 1.3,";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});Ratio to central case")
hFrame.GetXaxis().SetTitleOffset(1.1)
hFrame.GetYaxis().SetTitleOffset(1.3)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetLabelSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.05)
ROOT.gPad.SetTickx(1)
ROOT.gPad.SetTicky(1)

ROOT.gStyle.SetLineStyleString(11,"50 15")
lineAtOne = ROOT.TLine(0, 1, 24, 1)
lineAtOne.SetLineStyle(2)
lineAtOne.SetLineColor(ROOT.kGray+1)
lineAtOne.SetLineWidth(2)
lineAtOne.Draw("same")

hRatioPtSmearing = hPtSmear.Clone("hRatioPtSmearing")
hRatioPtSmearing.Divide(hCentral)
hRatioPtSmearing.Draw("same")

hRatioWorseReso = hWorseReso.Clone("hRatioWorseReso")
hRatioWorseReso.Divide(hCentral)
hRatioWorseReso.Draw("same")

leg = ROOT.TLegend(0.25, 0.15, 0.6, 0.3)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.05)
leg.AddEntry(hCentral, "Central", "lp")
leg.AddEntry(hPtSmear, "2#times #it{p}_{T} smearing", "lp")
leg.AddEntry(hWorseReso, "10% worse DCA resolution", "lp")
leg.Draw()

ROOT.gPad.RedrawAxis()

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Systematics/TrackTuner/TrackTunerSyst.pdf")
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Systematics/TrackTuner/TrackTunerSyst.png")



