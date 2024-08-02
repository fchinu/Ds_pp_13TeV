import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette

colors, ROOTcols = get_discrete_matplotlib_palette('tab10')
set_matplotlib_palette("coolwarm")

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

# Get Ds histos
infileDs = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root")

# Correlation matrix
hCov = infileDs.Get("hCorrMatrixCutSets_pt15_20")
hCov.SetDirectory(0)

# Efficiency
hEffPropmtDs = infileDs.Get("hEffPromptVsCut_pt15_20")
hEffPropmtDs.SetDirectory(0)
hEffNonPromptDs = infileDs.Get("hEffNonPromptVsCut_pt15_20")
hEffNonPromptDs.SetDirectory(0)

# Prompt fraction
hPromptFracDs = infileDs.Get("hFracPromptVsCut_pt15_20")
hPromptFracDs.SetDirectory(0)
hNonPromptFracDs = infileDs.Get("hFracNonPromptVsCut_pt15_20")
hNonPromptFracDs.SetDirectory(0)

# Raw yields
hRawDs = infileDs.Get("hRawYieldVsCut_pt15_20")
hRawDs.SetDirectory(0)
hRawPromptDs = infileDs.Get("hRawYieldPromptVsCut_pt15_20")
hRawPromptDs.SetDirectory(0)
hRawNonPromptDs = infileDs.Get("hRawYieldNonPromptVsCut_pt15_20")
hRawNonPromptDs.SetDirectory(0)
hRawSumDs = infileDs.Get("hRawYieldSumVsCut_pt15_20")
hRawSumDs.SetDirectory(0)


infileDs.Close()

# Set style

# Efficiency
hEffPropmtDs.SetLineColor(colors[3])
hEffPropmtDs.SetLineWidth(2)
hEffPropmtDs.SetMarkerColor(colors[3])
hEffPropmtDs.SetMarkerStyle(ROOT.kFullCircle)
hEffPropmtDs.SetMarkerSize(1)

hEffNonPromptDs.SetLineColor(colors[0])
hEffNonPromptDs.SetLineWidth(2)
hEffNonPromptDs.SetMarkerColor(colors[0])
hEffNonPromptDs.SetMarkerStyle(ROOT.kFullCircle)
hEffNonPromptDs.SetMarkerSize(1)

# Prompt fraction
hPromptFracDs.SetLineColor(colors[3])
hPromptFracDs.SetLineWidth(2)
hPromptFracDs.SetMarkerColor(colors[3])
hPromptFracDs.SetMarkerStyle(ROOT.kFullCircle)
hPromptFracDs.SetMarkerSize(1)

hNonPromptFracDs.SetLineColor(colors[0])
hNonPromptFracDs.SetLineWidth(2)
hNonPromptFracDs.SetMarkerColor(colors[0])
hNonPromptFracDs.SetMarkerStyle(ROOT.kFullCircle)
hNonPromptFracDs.SetMarkerSize(1)

# Raw yields
hRawDs.SetLineColor(ROOT.kBlack)
hRawDs.SetMarkerColor(ROOT.kBlack) 
hRawDs.SetMarkerStyle(ROOT.kFullCircle)
hRawDs.SetMarkerSize(1)
hRawPromptDs.SetLineColor(colors[3])
hRawPromptDs.SetLineWidth(2)
hRawPromptDs.SetMarkerColor(colors[3])
hRawPromptDs.SetMarkerStyle(ROOT.kFullCircle)
hRawPromptDs.SetMarkerSize(1)
hRawPromptDs.SetFillStyle(1001)
hRawPromptDs.SetFillColorAlpha(colors[3],0.5)
hRawNonPromptDs.SetLineColor(colors[0])
hRawNonPromptDs.SetLineWidth(2)
hRawNonPromptDs.SetMarkerColor(colors[0])
hRawNonPromptDs.SetMarkerStyle(ROOT.kFullCircle)
hRawNonPromptDs.SetMarkerSize(1)
hRawNonPromptDs.SetFillStyle(1001)
hRawNonPromptDs.SetFillColorAlpha(colors[0],0.5)
hRawSumDs.SetLineColor(colors[2])
hRawSumDs.SetLineWidth(2)
hRawSumDs.SetMarkerColor(colors[2])
hRawSumDs.SetMarkerStyle(ROOT.kFullDiamond)
hRawSumDs.SetMarkerSize(1)


c = ROOT.TCanvas("c","c",800,600)
c.Divide(2,2)

c.cd(1)

ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.11)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetRightMargin(0.15)

hFrame = ROOT.gPad.DrawFrame(-0.5,-0.5,16.5,16.5,";BDT selection;BDT selection;#rho")
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetXaxis().SetTitleOffset(1.)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.045)
hFrame.GetXaxis().SetLabelSize(0.045)

#hCov.GetZaxis().SetRangeUser(0,1)
hCov.GetZaxis().SetTitle("#rho")
hCov.GetZaxis().SetTitleSize(0.05)
hCov.GetZaxis().SetTitleOffset(0.6)
hCov.Draw("colz, same")

ROOT.gPad.RedrawAxis()

c.cd(2)

ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.11)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetRightMargin(0.05)

hFrame = ROOT.gPad.DrawFrame(-0.5,2.e-5,16.5,1.,";BDT selection;Acceptance #times Efficiency")
hFrame.GetYaxis().SetTitleOffset(1.25)
hFrame.GetXaxis().SetTitleOffset(1.)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.045)
hFrame.GetXaxis().SetLabelSize(0.045)
ROOT.gPad.SetLogy()

hEffPropmtDs.Draw("same")
hEffNonPromptDs.Draw("same")

legEff = ROOT.TLegend(0.15,0.15,0.45,0.3)
legEff.SetBorderSize(0)
legEff.SetFillStyle(0)
legEff.SetTextSize(0.04)
legEff.AddEntry(hEffPropmtDs,"Prompt D_{#lower[-0.2]{s}}^{+}","lp")
legEff.AddEntry(hEffNonPromptDs,"Non-prompt D_{#lower[-0.2]{s}}^{+}","lp")
legEff.Draw("same")

thesisText = ROOT.TLatex(0.5, 0.865, "This Thesis")
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.07)
thesisText.Draw()

ppText = ROOT.TLatex(0.5, 0.8, "pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
ppText.SetNDC()
ppText.SetTextFont(42)
ppText.SetTextSize(0.05)
ppText.Draw()

DecayText = ROOT.TLatex(0.5, 0.72, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.5, 0.67, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ptText = ROOT.TLatex(0.5, 0.6, '1.5 < #it{p}_{T} < 2.0 Ge#kern[-0.03]{V}/#it{c}')
ptText.SetNDC()
ptText.SetTextFont(42)
ptText.SetTextSize(0.05)
ptText.Draw("same")

ROOT.gPad.RedrawAxis()

c.cd(3)

ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.11)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetRightMargin(0.05)

hFrame = ROOT.gPad.DrawFrame(-0.5,0,16.5,1,";BDT selection;Fraction")
hFrame.GetYaxis().SetTitleOffset(1.25)
hFrame.GetXaxis().SetTitleOffset(1.)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.045)
hFrame.GetXaxis().SetLabelSize(0.045)

hPromptFracDs.Draw("same")
hNonPromptFracDs.Draw("same")

legFrac = ROOT.TLegend(0.4,0.75,0.7,0.9)
legFrac.SetBorderSize(0)
legFrac.SetFillStyle(0)
legFrac.SetTextSize(0.04)
legFrac.AddEntry(hPromptFracDs,"Prompt D_{#lower[-0.2]{s}}^{+}","lp")
legFrac.AddEntry(hNonPromptFracDs,"Non-prompt D_{#lower[-0.2]{s}}^{+}","lp")
legFrac.Draw("same")

ROOT.gPad.RedrawAxis()

c.cd(4)

ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.11)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetRightMargin(0.05)

hFrame = ROOT.gPad.DrawFrame(-0.5,0,16.5,5000,";BDT selection;Raw yield")
hFrame.GetYaxis().SetTitleOffset(1.25)
hFrame.GetXaxis().SetTitleOffset(1.)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetLabelSize(0.045)
hFrame.GetXaxis().SetLabelSize(0.045)

hRawDs.Draw("same")
hRawPromptDs.Draw("hist, same")
hRawNonPromptDs.Draw("hist, same")
hRawSumDs.Draw("hist, same")

chi2text = ROOT.TLatex(0.515,0.85,"#chi^{2}/ndf = 10.33/15")
chi2text.SetNDC()
chi2text.SetTextFont(42)
chi2text.SetTextSize(0.041)
chi2text.Draw("same")

legRY = ROOT.TLegend(0.5,0.6,0.8,0.84)
legRY.SetBorderSize(0)
legRY.SetFillStyle(0)
legRY.SetTextSize(0.04)
legRY.AddEntry(hRawDs,"Data","lp")
legRY.AddEntry(hRawPromptDs,"Prompt D_{#lower[-0.2]{s}}^{+}","f")
legRY.AddEntry(hRawNonPromptDs,"Non-prompt D_{#lower[-0.2]{s}}^{+}","f")
legRY.AddEntry(hRawSumDs,"Prompt + Non-prompt D_{#lower[-0.2]{s}}^{+}","l")
legRY.Draw("same")


ROOT.gPad.RedrawAxis()

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/fprompt/DsPromptFrac.pdf")
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/fprompt/DsPromptFrac.png")



