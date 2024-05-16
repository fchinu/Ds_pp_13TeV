import ROOT 
import numpy as np

TransparentBlue = ROOT.TColor.GetFreeColorIndex()
cTransparentBlue = ROOT.TColor(TransparentBlue, 79./255, 123./255, 161./255, 'TransparentBlue', 1.0)

# LHCb comparison

ratioFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/DsOverDplus_LHC24d3a.root'
systFileNameRawYields = '/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/SystUncertainties.root'
systFileNameBDT = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/output/BDT_syst.root"
systFileNameFD = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/FDSystematic.root"
systFileNameTopo = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/SystUncertainty.root"
LHCbFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/ds_over_dplus_lhcb.root'

ratioFile = ROOT.TFile(ratioFileName)
hRatio = ratioFile.Get('hRatio')
hRatio.SetDirectory(0)
ratioFile.Close()

systFileRawYields = ROOT.TFile(systFileNameRawYields)
hSystRawYields = systFileRawYields.Get('hSyst')
hSystRawYields.SetDirectory(0)
systFileRawYields.Close()

systFileBDT = ROOT.TFile(systFileNameBDT)
hSystBDT = systFileBDT.Get('hSyst')
hSystBDT.SetDirectory(0)
systFileBDT.Close()

systFileFD = ROOT.TFile(systFileNameFD)
hSystFD = systFileFD.Get('hSyst')
hSystFD.SetDirectory(0)
systFileFD.Close()

systFileTopo = ROOT.TFile(systFileNameTopo)
hSystTopo = systFileTopo.Get('hSyst')
hSystTopo.SetDirectory(0)
systFileTopo.Close()

gSyst = ROOT.TGraphErrors(hSystRawYields)

for iPt in range(gSyst.GetN()):
    unc = np.sqrt(hSystRawYields.GetBinContent(iPt+1)**2 + hSystBDT.GetBinContent(iPt+1)**2 + hSystFD.GetBinContent(iPt+1)**2 + hSystTopo.GetBinContent(iPt+1)**2)
    gSyst.SetPointError(iPt, 0.2, unc*hRatio.GetBinContent(iPt+1))
    gSyst.SetPointY(iPt, hRatio.GetBinContent(iPt+1))

LHCbFile = ROOT.TFile(LHCbFileName)
hLHCb_20_25 = LHCbFile.Get('h_ds_over_dp_rap0')
hLHCb_20_25.SetDirectory(0)
hLHCb_25_30 = LHCbFile.Get('h_ds_over_dp_rap1')
hLHCb_25_30.SetDirectory(0)
gLHCb_20_25 = LHCbFile.Get('gsys_ds_over_dp_rap0')
gLHCb_25_30 = LHCbFile.Get('gsys_ds_over_dp_rap1')
LHCbFile.Close()

# Set style
hRatio.SetLineColor(ROOT.kRed)
hRatio.SetMarkerColor(ROOT.kRed)
hRatio.SetMarkerStyle(ROOT.kFullCircle)
hRatio.SetMarkerSize(1.5)
hRatio.SetLineWidth(2)

hRatioBorder = hRatio.Clone('hRatioBorder')
hRatioBorder.SetLineColor(ROOT.kRed)
hRatioBorder.SetMarkerColor(ROOT.kBlack)
hRatioBorder.SetMarkerStyle(53)
hRatioBorder.SetLineWidth(3)

gSyst.SetLineColor(ROOT.kRed)
gSyst.SetMarkerColor(ROOT.kRed)
gSyst.SetMarkerStyle(ROOT.kFullCircle)
gSyst.SetMarkerSize(1.5)
gSyst.SetLineWidth(2)
gSyst.SetFillStyle(0)

hLHCb_20_25.SetLineColor(ROOT.kAzure+7)
hLHCb_20_25.SetMarkerColor(ROOT.kAzure+7)
hLHCb_20_25.SetMarkerStyle(ROOT.kFullSquare)
hLHCb_20_25.SetMarkerSize(1.5)
hLHCb_20_25.SetLineWidth(2)

hLHCb_20_25Border = hLHCb_20_25.Clone('hLHCb_20_25Border')
hLHCb_20_25Border.SetLineColor(ROOT.kAzure+7)
hLHCb_20_25Border.SetMarkerColor(ROOT.kBlack)
hLHCb_20_25Border.SetMarkerStyle(54)
hLHCb_20_25Border.SetLineWidth(3)

gLHCb_20_25.SetLineColor(ROOT.kAzure+7)
gLHCb_20_25.SetFillStyle(0)
gLHCb_20_25.SetLineWidth(2)
for iPt in range(gLHCb_20_25.GetN()):
    gLHCb_20_25.SetPointEXhigh(iPt, 0.2)
    gLHCb_20_25.SetPointEXlow(iPt, 0.2)

hLHCb_25_30.SetLineColor(ROOT.kGreen-2)
hLHCb_25_30.SetMarkerColor(ROOT.kGreen-2)
hLHCb_25_30.SetMarkerStyle(ROOT.kFullCross)
hLHCb_25_30.SetMarkerSize(1.5)
hLHCb_25_30.SetLineWidth(2)

hLHCb_25_30Border = hLHCb_25_30.Clone('hLHCb_25_30Border')
hLHCb_25_30Border.SetLineColor(ROOT.kGreen-2)
hLHCb_25_30Border.SetMarkerColor(ROOT.kBlack)
hLHCb_25_30Border.SetMarkerStyle(57)
hLHCb_25_30Border.SetLineWidth(3)

gLHCb_25_30.SetLineColor(ROOT.kGreen-2)
gLHCb_25_30.SetFillStyle(0)
gLHCb_25_30.SetLineWidth(2)
for iPt in range(gLHCb_25_30.GetN()):
    gLHCb_25_30.SetPointEXhigh(iPt, 0.2)
    gLHCb_25_30.SetPointEXlow(iPt, 0.2)

canvas = ROOT.TCanvas('canvas', 'canvas', 800, 800)
ROOT.gPad.SetLeftMargin(0.1)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.1)

hFrame = canvas.DrawFrame(0, 0, 24, 1., ';#it{p}_{T} (GeV/c);D_{s}^{+}/D^{+}')
hFrame.GetXaxis().SetTitleOffset(1.05)
hFrame.GetYaxis().SetTitleOffset(1.05)
hFrame.GetXaxis().SetTitleSize(0.043)
hFrame.GetYaxis().SetTitleSize(0.043)
gLHCb_20_25.Draw('5 same')
gLHCb_25_30.Draw('5 same')
gSyst.Draw('5 same')
hLHCb_20_25.Draw('p same')
hLHCb_25_30.Draw('p same')
hLHCb_20_25Border.Draw('p same')
hLHCb_25_30Border.Draw('p same')
hRatio.Draw('p same')
hRatioBorder.Draw('p same')
canvas.Update()

textALICE = ROOT.TLatex(0.15, 0.9, 'ALICE Preliminary')
textALICE.SetNDC()
textALICE.SetTextFont(42)
textALICE.SetTextSize(0.05)
textALICE.Draw()

textpp = ROOT.TLatex(0.155, 0.86, 'pp collisions')
textpp.SetNDC()
textpp.SetTextFont(42)
textpp.SetTextSize(0.04)
textpp.Draw()

texty = ROOT.TLatex(0.155, 0.81, '#sqrt{s} = 13.6 TeV')
texty.SetNDC()
texty.SetTextFont(42)
texty.SetTextSize(0.04)
texty.Draw()

textDsDecay = ROOT.TLatex(0.155, 0.76, 'D_{s}^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+} ')
textDsDecay.SetNDC()
textDsDecay.SetTextFont(42)
textDsDecay.SetTextSize(0.035)
textDsDecay.Draw()

textDsDecayConj = ROOT.TLatex(0.155, 0.72, 'and charge conj.')
textDsDecayConj.SetNDC()
textDsDecayConj.SetTextFont(42)
textDsDecayConj.SetTextSize(0.035)
textDsDecayConj.Draw()

textEnergy = ROOT.TLatex(0.155, 0.85, '#sqrt{s} = 13.6 TeV')
textEnergy.SetNDC()
textEnergy.SetTextFont(42)
textEnergy.SetTextSize(0.035)
#textEnergy.Draw()

textDecay = ROOT.TLatex(0.6, 0.9, 'D^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+} ')
textDecay.SetNDC()
textDecay.SetTextFont(42)
textDecay.SetTextSize(0.035)
textDecay.Draw()

textDecayConj = ROOT.TLatex(0.6, 0.86, 'and charge conj.')
textDecayConj.SetNDC()
textDecayConj.SetTextFont(42)
textDecayConj.SetTextSize(0.035)
textDecayConj.Draw()

textBR = ROOT.TLatex(0.6, 0.82, 'BR unc. (not shown): ^{#lower[5]{+3.7}}_{#lower[-0.2]{#font[122]{-}4.0}}#it{%} ')
textBR.SetNDC()
textBR.SetTextFont(42)
textBR.SetTextSize(0.03)
textBR.Draw()

Thislegend = ROOT.TLegend(0.6, 0.77, 0.8, 0.8)
Thislegend.SetBorderSize(0)
Thislegend.SetTextSize(0.027)
Thislegend.SetFillStyle(0)
Thislegend.AddEntry(hRatio, '|#it{y}|<0.5', 'lp')
Thislegend.Draw()

ThislegendBorder = ROOT.TLegend(0.6, 0.77, 0.8, 0.8)
ThislegendBorder.SetBorderSize(0)
ThislegendBorder.SetTextSize(0.027)
ThislegendBorder.SetFillStyle(0)
ThislegendBorder.AddEntry(hRatioBorder, '|#it{y}|<0.5', 'lp')
ThislegendBorder.Draw()

text = ROOT.TLatex(0.25, 0.35, 'D^{+}#rightarrow #pi^{+}K^{#font[122]{-}}#pi^{+}')
text.SetNDC()
text.SetTextFont(42)
text.SetTextSize(0.035)
text.Draw()

textConj = ROOT.TLatex(0.25, 0.31, 'and charge conj.')
textConj.SetNDC()
textConj.SetTextFont(42)
textConj.SetTextSize(0.035)
textConj.Draw()

textBRGoldenDplus = ROOT.TLatex(0.25, 0.27, 'BR unc. (not shown): 3.2#it{%}')
textBRGoldenDplus.SetNDC()
textBRGoldenDplus.SetTextFont(42)
textBRGoldenDplus.SetTextSize(0.03)
textBRGoldenDplus.Draw()

LegLHCbText = ROOT.TLatex(0.25, 0.23, 'LHCb, #sqrt{s} = 13 TeV, JHEP03(2016)159')
LegLHCbText.SetNDC()
LegLHCbText.SetTextFont(42)
LegLHCbText.SetTextSize(0.03)
LegLHCbText.Draw()

legLHCb = ROOT.TLegend(0.25, 0.13, 0.45, 0.21)
legLHCb.SetBorderSize(0)
legLHCb.SetTextSize(0.027)
legLHCb.SetFillStyle(0)
legLHCb.AddEntry(hLHCb_20_25, '2.0 < #it{y} < 2.5', 'lp')
legLHCb.AddEntry(hLHCb_25_30, '2.5 < #it{y} < 3.0', 'lp')
legLHCb.Draw()

legLHCbBorder = ROOT.TLegend(0.25, 0.13, 0.45, 0.21)
legLHCbBorder.SetBorderSize(0)
legLHCbBorder.SetTextSize(0.027)
legLHCbBorder.SetFillStyle(0)
legLHCbBorder.AddEntry(hLHCb_20_25Border, '2.0 < #it{y} < 2.5', 'lp')
legLHCbBorder.AddEntry(hLHCb_25_30Border, '2.5 < #it{y} < 3.0', 'lp')
legLHCbBorder.Draw()

canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonLHCb.png')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonLHCb.pdf')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonLHCb.eps')
