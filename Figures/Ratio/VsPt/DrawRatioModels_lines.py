import ROOT 
import numpy as np


TransparentBlue = ROOT.TColor.GetFreeColorIndex()
cTransparentBlue = ROOT.TColor(TransparentBlue, 79./255, 123./255, 161./255, 'TransparentBlue', 1.0)

# ALICE comparison

ratioFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/DsOverDplus_LHC24d3a.root'
systFileNameRawYields = '/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/SystUncertainties.root'
systFileNameBDT = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/output/BDT_syst.root"
systFileNameFD = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/FDSystematic.root"
systFileNameTopo = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/SystUncertainty.root"
PythiaFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/Models/PYTHIA8_Ds_over_Dplus_pp13dot6TeV.root"
PowlangFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/Models/POWLANG-pp-13-TeV.root"

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
    print(f"{hSystRawYields.GetBinLowEdge(iPt+1)} < pT < {hSystRawYields.GetBinLowEdge(iPt+1)+hSystRawYields.GetBinWidth(iPt+1)}: ", unc)
    gSyst.SetPointError(iPt, 0.2, unc*hRatio.GetBinContent(iPt+1))
    gSyst.SetPointY(iPt, hRatio.GetBinContent(iPt+1))


ModelsFile = ROOT.TFile(PythiaFileName)
hMonash = ModelsFile.Get('hist_ds_over_dplus_prompt_Monash_reb')
hMonash.SetDirectory(0)
hMonash.SetLineColorAlpha(ROOT.kAzure-3,0.6)
hMonash.SetFillColorAlpha(ROOT.kAzure-3,0.5)

hCRMode0 = ModelsFile.Get('hist_ds_over_dplus_prompt_CRMode0_reb')
hCRMode0.SetDirectory(0)
#hCRMode0.SetLineColor(ROOT.kTeal-5)
hCRMode0.SetFillColorAlpha(ROOT.kTeal-5,0.5)
hCRMode0.SetLineColorAlpha(ROOT.kTeal-5,0.5)


hCRMode2 = ModelsFile.Get('hist_ds_over_dplus_prompt_CRMode2_reb')
hCRMode2.SetDirectory(0)
#hCRMode2.SetLineColor(ROOT.kViolet-4)
hCRMode2.SetFillColorAlpha(ROOT.kViolet-4,0.6)
hCRMode2.SetLineColorAlpha(ROOT.kViolet-4,0.6)

hCRMode3 = ModelsFile.Get('hist_ds_over_dplus_prompt_CRMode3_reb')
hCRMode3.SetDirectory(0)
#hCRMode3.SetLineColor(ROOT.kOrange-3)
hCRMode3.SetFillColorAlpha(ROOT.kOrange-3,0.5)
hCRMode3.SetLineColorAlpha(ROOT.kOrange-3,0.5)

ModelsFile.Close()

hMonash = ROOT.TGraphErrors(hMonash)
hMonash.SetLineWidth(2)
hCRMode0 = ROOT.TGraphErrors(hCRMode0)
hCRMode0.SetLineWidth(2)
hCRMode2 = ROOT.TGraphErrors(hCRMode2)
hCRMode2.SetLineWidth(2)
hCRMode3 = ROOT.TGraphErrors(hCRMode3)
hCRMode3.SetLineWidth(2)

ROOT.gStyle.SetLineStyleString(11,"40 15 10 15")
PowlangFile = ROOT.TFile(PowlangFileName)
hHTL = PowlangFile.Get('hHTL')
hHTL.SetDirectory(0)
hHTL.SetLineColorAlpha(ROOT.kAzure-3,0.6)
hHTL.SetLineStyle(2)
hHTL.SetLineWidth(2)
hLQCD = PowlangFile.Get('hlQCD')
hLQCD.SetDirectory(0)
hLQCD.SetLineColor(ROOT.kTeal-5)
hLQCD.SetLineStyle(11)
hLQCD.SetLineWidth(2)
PowlangFile.Close()


#for iPt in range(hMonash.GetNbinsX()):
#    hMonash.SetBinError(iPt+1, 0)
#    hCRMode0.SetBinError(iPt+1, 0)
#    hCRMode2.SetBinError(iPt+1, 0)
#    hCRMode3.SetBinError(iPt+1, 0)



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

canvas = ROOT.TCanvas('canvas', 'canvas', 800, 800)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetRightMargin(0.04)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.1)

hFrame = canvas.DrawFrame(0, 0, 24, 1.3, ';#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});D_{s}^{+}/D^{+}')
hFrame.GetXaxis().SetTitleOffset(1.05)
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetXaxis().SetTitleSize(0.043)
hFrame.GetYaxis().SetTitleSize(0.043)
hFrame.GetXaxis().SetLabelSize(0.04)
hFrame.GetYaxis().SetLabelSize(0.04)
ROOT.gPad.SetTickx(1)
ROOT.gPad.SetTicky(1)
hHTL.Draw('L same')
hLQCD.Draw('L same')
hMonash.Draw('LZX same')
hCRMode0.Draw('LZX same')
hCRMode2.Draw('LZX same')
hCRMode3.Draw('LZX same')
gSyst.Draw('5 same')
hRatio.Draw('p same')
hRatioBorder.Draw('p same')
canvas.Update()

textALICE = ROOT.TLatex(0.15, 0.88, 'ALICE Preliminary')
textALICE.SetNDC()
textALICE.SetTextFont(42)
textALICE.SetTextSize(0.05)
textALICE.Draw()

textpp = ROOT.TLatex(0.15, 0.83, 'pp, #sqrt{s} = 13.6 Te#kern[-0.05]{V}')
textpp.SetNDC()
textpp.SetTextFont(42)
textpp.SetTextSize(0.04)
textpp.Draw()

texty = ROOT.TLatex(0.15, 0.78, '|#it{y}|#kern[0.4]{<}#kern[0.2]{0.5}')
texty.SetNDC()
texty.SetTextFont(42)
texty.SetTextSize(0.04)
texty.Draw()

textDsDecay = ROOT.TLatex(0.15, 0.71, 'D_{s}^{+}, D^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+} ')
textDsDecay.SetNDC()
textDsDecay.SetTextFont(42)
textDsDecay.SetTextSize(0.035)
textDsDecay.Draw()

textDsDecayConj = ROOT.TLatex(0.15, 0.67, 'and charge conj.')
textDsDecayConj.SetNDC()
textDsDecayConj.SetTextFont(42)
textDsDecayConj.SetTextSize(0.035)
textDsDecayConj.Draw()

#textDecay = ROOT.TLatex(0.155, 0.84, 'D^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+} ')
#textDecay.SetNDC()
#textDecay.SetTextFont(42)
#textDecay.SetTextSize(0.035)
#textDecay.Draw()
#
#textDecayConj = ROOT.TLatex(0.155, 0.8, 'and charge conj.')
#textDecayConj.SetNDC()
#textDecayConj.SetTextFont(42)
#textDecayConj.SetTextSize(0.035)
#textDecayConj.Draw()
#
textBR = ROOT.TLatex(0.15, 0.15, 'BR unc. (not shown): ^{#lower[5]{+3.7}}_{#lower[-0.2]{#font[122]{-}4.0}}#it{%} ')
textBR.SetNDC()
textBR.SetTextFont(42)
textBR.SetTextSize(0.03)
textBR.Draw()

Thislegend = ROOT.TLegend(0.6, 0.85, 0.9, 0.9)
Thislegend.SetBorderSize(0)
Thislegend.SetTextSize(0.03)
Thislegend.SetFillStyle(0)
Thislegend.AddEntry(hRatio, 'Data', 'lp')
Thislegend.Draw()

ThislegendBorder = ROOT.TLegend(0.6, 0.85, 0.9, 0.9)
ThislegendBorder.SetBorderSize(0)
ThislegendBorder.SetTextSize(0.03)
ThislegendBorder.SetFillStyle(0)
ThislegendBorder.AddEntry(hRatioBorder, 'Data', 'lp')
ThislegendBorder.Draw()

legendPythia = ROOT.TLegend(0.6, 0.57, 0.9, 0.83)
legendPythia.SetBorderSize(0)
legendPythia.SetTextSize(0.03)
legendPythia.SetFillStyle(0)
legendPythia.SetHeader('PYTHIA 8')
legendPythia.AddEntry(hMonash, 'Monash', 'l')
legendPythia.AddEntry(hCRMode0, 'CR - Mode 0', 'l')
legendPythia.AddEntry(hCRMode2, 'CR - Mode 2', 'l')
legendPythia.AddEntry(hCRMode3, 'CR - Mode 3', 'l')
legendPythia.Draw()


legendPowlang = ROOT.TLegend(0.6, 0.14, 0.9, 0.27)
legendPowlang.SetBorderSize(0)
legendPowlang.SetTextSize(0.03)
legendPowlang.SetFillStyle(0)
legendPowlang.SetHeader('POWLANG')
legendPowlang.AddEntry(hHTL, 'HTL', 'l')
legendPowlang.AddEntry(hLQCD, 'lQCD', 'l')
legendPowlang.Draw()
#
#legendBorder = ROOT.TLegend(0.25, 0.14, 0.4, 0.26)
#legendBorder.SetBorderSize(0)
#legendBorder.SetTextSize(0.027)
#legendBorder.SetFillStyle(0)
#legendBorder.SetHeader('PYTHIA 8')
#legendBorder.AddEntry(hMonash, 'Monash', 'f')
#legendBorder.AddEntry(hCRMode0, 'CR - Mode 0', 'f')
#legendBorder.AddEntry(hCRMode2, 'CR - Mode 2', 'f')
#legendBorder.AddEntry(hCRMode3, 'CR - Mode 3', 'f')
#legendBorder.Draw()

#text = ROOT.TLatex(0.25, 0.35, 'D^{+}#rightarrow #pi^{+}K^{#font[122]{-}}#pi^{+}')
#text.SetNDC()
#text.SetTextFont(42)
#text.SetTextSize(0.035)
#text.Draw()
#
#textConj = ROOT.TLatex(0.25, 0.31, 'and charge conj.')
#textConj.SetNDC()
#textConj.SetTextFont(42)
#textConj.SetTextSize(0.035)
#textConj.Draw()
#
#textBRGoldenDplus = ROOT.TLatex(0.25, 0.27, 'BR unc. (not shown): 3.2#kern[-0.5]{#it{%}}')
#textBRGoldenDplus.SetNDC()
#textBRGoldenDplus.SetTextFont(42)
#textBRGoldenDplus.SetTextSize(0.03)
#textBRGoldenDplus.Draw()




canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonModels_lines.png')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonModels_lines.pdf')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonModels_lines.eps')
