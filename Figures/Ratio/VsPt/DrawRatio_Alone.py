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

hFrame = canvas.DrawFrame(0, 0, 24, 1, ';#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});D_{s}^{+}/D^{+}')
hFrame.GetXaxis().SetTitleOffset(1.05)
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetXaxis().SetTitleSize(0.043)
hFrame.GetYaxis().SetTitleSize(0.043)
hFrame.GetXaxis().SetLabelSize(0.04)
hFrame.GetYaxis().SetLabelSize(0.04)
ROOT.gPad.SetTickx(1)
ROOT.gPad.SetTicky(1)
gSyst.Draw('5 same')
hRatio.Draw('p same')
hRatioBorder.Draw('p same')
canvas.Update()

textALICE = ROOT.TLatex(0.15, 0.88, 'This Thesis')
textALICE.SetNDC()
textALICE.SetTextFont(42)
textALICE.SetTextSize(0.05)
textALICE.Draw()

textpp = ROOT.TLatex(0.155, 0.84, 'pp collisions')
textpp.SetNDC()
textpp.SetTextFont(42)
textpp.SetTextSize(0.04)
textpp.Draw()

texty = ROOT.TLatex(0.155, 0.79, '|#it{y}|#kern[0.4]{<}#kern[0.2]{0.5}')
texty.SetNDC()
texty.SetTextFont(42)
texty.SetTextSize(0.04)
texty.Draw()

textDsDecay = ROOT.TLatex(0.155, 0.74, 'D_{s}^{+}, D^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+} ')
textDsDecay.SetNDC()
textDsDecay.SetTextFont(42)
textDsDecay.SetTextSize(0.035)
textDsDecay.Draw()

textDsDecayConj = ROOT.TLatex(0.155, 0.7, 'and charge conj.')
textDsDecayConj.SetNDC()
textDsDecayConj.SetTextFont(42)
textDsDecayConj.SetTextSize(0.035)
textDsDecayConj.Draw()


textBR = ROOT.TLatex(0.6, 0.8, 'BR unc. (not shown): ^{#lower[5]{+3.7}}_{#lower[-0.2]{#font[122]{-}4.0}}#it{%} ')
textBR.SetNDC()
textBR.SetTextFont(42)
textBR.SetTextSize(0.03)
textBR.Draw()

ROOT.gPad.RedrawAxis()

canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplus.png')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplus.pdf')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplus.eps')
