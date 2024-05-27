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
Alice5TevFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/RatioComparisonVsEnergy.root'
Alice7TevFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/RatioComparisonVsEnergy.root'
Alice13TevFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/RatioComparisonVsEnergy.root'

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


Alice5TevFile = ROOT.TFile(Alice5TevFileName)
hAlice5Tev = Alice5TevFile.Get('hDsOverDplus_5TeV')
hAlice5Tev.SetDirectory(0)
gAlice5Tev = Alice5TevFile.Get('gDsOverDplus_5TeV')
Alice5TevFile.Close()

Alice7TevFile = ROOT.TFile(Alice7TevFileName)
hAlice7Tev = Alice7TevFile.Get('hDsOverDplus_7TeV')
hAlice7Tev.SetDirectory(0)
gAlice7Tev = Alice7TevFile.Get('gDsOverDplus_7TeV')
Alice7TevFile.Close()

Alice13TevFile = ROOT.TFile(Alice13TevFileName)
hAlice13Tev = Alice13TevFile.Get('hDsOverDplus')
hAlice13Tev.SetDirectory(0)
gAlice13Tev = Alice13TevFile.Get('gDsOverDplus')
Alice13TevFile.Close()

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

hAlice5Tev.SetLineColor(ROOT.kAzure+7)
hAlice5Tev.SetMarkerColor(ROOT.kAzure+7)
hAlice5Tev.SetMarkerStyle(ROOT.kFullSquare)
hAlice5Tev.SetMarkerSize(1.5)
hAlice5Tev.SetLineWidth(2)

hAlice5TevBorder = hAlice5Tev.Clone('hAlice5TevBorder')
hAlice5TevBorder.SetLineColor(ROOT.kAzure+7)
hAlice5TevBorder.SetMarkerColor(ROOT.kBlack)
hAlice5TevBorder.SetMarkerStyle(54)
hAlice5TevBorder.SetLineWidth(3)

gAlice5Tev.SetLineColor(ROOT.kAzure+7)
gAlice5Tev.SetFillStyle(0)
gAlice5Tev.SetLineWidth(2)
for iPt in range(gAlice5Tev.GetN()):
    gAlice5Tev.SetPointEXhigh(iPt, 0.2)
    gAlice5Tev.SetPointEXlow(iPt, 0.2)

hAlice7Tev.SetLineColor(ROOT.kGreen-2)
hAlice7Tev.SetMarkerColor(ROOT.kGreen-2)
hAlice7Tev.SetMarkerStyle(ROOT.kFullCross)
hAlice7Tev.SetMarkerSize(1.7)
hAlice7Tev.SetLineWidth(2)

hAlice7TevBorder = hAlice7Tev.Clone('hAlice7TevBorder')
hAlice7TevBorder.SetLineColor(ROOT.kGreen-2)
hAlice7TevBorder.SetMarkerColor(ROOT.kBlack)
hAlice7TevBorder.SetMarkerStyle(57)
hAlice7TevBorder.SetLineWidth(3)

gAlice7Tev.SetLineColor(ROOT.kGreen-2)
gAlice7Tev.SetFillStyle(0)
gAlice7Tev.SetLineWidth(2)
for iPt in range(gAlice7Tev.GetN()):
    gAlice7Tev.SetPointEXhigh(iPt, 0.2)
    gAlice7Tev.SetPointEXlow(iPt, 0.2)

hAlice13Tev.SetLineColor(ROOT.kOrange-3)
hAlice13Tev.SetMarkerColor(ROOT.kOrange-3)
hAlice13Tev.SetMarkerStyle(ROOT.kFullDiamond)
hAlice13Tev.SetMarkerSize(2)
hAlice13Tev.SetLineWidth(2)

hAlice13TevBorder = hAlice13Tev.Clone('hAlice13TevBorder')
hAlice13TevBorder.SetLineColor(ROOT.kOrange-3)
hAlice13TevBorder.SetMarkerColor(ROOT.kBlack)
hAlice13TevBorder.SetMarkerStyle(56)
hAlice13TevBorder.SetLineWidth(3)

gAlice13Tev.SetLineColor(ROOT.kOrange-3)
gAlice13Tev.SetFillStyle(0)
gAlice13Tev.SetLineWidth(2)
for iPt in range(gAlice13Tev.GetN()):
    gAlice13Tev.SetPointEXhigh(iPt, 0.2)
    gAlice13Tev.SetPointEXlow(iPt, 0.2)

canvas = ROOT.TCanvas('canvas', 'canvas', 800, 800)
ROOT.gPad.SetLeftMargin(0.12)
ROOT.gPad.SetRightMargin(0.04)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.1)

hFrame = canvas.DrawFrame(0, 0, 36, 1.05, ';#it{p}_{T} (GeV/#it{c});D_{s}^{+}/D^{+}')
hFrame.GetXaxis().SetTitleOffset(1.05)
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetXaxis().SetTitleSize(0.043)
hFrame.GetYaxis().SetTitleSize(0.043)
hFrame.GetXaxis().SetLabelSize(0.04)
hFrame.GetYaxis().SetLabelSize(0.04)
ROOT.gPad.SetTickx(1)
ROOT.gPad.SetTicky(1)
gAlice5Tev.Draw('5 same')
gAlice7Tev.Draw('5 same')
gAlice13Tev.Draw('5 same')
gSyst.Draw('5 same')
hAlice5Tev.Draw('p same')
hAlice5TevBorder.Draw('p same')
hAlice7Tev.Draw('p same')
hAlice7TevBorder.Draw('p same')
hAlice13Tev.Draw('p same')
hAlice13TevBorder.Draw('p same')
hRatio.Draw('p same')
hRatioBorder.Draw('p same')
canvas.Update()

textALICE = ROOT.TLatex(0.15, 0.88, 'ALICE Preliminary')
textALICE.SetNDC()
textALICE.SetTextFont(42)
textALICE.SetTextSize(0.05)
textALICE.Draw()

textpp = ROOT.TLatex(0.6, 0.88, 'pp collisions')
textpp.SetNDC()
textpp.SetTextFont(42)
textpp.SetTextSize(0.04)
textpp.Draw()

texty = ROOT.TLatex(0.6, 0.84, '|#it{y}|<0.5')
texty.SetNDC()
texty.SetTextFont(42)
texty.SetTextSize(0.04)
texty.Draw()

textDsDecay = ROOT.TLatex(0.155, 0.84, 'D_{s}^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
textDsDecay.SetNDC()
textDsDecay.SetTextFont(42)
textDsDecay.SetTextSize(0.035)
textDsDecay.Draw()

textDsDecayConj = ROOT.TLatex(0.155, 0.8, 'and charge conj.')
textDsDecayConj.SetNDC()
textDsDecayConj.SetTextFont(42)
textDsDecayConj.SetTextSize(0.035)
textDsDecayConj.Draw()

textDecay = ROOT.TLatex(0.155, 0.75, 'D^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
textDecay.SetNDC()
textDecay.SetTextFont(42)
textDecay.SetTextSize(0.035)
textDecay.Draw()

textDecayConj = ROOT.TLatex(0.155, 0.71, 'and charge conj.')
textDecayConj.SetNDC()
textDecayConj.SetTextFont(42)
textDecayConj.SetTextSize(0.035)
textDecayConj.Draw()

textBR = ROOT.TLatex(0.155, 0.67, 'BR unc. (not shown): ^{#lower[5]{+3.7}}_{#lower[-0.2]{#font[122]{-}4.0}}#it{%} ')
textBR.SetNDC()
textBR.SetTextFont(42)
textBR.SetTextSize(0.03)
textBR.Draw()

Thislegend = ROOT.TLegend(0.155, 0.61, 0.355, 0.65)
Thislegend.SetBorderSize(0)
Thislegend.SetTextSize(0.027)
Thislegend.SetFillStyle(0)
Thislegend.AddEntry(hRatio, '#sqrt{s} = 13.6 Te#kern[-0.1]{V}', 'pl')
Thislegend.Draw()

ThislegendBorder = ROOT.TLegend(0.155, 0.61, 0.355, 0.65)
ThislegendBorder.SetBorderSize(0)
ThislegendBorder.SetTextSize(0.027)
ThislegendBorder.SetFillStyle(0)
ThislegendBorder.AddEntry(hRatioBorder, '#sqrt{s} = 13.6 Te#kern[-0.1]{V}', 'pl')
ThislegendBorder.Draw()

text = ROOT.TLatex(0.3, 0.35, 'D^{+}#rightarrow #pi^{+}K^{#font[122]{-}}#pi^{+}')
text.SetNDC()
text.SetTextFont(42)
text.SetTextSize(0.035)
text.Draw()

textConj = ROOT.TLatex(0.3, 0.31, 'and charge conj.')
textConj.SetNDC()
textConj.SetTextFont(42)
textConj.SetTextSize(0.035)
textConj.Draw()

textBRGoldenDplus = ROOT.TLatex(0.3, 0.27, 'BR unc. (not shown): 3.2#it{%}')
textBRGoldenDplus.SetNDC()
textBRGoldenDplus.SetTextFont(42)
textBRGoldenDplus.SetTextSize(0.03)
textBRGoldenDplus.Draw()

legend = ROOT.TLegend(0.3, 0.14, 0.5, 0.26)
legend.SetBorderSize(0)
legend.SetTextSize(0.027)
legend.SetFillStyle(0)
legend.AddEntry(hAlice13Tev, '#sqrt{s} = 13 Te#kern[-0.1]{V}, JHEP12(2023)086', 'pl')
legend.AddEntry(hAlice7Tev, '#sqrt{s} = 7 Te#kern[-0.1]{V}, EPJC - s10052-019-6873-6', 'pl')
legend.AddEntry(hAlice5Tev, '#sqrt{s} = 5 Te#kern[-0.1]{V}, JHEP05(2021)220', 'pl')
legend.Draw()

legendBorder = ROOT.TLegend(0.3, 0.14, 0.5, 0.26)
legendBorder.SetBorderSize(0)
legendBorder.SetTextSize(0.027)
legendBorder.SetFillStyle(0)
legendBorder.AddEntry(hAlice13TevBorder, '#sqrt{s} = 13 Te#kern[-0.1]{V}', 'pl')
legendBorder.AddEntry(hAlice7TevBorder, '#sqrt{s} = 7 Te#kern[-0.1]{V}', 'pl')
legendBorder.AddEntry(hAlice5TevBorder, '#sqrt{s} = 5 Te#kern[-0.1]{V}', 'pl')
legendBorder.Draw()



canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonALICE_andcc.png')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonALICE_andcc.pdf')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonALICE_andcc.eps')
