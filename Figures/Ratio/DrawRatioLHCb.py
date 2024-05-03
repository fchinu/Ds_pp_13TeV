import ROOT 


TransparentBlue = ROOT.TColor.GetFreeColorIndex()
cTransparentBlue = ROOT.TColor(TransparentBlue, 79./255, 123./255, 161./255, 'TransparentBlue', 1.0)

# LHCb comparison

ratioFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/DsOverDplus_newMC.root'
systFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/SystUncertainties.root'
LHCbFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/ds_over_dplus_lhcb.root'

ratioFile = ROOT.TFile(ratioFileName)
hRatio = ratioFile.Get('hRatio')
hRatio.SetDirectory(0)
ratioFile.Close()

systFile = ROOT.TFile(systFileName)
hSystRel = systFile.Get('hSystRel')
hSystRel.SetDirectory(0)
systFile.Close()

gSyst = ROOT.TGraphErrors(hSystRel)

for iPt in range(gSyst.GetN()):
    gSyst.SetPointError(iPt, 0.2, gSyst.GetPointY(iPt)*hRatio.GetBinContent(iPt+1))
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

hFrame = canvas.DrawFrame(0, 0, 24, 0.7, ';#it{p}_{T} (GeV/c);D_{s}^{+}/D^{+}')
hFrame.GetXaxis().SetTitleOffset(1.2)
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

textALICE = ROOT.TLatex(0.15, 0.9, 'This Work')
textALICE.SetNDC()
textALICE.SetTextFont(42)
textALICE.SetTextSize(0.05)
textALICE.Draw()

textpp = ROOT.TLatex(0.155, 0.85, 'pp collisions')
textpp.SetNDC()
textpp.SetTextFont(42)
textpp.SetTextSize(0.04)
textpp.Draw()

text = ROOT.TLatex(0.2, 0.4, 'Prompt D_{s}^{+}, D^{+}')
text.SetNDC()
text.SetTextFont(42)
text.SetTextSize(0.04)
text.Draw()

LegAliceText = ROOT.TLatex(0.2, 0.35, 'ALICE #font[122]{-} Run 3, #sqrt{s} = 13.6 TeV')
LegAliceText.SetNDC()
LegAliceText.SetTextFont(42)
LegAliceText.SetTextSize(0.035)
LegAliceText.Draw()

LegLHCbText = ROOT.TLatex(0.2, 0.25, '#splitline{LHCb #font[122]{-} Run 2, #sqrt{s} = 13 TeV}{JHEP03(2016)159}')
LegLHCbText.SetNDC()
LegLHCbText.SetTextFont(42)
LegLHCbText.SetTextSize(0.035)
LegLHCbText.Draw()

legAlice = ROOT.TLegend(0.65, 0.335, 0.9, 0.385)
legAlice.SetBorderSize(0)
legAlice.SetTextSize(0.027)
legAlice.SetFillStyle(0)
legAlice.AddEntry(hRatio, '|#it{y}|<0.5', 'lp')
legAlice.Draw()

legAliceBorder = ROOT.TLegend(0.65, 0.335, 0.9, 0.385)
legAliceBorder.SetBorderSize(0)
legAliceBorder.SetTextSize(0.027)
legAliceBorder.SetFillStyle(0)
legAliceBorder.AddEntry(hRatioBorder, '|#it{y}|<0.5', 'lp')
legAliceBorder.Draw()

legLHCb = ROOT.TLegend(0.65, 0.22, 0.9, 0.32)
legLHCb.SetBorderSize(0)
legLHCb.SetTextSize(0.027)
legLHCb.SetFillStyle(0)
legLHCb.AddEntry(hLHCb_20_25, '2 < #it{y} < 2.5', 'lp')
legLHCb.AddEntry(hLHCb_25_30, '2.5 < #it{y} < 3', 'lp')
legLHCb.Draw()

legLHCbBorder = ROOT.TLegend(0.65, 0.22, 0.9, 0.32)
legLHCbBorder.SetBorderSize(0)
legLHCbBorder.SetTextSize(0.027)
legLHCbBorder.SetFillStyle(0)
legLHCbBorder.AddEntry(hLHCb_20_25Border, '2 < #it{y} < 2.5', 'lp')
legLHCbBorder.AddEntry(hLHCb_25_30Border, '2.5 < #it{y} < 3', 'lp')
legLHCbBorder.Draw()

canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonLHCb.png')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonLHCb.pdf')
