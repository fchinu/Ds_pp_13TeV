import ROOT 


TransparentBlue = ROOT.TColor.GetFreeColorIndex()
cTransparentBlue = ROOT.TColor(TransparentBlue, 79./255, 123./255, 161./255, 'TransparentBlue', 1.0)

ratioFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/DsOverDplus_newMC.root'
systFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/SystUncertainties.root'
LHCbFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/ds_over_dplus_lhcb.root'
ALICEFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/v1Figure7.root'

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

ALICEFile = ROOT.TFile(ALICEFileName)
hAlice = ALICEFile.Get('histoStat')
hAlice.SetDirectory(0)
gAlice = ALICEFile.Get('syst')
ALICEFile.Close()

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

hAlice.SetLineColor(ROOT.kOrange-3)
hAlice.SetMarkerColor(ROOT.kOrange-3)
hAlice.SetMarkerStyle(ROOT.kFullDiamond)
hAlice.SetMarkerSize(1.5)
hAlice.SetLineWidth(2)

hAliceBorder = hAlice.Clone('hAliceBorder')
hAliceBorder.SetLineColor(ROOT.kOrange-3)
hAliceBorder.SetMarkerColor(ROOT.kBlack)
hAliceBorder.SetMarkerStyle(56)
hAliceBorder.SetLineWidth(3)

gAlice.SetLineColor(ROOT.kOrange-3)
gAlice.SetFillStyle(0)
gAlice.SetLineWidth(2)
for iPt in range(gAlice.GetN()):
    gAlice.SetPointEXhigh(iPt, 0.2)
    gAlice.SetPointEXlow(iPt, 0.2)

canvas = ROOT.TCanvas('canvas', 'canvas', 800, 800)
ROOT.gPad.SetLeftMargin(0.1)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.1)

hFrame = canvas.DrawFrame(0, 0, 24, 0.7, ';#it{p}_{T} (GeV/c);D_{s}^{+}/D^{+}')
hFrame.GetXaxis().SetTitleOffset(1.2)
gLHCb_20_25.Draw('5 same')
gLHCb_25_30.Draw('5 same')
gAlice.Draw('5 same')
gSyst.Draw('5 same')
hLHCb_20_25.Draw('p same')
hLHCb_25_30.Draw('p same')
hAlice.Draw('p same')
hLHCb_20_25Border.Draw('p same')
hLHCb_25_30Border.Draw('p same')
hAliceBorder.Draw('p same')
hRatio.Draw('p same')
hRatioBorder.Draw('p same')
canvas.Update()

legend = ROOT.TLegend(0.5, 0.2, 0.9, 0.4)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.AddEntry(hRatio, 'ALICE #font[122]{-} Run 3, |#it{y}|<0.5', 'lp')
legend.AddEntry(hLHCb_20_25, 'LHCb #font[122]{-} Run 2, 2 < #it{y} < 2.5', 'lp')
legend.AddEntry(hLHCb_25_30, 'LHCb #font[122]{-} Run 2, 2.5 < #it{y} < 3', 'lp')
legend.AddEntry(hAlice, 'ALICE #font[122]{-} Run 2, |#it{y}|<0.5', 'lp')
legend.Draw()

legendBorders = ROOT.TLegend(0.5, 0.2, 0.9, 0.4)
legendBorders.SetBorderSize(0)
legendBorders.SetFillStyle(0)
legendBorders.AddEntry(hRatioBorder, 'ALICE #font[122]{-} Run 3, |#it{y}|<0.5', 'lp')
legendBorders.AddEntry(hLHCb_20_25Border, 'LHCb #font[122]{-} Run 2, 2 < #it{y} < 2.5', 'lp')
legendBorders.AddEntry(hLHCb_25_30Border, 'LHCb #font[122]{-} Run 2, 2.5 < #it{y} < 3', 'lp')
legendBorders.AddEntry(hAliceBorder, 'ALICE #font[122]{-} Run 2, |#it{y}|<0.5', 'lp')
legendBorders.Draw()

canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparison.png')
