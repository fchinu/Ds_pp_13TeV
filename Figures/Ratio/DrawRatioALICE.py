import ROOT 


TransparentBlue = ROOT.TColor.GetFreeColorIndex()
cTransparentBlue = ROOT.TColor(TransparentBlue, 79./255, 123./255, 161./255, 'TransparentBlue', 1.0)

# ALICE comparison

ratioFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/DsOverDplus_newMC.root'
systFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/SystUncertainties.root'
Alice5TevFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/RatioComparisonVsEnergy.root'
Alice7TevFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/RatioComparisonVsEnergy.root'
Alice13TevFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Ratios/pp_data/RatioComparisonVsEnergy.root'

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
hAlice7Tev.SetMarkerSize(1.5)
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
hAlice13Tev.SetMarkerSize(1.5)
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
ROOT.gPad.SetLeftMargin(0.1)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.1)

hFrame = canvas.DrawFrame(0, 0, 24, 1, ';#it{p}_{T} (GeV/c);D_{s}^{+}/D^{+}')
hFrame.GetXaxis().SetTitleOffset(1.2)
gAlice5Tev.Draw('5 same')
gAlice7Tev.Draw('5 same')
gAlice13Tev.Draw('5 same')
gSyst.Draw('5 same')
hAlice5Tev.Draw('p same')
hAlice7Tev.Draw('p same')
hAlice13Tev.Draw('p same')
hAlice5TevBorder.Draw('p same')
hAlice7TevBorder.Draw('p same')
hAlice13TevBorder.Draw('p same')
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

texty = ROOT.TLatex(0.155, 0.8, '|#it{y}|<0.5')
texty.SetNDC()
texty.SetTextFont(42)
texty.SetTextSize(0.04)
texty.Draw()

text = ROOT.TLatex(0.3, 0.35, 'Prompt D_{s}^{+}, D^{+}')
text.SetNDC()
text.SetTextFont(42)
text.SetTextSize(0.04)
text.Draw()

legend = ROOT.TLegend(0.3, 0.14, 0.5, 0.34)
legend.SetBorderSize(0)
legend.SetTextSize(0.027)
legend.SetFillStyle(0)
legend.AddEntry(hRatio, '#sqrt{s} = 13.6 TeV', 'lp')
legend.AddEntry(hAlice13Tev, '#sqrt{s} = 13 TeV, JHEP12(2023)086', 'lp')
legend.AddEntry(hAlice7Tev, '#sqrt{s} = 7 TeV, EPJC - s10052-019-6873-6', 'lp')
legend.AddEntry(hAlice5Tev, '#sqrt{s} = 5 TeV, JHEP05(2021)220', 'lp')
legend.Draw()

legendBorder = ROOT.TLegend(0.3, 0.14, 0.5, 0.34)
legendBorder.SetBorderSize(0)
legendBorder.SetTextSize(0.027)
legendBorder.SetFillStyle(0)
legendBorder.AddEntry(hRatioBorder, '#sqrt{s} = 13.6 TeV', 'lp')
legendBorder.AddEntry(hAlice13TevBorder, '#sqrt{s} = 13 TeV', 'lp')
legendBorder.AddEntry(hAlice7TevBorder, '#sqrt{s} = 7 TeV', 'lp')
legendBorder.AddEntry(hAlice5TevBorder, '#sqrt{s} = 5 TeV', 'lp')
legendBorder.Draw()



canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonALICE.png')
canvas.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/DsOverDplusComparisonALICE.pdf')
