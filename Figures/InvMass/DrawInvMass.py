import ROOT 
import numpy as np

def get_function_fromhist(hist):
    """

    """
    spline = ROOT.TSpline3(hist)
    spline.SetLineWidth(2)

    return spline


def SetStyle(graypalette) :
  
  ROOT.gStyle.Reset("Plain")
  ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetOptStat(0)
  if(graypalette):
    ROOT.gStyle.SetPalette(8,0)
  else:
    ROOT.gStyle.SetPalette(1)
  ROOT.gStyle.SetCanvasColor(10)
  ROOT.gStyle.SetCanvasBorderMode(0)
  ROOT.gStyle.SetFrameLineWidth(1)
  ROOT.gStyle.SetFrameFillColor(ROOT.kWhite)
  ROOT.gStyle.SetPadColor(10)
  ROOT.gStyle.SetPadTickX(1)
  ROOT.gStyle.SetPadTickY(1)
  ROOT.gStyle.SetPadBottomMargin(0.14)
  ROOT.gStyle.SetPadLeftMargin(0.13)
  ROOT.gStyle.SetPadRightMargin(0.035)
  ROOT.gStyle.SetPadTopMargin(0.035)
  ROOT.gStyle.SetHistLineWidth(1)
  ROOT.gStyle.SetHistLineColor(ROOT.kRed)
  ROOT.gStyle.SetFuncWidth(2)
  ROOT.gStyle.SetFuncColor(ROOT.kGreen)
  ROOT.gStyle.SetLineWidth(2)
  ROOT.gStyle.SetLabelSize(0.045,"xyz")
  ROOT.gStyle.SetLabelOffset(0.01,"y")
  ROOT.gStyle.SetLabelOffset(0.01,"x")
  ROOT.gStyle.SetLabelColor(ROOT.kBlack,"xyz")
  ROOT.gStyle.SetTitleSize(0.05,"xyz")
  ROOT.gStyle.SetTitleOffset(1.2,"y")
  ROOT.gStyle.SetTitleOffset(1.2,"x")
  ROOT.gStyle.SetTitleFillColor(ROOT.kWhite)
  ROOT.gStyle.SetTextSizePixels(26)
  ROOT.gStyle.SetTextFont(42)
  #  ROOT.gStyle.SetTickLength(0.04,"X")  ROOT.gStyle.SetTickLength(0.04,"Y") 

  ROOT.gStyle.SetLegendBorderSize(0)
  ROOT.gStyle.SetLegendFillColor(ROOT.kWhite)
  #  ROOT.gStyle.SetFillColor(kWhite)
  ROOT.gStyle.SetLegendFont(42)


cFillGreen = ROOT.TColor.GetFreeColorIndex()
FillGreen = ROOT.TColor(cFillGreen, 164/255., 223/255., 181/255., 'kFillGreen', 1.0)

cFillBlue = ROOT.TColor.GetFreeColorIndex()
FillBlue = ROOT.TColor(cFillBlue, 146/255., 190/255., 198/255., 'kFillBlue', 1.0)

# Open the file
InputFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/Data/fits/fits.root")
hData2_2p5 = InputFile.Get("hdata2.0_2.5")
hData2_2p5.SetDirectory(0)
hTotal_func2_2p5 = InputFile.Get("total_func2.0_2.5")
hTotal_func2_2p5.SetDirectory(0)
hBkg0_2_2p5 = InputFile.Get("bkg_02.0_2.5")
hBkg0_2_2p5.SetDirectory(0)
hBkg1_2_2p5 = InputFile.Get("bkg_12.0_2.5")
hBkg1_2_2p5.SetDirectory(0)
hSignal0_2_2p5 = InputFile.Get("signal_02.0_2.5")
hSignal0_2_2p5.SetDirectory(0)
hSignal1_2_2p5 = InputFile.Get("signal_12.0_2.5")
hSignal1_2_2p5.SetDirectory(0)
InputFile.Close()

centralFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/Data/RawYields_Data.root")
hRawYields = centralFile.Get("hRawYields")
hRawYields.SetDirectory(0)
hRawYieldsSecondPeak = centralFile.Get("hRawYieldsSecondPeak")
hRawYieldsSecondPeak.SetDirectory(0)
hRawYieldsMean = centralFile.Get("hRawYieldsMean")
hRawYieldsMean.SetDirectory(0)
hRawYieldsMeanSecondPeak = centralFile.Get("hRawYieldsMeanSecondPeak")
hRawYieldsMeanSecondPeak.SetDirectory(0)
hRawYieldsSigma = centralFile.Get("hRawYieldsSigma")
hRawYieldsSigma.SetDirectory(0)
hRawYieldsSigmaSecondPeak = centralFile.Get("hRawYieldsSigmaSecondPeak")
hRawYieldsSigmaSecondPeak.SetDirectory(0)
hSignificance = centralFile.Get("hRawYieldsSignificance")
hSignificance.SetDirectory(0)
hSignificanceSecondPeak = centralFile.Get("hRawYieldsSignificanceSecondPeak")
hSignificanceSecondPeak.SetDirectory(0)
hBkg = centralFile.Get("hRawYieldsBkg")
hBkg.SetDirectory(0)
hBkgSecondPeak = centralFile.Get("hRawYieldsBkgSecondPeak")
hBkgSecondPeak.SetDirectory(0)
centralFile.Close()

# Set style
SetStyle(0)

hData2_2p5.SetMarkerStyle(ROOT.kFullCircle)
hData2_2p5.SetMarkerSize(0.9)
hData2_2p5.SetMarkerColor(ROOT.kBlack)
hData2_2p5.SetLineColor(ROOT.kBlack)
hData2_2p5.SetLineWidth(2)

hTotal_func2_2p5 = get_function_fromhist(hTotal_func2_2p5)
hTotal_func2_2p5.SetLineColor(ROOT.kAzure-3)
hTotal_func2_2p5.SetLineWidth(3)

hBkg0_2_2p5 = get_function_fromhist(hBkg0_2_2p5)
hBkg0_2_2p5.SetLineColor(ROOT.kViolet-4)
hBkg0_2_2p5.SetLineWidth(2)
hBkg0_2_2p5.SetLineStyle(2)

hBkg1_2_2p5 = get_function_fromhist(hBkg1_2_2p5)
hBkg1_2_2p5.SetLineColor(ROOT.kRed)
hBkg1_2_2p5.SetLineWidth(2)

hSignal0_2_2p5.SetLineColor(ROOT.kBlack)
hSignal0_2_2p5.SetLineWidth(2)
hSignal0_2_2p5.SetFillColor(cFillBlue)

hSignal1_2_2p5.SetLineColor(ROOT.kBlack)
hSignal1_2_2p5.SetLineWidth(2)
hSignal1_2_2p5.SetFillColor(cFillGreen)

# Add text

ALICEtext = ROOT.TLatex(0.19, 0.87, "ALICE Preliminary")
ALICEtext.SetNDC()
ALICEtext.SetTextSize(0.06)
ALICEtext.SetTextFont(42)

ppEnergytext = ROOT.TLatex(0.19, 0.81, "pp, #sqrt{s} = 13.6 TeV")
ppEnergytext.SetNDC()
ppEnergytext.SetTextSize(0.05)
ppEnergytext.SetTextFont(42)

textDsDecay = ROOT.TLatex(0.19, 0.73, "D_{s}^{+}, D^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+} ")
textDsDecay.SetNDC()
textDsDecay.SetTextSize(0.05)
textDsDecay.SetTextFont(42)

textDsDecayConj = ROOT.TLatex(0.19, 0.68, 'and charge conj.')
textDsDecayConj.SetNDC()
textDsDecayConj.SetTextFont(42)
textDsDecayConj.SetTextSize(0.05)

textPt = ROOT.TLatex(0.19, 0.605, '2.0 < #it{p}_{T} < 2.5 GeV/#it{c}')
textPt.SetNDC()
textPt.SetTextFont(42)
textPt.SetTextSize(0.05)

canvas_2_2p5 = ROOT.TCanvas("canvas_2_2p5", "canvas_2_2p5", 700, 800)

padFit = ROOT.TPad("padFit", "padFit", 0.0, 0.3, 1.0, 1.0)
padFit.SetRightMargin(0.05)
padFit.SetTopMargin(0.05)
padFit.SetLeftMargin(0.14)
padFit.SetBottomMargin(0.)
padFit.Draw()

hFrame = padFit.cd().DrawFrame(1.75, 0.1, 2.1, 3000, ";#it{M}(KK#pi) (Ge#kern[-0.05]{V}/#it{c}^{2});Counts per 2.0 Me#kern[-0.05]{V}/#it{c}^{2}")
hFrame.GetXaxis().SetTitleOffset(1.)
hFrame.GetYaxis().SetTitleOffset(1.4)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetLabelSize(0.045)
hFrame.GetYaxis().SetLabelSize(0.045)

legend = ROOT.TLegend(0.585, 0.6, 0.885, 0.92)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.043)
legend.AddEntry(hData2_2p5, "Data", "pe")
legend.AddEntry(hSignal0_2_2p5, "D_{s}^{+} signal", "f")
legend.AddEntry(hSignal1_2_2p5, "D^{+} signal", "f")
legend.AddEntry(hBkg1_2_2p5, "Combinatorial bkg.", "l")
legend.AddEntry(hBkg0_2_2p5, "D^{+}#rightarrow#pi^{+}K^{#font[122]{-}}#pi^{+} bkg.", "l")
legend.AddEntry(hTotal_func2_2p5, "Total fit function", "l")


hSignal0_2_2p5.Draw("same")
hSignal1_2_2p5.Draw("same")
hBkg0_2_2p5.Draw("same")
hBkg1_2_2p5.Draw("same")
hTotal_func2_2p5.Draw("same")
hData2_2p5.Draw("pe same")
ALICEtext.Draw("same")
ppEnergytext.Draw("same")
textPt.Draw("same")
textDsDecay.Draw("same")
textDsDecayConj.Draw("same")
legend.Draw("same")

ROOT.gPad.RedrawAxis()

canvas_2_2p5.cd()

padRes = ROOT.TPad("padRes", "padRes", 0.0, 0.0, 1.0, 0.3)
padRes.SetRightMargin(0.05)
padRes.SetTopMargin(0.)
padRes.SetLeftMargin(0.14)
padRes.SetBottomMargin(0.3)
padRes.Draw()

hFrameRes = padRes.cd().DrawFrame(1.75, -300., 2.1, 1750, ";#it{M}(KK#pi) (Ge#kern[-0.05]{V}/#it{c}^{2});Residuals")
hFrameRes.GetXaxis().SetTitleOffset(1)
hFrameRes.GetYaxis().SetTitleOffset(0.57)
hFrameRes.GetXaxis().SetTitleSize(0.125)
hFrameRes.GetYaxis().SetTitleSize(0.125)
hFrameRes.GetXaxis().SetLabelSize(0.10)
hFrameRes.GetYaxis().SetLabelSize(0.10)
hFrameRes.GetYaxis().SetLabelOffset(0.01)
hFrameRes.GetYaxis().SetNdivisions(505)


hRes = hData2_2p5.Clone("hRes")

for iPt in range(hData2_2p5.GetNbinsX()):
  hRes.SetBinContent(iPt+1, hData2_2p5.GetBinContent(iPt+1) - hBkg0_2_2p5.Eval(hData2_2p5.GetBinCenter(iPt+1)) - hBkg1_2_2p5.Eval(hData2_2p5.GetBinCenter(iPt+1)))

# Set Style
hRes.SetMarkerStyle(ROOT.kFullCircle)
hRes.SetMarkerSize(0.9)
hRes.SetMarkerColor(ROOT.kBlack)
hRes.SetLineColor(ROOT.kBlack)
hRes.SetLineWidth(2)

DsmesonText = ROOT.TLatex(0.68, 0.87, "D_{s}^{+} meson")
DsmesonText.SetNDC()
DsmesonText.SetTextSize(0.1)
DsMeanText = ROOT.TLatex(0.68, 0.77, f"#mu = ({hRawYieldsMean.GetBinContent(hRawYieldsMean.FindBin(2.01))*1000:.1f} #pm {hRawYieldsMean.GetBinError(hRawYieldsMean.FindBin(2.01))*1000:.1f}) Me#kern[-0.05]{{V}}/#it{{c}}^{{2}}")
DsMeanText.SetNDC()
DsMeanText.SetTextSize(0.07)
DsSigmaText = ROOT.TLatex(0.68, 0.68, f"#sigma = ({hRawYieldsSigma.GetBinContent(hRawYieldsSigma.FindBin(2.01))*1000:.1f} #pm {hRawYieldsSigma.GetBinError(hRawYieldsSigma.FindBin(2.01))*1000:.1f}) Me#kern[-0.05]{{V}}/#it{{c}}^{{2}}")
DsSigmaText.SetNDC()
DsSigmaText.SetTextSize(0.07)
DsSignalText = ROOT.TLatex(0.68, 0.59, f"#it{{S}} = {hRawYields.GetBinContent(hRawYields.FindBin(2.01)):.0f} #pm {hRawYields.GetBinError(hRawYields.FindBin(2.01)):.0f}")
DsSignalText.SetNDC()
DsSignalText.SetTextSize(0.07)

DplusmesonText = ROOT.TLatex(0.17, 0.87, "D^{+} meson")
DplusmesonText.SetNDC()
DplusmesonText.SetTextSize(0.1)
DplusMeanText = ROOT.TLatex(0.17, 0.77, f"#mu = ({hRawYieldsMeanSecondPeak.GetBinContent(hRawYieldsMeanSecondPeak.FindBin(2.01))*1000:.1f} #pm {hRawYieldsMeanSecondPeak.GetBinError(hRawYieldsMeanSecondPeak.FindBin(2.01))*1000:.1f}) Me#kern[-0.05]{{V}}/#it{{c}}^{{2}}")
DplusMeanText.SetNDC()
DplusMeanText.SetTextSize(0.07)
DplusSigmaText = ROOT.TLatex(0.17, 0.68, f"#sigma = ({hRawYieldsSigmaSecondPeak.GetBinContent(hRawYieldsSigmaSecondPeak.FindBin(2.01))*1000:.1f} #pm {hRawYieldsSigmaSecondPeak.GetBinError(hRawYieldsSigmaSecondPeak.FindBin(2.01))*1000:.1f}) Me#kern[-0.05]{{V}}/#it{{c}}^{{2}}")
DplusSigmaText.SetNDC()
DplusSigmaText.SetTextSize(0.07)
DplusSignalText = ROOT.TLatex(0.17, 0.59, f"#it{{S}} = {hRawYieldsSecondPeak.GetBinContent(hRawYieldsSecondPeak.FindBin(2.01)):.0f} #pm {hRawYieldsSecondPeak.GetBinError(hRawYieldsSecondPeak.FindBin(2.01)):.0f}")
DplusSignalText.SetNDC()
DplusSignalText.SetTextSize(0.07)

hSignal0_2_2p5.Draw("same")
hSignal1_2_2p5.Draw("same")
hRes.Draw("pe same")
DsmesonText.Draw("same")
DsMeanText.Draw("same")
DsSigmaText.Draw("same")
DsSignalText.Draw("same")

DplusmesonText.Draw("same")
DplusMeanText.Draw("same")
DplusSigmaText.Draw("same")
DplusSignalText.Draw("same")


canvas_2_2p5.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/InvMassFit_2_2p5.png")
canvas_2_2p5.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/InvMassFit_2_2p5.pdf")
canvas_2_2p5.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/InvMassFit_2_2p5.eps")









# Open the file
InputFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/Data/fits/fits.root")
hData5_5p5 = InputFile.Get("hdata5.0_5.5")
hData5_5p5.SetDirectory(0)
hTotal_func5_5p5 = InputFile.Get("total_func5.0_5.5")
hTotal_func5_5p5.SetDirectory(0)
hBkg0_5_5p5 = InputFile.Get("bkg_05.0_5.5")
hBkg0_5_5p5.SetDirectory(0)
hBkg1_5_5p5 = InputFile.Get("bkg_15.0_5.5")
hBkg1_5_5p5.SetDirectory(0)
hSignal0_5_5p5 = InputFile.Get("signal_05.0_5.5")
hSignal0_5_5p5.SetDirectory(0)
hSignal1_5_5p5 = InputFile.Get("signal_15.0_5.5")
hSignal1_5_5p5.SetDirectory(0)
InputFile.Close()

# Set style

hData5_5p5.SetMarkerStyle(ROOT.kFullCircle)
hData5_5p5.SetMarkerSize(0.9)
hData5_5p5.SetMarkerColor(ROOT.kBlack)
hData5_5p5.SetLineColor(ROOT.kBlack)
hData5_5p5.SetLineWidth(2)

hTotal_func5_5p5 = get_function_fromhist(hTotal_func5_5p5)
hTotal_func5_5p5.SetLineColor(ROOT.kAzure-3)
hTotal_func5_5p5.SetLineWidth(3)

hBkg0_5_5p5 = get_function_fromhist(hBkg0_5_5p5)
hBkg0_5_5p5.SetLineColor(ROOT.kViolet-4)
hBkg0_5_5p5.SetLineWidth(2)
hBkg0_5_5p5.SetLineStyle(2)

hBkg1_5_5p5 = get_function_fromhist(hBkg1_5_5p5)
hBkg1_5_5p5.SetLineColor(ROOT.kRed)
hBkg1_5_5p5.SetLineWidth(2)

hSignal0_5_5p5.SetLineColor(ROOT.kBlack)
hSignal0_5_5p5.SetLineWidth(2)
hSignal0_5_5p5.SetFillColor(cFillBlue)

hSignal1_5_5p5.SetLineColor(ROOT.kBlack)
hSignal1_5_5p5.SetLineWidth(2)
hSignal1_5_5p5.SetFillColor(cFillGreen)

# Add text

ALICEtext = ROOT.TLatex(0.19, 0.87, "ALICE Preliminary")
ALICEtext.SetNDC()
ALICEtext.SetTextSize(0.06)
ALICEtext.SetTextFont(42)

ppEnergytext = ROOT.TLatex(0.19, 0.81, "pp, #sqrt{s} = 13.6 Te#kern[-0.05]{V}")
ppEnergytext.SetNDC()
ppEnergytext.SetTextSize(0.05)
ppEnergytext.SetTextFont(42)

textDsDecay = ROOT.TLatex(0.19, 0.73, "D_{s}^{+}, D^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+} ")
textDsDecay.SetNDC()
textDsDecay.SetTextSize(0.05)
textDsDecay.SetTextFont(42)

textDsDecayConj = ROOT.TLatex(0.19, 0.68, 'and charge conj.')
textDsDecayConj.SetNDC()
textDsDecayConj.SetTextFont(42)
textDsDecayConj.SetTextSize(0.05)

textPt = ROOT.TLatex(0.19, 0.605, '5.0 < #it{p}_{T} < 5.5 Ge#kern[-0.05]{V}/#it{c}')
textPt.SetNDC()
textPt.SetTextFont(42)
textPt.SetTextSize(0.05)

canvas_5_5p5 = ROOT.TCanvas("canvas_5_5p5", "canvas_5_5p5", 700, 800)

padFit = ROOT.TPad("padFit", "padFit", 0.0, 0.3, 1.0, 1.0)
padFit.SetRightMargin(0.05)
padFit.SetTopMargin(0.05)
padFit.SetLeftMargin(0.14)
padFit.SetBottomMargin(0.)
padFit.Draw()

hFrame = padFit.cd().DrawFrame(1.75, 0.1, 2.15, 1250, ";#it{M}(KK#pi) (Ge#kern[-0.05]{V}/#it{c}^{2});Counts per 2.0 Me#kern[-0.05]{V}/#it{c}^{2}")
hFrame.GetXaxis().SetTitleOffset(1.)
hFrame.GetYaxis().SetTitleOffset(1.4)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetLabelSize(0.045)
hFrame.GetYaxis().SetLabelSize(0.045)

legend = ROOT.TLegend(0.585, 0.6, 0.885, 0.92)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.043)
legend.AddEntry(hData5_5p5, "Data", "pe")
legend.AddEntry(hSignal0_5_5p5, "D_{s}^{+} signal", "f")
legend.AddEntry(hSignal1_5_5p5, "D^{+} signal", "f")
legend.AddEntry(hBkg1_5_5p5, "Combinatorial bkg.", "l")
legend.AddEntry(hBkg0_5_5p5, "D^{+}#rightarrow#pi^{+}K^{#font[122]{-}}#pi^{+} bkg.", "l")
legend.AddEntry(hTotal_func5_5p5, "Total fit function", "l")


hSignal0_5_5p5.Draw("same")
hSignal1_5_5p5.Draw("same")
hBkg0_5_5p5.Draw("same")
hBkg1_5_5p5.Draw("same")
hTotal_func5_5p5.Draw("same")
hData5_5p5.Draw("pe same")
ALICEtext.Draw("same")
ppEnergytext.Draw("same")
textPt.Draw("same")
textDsDecay.Draw("same")
textDsDecayConj.Draw("same")
legend.Draw("same")

ROOT.gPad.RedrawAxis()

canvas_5_5p5.cd()

padRes = ROOT.TPad("padRes", "padRes", 0.0, 0.0, 1.0, 0.3)
padRes.SetRightMargin(0.05)
padRes.SetTopMargin(0.)
padRes.SetLeftMargin(0.14)
padRes.SetBottomMargin(0.3)
padRes.Draw()

hFrameRes = padRes.cd().DrawFrame(1.75, -150., 2.15, 650, ";#it{M}(KK#pi) (Ge#kern[-0.05]{V}/#it{c}^{2});Residuals")
hFrameRes.GetXaxis().SetTitleOffset(1)
hFrameRes.GetYaxis().SetTitleOffset(0.57)
hFrameRes.GetXaxis().SetTitleSize(0.125)
hFrameRes.GetYaxis().SetTitleSize(0.125)
hFrameRes.GetXaxis().SetLabelSize(0.10)
hFrameRes.GetYaxis().SetLabelSize(0.10)
hFrameRes.GetYaxis().SetLabelOffset(0.01)
hFrameRes.GetYaxis().SetNdivisions(505)


hRes = hData5_5p5.Clone("hRes")

for iPt in range(hData5_5p5.GetNbinsX()):
  hRes.SetBinContent(iPt+1, hData5_5p5.GetBinContent(iPt+1) - hBkg0_5_5p5.Eval(hData5_5p5.GetBinCenter(iPt+1)) - hBkg1_5_5p5.Eval(hData5_5p5.GetBinCenter(iPt+1)))

# Set Style
hRes.SetMarkerStyle(ROOT.kFullCircle)
hRes.SetMarkerSize(0.9)
hRes.SetMarkerColor(ROOT.kBlack)
hRes.SetLineColor(ROOT.kBlack)
hRes.SetLineWidth(2)

DsmesonText = ROOT.TLatex(0.68, 0.87, "D_{s}^{+} meson")
DsmesonText.SetNDC()
DsmesonText.SetTextSize(0.1)
DsMeanText = ROOT.TLatex(0.68, 0.77, f"#mu = ({hRawYieldsMean.GetBinContent(hRawYieldsMean.FindBin(5.01))*1000:.1f} #pm {hRawYieldsMean.GetBinError(hRawYieldsMean.FindBin(5.01))*1000:.1f}) Me#kern[-0.05]{{V}}/#it{{c}}^{{2}}")
DsMeanText.SetNDC()
DsMeanText.SetTextSize(0.07)
DsSigmaText = ROOT.TLatex(0.68, 0.68, f"#sigma = ({hRawYieldsSigma.GetBinContent(hRawYieldsSigma.FindBin(5.01))*1000:.1f} #pm {hRawYieldsSigma.GetBinError(hRawYieldsSigma.FindBin(5.01))*1000:.1f}) Me#kern[-0.05]{{V}}/#it{{c}}^{{2}}")
DsSigmaText.SetNDC()
DsSigmaText.SetTextSize(0.07)
DsSignalText = ROOT.TLatex(0.68, 0.59, f"#it{{S}} = {hRawYields.GetBinContent(hRawYields.FindBin(5.01)):.0f} #pm {hRawYields.GetBinError(hRawYields.FindBin(5.01)):.0f}")
DsSignalText.SetNDC()
DsSignalText.SetTextSize(0.07)

DplusmesonText = ROOT.TLatex(0.17, 0.87, "D^{+} meson")
DplusmesonText.SetNDC()
DplusmesonText.SetTextSize(0.1)
DplusMeanText = ROOT.TLatex(0.17, 0.77, f"#mu = ({hRawYieldsMeanSecondPeak.GetBinContent(hRawYieldsMeanSecondPeak.FindBin(5.01))*1000:.1f} #pm {hRawYieldsMeanSecondPeak.GetBinError(hRawYieldsMeanSecondPeak.FindBin(5.01))*1000:.1f}) Me#kern[-0.05]{{V}}/#it{{c}}^{{2}}")
DplusMeanText.SetNDC()
DplusMeanText.SetTextSize(0.07)
DplusSigmaText = ROOT.TLatex(0.17, 0.68, f"#sigma = ({hRawYieldsSigmaSecondPeak.GetBinContent(hRawYieldsSigmaSecondPeak.FindBin(5.01))*1000:.1f} #pm {hRawYieldsSigmaSecondPeak.GetBinError(hRawYieldsSigmaSecondPeak.FindBin(5.01))*1000:.1f}) Me#kern[-0.05]{{V}}/#it{{c}}^{{2}}")
DplusSigmaText.SetNDC()
DplusSigmaText.SetTextSize(0.07)
DplusSignalText = ROOT.TLatex(0.17, 0.59, f"#it{{S}} = {hRawYieldsSecondPeak.GetBinContent(hRawYieldsSecondPeak.FindBin(5.01)):.0f} #pm {hRawYieldsSecondPeak.GetBinError(hRawYieldsSecondPeak.FindBin(5.01)):.0f}")
DplusSignalText.SetNDC()
DplusSignalText.SetTextSize(0.07)

hSignal0_5_5p5.Draw("same")
hSignal1_5_5p5.Draw("same")
hRes.Draw("pe same")
DsmesonText.Draw("same")
DsMeanText.Draw("same")
DsSigmaText.Draw("same")
DsSignalText.Draw("same")

DplusmesonText.Draw("same")
DplusMeanText.Draw("same")
DplusSigmaText.Draw("same")
DplusSignalText.Draw("same")


canvas_5_5p5.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/InvMassFit_5_5p5.png")
canvas_5_5p5.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/InvMassFit_5_5p5.pdf")
canvas_5_5p5.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/InvMassFit_5_5p5.eps")


