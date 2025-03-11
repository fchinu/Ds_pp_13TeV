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

ptMins = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 12] 
ptMaxs = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 12, 24]

hData = []
hTotal_func = []
hBkg0 = []
hBkg1 = []
hSignal0 = []
hSignal1 = []

# Open the file
InputFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/fits/fits.root")

for iPt in range(len(ptMins)):
  hData.append(InputFile.Get(f"hdata{ptMins[iPt]:.1f}_{ptMaxs[iPt]:.1f}"))
  hData[-1].SetDirectory(0)
  hTotal_func.append(InputFile.Get(f"total_func{ptMins[iPt]:.1f}_{ptMaxs[iPt]:.1f}"))
  hTotal_func[-1].SetDirectory(0)
  hBkg0.append(InputFile.Get(f"bkg_0{ptMins[iPt]:.1f}_{ptMaxs[iPt]:.1f}"))
  hBkg0[-1].SetDirectory(0)
  hBkg1.append(InputFile.Get(f"bkg_1{ptMins[iPt]:.1f}_{ptMaxs[iPt]:.1f}"))
  hBkg1[-1].SetDirectory(0)
  hSignal0.append(InputFile.Get(f"signal_0{ptMins[iPt]:.1f}_{ptMaxs[iPt]:.1f}"))
  hSignal0[-1].SetDirectory(0)
  hSignal1.append(InputFile.Get(f"signal_1{ptMins[iPt]:.1f}_{ptMaxs[iPt]:.1f}"))
  hSignal1[-1].SetDirectory(0)
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

for iPt in range(len(ptMins)):
  hData[iPt].SetMarkerStyle(ROOT.kFullCircle)
  hData[iPt].SetMarkerSize(0.9)
  hData[iPt].SetMarkerColor(ROOT.kBlack)
  hData[iPt].SetLineColor(ROOT.kBlack)
  hData[iPt].SetLineWidth(2)

  hTotal_func[iPt] = get_function_fromhist(hTotal_func[iPt])
  hTotal_func[iPt].SetLineColor(ROOT.kAzure-3)
  hTotal_func[iPt].SetLineWidth(3)

  hBkg0[iPt] = get_function_fromhist(hBkg0[iPt])
  hBkg0[iPt].SetLineColor(ROOT.kViolet-4)
  hBkg0[iPt].SetLineWidth(2)
  hBkg0[iPt].SetLineStyle(2)

  hBkg1[iPt] = get_function_fromhist(hBkg1[iPt])
  hBkg1[iPt].SetLineColor(ROOT.kRed)
  hBkg1[iPt].SetLineWidth(2)

  hSignal0[iPt].SetLineColor(ROOT.kBlack)
  hSignal0[iPt].SetLineWidth(2)
  hSignal0[iPt].SetFillColor(cFillBlue)

  hSignal1[iPt].SetLineColor(ROOT.kBlack)
  hSignal1[iPt].SetLineWidth(2)
  hSignal1[iPt].SetFillColor(cFillGreen)

# Add text

ALICEtext = ROOT.TLatex(0.17, 0.87, "This Thesis")
ALICEtext.SetNDC()
ALICEtext.SetTextSize(0.06)
ALICEtext.SetTextFont(42)

ppEnergytext = ROOT.TLatex(0.17, 0.81, "pp, #sqrt{#it{s}} = 13.6 Te#kern[-0.05]{V}")
ppEnergytext.SetNDC()
ppEnergytext.SetTextSize(0.05)
ppEnergytext.SetTextFont(42)

textDsDecay = ROOT.TLatex(0.17, 0.73, "D_{s}^{+}, D^{+}#rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+} ")
textDsDecay.SetNDC()
textDsDecay.SetTextSize(0.05)
textDsDecay.SetTextFont(42)

textDsDecayConj = ROOT.TLatex(0.17, 0.68, 'and charge conj.')
textDsDecayConj.SetNDC()
textDsDecayConj.SetTextFont(42)
textDsDecayConj.SetTextSize(0.05)

textPt = []
for iPt in range(len(ptMins)):
  textPt.append(ROOT.TLatex(0.17, 0.605, f'{ptMins[iPt]:.1f} < #it{{p}}_{{T}} < {ptMaxs[iPt]:.1f} Ge#kern[-0.05]{{V}}/#it{{c}}'))
  textPt[-1].SetNDC()
  textPt[-1].SetTextFont(42)
  textPt[-1].SetTextSize(0.05)

legend = ROOT.TLegend(0.565, 0.6, 0.865, 0.92)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.043)
legend.AddEntry(hData[0], "Data", "pe")
legend.AddEntry(hSignal0[0], "D_{s}^{+} signal", "f")
legend.AddEntry(hSignal1[0], "D^{+} signal", "f")
legend.AddEntry(hBkg1[0], "Combinatorial bkg.", "l")
legend.AddEntry(hBkg0[0], "D^{+}#rightarrow#pi^{+}K^{#font[122]{-}}#pi^{+} bkg.", "l")
legend.AddEntry(hTotal_func[0], "Total fit function", "l")

canvas = []
for iPt in range(len(ptMins)):
  canvas.append(ROOT.TCanvas(f"canvas_{ptMins[iPt]}_{ptMaxs[iPt]}", f"canvas_{ptMins[iPt]}_{ptMaxs[iPt]}", 800, 700))

  hFrame = canvas[-1].DrawFrame(hTotal_func[iPt].GetXmin(), 0., hData[iPt].GetBinLowEdge(hData[iPt].GetNbinsX())+hData[iPt].GetBinWidth(hData[iPt].GetNbinsX()), hData[iPt].GetMaximum()*2, f";#it{{M}}(KK#pi) (Ge#kern[-0.05]{{V}}/#it{{c}}^{{2}});Counts per {hData[iPt].GetBinWidth(1)*1000:.0f} Me#kern[-0.05]{{V}}/#it{{c}}^{{2}}")
  hFrame.GetXaxis().SetTitleOffset(1.)
  hFrame.GetYaxis().SetTitleOffset(1.4)
  hFrame.GetXaxis().SetTitleSize(0.05)
  hFrame.GetYaxis().SetTitleSize(0.05)
  hFrame.GetXaxis().SetLabelSize(0.045)
  hFrame.GetYaxis().SetLabelSize(0.045)

  hSignal0[iPt].Draw("same")
  hSignal1[iPt].Draw("same")
  hBkg0[iPt].Draw("same")
  hBkg1[iPt].Draw("same")
  hTotal_func[iPt].Draw("same")
  hData[iPt].Draw("pe same")
  ALICEtext.Draw("same")
  ppEnergytext.Draw("same")
  textPt[iPt].Draw("same")
  textDsDecay.Draw("same")
  textDsDecayConj.Draw("same")
  legend.Draw("same")

  ROOT.gPad.RedrawAxis()


  canvas[-1].SaveAs(f"/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/AllPtBins/pt_{ptMins[iPt]:.1f}_{ptMaxs[iPt]:.1f}.png")
  canvas[-1].SaveAs(f"/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/AllPtBins/pt_{ptMins[iPt]:.1f}_{ptMaxs[iPt]:.1f}.pdf")
  canvas[-1].SaveAs(f"/home/fchinu/Run3/Ds_pp_13TeV/Figures/InvMass/AllPtBins/pt_{ptMins[iPt]:.1f}_{ptMaxs[iPt]:.1f}.eps")





