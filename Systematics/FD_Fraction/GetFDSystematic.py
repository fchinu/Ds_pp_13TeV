import ROOT 
import numpy as np
import itertools

Ds_FullFileName = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root"
Dplus_FullFileName = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root"
Ds_CentralFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/Central/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root"
Dplus_CentralFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/Central/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root"
Ds_FirstFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/First/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root"
Dplus_FirstFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/First/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root"
Ds_LastFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/Last/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root"
Dplus_LastFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/Last/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root"

ROOT.gStyle.SetPalette(ROOT.kRainbow)

AssignedSyst = [0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03]

# Open the files and get the fractions
Ds_FullFile = ROOT.TFile(Ds_FullFileName)
hDsFull = Ds_FullFile.Get("hRawFracPrompt")
hDsFull.SetDirectory(0)
Ds_FullFile.Close()

Dplus_FullFile = ROOT.TFile(Dplus_FullFileName)
hDplusFull = Dplus_FullFile.Get("hRawFracPrompt")
hDplusFull.SetDirectory(0)
Dplus_FullFile.Close()

Ds_CentralFile = ROOT.TFile(Ds_CentralFileName)
hDsCentral = Ds_CentralFile.Get("hRawFracPrompt")
hDsCentral.SetDirectory(0)
Ds_CentralFile.Close()

Dplus_CentralFile = ROOT.TFile(Dplus_CentralFileName)
hDplusCentral = Dplus_CentralFile.Get("hRawFracPrompt")
hDplusCentral.SetDirectory(0)
Dplus_CentralFile.Close()

Ds_FirstFile = ROOT.TFile(Ds_FirstFileName)
hDsFirst = Ds_FirstFile.Get("hRawFracPrompt")
hDsFirst.SetDirectory(0)
Ds_FirstFile.Close()

Dplus_FirstFile = ROOT.TFile(Dplus_FirstFileName)
hDplusFirst = Dplus_FirstFile.Get("hRawFracPrompt")
hDplusFirst.SetDirectory(0)
Dplus_FirstFile.Close()

Ds_LastFile = ROOT.TFile(Ds_LastFileName)
hDsLast = Ds_LastFile.Get("hRawFracPrompt")
hDsLast.SetDirectory(0)
Ds_LastFile.Close()

Dplus_LastFile = ROOT.TFile(Dplus_LastFileName)
hDplusLast = Dplus_LastFile.Get("hRawFracPrompt")
hDplusLast.SetDirectory(0)
Dplus_LastFile.Close()

# Set names
hDsFull.SetName("hDsFull")
hDplusFull.SetName("hDplusFull")
hDsCentral.SetName("hDsCentral")
hDplusCentral.SetName("hDplusCentral")
hDsFirst.SetName("hDsFirst")
hDplusFirst.SetName("hDplusFirst")
hDsLast.SetName("hDsLast")
hDplusLast.SetName("hDplusLast")

# Get the systematic errors
hRMS = hDsFull.Clone("hRMS")
hSyst = hDsFull.Clone("hSyst")
hCombinations = []
for i in range(16):
    hComb = hDsFull.Clone(f"hComb{i}")
    hCombinations.append(hComb)

for iPt in range(1, hRMS.GetNbinsX()+1):
    fracDs = [hDsFull.GetBinContent(iPt), hDsCentral.GetBinContent(iPt), hDsFirst.GetBinContent(iPt), hDsLast.GetBinContent(iPt)]
    fracDplus = [hDplusFull.GetBinContent(iPt), hDplusCentral.GetBinContent(iPt), hDplusFirst.GetBinContent(iPt), hDplusLast.GetBinContent(iPt)]
    combs = itertools.product(fracDs, fracDplus)
    ratios = []
    for iComb, comb in enumerate(combs):
        ratios.append(comb[0]/comb[1])
        hCombinations[iComb].SetBinContent(iPt, comb[0]/comb[1])
        hCombinations[iComb].SetBinError(iPt, 1.e-8)
    syst = np.std(ratios)
    hRMS.SetBinContent(iPt, syst)
    hRMS.SetBinError(iPt, 0)
    print(f"{hRMS.GetBinLowEdge(iPt)} < pT < {hRMS.GetBinLowEdge(iPt)+hRMS.GetBinWidth(iPt)}: ", syst)
    hSyst.SetBinContent(iPt, AssignedSyst[iPt-1])
    hSyst.SetBinError(iPt, 0)

# Save the histograms
outputFile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/FDSystematic.root", "RECREATE")

canvasSyst = ROOT.TCanvas("canvasSyst", "canvasSyst", 800, 600)
canvasSyst.DrawFrame(0, 0, 24, 0.2, ";#it{p}_{T} (GeV/#it{c});Syst")
hRMS.Draw("hist same")
canvasSyst.Write()
hRMS.Write("hRMS")
hSyst.Write("hSyst")


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kRainbow)
canvasDs = ROOT.TCanvas("canvasDs", "canvasDs", 800, 600)#
canvasDs.DrawFrame(0, 0.5, 24, 1.5, ";#it{p}_{T} (GeV/#it{c});#it{f}_{D_{s}^{+}}")
hDsFull.Draw("pe PLC PMC same")
hDsCentral.Draw("pe PLC PMC same")
hDsFirst.Draw("pe PLC PMC same")
hDsLast.Draw("pe PLC PMC same")

legendDs = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
legendDs.AddEntry(hDsFull, "Full cut set", "pe")
legendDs.AddEntry(hDsCentral, "Central 13 cuts", "pe")
legendDs.AddEntry(hDsFirst, "First 13 cuts", "pe")
legendDs.AddEntry(hDsLast, "Last 13 cuts", "pe")
legendDs.SetBorderSize(0)
legendDs.SetFillStyle(0)
legendDs.Draw()
canvasDs.Write()

canvasDplus = ROOT.TCanvas("canvasDplus", "canvasDplus", 800, 600)
canvasDplus.DrawFrame(0, 0.5, 24, 1.5, ";#it{p}_{T} (GeV/#it{c});#it{f}_{D^{+}}")
hDplusFull.Draw("pe PLC PMC same")
hDplusCentral.Draw("pe PLC PMC same")
hDplusFirst.Draw("pe PLC PMC same")
hDplusLast.Draw("pe PLC PMC same")
legendDplus = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
legendDplus.AddEntry(hDplusFull, "Full cut set", "pe")
legendDplus.AddEntry(hDplusCentral, "Central 13 cuts", "pe")
legendDplus.AddEntry(hDplusFirst, "First 13 cuts", "pe")
legendDplus.AddEntry(hDplusLast, "Last 13 cuts", "pe")
legendDplus.SetBorderSize(0)
legendDplus.SetFillStyle(0)
legendDplus.Draw()
canvasDplus.Write()
ROOT.gStyle.SetErrorX(1)
canvasCombinations = ROOT.TCanvas("canvasCombinations", "canvasCombinations", 800, 600)
canvasCombinations.DrawFrame(0, 0.5, 24, 1.5, ";#it{p}_{T} (GeV/#it{c});#it{f}_{D_{s}^{+}}/#it{f}_{D^{+}}")
for hComb in hCombinations:
    hComb.Draw("pe PLC PMC same")
canvasCombinations.Write()

outputFile.Close()

canvasCombinations.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/test.png")
