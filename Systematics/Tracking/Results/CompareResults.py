import ROOT
import numpy as np
import pandas as pd
import sys
sys.path.append("/home/fchinu/DmesonAnalysis/utils")
from StyleFormatter import SetGlobalStyle

ColsToUse = [ROOT.kRed,ROOT.kAzure+3,ROOT.kTeal-7]
SetGlobalStyle(padrightmargin=0.05, padleftmargin=0.15, padbottommargin=0.13, padtopmargin=0.05, opttitle=0, titleoffsetx=1.2, titleoffsety=1.5)

AssignedSyst = [0.08, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02]

#################
# Efficiency
#################

EffCentralFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/Efficiencies_LHC24d3a.root"
EffPtSmearingFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/ptSmearing/Efficiency/Efficiency.root"
EffWorseResoPullsFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/worseResoPullsData1p1/Efficiency/Efficiency.root"

# Load the histograms

histoNames = ["Eff_DsPrompt", "Eff_DplusPrompt", "Eff_DsFD", "Eff_DplusFD"]

EffsCentral = []
EffsPtSmearing = []
EffsWorseResoPulls = []

EffCentralFile = ROOT.TFile.Open(EffCentralFileName)
EffPtSmearingFile = ROOT.TFile.Open(EffPtSmearingFileName)
EffWorseResoPullsFile = ROOT.TFile.Open(EffWorseResoPullsFileName)

for histoName in histoNames:
    EffsCentral.append(EffCentralFile.Get(histoName))
    EffsPtSmearing.append(EffPtSmearingFile.Get(histoName))
    EffsWorseResoPulls.append(EffWorseResoPullsFile.Get(histoName))
    EffsCentral[-1].SetDirectory(0)
    EffsCentral[-1].SetMarkerStyle(ROOT.kFullCircle)
    EffsCentral[-1].SetMarkerColor(ColsToUse[0])
    EffsCentral[-1].SetLineColor(ColsToUse[0])
    EffsPtSmearing[-1].SetDirectory(0)
    EffsPtSmearing[-1].SetMarkerStyle(ROOT.kFullCircle)
    EffsPtSmearing[-1].SetMarkerColor(ColsToUse[1])
    EffsPtSmearing[-1].SetLineColor(ColsToUse[1])
    EffsWorseResoPulls[-1].SetDirectory(0)
    EffsWorseResoPulls[-1].SetMarkerStyle(ROOT.kFullCircle)
    EffsWorseResoPulls[-1].SetMarkerColor(ColsToUse[2])
    EffsWorseResoPulls[-1].SetLineColor(ColsToUse[2])

EffCentralFile.Close()
EffPtSmearingFile.Close()
EffWorseResoPullsFile.Close()

# Evaluate ratios to central case

EffRatiosPtSmearing = []
EffRatiosWorseResoPulls = []

for i in range(len(histoNames)):
    EffRatiosPtSmearing.append(EffsPtSmearing[i].Clone("EffRatiosPtSmearing_" + histoNames[i]))
    EffRatiosWorseResoPulls.append(EffsWorseResoPulls[i].Clone("EffRatiosWorseResoPulls_" + histoNames[i]))
    EffRatiosPtSmearing[-1].Divide(EffsCentral[i])
    EffRatiosWorseResoPulls[-1].Divide(EffsCentral[i])

# Plot the histograms

cEffDsPrompt = ROOT.TCanvas("cEffDsPrompt", "cEffDsPrompt", 800, 600)
cEffDsPrompt.Divide(2)
cEffDsPrompt.cd(1).DrawFrame(0., 0.0001, 24, 1.0, ";#it{p}_{T} (GeV/c);Prompt D_{s}^{+} Efficiency")
ROOT.gPad.SetLogy()

EffsCentral[0].Draw("e same")
EffsPtSmearing[0].Draw("e same")
EffsWorseResoPulls[0].Draw("e same")

legEffDsPrompt = ROOT.TLegend(0.6, 0.3, 0.9, 0.45)
legEffDsPrompt.SetTextSize(0.035)
legEffDsPrompt.AddEntry(EffsCentral[0], "Central", "PL")
legEffDsPrompt.AddEntry(EffsPtSmearing[0], "Pt smearing", "PL")
legEffDsPrompt.AddEntry(EffsWorseResoPulls[0], "Worse reso pulls", "PL")
legEffDsPrompt.SetBorderSize(0)
legEffDsPrompt.SetFillStyle(0)
legEffDsPrompt.Draw()

cEffDsPrompt.cd(2).DrawFrame(0., 0.75, 24, 1.25, ";#it{p}_{T} (GeV/c);Ratio")

EffRatiosPtSmearing[0].Draw("e same")
EffRatiosWorseResoPulls[0].Draw("e same")

cEffDsPrompt.Modified()
cEffDsPrompt.Update()

cEffDplusPrompt = ROOT.TCanvas("cEffDplusPrompt", "cEffDplusPrompt", 800, 600)
cEffDplusPrompt.Divide(2)
cEffDplusPrompt.cd(1).DrawFrame(0., 0.0001, 24, 1.0, ";#it{p}_{T} (GeV/c);Prompt D^{+} Efficiency")
ROOT.gPad.SetLogy()

EffsCentral[1].Draw("e same")
EffsPtSmearing[1].Draw("e same")
EffsWorseResoPulls[1].Draw("e same")

legEffDplusPrompt = ROOT.TLegend(0.6, 0.3, 0.9, 0.45)
legEffDplusPrompt.SetTextSize(0.035)
legEffDplusPrompt.AddEntry(EffsCentral[1], "Central", "PL")
legEffDplusPrompt.AddEntry(EffsPtSmearing[1], "Pt smearing", "PL")
legEffDplusPrompt.AddEntry(EffsWorseResoPulls[1], "Worse reso pulls", "PL")
legEffDplusPrompt.SetBorderSize(0)
legEffDplusPrompt.SetFillStyle(0)
legEffDplusPrompt.Draw()

cEffDplusPrompt.cd(2).DrawFrame(0., 0.75, 24, 1.25, ";#it{p}_{T} (GeV/c);Ratio")

EffRatiosPtSmearing[1].Draw("e same")
EffRatiosWorseResoPulls[1].Draw("e same")

cEffDsFD = ROOT.TCanvas("cEffDsFD", "cEffDsFD", 800, 600)
cEffDsFD.Divide(2)
cEffDsFD.cd(1).DrawFrame(0., 0.0001, 24, 1.0, ";#it{p}_{T} (GeV/c);FD D_{s}^{+} Efficiency")
ROOT.gPad.SetLogy()

EffsCentral[2].Draw("e same")
EffsPtSmearing[2].Draw("e same")
EffsWorseResoPulls[2].Draw("e same")

legEffDsFD = ROOT.TLegend(0.6, 0.3, 0.9, 0.45)
legEffDsFD.SetTextSize(0.035)
legEffDsFD.AddEntry(EffsCentral[2], "Central", "PL")
legEffDsFD.AddEntry(EffsPtSmearing[2], "Pt smearing", "PL")
legEffDsFD.AddEntry(EffsWorseResoPulls[2], "Worse reso pulls", "PL")
legEffDsFD.SetBorderSize(0)
legEffDsFD.SetFillStyle(0)
legEffDsFD.Draw()

cEffDsFD.cd(2).DrawFrame(0., 0.75, 24, 1.25, ";#it{p}_{T} (GeV/c);Ratio")

EffRatiosPtSmearing[2].Draw("e same")
EffRatiosWorseResoPulls[2].Draw("e same")

cEffDplusFD = ROOT.TCanvas("cEffDplusFD", "cEffDplusFD", 800, 600)
cEffDplusFD.Divide(2)
cEffDplusFD.cd(1).DrawFrame(0., 0.0001, 24, 1.0, ";#it{p}_{T} (GeV/c);FD D^{+} Efficiency")
ROOT.gPad.SetLogy()

EffsCentral[3].Draw("e same")
EffsPtSmearing[3].Draw("e same")
EffsWorseResoPulls[3].Draw("e same")

legEffDplusFD = ROOT.TLegend(0.6, 0.3, 0.9, 0.45)
legEffDplusFD.SetTextSize(0.035)
legEffDplusFD.AddEntry(EffsCentral[3], "Central", "PL")
legEffDplusFD.AddEntry(EffsPtSmearing[3], "Pt smearing", "PL")
legEffDplusFD.AddEntry(EffsWorseResoPulls[3], "Worse reso pulls", "PL")
legEffDplusFD.SetBorderSize(0)
legEffDplusFD.SetFillStyle(0)
legEffDplusFD.Draw()

cEffDplusFD.cd(2).DrawFrame(0., 0.75, 24, 1.25, ";#it{p}_{T} (GeV/c);Ratio")

EffRatiosPtSmearing[3].Draw("e same")
EffRatiosWorseResoPulls[3].Draw("e same")

#################
# FD fraction
#################

FDsCentralFileNameDs = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root"
FDsCentralFileNameDplus = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root"
FDPtSmearingFileNameDs = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/ptSmearing/FD_Fraction/Ds/CutVarDs_pp13TeV_MB_newMC.root"
FDPtSmearingFileNameDplus = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/ptSmearing/FD_Fraction/Dplus/CutVarDs_pp13TeV_MB_newMC.root"
FDWorseResoPullsFileNameDs = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/worseResoPullsData1p1/FD_Fraction/Ds/CutVarDs_pp13TeV_MB_newMC.root"
FDWorseResoPullsFileNameDplus = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/worseResoPullsData1p1/FD_Fraction/Dplus/CutVarDs_pp13TeV_MB_newMC.root"

# Load the histograms

FDsCentralFileDs = ROOT.TFile.Open(FDsCentralFileNameDs)
hFDDsCentral = FDsCentralFileDs.Get("hRawFracPrompt")
hFDDsCentral.SetDirectory(0)
hFDDsCentral.SetMarkerStyle(ROOT.kFullCircle)
hFDDsCentral.SetMarkerColor(ColsToUse[0])
hFDDsCentral.SetLineColor(ColsToUse[0])
FDsCentralFileDs.Close()
FDsCentralFileDplus = ROOT.TFile.Open(FDsCentralFileNameDplus)
hFDDplusCentral = FDsCentralFileDplus.Get("hRawFracPrompt")
hFDDplusCentral.SetDirectory(0)
hFDDplusCentral.SetMarkerStyle(ROOT.kFullCircle)
hFDDplusCentral.SetMarkerColor(ColsToUse[0])
hFDDplusCentral.SetLineColor(ColsToUse[0])
FDsCentralFileDplus.Close()
FDsPtSmearingFileDs = ROOT.TFile.Open(FDPtSmearingFileNameDs)
hFDDsPtSmearing = FDsPtSmearingFileDs.Get("hRawFracPrompt")
hFDDsPtSmearing.SetDirectory(0)
hFDDsPtSmearing.SetMarkerStyle(ROOT.kFullCircle)
hFDDsPtSmearing.SetMarkerColor(ColsToUse[1])
hFDDsPtSmearing.SetLineColor(ColsToUse[1])
FDsPtSmearingFileDs.Close()
FDsPtSmearingFileDplus = ROOT.TFile.Open(FDPtSmearingFileNameDplus)
hFDDplusPtSmearing = FDsPtSmearingFileDplus.Get("hRawFracPrompt")
hFDDplusPtSmearing.SetDirectory(0)
hFDDplusPtSmearing.SetMarkerStyle(ROOT.kFullCircle)
hFDDplusPtSmearing.SetMarkerColor(ColsToUse[1])
hFDDplusPtSmearing.SetLineColor(ColsToUse[1])
FDsPtSmearingFileDplus.Close()
FDsWorseResoPullsFileDs = ROOT.TFile.Open(FDWorseResoPullsFileNameDs)
hFDDsWorseResoPulls = FDsWorseResoPullsFileDs.Get("hRawFracPrompt")
hFDDsWorseResoPulls.SetDirectory(0)
hFDDsWorseResoPulls.SetMarkerStyle(ROOT.kFullCircle)
hFDDsWorseResoPulls.SetMarkerColor(ColsToUse[2])
hFDDsWorseResoPulls.SetLineColor(ColsToUse[2])
FDsWorseResoPullsFileDs.Close()
FDsWorseResoPullsFileDplus = ROOT.TFile.Open(FDWorseResoPullsFileNameDplus)
hFDDplusWorseResoPulls = FDsWorseResoPullsFileDplus.Get("hRawFracPrompt")
hFDDplusWorseResoPulls.SetDirectory(0)
hFDDplusWorseResoPulls.SetMarkerStyle(ROOT.kFullCircle)
hFDDplusWorseResoPulls.SetMarkerColor(ColsToUse[2])
hFDDplusWorseResoPulls.SetLineColor(ColsToUse[2])
FDsWorseResoPullsFileDplus.Close()

# Evaluate ratios to central case

hFDDsRatiosPtSmearing = hFDDsPtSmearing.Clone("hFDDsRatiosPtSmearing")
hFDDsRatiosPtSmearing.Divide(hFDDsCentral)
hFDDplusRatiosPtSmearing = hFDDplusPtSmearing.Clone("hFDDplusRatiosPtSmearing")
hFDDplusRatiosPtSmearing.Divide(hFDDplusCentral)

hFDDsRatiosWorseResoPulls = hFDDsWorseResoPulls.Clone("hFDDsRatiosWorseResoPulls")
hFDDsRatiosWorseResoPulls.Divide(hFDDsCentral)
hFDDplusRatiosWorseResoPulls = hFDDplusWorseResoPulls.Clone("hFDDplusRatiosWorseResoPulls")
hFDDplusRatiosWorseResoPulls.Divide(hFDDplusCentral)

# Plot the histograms

cFDDs = ROOT.TCanvas("cFDDs", "cFDDs", 800, 600)
cFDDs.Divide(2)
cFDDs.cd(1).DrawFrame(0., 0.75, 24, 1.05, ";#it{p}_{T} (GeV/c);D_{s}^{+} prompt fraction")

hFDDsCentral.Draw("e same")
hFDDsPtSmearing.Draw("e same")
hFDDsWorseResoPulls.Draw("e same")

legFDDs = ROOT.TLegend(0.3, 0.3, 0.6, 0.45)
legFDDs.SetTextSize(0.035)
legFDDs.AddEntry(hFDDsCentral, "Central", "PL")
legFDDs.AddEntry(hFDDsPtSmearing, "Pt smearing", "PL")
legFDDs.AddEntry(hFDDsWorseResoPulls, "Worse reso pulls", "PL")
legFDDs.SetBorderSize(0)
legFDDs.SetFillStyle(0)
legFDDs.Draw()

cFDDs.cd(2).DrawFrame(0., 0.9, 24, 1.05, ";#it{p}_{T} (GeV/c);Ratio")

hFDDsRatiosPtSmearing.Draw("e same")
hFDDsRatiosWorseResoPulls.Draw("e same")

cFDDplus = ROOT.TCanvas("cFDDplus", "cFDDplus", 800, 600)
cFDDplus.Divide(2)
cFDDplus.cd(1).DrawFrame(0., 0.75, 24, 1.05, ";#it{p}_{T} (GeV/c);D^{+} prompt fraction")

hFDDplusCentral.Draw("e same")
hFDDplusPtSmearing.Draw("e same")
hFDDplusWorseResoPulls.Draw("e same")

legFDDplus = ROOT.TLegend(0.3, 0.3, 0.6, 0.45)
legFDDplus.SetTextSize(0.035)
legFDDplus.AddEntry(hFDDplusCentral, "Central", "PL")
legFDDplus.AddEntry(hFDDplusPtSmearing, "Pt smearing", "PL")
legFDDplus.AddEntry(hFDDplusWorseResoPulls, "Worse reso pulls", "PL")
legFDDplus.SetBorderSize(0)
legFDDplus.SetFillStyle(0)
legFDDplus.Draw()

cFDDplus.cd(2).DrawFrame(0., 0.9, 24, 1.05, ";#it{p}_{T} (GeV/c);Ratio")

hFDDplusRatiosPtSmearing.Draw("e same")
hFDDplusRatiosWorseResoPulls.Draw("e same")

#################
# Ratios
#################

RatioCentralFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Ratios/DsOverDplus_LHC24d3a.root"
RatioPtSmearingFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/ptSmearing/Ratio/DsOverDplus_LHC24d3a.root"
RatioWorseResoPullsFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/worseResoPullsData1p1/Ratio/DsOverDplus_LHC24d3a.root"

# Load the histograms

RatioCentralFile = ROOT.TFile.Open(RatioCentralFileName)
hRatioCentral = RatioCentralFile.Get("hRatio")
hRatioCentral.SetDirectory(0)
hRatioCentral.SetMarkerStyle(ROOT.kFullCircle)
hRatioCentral.SetMarkerColor(ColsToUse[0])
hRatioCentral.SetLineColor(ColsToUse[0])
RatioCentralFile.Close()
RatioPtSmearingFile = ROOT.TFile.Open(RatioPtSmearingFileName)
hRatioPtSmearing = RatioPtSmearingFile.Get("hRatio")
hRatioPtSmearing.SetDirectory(0)
hRatioPtSmearing.SetMarkerStyle(ROOT.kFullCircle)
hRatioPtSmearing.SetMarkerColor(ColsToUse[1])
hRatioPtSmearing.SetLineColor(ColsToUse[1])
RatioPtSmearingFile.Close()
RatioWorseResoPullsFile = ROOT.TFile.Open(RatioWorseResoPullsFileName)
hRatioWorseResoPulls = RatioWorseResoPullsFile.Get("hRatio")
hRatioWorseResoPulls.SetDirectory(0)
hRatioWorseResoPulls.SetMarkerStyle(ROOT.kFullCircle)
hRatioWorseResoPulls.SetMarkerColor(ColsToUse[2])
hRatioWorseResoPulls.SetLineColor(ColsToUse[2])
RatioWorseResoPullsFile.Close()

# Plot the histograms

cRatio = ROOT.TCanvas("cRatio", "cRatio", 800, 600)
cRatio.Divide(2)
cRatio.cd(1).DrawFrame(0., 0., 24, 1., ";#it{p}_{T} (GeV/c);D_{s}^{+}/D^{+} ratio")

hRatioCentral.Draw("e same")
hRatioPtSmearing.Draw("e same")
hRatioWorseResoPulls.Draw("e same")

legRatio = ROOT.TLegend(0.5, 0.2, 0.8, 0.35)
legRatio.SetTextSize(0.035)
legRatio.AddEntry(hRatioCentral, "Central", "PL")
legRatio.AddEntry(hRatioPtSmearing, "Pt smearing", "PL")
legRatio.AddEntry(hRatioWorseResoPulls, "Worse reso pulls", "PL")
legRatio.SetBorderSize(0)
legRatio.SetFillStyle(0)
legRatio.Draw()

cRatio.cd(2).DrawFrame(0., 0.7, 24, 1.4, ";#it{p}_{T} (GeV/c);Ratio")

hRatiosRatioPtSmearing = hRatioPtSmearing.Clone("hRatiosRatioPtSmearing")
hRatiosRatioWorseResoPulls = hRatioWorseResoPulls.Clone("hRatiosRatioWorseResoPulls")

hRatiosRatioPtSmearing.Divide(hRatioCentral)
hRatiosRatioWorseResoPulls.Divide(hRatioCentral)

hRatiosRatioPtSmearing.Draw("e same")
hRatiosRatioWorseResoPulls.Draw("e same")

#################
# Save the histograms
#################

outputFile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/CompareResults.root", "RECREATE")

cEffDsPrompt.Write()
cEffDplusPrompt.Write()
cEffDsFD.Write()
cEffDplusFD.Write()
cFDDs.Write()
cFDDplus.Write()
cRatio.Write()

outputFile.Close()


cEffDsPrompt.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/EffDsPrompt.png")
cEffDplusPrompt.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/EffDplusPrompt.png")
cEffDsFD.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/EffDsFD.png")
cEffDplusFD.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/EffDplusFD.png")
cFDDs.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/FDDs.png")
cFDDplus.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/FDDplus.png")
cRatio.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/Ratio.png")


#################
# Save systematic uncertainty
#################

hSemidisp = hRatioCentral.Clone("hSemidisp")
hSyst = hRatioCentral.Clone("hSyst")
for iPt in range(hRatioCentral.GetNbinsX()):
    maxRatio = max(hRatioPtSmearing.GetBinContent(iPt+1), hRatioWorseResoPulls.GetBinContent(iPt+1), hRatioCentral.GetBinContent(iPt+1))
    minRatio = min(hRatioPtSmearing.GetBinContent(iPt+1), hRatioWorseResoPulls.GetBinContent(iPt+1), hRatioCentral.GetBinContent(iPt+1))
    syst = 0.5*(maxRatio - minRatio)/(hRatioCentral.GetBinContent(iPt+1))
    hSemidisp.SetBinContent(iPt+1, syst)
    hSemidisp.SetBinError(iPt+1, 0.)
    print(f"{hRatioCentral.GetBinLowEdge(iPt+1):.3f} < pT < {hRatioCentral.GetBinLowEdge(iPt+2):.3f}: {syst}")
    hSyst.SetBinContent(iPt+1, AssignedSyst[iPt])
    hSyst.SetBinError(iPt+1, 0.)

outputFileSyst = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/Results/SystUncertainty.root", "RECREATE")
hSemidisp.Write()
hSyst.Write()
outputFileSyst.Close()