import ROOT 
import numpy as np

Ds_CentralFileName = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_newMC.root"
Dplus_CentralFileName = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB_newMC.root"
Ds_CentralVarFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/Central/Ds/CutVarDs_pp13TeV_MB_newMC.root"
Dplus_CentralVarFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/Central/Dplus/CutVarDplus_pp13TeV_MB_newMC.root"
Ds_FirstFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/First/Ds/CutVarDs_pp13TeV_MB_newMC.root"
Dplus_FirstFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/First/Dplus/CutVarDplus_pp13TeV_MB_newMC.root"
Ds_LastFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/Last/Ds/CutVarDs_pp13TeV_MB_newMC.root"
Dplus_LastFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/Last/Dplus/CutVarDplus_pp13TeV_MB_newMC.root"

# Open the files and get the fractions
Ds_CentralFile = ROOT.TFile(Ds_CentralFileName)
hDsCentral = Ds_CentralFile.Get("hRawFracPrompt")
hDsCentral.SetDirectory(0)
Ds_CentralFile.Close()

Dplus_CentralFile = ROOT.TFile(Dplus_CentralFileName)
hDplusCentral = Dplus_CentralFile.Get("hRawFracPrompt")
hDplusCentral.SetDirectory(0)
Dplus_CentralFile.Close()

Ds_CentralVarFile = ROOT.TFile(Ds_CentralVarFileName)
hDsCentralVar = Ds_CentralVarFile.Get("hRawFracPrompt")
hDsCentralVar.SetDirectory(0)
Ds_CentralVarFile.Close()

Dplus_CentralVarFile = ROOT.TFile(Dplus_CentralVarFileName)
hDplusCentralVar = Dplus_CentralVarFile.Get("hRawFracPrompt")
hDplusCentralVar.SetDirectory(0)
Dplus_CentralVarFile.Close()

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
hDsCentral.SetName("hDsCentral")
hDplusCentral.SetName("hDplusCentral")
hDsCentralVar.SetName("hDsCentralVar")
hDplusCentralVar.SetName("hDplusCentralVar")
hDsFirst.SetName("hDsFirst")
hDplusFirst.SetName("hDplusFirst")
hDsLast.SetName("hDsLast")
hDplusLast.SetName("hDplusLast")

# Get the ratios
hRatioCentral = hDsCentral.Clone("hRatioCentral")
hRatioCentral.Divide(hDplusCentral)
hRatioCentralVar = hDsCentralVar.Clone("hRatioCentralVar")
hRatioCentralVar.Divide(hDplusCentralVar)
hRatioFirst = hDsFirst.Clone("hRatioFirst")
hRatioFirst.Divide(hDplusFirst)
hRatioLast = hDsLast.Clone("hRatioLast")
hRatioLast.Divide(hDplusLast)

# Get the systematic errors
hSys = hRatioCentral.Clone("hSys")
for iPt in range(1, hSys.GetNbinsX()+1):
    syst = np.std([hRatioCentral.GetBinContent(iPt), hRatioCentralVar.GetBinContent(iPt), hRatioFirst.GetBinContent(iPt), hRatioLast.GetBinContent(iPt)])
    hSys.SetBinContent(iPt, syst)
    print(f"{hSys.GetBinLowEdge(iPt)} < pT < {hSys.GetBinLowEdge(iPt)+hSys.GetBinWidth(iPt)}: ", syst)

# Save the histograms
outputFile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/FDSystematic.root", "RECREATE")

canvasSyst = ROOT.TCanvas("canvasSyst", "canvasSyst", 800, 600)
canvasSyst.DrawFrame(0, 0, 24, 0.2, ";#it{p}_{T} (GeV/#it{c});Syst")
hSys.Draw("p same")
canvasSyst.Write()


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)
canvasDs = ROOT.TCanvas("canvasDs", "canvasDs", 800, 600)#
canvasDs.DrawFrame(0, 0.5, 24, 1.5, ";#it{p}_{T} (GeV/#it{c});#it{f}_{D_{s}^{+}}")
hDsCentral.Draw("pe PLC PMC same")
hDsCentralVar.Draw("pe PLC PMC same")
hDsFirst.Draw("pe PLC PMC same")
hDsLast.Draw("pe PLC PMC same")

legendDs = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
legendDs.AddEntry(hDsCentral, "Full cut set", "pe")
legendDs.AddEntry(hDsCentralVar, "Central 9 cuts", "pe")
legendDs.AddEntry(hDsFirst, "First 9 cuts", "pe")
legendDs.AddEntry(hDsLast, "Last 9 cuts", "pe")
legendDs.SetBorderSize(0)
legendDs.SetFillStyle(0)
legendDs.Draw()
canvasDs.Write()

canvasDplus = ROOT.TCanvas("canvasDplus", "canvasDplus", 800, 600)
canvasDplus.DrawFrame(0, 0.5, 24, 1.5, ";#it{p}_{T} (GeV/#it{c});#it{f}_{D^{+}}")
hDplusCentral.Draw("pe PLC PMC same")
hDplusCentralVar.Draw("pe PLC PMC same")
hDplusFirst.Draw("pe PLC PMC same")
hDplusLast.Draw("pe PLC PMC same")
legendDplus = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
legendDplus.AddEntry(hDplusCentral, "Full cut set", "pe")
legendDplus.AddEntry(hDplusCentralVar, "Central 9 cuts", "pe")
legendDplus.AddEntry(hDplusFirst, "First 9 cuts", "pe")
legendDplus.AddEntry(hDplusLast, "Last 9 cuts", "pe")
legendDplus.SetBorderSize(0)
legendDplus.SetFillStyle(0)
legendDplus.Draw()
canvasDplus.Write()

outputFile.Close()
