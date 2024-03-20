from ROOT import TFile, TH1F, TCanvas, TLegend, kRed, kAzure, kTeal, kOrange
import os
import sys
sys.path.append("/home/fchinu/DmesonAnalysis/utils")
from StyleFormatter import SetObjectStyle

EffFolder = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/Efficiency"

file_list = os.listdir(EffFolder)
file_list.sort()
histosDsPrompt = []
histosDplusPrompt = []
histosDsFD = []
histosDplusFD = []

for idx, file_name in enumerate(file_list):
    print(file_name)
    file_path = os.path.join(EffFolder, file_name)
    if idx == 0:    # Open the first file to get the binning
        EffFile = TFile.Open(file_path)
        histo0 = EffFile.Get("Eff_DsPrompt")
        histo0.SetDirectory(0)
        EffFile.Close()
        for iPt in range(histo0.GetNbinsX()):
            histosDsPrompt.append(TH1F(f"Eff_DsPrompt_{histo0.GetBinLowEdge(iPt+1)}_{histo0.GetBinLowEdge(iPt+1)+histo0.GetBinWidth(iPt+1)}", f"Efficiency; Cutset; Efficiency",len(file_list), 0,len(file_list)))
            histosDplusPrompt.append(TH1F(f"Eff_DplusPrompt_{histo0.GetBinLowEdge(iPt+1)}_{histo0.GetBinLowEdge(iPt+1)+histo0.GetBinWidth(iPt+1)}", f"Efficiency; Cutset; Efficiency",len(file_list), 0,len(file_list)))
            histosDsFD.append(TH1F(f"Eff_DsFD_{histo0.GetBinLowEdge(iPt+1)}_{histo0.GetBinLowEdge(iPt+1)+histo0.GetBinWidth(iPt+1)}", f"Efficiency; Cutset; Efficiency",len(file_list), 0,len(file_list)))
            histosDplusFD.append(TH1F(f"Eff_DplusFD_{histo0.GetBinLowEdge(iPt+1)}_{histo0.GetBinLowEdge(iPt+1)+histo0.GetBinWidth(iPt+1)}", f"Efficiency; Cutset; Efficiency",len(file_list), 0,len(file_list)))

    for iPt in range(histo0.GetNbinsX()):
        EffFile = TFile.Open(file_path)
        histosDsPrompt[iPt].SetBinContent(idx+1, EffFile.Get("Eff_DsPrompt").GetBinContent(iPt+1))
        histosDsPrompt[iPt].SetBinError(idx+1, EffFile.Get("Eff_DsPrompt").GetBinError(iPt+1))
        histosDplusPrompt[iPt].SetBinContent(idx+1, EffFile.Get("Eff_DplusPrompt").GetBinContent(iPt+1))
        histosDplusPrompt[iPt].SetBinError(idx+1, EffFile.Get("Eff_DplusPrompt").GetBinError(iPt+1))
        histosDsFD[iPt].SetBinContent(idx+1, EffFile.Get("Eff_DsFD").GetBinContent(iPt+1))
        histosDsFD[iPt].SetBinError(idx+1, EffFile.Get("Eff_DsFD").GetBinError(iPt+1))
        histosDplusFD[iPt].SetBinContent(idx+1, EffFile.Get("Eff_DplusFD").GetBinContent(iPt+1))
        histosDplusFD[iPt].SetBinError(idx+1, EffFile.Get("Eff_DplusFD").GetBinError(iPt+1))


outputFile = TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/Efficiency.root", "RECREATE")
for iPt in range(histo0.GetNbinsX()):   # First save the canvases so they appear first in the file
    c = TCanvas()
    c.SetLogy()
    c.DrawFrame(0,1.e-6,len(file_list),1, ";Cutset;Efficiency")
    SetObjectStyle(histosDsPrompt[iPt], color=kRed, markerstyle=20, markersize=1.5, linewidth=1)
    histosDsPrompt[iPt].Draw("same")
    SetObjectStyle(histosDplusPrompt[iPt], color=kAzure+3, markerstyle=20, markersize=1.5, linewidth=1)
    histosDplusPrompt[iPt].Draw("same")
    SetObjectStyle(histosDsFD[iPt], color=kTeal-7, markerstyle=20, markersize=1.5, linewidth=1)
    histosDsFD[iPt].Draw("same")
    SetObjectStyle(histosDplusFD[iPt], color=kOrange-3, markerstyle=20, markersize=1.5, linewidth=1)
    histosDplusFD[iPt].Draw("same")

    leg = TLegend(0.6,0.7,0.9,0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(histosDsPrompt[iPt], "D_{s}^{+} Prompt", "lp")
    leg.AddEntry(histosDplusPrompt[iPt], "D^{+} Prompt", "lp")
    leg.AddEntry(histosDsFD[iPt], "D_{s}^{+} FD", "lp")
    leg.AddEntry(histosDplusFD[iPt], "D^{+} FD", "lp")
    leg.Draw()

    c.Write(f"c_{histo0.GetBinLowEdge(iPt+1)}_{histo0.GetBinLowEdge(iPt+1)+histo0.GetBinWidth(iPt+1)}")

for iPt in range(histo0.GetNbinsX()):
    histosDsPrompt[iPt].Write()
    histosDplusPrompt[iPt].Write()
    histosDsFD[iPt].Write()
    histosDplusFD[iPt].Write()

outputFile.Close()
