from ROOT import TFile, TCanvas, TH1F
import pandas as pd
import numpy as np

RawYieldsFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Tests/Compare_w_wo_Normalized_DL/w_norm_DL/Projections_RawYields/RawYields_Data.root"
EffFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Tests/Compare_w_wo_Normalized_DL/w_norm_DL/Efficiency/Efficiency.root"
outputFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Tests/Compare_w_wo_Normalized_DL/w_norm_DL/Ratio/DsOverDplus.root"

RawYieldsFile = TFile.Open(RawYieldsFileName)
EffFile = TFile.Open(EffFileName)
outputFile = TFile.Open(outputFileName, "RECREATE")

DsRawYields = RawYieldsFile.Get("hRawYields")
DsEff = EffFile.Get("Eff_DsPrompt")
DplusRawYields = RawYieldsFile.Get("hRawYieldsSecondPeak")
DplusEff = EffFile.Get("Eff_DplusPrompt")

CorrectedDsYields = DsRawYields.Clone("hCorrectedDsYields")
CorrectedDsYields.Divide(DsEff)
CorrectedDsYields.Scale(1/2.21e-2)  #BR https://pdg.lbl.gov/2023/listings/rpp2023-list-Ds-plus-minus.pdf

CorrectedDplusYields = DplusRawYields.Clone("hCorrectedDplusYields")
CorrectedDplusYields.Divide(DplusEff)
CorrectedDplusYields.Scale(1/2.69e-3)  #BR https://pdg.lbl.gov/2023/listings/rpp2023-list-D-plus-minus.pdf

Ratio = CorrectedDsYields.Clone("hRatio")
Ratio.SetTitle(";p_{T} (GeV/c);D_{s}^{+}/D^{+} Ratio")
Ratio.Divide(CorrectedDplusYields)

UncorrectedRatio = DsRawYields.Clone("hUncorrectedRatio")
UncorrectedRatio.SetTitle(";p_{T} (GeV/c);D_{s}^{+}/D^{+} Uncorrected ratio")
UncorrectedRatio.Divide(DplusRawYields)
UncorrectedRatio.Scale(2.69e-3/2.21e-2)

outputFile.cd()
Ratio.Write()
UncorrectedRatio.Write()
outputFile.Close()