from ROOT import TFile, TCanvas, TH1F
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Calculate Ds+/D+ ratio')
parser.add_argument('-r', '--raw_yields_file', metavar='raw_yields_file', type=str, help='Path to the raw yields file', default="/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/Data/RawYields_DataFlareFly.root")
parser.add_argument('-e', '--efficiency_file', type=str, help='Path to the efficiency file', default="/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/Efficiencies.root")
parser.add_argument('-o', '--output_file', type=str, help='Path to the output file', default="/home/fchinu/Run3/Ds_pp_13TeV/Ratios/DsOverDplusFlareFly.root")

args = parser.parse_args()

RawYieldsFileName = args.raw_yields_file
EffFileName = args.efficiency_file
outputFileName = args.output_file

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