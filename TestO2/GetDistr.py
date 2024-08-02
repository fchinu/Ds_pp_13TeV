import ROOT
import sys
sys.path.append("Run3/ThesisUtils")
from DfUtils import read_parquet_in_batches

# Open the file
df = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt2_2.5/Data_pT_2_2.5_ModelApplied.parquet.gzip")

hM = ROOT.TH1F("hM", "hM", 400, 1.7, 2.1)
hM.SetDirectory(0)
for mass in df['fM']:
    hM.Fill(mass)

hBkg = ROOT.TH1F("Bkg_score", "Bkg_score", 1000, 0, 1)
hBkg.SetDirectory(0)
for mass in df['ML_output_Bkg']:
    hBkg.Fill(mass)

hPrompt = ROOT.TH1F("Prompt_score", "Prompt_score", 100, 0, 1)
hPrompt.SetDirectory(0)
for mass in df['ML_output_Prompt']:
    hPrompt.Fill(mass)

outputFile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/TestO2/massTree.root", "RECREATE")
hM.Write()
hBkg.Write()
hPrompt.Write()
outputFile.Close()
