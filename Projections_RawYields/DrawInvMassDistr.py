import pandas as pd
from ROOT import TFile, TH1D, TF1
import numpy as np

# Input Parquet file containing columns: fM, fPt, ML_output
parquet_file_path = "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt2_4/Data_pT_2_4_ModelApplied.parquet.gzip"

# Output ROOT file
output_root_file = TFile("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/Inv_mass_2_4.root", "RECREATE")

# Specify the pt interval and ML thresholds
pt_min, pt_max = 2, 4  # Choose your desired pt interval
ml_thresholds = np.linspace(0.9, 0.999, 99).tolist() # Define different ML thresholds
#ml_thresholds = [0.9] # Define different ML thresholds


# Create histograms for the specified pt interval and ML thresholds
histograms = {}
fits = {}

for ml_threshold in ml_thresholds:
    ml_key = f"ML_{ml_threshold}"
    histogram_name = f"invariant_mass_{pt_min}_{pt_max}_{ml_key}"
    histogram_title = f"Invariant Mass Distribution ({pt_min} <= fPt < {pt_max} GeV/c) - ML > {ml_threshold}"

    histogram = TH1D(histogram_name, histogram_title, 200, 1.75, 2.2)  # Customize the binning as needed
    histograms[ml_key] = histogram

# Read the Parquet file and fill histograms
df = pd.read_parquet(parquet_file_path)

for ml_threshold in ml_thresholds:
    ml_key = f"ML_{ml_threshold}"
    dfSel = df.query(f"ML_output > {ml_threshold}")
    histograms[ml_key].FillN(len(dfSel),np.asarray(dfSel["fM"],"d"),np.asarray([1]*len(dfSel),"d"))
    del dfSel

    fits[ml_key] = TF1(f"fit_{ml_key}", "expo + gaus(2) + [5]", 1.9, 2.05)
    fits[ml_key].SetParameters(1,1,1,20,1.96,0.2)
    fits[ml_key].SetParLimits(3,1.9, 2.05)
    histograms[ml_key].Fit(fits[ml_key], "R")

# Save histograms to the output ROOT file
output_root_file.cd()
for histogram in histograms.values():
    histogram.Write()

# Close the output ROOT file
output_root_file.Close()

print("Histograms have been created and saved to output.root for the specified pt interval.")
