import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ROOT import TFile, THnSparseF

# Open the file
file = TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train208051/LHC24d3a_AnalysisResults.root")
# Get the THnSparse
sparse = file.Get("hf-task-ds/MC/Ds/Prompt/hSparseMass")

# Get the number of bins in each dimension
nBins = [sparse.GetAxis(i).GetNbins() for i in range(sparse.GetNdimensions())]
# Get the bin limits
binLimits = [sparse.GetAxis(i).GetXbins() for i in range(sparse.GetNdimensions())]

# Project the THnSparse into a TH1D for each dimension
projections = [sparse.Projection(i) for i in range(sparse.GetNdimensions())]

outfilename = "/home/fchinu/Run3/Ds_pp_13TeV/TestO2/StudySparse.root"
outfile = TFile(outfilename, "RECREATE")
for i, proj in enumerate(projections):
    proj.Write()
outfile.Close()