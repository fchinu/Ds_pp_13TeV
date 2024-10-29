import pandas as pd
import numpy as np
import ROOT

df = pd.read_csv('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/Models/Ds_Dplus_coalfragm.dat', sep=' ', header=None, names=['pT', 'DsOverDplus'])

gCatania = ROOT.TGraphErrors(len(df))
for i, row in df.iterrows():
    gCatania.SetPoint(i, row['pT'], row['DsOverDplus'])
    gCatania.SetPointError(i, 0.1, 0.0)

outFile = ROOT.TFile('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/Models/CATANIA.root', 'RECREATE')
gCatania.Write("gCatania")
outFile.Close()
