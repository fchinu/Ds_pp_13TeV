import pandas as pd
import numpy as np
import ROOT

df = pd.read_csv('Catania_Ds_Dplus_pp13dot6TeV_coalfragm.dat', sep=' ', header=None, names=['pT', 'DsOverDplus'])

gCatania = ROOT.TGraphErrors(len(df))
for i, row in df.iterrows():
    gCatania.SetPoint(i, row['pT'], row['DsOverDplus'])
    gCatania.SetPointError(i, 0.1, 0.0)

outFile = ROOT.TFile('CATANIA.root', 'RECREATE')
gCatania.Write("gCatania")
outFile.Close()
