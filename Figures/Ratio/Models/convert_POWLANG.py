import pandas as pd
import numpy as np
import ROOT 

HTLFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/Models/POWLANG-pp-13-TeV-HTL.txt"
lQCDFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/Models/POWLANG-pp-13-TeV-lQCD.txt"

dfHTL = pd.read_csv(HTLFileName, sep='\t')
dflQCD = pd.read_csv(lQCDFileName, sep='\t')

print(dflQCD.columns)

dfHTL["Ds/Dp"] = dfHTL["Ds/D0"]/dfHTL["D+/D0"]
dflQCD["Ds/Dp"] = dflQCD["Ds/D0"]/dflQCD["D+/D0"]

ptEdges = dfHTL["PtMin"].to_list()
ptEdges.append(dfHTL["PtMax"].iloc[-1])

hHTL = ROOT.TH1F("hHTL", "hHTL", len(ptEdges)-1, np.asarray(ptEdges,"d"))
hlQCD = ROOT.TH1F("hlQCD", "hlQCD", len(ptEdges)-1, np.asarray(ptEdges,"d"))

for i in range(len(ptEdges)-1):
    hHTL.SetBinContent(i+1, dfHTL["Ds/Dp"].iloc[i])
    hHTL.SetBinError(i+1, 0)
    hlQCD.SetBinContent(i+1, dflQCD["Ds/Dp"].iloc[i])
    hlQCD.SetBinError(i+1, 0)
    
outFile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/Models/POWLANG-pp-13-TeV.root", "RECREATE")
hHTL.Write()
hlQCD.Write()
outFile.Close()