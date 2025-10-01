import pandas as pd
import ROOT

FileNamePP="TAMU_pp5dot02TeV_DsOverD0.txt"
FileNamePbPb010="TAMU_PbPb010_5dot02TeV_DsOverD0.txt"
FileNamePbPb3050="TAMU_PbPb3050_5dot02TeV_DsOverD0.txt"

dfpp = pd.read_csv(FileNamePP, sep=' ', header=None, names=['pT', 'DsOverD0Min', 'DsOverD0Max'])
dfPbPb010 = pd.read_csv(FileNamePbPb010, sep=' ', header=None, names=['pT', 'DsOverD0Min', 'DsOverD0Max'])
dfPbPb3050 = pd.read_csv(FileNamePbPb3050, sep=' ', header=None, names=['pT', 'DsOverD0Min', 'DsOverD0Max'])

DpOverD0 = 0.4391 # from https://arxiv.org/abs/1902.08889 (RQM, 170 MeV), assuming independent of everything (pt, cent, etc)

gTAMU_pp = ROOT.TGraphAsymmErrors(len(dfpp))
gTAMU_pp.SetName("gTAMU_pp")
gTAMU_PbPb010 = ROOT.TGraphAsymmErrors(len(dfPbPb010))
gTAMU_PbPb010.SetName("gTAMU_PbPb010")
gTAMU_PbPb3050 = ROOT.TGraphAsymmErrors(len(dfPbPb3050))
gTAMU_PbPb3050.SetName("gTAMU_PbPb3050")

for ipt, (pt, ratioMin, ratioMax) in enumerate(
    zip(dfpp["pT"].to_numpy(), dfpp["DsOverD0Min"].to_numpy(), dfpp["DsOverD0Max"].to_numpy())):
    ratio = (ratioMin + ratioMax) / 2
    if ratioMin < ratioMax:
        ratioUncLow = ratio - ratioMin
        ratioUncHigh = ratioMax - ratio
    else:
        ratioUncLow = ratioMax - ratio
        ratioUncHigh = ratio - ratioMin
    gTAMU_pp.SetPoint(ipt, pt, ratio / DpOverD0)
    gTAMU_pp.SetPointError(ipt, 0., 0., ratioUncLow, ratioUncHigh)

for ipt, (pt, ratioMin, ratioMax) in enumerate(
    zip(dfPbPb010["pT"].to_numpy(), dfPbPb010["DsOverD0Min"].to_numpy(), dfPbPb010["DsOverD0Max"].to_numpy())):
    ratio = (ratioMin + ratioMax) / 2
    if ratioMin < ratioMax:
        ratioUncLow = ratio - ratioMin
        ratioUncHigh = ratioMax - ratio
    else:
        ratioUncLow = ratioMax - ratio
        ratioUncHigh = ratio - ratioMin
    gTAMU_PbPb010.SetPoint(ipt, pt, ratio / DpOverD0)
    gTAMU_PbPb010.SetPointError(ipt, 0., 0., ratioUncLow, ratioUncHigh)

for ipt, (pt, ratioMin, ratioMax) in enumerate(
    zip(dfPbPb3050["pT"].to_numpy(), dfPbPb3050["DsOverD0Min"].to_numpy(), dfPbPb3050["DsOverD0Max"].to_numpy())):
    ratio = (ratioMin + ratioMax) / 2
    if ratioMin < ratioMax:
        ratioUncLow = ratio - ratioMin
        ratioUncHigh = ratioMax - ratio
    else:
        ratioUncLow = ratioMax - ratio
        ratioUncHigh = ratio - ratioMin
    gTAMU_PbPb3050.SetPoint(ipt, pt, ratio / DpOverD0)
    gTAMU_PbPb3050.SetPointError(ipt, 0., 0., ratioUncLow, ratioUncHigh)

outFile = ROOT.TFile("TAMU_allSyst_5dot02TeV_DsOverDp.root", "recreate")
gTAMU_pp.Write()
gTAMU_PbPb010.Write()
gTAMU_PbPb3050.Write()
outFile.Close()
