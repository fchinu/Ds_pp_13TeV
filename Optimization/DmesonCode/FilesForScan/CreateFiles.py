from ROOT import TFile, TCanvas, TMath, TObject, TH1D, TKey, TIter, TDatabasePDG, gStyle, kDarkBodyRadiator
import pandas as pd
import numpy as np

def CreateAccFiles(outputname):
    outfile = TFile.Open(outputname,"RECREATE")
    hPtGenAcc = TH1D("hPtGenAcc","hPtGenAcc",1,0,50)
    hPtGenLimAcc = TH1D("hPtGenLimAcc","hPtGenLimAcc",1,0,50)
    hPtGenAcc.SetBinContent(1,1)
    hPtGenLimAcc.SetBinContent(1,1)
    hPtGenAcc.Write()
    hPtGenLimAcc.Write()
    outfile.Close()

def CreatePreselEffFile(outputname):
    outfile = TFile.Open(outputname,"RECREATE")
    hPtGenAcc = TH1D("hEffPrompt","hEffPrompt",1,0,50)
    hPtGenLimAcc = TH1D("hEffFD","hEffFD",1,0,50)
    hPtGenAcc.SetBinContent(1,1)
    hPtGenLimAcc.SetBinContent(1,1)
    hPtGenAcc.Write()
    hPtGenLimAcc.Write()
    outfile.Close()

if __name__ == "__main__":
    CreateAccFiles("/home/fchinu/Run3/Ds_pp_13TeV/Optimization/FilesForScan/AcceptanceAtOneFile.root")
    CreatePreselEffFile("/home/fchinu/Run3/Ds_pp_13TeV/Optimization/FilesForScan/PreselEffAtOneFile.root")
