import ROOT
from ROOT import AliHFInvMassFitter, AliVertexingHFUtils
from ROOT import TH1D, TH2F, TFile, TCanvas, TMath, TObject, TH1F, TKey, TIter, TDatabasePDG, gStyle, kDarkBodyRadiator
import pandas as pd
import numpy as np
import ctypes
import math
import itertools

inputFileNames = [  "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt2_4/Data_pT_2_4_ModelApplied.parquet.gzip",
                    "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt4_6/Data_pT_4_6_ModelApplied.parquet.gzip",
                    "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt6_8/Data_pT_6_8_ModelApplied.parquet.gzip",
                    "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt8_12/Data_pT_8_12_ModelApplied.parquet.gzip",
                    "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt12_50/Data_pT_12_50_ModelApplied.parquet.gzip"
                    ]
outputFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Optimization/Significance_ptDiff.root"
pTmin = [2.,4.,6.,8.,12.]
pTmax = [4.,6.,8.,12.,50.]
BDT_prompt_mins = [0.2,0.2,0.1,0.2,0.2,0.2,0.2]
BDT_prompt_maxs = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
BDT_bkg_mins = [0.01,0.01,0.01,0.01,0.01,0.01,0.01]
BDT_bkg_maxs = [0.25,0.25,0.25,0.25,0.25,0.25,0.25]
stepBDT = 0.01
significanceLimits = [35.,35.,35.,35.,35.,35.,35.]
rebins = [4,4,4,6,8,8,8]

PtEdges = [[a, b] for a, b in zip(pTmin, pTmax)]
ptBins = [i for i in pTmin]
ptBins.append(pTmax[-1])
hSignificancepT = TH1D("hSignificancepT","hSignificancepT",len(pTmin),np.asarray(ptBins,"d"))
hSigma = TH1D("hSigma","hSigma",len(pTmin),np.asarray(ptBins,"d"))
hMean = TH1D("hMean","hMean",len(pTmin),np.asarray(ptBins,"d"))
hSignal = TH1D("hSignal","hSignal",len(pTmin),np.asarray(ptBins,"d"))
hSigmaSecPeak = TH1D("hSigmaSecPeak","hSigmaSecPeak",len(pTmin),np.asarray(ptBins,"d"))

massDplus = TDatabasePDG.Instance().GetParticle(411).Mass()
#Get Events
#hEventsFile = TFile(AnalysisResultsName)
#hEvents= hEventsFile.Get("event-selection-task/hColCounterAcc")
#nEv = hEvents.GetEntries()
#hEventsFile.Close()
#hEv = TH1D("hEv","hEv",1,0,1)
#hEv.SetBinContent(1,nEv)


#Start the loop
gStyle.SetPalette(kDarkBodyRadiator, 0)

outputFile = TFile.Open(outputFileName,"RECREATE")
for iPt, (minPt,maxPt) in enumerate(PtEdges):
    df = pd.read_parquet(inputFileNames[iPt], engine = "fastparquet")
    
    histos = []
    fitters = []
    significance = []
    significanceerr = []
    signal = []
    signalerr = []
    BDT_prompt_cuts = np.arange(BDT_prompt_mins[iPt],BDT_prompt_maxs[iPt],stepBDT)
    BDT_bkg_cuts = np.arange(BDT_bkg_mins[iPt],BDT_bkg_maxs[iPt],stepBDT)

    #Draw the significance scan
    hSignificanceScan = TH2F(f"hSignificanceScan_{minPt}_{maxPt}",f"hSignificanceScan_{minPt}_{maxPt};ML output bkg;ML output prompt",len(BDT_bkg_cuts),BDT_bkg_mins[iPt]-stepBDT/2,BDT_bkg_maxs[iPt]+stepBDT/2,len(BDT_prompt_cuts),BDT_prompt_mins[iPt]-stepBDT/2,BDT_prompt_maxs[iPt]+stepBDT/2)
    #Loop over the histograms
    for (BDT_prompt_cut, BDT_bkg_cut)  in itertools.product(BDT_prompt_cuts,BDT_bkg_cuts):
        dfSel = df.query(f"ML_output_Prompt > {BDT_prompt_cut} and ML_output_Bkg < {BDT_bkg_cut}")
        histos.append(TH1F(f"histo_{BDT_prompt_cut}_{BDT_bkg_cut}",f"histo_{BDT_prompt_cut}_{BDT_bkg_cut}",200,1.7,2.1))
        for mass in dfSel["fM"].to_numpy():
            histos[-1].Fill(mass)        
        histos[-1].Rebin(rebins[iPt])
        histos[-1].SetDirectory(0)
        histos[-1].SetTitle(";Mass (GeV/c^{2});"+f"Counts per {rebins[iPt]} " + "MeV/c^{2}")
        fitters.append(AliHFInvMassFitter(histos[-1], 1.7,2.1))
        fitters[-1].SetInitialGaussianMean(1.96)
        fitters[-1].SetInitialGaussianSigma(0.01)
        fitters[-1].IncludeSecondGausPeak(1.86, False, 0.02, False)
        signif, signiferr = ctypes.c_double(), ctypes.c_double()
        sgn, sgnerr = ctypes.c_double(), ctypes.c_double()
        if(fitters[-1].MassFitter(False)):            
            fitters[-1].Significance(3,signif,signiferr)
            fitters[-1].Signal(3,sgn,sgnerr)
        significance.append(signif.value)
        significanceerr.append(signiferr.value)
        signal.append(sgn.value)
        signalerr.append(sgnerr.value)
        hSignificanceScan.SetBinContent(hSignificanceScan.FindBin(BDT_bkg_cut,BDT_prompt_cut),signif.value)
    outputFile.cd()

    hSignificanceScan.SetDrawOption("colz")
    hSignificanceScan.Write()
    

    #Find the maximum significance
    maxSignificance = max(significance)
    maxSignificanceIndex = significance.index(maxSignificance)
    c = TCanvas(("canvas"+"_"+"Pt"+str(int(minPt))+"_"+str(int(maxPt))),"c",800,600)
    c.DrawFrame(1.7,0,2,1000,";Mass (GeV/c^{2});"+f"Counts per {1000*histos[maxSignificanceIndex].GetBinWidth(1)}" + "MeV/c^{2}")
    fitters[maxSignificanceIndex].DrawHere(c,3,2)
    c.Write()


    hSignificancepT.Fill((minPt+maxPt)/2,maxSignificance)
    hSignificancepT.SetBinError(iPt+1,significanceerr[maxSignificanceIndex])
    hSigma.Fill((minPt+maxPt)/2,fitters[maxSignificanceIndex].GetSigma())
    hSigma.SetBinError(iPt+1,fitters[maxSignificanceIndex].GetSigmaUncertainty())
    hSignal.Fill((minPt+maxPt)/2,signal[maxSignificanceIndex])
    hMean.Fill((minPt+maxPt)/2,fitters[maxSignificanceIndex].GetMean())
    hMean.SetBinError(iPt+1,fitters[maxSignificanceIndex].GetMeanUncertainty())
    if math.isnan(signalerr[maxSignificanceIndex]):
        hSignal.SetBinError(iPt+1,0)
    else:
        hSignal.SetBinError(iPt+1,signalerr[maxSignificanceIndex])
    fTotFunc = fitters[maxSignificanceIndex].GetMassFunc()
    hSigmaSecPeak.Fill((minPt+maxPt)/2,fTotFunc.GetParameter(fTotFunc.GetNpar()-1))
    hSigmaSecPeak.SetBinError(iPt+1,fTotFunc.GetParError(fTotFunc.GetNpar()-1))
    

#Write events and significance over sqrt(events) and signal over events
#hSignificanceOverNev = hSignificancepT.Clone(f"hSignificanceOverNev")
#hSignificanceOverNev.Scale(1./TMath.Sqrt(nEv))
#hSignalOverNev = hSignal.Clone(f"hSignalOverNev")
#hSignalOverNev.Scale(1./nEv)
#hSignificanceOverNev.Write()
#hSignalOverNev.Write()
#hEv.Write()
hSignal.Write()
hSigmaSecPeak.Write()

hSignificancepT.Write()
hSigma.Write()
hMean.Write()
outputFile.Close()
