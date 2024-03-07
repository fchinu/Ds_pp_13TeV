import ROOT
from ROOT import AliHFInvMassFitter, AliVertexingHFUtils
from ROOT import TH1D, TH2F, TFile, TCanvas, TMath, TObject, TH1F, TKey, TIter, TDatabasePDG, gStyle, kDarkBodyRadiator
import pandas as pd
import numpy as np
import ctypes
import math
import itertools
from hipe4ml import plot_utils
import matplotlib.pyplot as plt

inputFileNames = [             
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt0_2/Data_pT_0_2_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt2_2.5/Data_pT_2_2.5_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt2.5_3/Data_pT_2.5_3_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt3_3.5/Data_pT_3_3.5_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt3.5_4/Data_pT_3.5_4_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt4_4.5/Data_pT_4_4.5_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt4.5_5/Data_pT_4.5_5_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt5_5.5/Data_pT_5_5.5_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt5.5_6/Data_pT_5.5_6_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt6_8/Data_pT_6_8_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt8_12/Data_pT_8_12_ModelApplied.parquet.gzip",
            "/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt12_24/Data_pT_12_24_ModelApplied.parquet.gzip"
                    ]
outputFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Optimization/Significance_ptDiff.root"
pTmin = [0,2,2.5,3,3.5,4,4.5,5,5.5,6,8,12] # list 
pTmax = [2,2.5,3,3.5,4,4.5,5,5.5,6,8,12,24] # list 
BDT_prompt_mins = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
BDT_prompt_maxs = [0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]
stepPrompt = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
BDT_bkg_mins = [0.001,0.001,0.005,0.005,0.005,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]
BDT_bkg_maxs = [0.01,0.01,0.02,0.04,0.04,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.4]
stepBkg = [0.001, 0.001, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005]
significanceLimits = [50.,50.,50.,50.,50.,50.,50.,50.,50.,50.,50.,50.,50.]
rebins = [4,4,4,4,4,4,4,4,4,4,4,4,4]
fitFunction = [2,2,0,0,0,0,0,0,0,0,0,0,0]

saveAllFits = True

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
open("/home/fchinu/Run3/Ds_pp_13TeV/Optimization/BestCuts.txt", "w").close()
if saveAllFits:
    outputFileFits = TFile.Open(outputFileName.replace(".root","_fits.root"),"RECREATE")
outputFile = TFile.Open(outputFileName,"RECREATE")
for iPt, (minPt,maxPt) in enumerate(PtEdges):
    print("\033[91m" + f"{minPt} < pT < {maxPt}" + "\033[0m")
    if saveAllFits:
        outputFileFits.mkdir(f"{minPt}-{maxPt}")
    outputFile.cd()
    df = pd.read_parquet(inputFileNames[iPt], engine = "fastparquet")
    
    histos = []
    fitters = []
    significance = []
    significanceerr = []
    signal = []
    signalerr = []
    BDT_prompt_cuts = np.arange(BDT_prompt_mins[iPt],BDT_prompt_maxs[iPt],stepPrompt[iPt])
    BDT_bkg_cuts = np.arange(BDT_bkg_mins[iPt],BDT_bkg_maxs[iPt],stepBkg[iPt])

    #Draw the significance scan
    hSignificanceScan = TH2F(f"hSignificanceScan_{minPt}_{maxPt}",f"hSignificanceScan_{minPt}_{maxPt};ML output bkg;ML output prompt",len(BDT_bkg_cuts),BDT_bkg_mins[iPt]-stepBkg[iPt]/2,BDT_bkg_maxs[iPt]+stepBkg[iPt]/2,len(BDT_prompt_cuts),BDT_prompt_mins[iPt]-stepPrompt[iPt]/2,BDT_prompt_maxs[iPt]+stepPrompt[iPt]/2)
    #Loop over the histograms
    cuts = []
    for (BDT_prompt_cut, BDT_bkg_cut)  in itertools.product(BDT_prompt_cuts,BDT_bkg_cuts):
        cuts.append((BDT_prompt_cut, BDT_bkg_cut))
        dfSel = df.query(f"ML_output_Prompt > {BDT_prompt_cut} and ML_output_Bkg < {BDT_bkg_cut}")
        histos.append(TH1F(f"histo_{BDT_prompt_cut}_{BDT_bkg_cut}",f"histo_{BDT_prompt_cut}_{BDT_bkg_cut}",400,1.7,2.1))
        for mass in dfSel["fM"].to_numpy():
            histos[-1].Fill(mass)        
        histos[-1].Rebin(rebins[iPt])
        histos[-1].SetDirectory(0)
        histos[-1].SetTitle(";Mass (GeV/c^{2});"+f"Counts per {rebins[iPt]} " + "MeV/c^{2}")
        fitters.append(AliHFInvMassFitter(histos[-1], 1.7,2.1,fitFunction[iPt]))
        fitters[-1].SetInitialGaussianMean(1.96)
        fitters[-1].SetInitialGaussianSigma(0.01)
        fitters[-1].IncludeSecondGausPeak(1.865, False, 0.01, False)
        signif, signiferr = ctypes.c_double(), ctypes.c_double()
        sgn, sgnerr = ctypes.c_double(), ctypes.c_double()
        if(fitters[-1].MassFitter(False)):            
            fitters[-1].Significance(2,signif,signiferr)
            fitters[-1].Signal(2,sgn,sgnerr)
        significance.append(signif.value)
        significanceerr.append(signiferr.value)
        signal.append(sgn.value)
        signalerr.append(sgnerr.value)
        hSignificanceScan.SetBinContent(hSignificanceScan.FindBin(BDT_bkg_cut,BDT_prompt_cut),signif.value)
        if saveAllFits:
            outputFileFits.cd(f"{minPt}-{maxPt}")
            fitCanvas = TCanvas(f"Prompt > {BDT_prompt_cut} and Bkg < {BDT_bkg_cut}","c",800,600)
            fitCanvas.DrawFrame(1.7,0,2,1000,";Mass (GeV/c^{2});"+f"Counts per {1000*histos[-1].GetBinWidth(1)}" + "MeV/c^{2}")
            fitters[-1].DrawHere(fitCanvas,3,2)
            fitCanvas.Write(f"Prompt > {BDT_prompt_cut} and Bkg < {BDT_bkg_cut}")
            outputFile.cd()

    hSignificanceScan.SetDrawOption("colz")
    hSignificanceScan.Write()
    

    #Find the maximum significance
    maxSignificance = max(significance)
    maxSignificanceIndex = significance.index(maxSignificance)
    c = TCanvas(("canvas"+"_"+f"Pt{minPt:.1f}_{maxPt:.1f}"),"c",800,600)
    c.DrawFrame(1.7,0,2,1000,";Mass (GeV/c^{2});"+f"Counts per {1000*histos[maxSignificanceIndex].GetBinWidth(1)}" + "MeV/c^{2}")
    fitters[maxSignificanceIndex].DrawHere(c,3,2)
    c.Write()

    with open("/home/fchinu/Run3/Ds_pp_13TeV/Optimization/BestCuts.txt", "a") as text_file:
        text_file.write(str(cuts[maxSignificanceIndex][0]) + "\t" + str(cuts[maxSignificanceIndex][1]) + "\n")
    ##Draw distributions comparison of all candidates and selected candidates
    #BDT_prompt_max, BDT_bkg_max = list(itertools.product(BDT_prompt_cuts,BDT_bkg_cuts))[maxSignificanceIndex]
    #dfSel = df.query(f"ML_output_Prompt > {BDT_prompt_max} and ML_output_Bkg < {BDT_bkg_max}")
    #VarsToDraw = ['fCpa', 'fCpaXY', 'fDecayLength', 'fDecayLengthXY', 'fDeltaMassPhi', 'fImpactParameterXY', 'fAbsCos3PiK','fMaxNormalisedDeltaIP']
    #plot_utils.plot_distr([df, dfSel], column=VarsToDraw, bins=100, labels=['df','df_sel'],
    #                      alpha=0.3, log=True, grid=False, density=True)
    #plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    #plt.savefig(outputFileName.replace(".root","_Pt"+str(int(minPt))+"_"+str(int(maxPt))+"_maxSignifDistr.pdf"))
    #plt.close('all')


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
