import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ROOT import TFile, TGraph, TGraphErrors, TGraphAsymmErrors, TH1D, TH2D, TCanvas, gStyle, kGray, kRed, kBlue, kGreen, kMagenta, kOrange, kCyan, kYellow, kBlack, kWhite, kAzure
from hipe4ml import plot_utils
from hipe4ml.model_handler import ModelHandler
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from DfUtils import apply_model_in_batches
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
from particle import Particle


Variables = ["fCpa", "fCpaXY", "fDecayLength", "fDecayLengthXY", "fImpactParameterXY", "fAbsCos3PiK", "fChi2PCA",
                       "fNSigTpcPi0", "fNSigTpcKa0", "fNSigTpcPi1", "fNSigTpcKa1", "fNSigTpcPi2", "fNSigTpcKa2", "fPtProng0", "fPtProng1", "fPtProng2"] 

Ranges = {
            "fCpa": (0.95,1.),
            "fCpaXY": (0.95,1.),
            "fDecayLength": (0.,0.2),
            "fDecayLengthXY": (0.,0.2),
            "fImpactParameterXY": (-0.01,0.01),
            "fAbsCos3PiK": (0.,1.),
            "fChi2PCA": (0.,4.),
            "fNSigTpcPi0": (-10.,10.),
            "fNSigTpcKa0": (-10.,10.),
            "fNSigTpcPi1": (-10.,10.),
            "fNSigTpcKa1": (-10.,10.),
            "fNSigTpcPi2": (-10.,10.),
            "fNSigTpcKa2": (-10.,10.),
            "fPtProng0": (0., 5.),
            "fPtProng1": (0., 5.),
            "fPtProng2": (0., 5.)
            }

ptMins = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,8,12] 
ptMaxs = [1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,8,12,24]

ModelDirs = [
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt0_1.5/ModelHandler_pT_0_1.5.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt1.5_2/ModelHandler_pT_1.5_2.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt2_3/ModelHandler_pT_2_3.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt2_3/ModelHandler_pT_2_3.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt3_4/ModelHandler_pT_3_4.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt3_4/ModelHandler_pT_3_4.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt4_5/ModelHandler_pT_4_5.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt4_5/ModelHandler_pT_4_5.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt5_6/ModelHandler_pT_5_6.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt5_6/ModelHandler_pT_5_6.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt6_8/ModelHandler_pT_6_8.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt8_12/ModelHandler_pT_8_12.pickle",
                "/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt12_24/ModelHandler_pT_12_24.pickle"
            ] 

BkgScoreMin = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
BkgScoreMax = [0.1,0.1,0.3,0.45,0.25,0.90,0.3,0.60,0.55,0.7,0.60,0.75,0.75,0.70]
PromptScoreMin = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
PromptScoreMax = [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]


DataFile = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train191566/LHC22o.parquet"
PromptDsFile = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train194380/LHC23l2a_PromptDs.parquet"
NonPromptDsFile = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train194380/LHC23l2a_NonPromptDs.parquet"
PromptDplusFile = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train194380/LHC23l2a_PromptDplus.parquet"
NonPromptDplusFile = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train194380/LHC23l2a_NonPromptDplus.parquet"

PromptFractionDsFile = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/theory_driven/FD_fraction_Ds_newMC.root"
PromptFractionDplusFile = "/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/theory_driven/FD_fraction_Dplus_newMC.root"

TemplateFile = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/DplusForTemplateHistos_Train165702.root"

nSigma = 3 # Number of sigmas to consider for the signal region

infilePromptFracDs = TFile(PromptFractionDsFile)
gPromptDs = infilePromptFracDs.Get("gfraction")
infilePromptFracDs.Close()

infilePromptFracDplus = TFile(PromptFractionDplusFile)
gPromptDplus = infilePromptFracDplus.Get("gfraction")
infilePromptFracDplus.Close()   
VarsToSave = Variables + ["fPt", "fM"]

for iPt, (ptMin, ptMax, ModelDir) in enumerate(zip(ptMins, ptMaxs, ModelDirs)):
    ModelHandl = ModelHandler()
    ModelHandl.load_model_handler(ModelDir)
    DataDf = apply_model_in_batches(ModelHandl, VarsToSave, DataFile, f"fPt > {ptMin} and fPt < {ptMax}")
    DataDf = DataDf.query(f"{BkgScoreMin[iPt]} < ML_output_Bkg < {BkgScoreMax[iPt]} and {PromptScoreMin[iPt]} < ML_output_Prompt < {PromptScoreMax[iPt]}")
    PromptDsDf = apply_model_in_batches(ModelHandl, VarsToSave, PromptDsFile, f"fPt > {ptMin} and fPt < {ptMax}")
    PromptDsDf = PromptDsDf.query(f"{BkgScoreMin[iPt]} < ML_output_Bkg < {BkgScoreMax[iPt]} and {PromptScoreMin[iPt]} < ML_output_Prompt < {PromptScoreMax[iPt]}")
    NonPromptDsDf = apply_model_in_batches(ModelHandl, VarsToSave, NonPromptDsFile, f"fPt > {ptMin} and fPt < {ptMax}")
    NonPromptDsDf = NonPromptDsDf.query(f"{BkgScoreMin[iPt]} < ML_output_Bkg < {BkgScoreMax[iPt]} and {PromptScoreMin[iPt]} < ML_output_Prompt < {PromptScoreMax[iPt]}")
    PromptDplusDf = apply_model_in_batches(ModelHandl, VarsToSave, PromptDplusFile, f"fPt > {ptMin} and fPt < {ptMax}")
    PromptDplusDf = PromptDplusDf.query(f"{BkgScoreMin[iPt]} < ML_output_Bkg < {BkgScoreMax[iPt]} and {PromptScoreMin[iPt]} < ML_output_Prompt < {PromptScoreMax[iPt]}")
    NonPromptDplusDf = apply_model_in_batches(ModelHandl, VarsToSave, NonPromptDplusFile, f"fPt > {ptMin} and fPt < {ptMax}")
    NonPromptDplusDf = NonPromptDplusDf.query(f"{BkgScoreMin[iPt]} < ML_output_Bkg < {BkgScoreMax[iPt]} and {PromptScoreMin[iPt]} < ML_output_Prompt < {PromptScoreMax[iPt]}")

    promptFracDs = gPromptDs.GetY()[iPt]
    promptFracDplus = gPromptDplus.GetY()[iPt]

    # Fit the data
    data_hdl = DataHandler(data=DataDf, var_name="fM",
                            limits=[1.75,2.1], rebin=2)

    fitter = F2MassFitter(data_hdl, name_signal_pdf=["gaussian", "gaussian"],
                        name_background_pdf=["chebpol2"],
                        name=f"ds_pt{ptMin*10:.0f}_{ptMax*10:.0f}", chi2_loss=False,
                        verbosity=1, tol=1.e-1)

    # bkg initialisation
    fitter.set_background_initpar(0, "c0", 0.6)
    fitter.set_background_initpar(0, "c1", -0.2)
    fitter.set_background_initpar(0, "c2", 0.01)

    # signals initialisation
    fitter.set_particle_mass(0, pdg_id=431, limits=[Particle.from_pdgid(431).mass*0.99e-3,
                                                    Particle.from_pdgid(431).mass*1.01e-3])
    fitter.set_signal_initpar(0, "sigma", 0.008, limits=[0.002, 0.030])
    fitter.set_signal_initpar(0, "frac", 0.1, limits=[0., 1.])
    fitter.set_particle_mass(1, pdg_id=411,
                             limits=[Particle.from_pdgid(411).mass*0.99e-3,
                                     Particle.from_pdgid(411).mass*1.01e-3])
    fitter.set_signal_initpar(1, "sigma", 0.006, limits=[0.002, 0.030])
    fitter.set_signal_initpar(1, "frac", 0.1, limits=[0., 1.])

    fit_result = fitter.mass_zfit()

    meanDs = fitter.get_mass(0)[0]
    sigmaDs = fitter.get_sigma(0)[0]
    meanDplus = fitter.get_mass(1)[0]
    sigmaDplus = fitter.get_sigma(1)[0]

    bkgDs = fitter.get_background(0)[0]
    bkgDplus = fitter.get_background(1)[0]

    # Get signal and background distributions
    bkgDf = DataDf.query(f"fM > {meanDs + nSigma*sigmaDs} or fM < {meanDplus - nSigma*sigmaDplus}")
    peakDsDf = DataDf.query(f"{meanDs - nSigma*sigmaDs} < fM < {meanDs + nSigma*sigmaDs}")
    peakDplusDf = DataDf.query(f"{meanDplus - nSigma*sigmaDplus} < fM < {meanDplus + nSigma*sigmaDplus}")

    nRows = 4
    nCols = len(Variables)//nRows if len(Variables)%nRows == 0 else len(Variables)//nRows + 1

    figDs, axsDs = plt.subplots(nRows, nCols, figsize=(20, 20))   
    figDplus, axsDplus = plt.subplots(nRows, nCols, figsize=(20, 20))   
    figDsVsDplus, axsDsVsDplus = plt.subplots(nRows, nCols, figsize=(20, 20))
    for iVar, var in enumerate(Variables):
        range = Ranges[var]
        peakDsHist, edges = np.histogram(peakDsDf[var].values, bins=100, range=range)
        peakDplusHist, _ = np.histogram(peakDplusDf[var].values, bins=100, range=range)
        bkgHist, _ = np.histogram(bkgDf[var].values, bins=100, range=range)
        bkgHist = bkgHist/len(bkgDf)
        signalDsHist = peakDsHist - bkgHist*bkgDs
        signalDsHist = signalDsHist/signalDsHist.sum()
        signalDplusHist = peakDplusHist - bkgHist*bkgDplus
        signalDplusHist = signalDplusHist/signalDplusHist.sum()

        PromptDsHist, _ = np.histogram(PromptDsDf[var].values, bins=100, range=range)
        PromptDsHist = PromptDsHist/len(PromptDsDf)
        NonPromptDsHist, _ = np.histogram(NonPromptDsDf[var].values, bins=100, range=range)
        NonPromptDsHist = NonPromptDsHist/len(NonPromptDsDf)

        PromptDplusHist, _ = np.histogram(PromptDplusDf[var].values, bins=100, range=range)
        PromptDplusHist = PromptDplusHist/len(PromptDplusDf)
        NonPromptDplusHist, _ = np.histogram(NonPromptDplusDf[var].values, bins=100, range=range)
        NonePromptDplusHist = NonPromptDplusHist/len(NonPromptDplusDf)

        McHybridDs = PromptDsHist + NonPromptDsHist*(1/promptFracDs - 1)
        McHibridDs = McHybridDs/McHybridDs
        McHybridDplus = PromptDplusHist + NonPromptDplusHist*(1/promptFracDplus - 1)
        McHybridDplus = McHybridDplus/McHybridDplus.sum()

        # Plot the signal and background distributions
        axsDs[iVar//nRows, iVar%nRows].bar(edges[:-1], signalDsHist, width=np.diff(edges), align='edge', alpha=0.5, label="Signal Ds")
        axsDs[iVar//nRows, iVar%nRows].bar(edges[:-1], McHybridDs, width=np.diff(edges), align='edge', alpha=0.5, label="MC Hybrid Ds")
        axsDs[iVar//nRows, iVar%nRows].set_title(var)
        axsDs[iVar//nRows, iVar%nRows].legend()

        axsDplus[iVar//nRows, iVar%nRows].bar(edges[:-1], signalDplusHist, width=np.diff(edges), align='edge', alpha=0.5, label="Signal Dplus")
        axsDplus[iVar//nRows, iVar%nRows].bar(edges[:-1], McHybridDplus, width=np.diff(edges), align='edge', alpha=0.5, label="MC Hybrid Dplus")
        axsDplus[iVar//nRows, iVar%nRows].set_title(var)
        axsDplus[iVar//nRows, iVar%nRows].legend()
       
        axsDsVsDplus[iVar//nRows, iVar%nRows].bar(edges[:-1], signalDsHist, width=np.diff(edges), align='edge', alpha=0.5, label="Signal Ds")
        axsDsVsDplus[iVar//nRows, iVar%nRows].bar(edges[:-1], signalDplusHist, width=np.diff(edges), align='edge', alpha=0.5, label="Signal Dplus")
        axsDsVsDplus[iVar//nRows, iVar%nRows].set_title(var)
        axsDsVsDplus[iVar//nRows, iVar%nRows].legend()
    plt.tight_layout()
    figDs.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Tests/CheckSignalDistribution/SignalDistributionsDs_{ptMin}_{ptMax}.png")
    figDplus.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Tests/CheckSignalDistribution/SignalDistributionsDplus_{ptMin}_{ptMax}.png")
    figDsVsDplus.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Tests/CheckSignalDistribution/DsVsDplus/SignalDistributionsDsVsDplus_{ptMin}_{ptMax}.png")

        


