import ROOT 
import numpy as np
import pandas as pd
import yaml
from hipe4ml.model_handler import ModelHandler
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from DfUtils import apply_model_in_batches

ptMins = [0.5, 1., 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 12] 
ptMaxs = [1, 1.5, 2., 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 12, 24]
ptEdges = ptMins + [ptMaxs[-1]]

infileFONLL = ROOT.TFile.Open("/home/fchinu/DmesonAnalysis/models/fonll/feeddown/DmesonLcPredictions_13TeV_y05_FFee_BRPDG_SepContr_PDG2020.root")
infilePYTHIA = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/Models/PYTHIA8_Ds_over_Dplus_pp13dot6TeV.root")
infileAnalysisResults = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train213175/LHC24d3a_AnalysisResults.root")
infileRecDs = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train213175/LHC24d3a_PromptDs.parquet"
infileRecDplus = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train213175/LHC24d3a_PromptDplus.parquet"
infileCutset = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/cutset_pp13TeV_binary.yml"
infileConfig = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/PtShape/Config_Efficiency_newestMC.yaml"

with open(infileConfig, 'r') as f:
    config = yaml.safe_load(f)

hFONLLDs = infileFONLL.Get("hDsPhipitoKkpipred_central")
hFONLLDplus = infileFONLL.Get("hDpluskpipipred_central")

# Remove the BR
hFONLLDs.Scale(1/0.0221)
hFONLLDplus.Scale(1/0.0938)

# Set Range and rebin
hFONLLDs = hFONLLDs.Rebin(len(np.arange(0.5,24,0.1)), "hFONLLDs", np.asarray(np.arange(0.5,24.1,0.1), "d"))
hFONLLDs.SetDirectory(0)
hFONLLDplus = hFONLLDplus.Rebin(len(np.arange(0.5,24,0.1)), "hFONLLDplus", np.asarray(np.arange(0.5,24.1,0.1), "d"))
hFONLLDplus.SetDirectory(0)

infileFONLL.Close()

hPYTHIADs = infilePYTHIA.Get("hist_ds_prompt_CRMode2")
hPYTHIADplus = infilePYTHIA.Get("hist_dplus_prompt_CRMode2")

# Set Range
hPYTHIADs = hPYTHIADs.Rebin(len(np.arange(0.5,24,0.1)), "hPYTHIADs", np.asarray(np.arange(0.5,24.1,0.1), "d"))
hPYTHIADs.SetDirectory(0)
hPYTHIADplus = hPYTHIADplus.Rebin(len(np.arange(0.5,24,0.1)), "hPYTHIADplus", np.asarray(np.arange(0.5,24.1,0.1), "d"))
hPYTHIADplus.SetDirectory(0)

infilePYTHIA.Close()

hGenDs = infileAnalysisResults.Get("hf-task-ds/MC/Ds/Prompt/hPtGen")
hGenDplus = infileAnalysisResults.Get("hf-task-ds/MC/Dplus/Prompt/hPtGen")

hGenDs.SetName("hGenDs")
hGenDs.SetDirectory(0)
hGenDplus.SetName("hGenDplus")
hGenDplus.SetDirectory(0)

infileAnalysisResults.Close()

# SetRange
hGenDs = hGenDs.Rebin(len(np.arange(0.5,24,0.1)), "hGenDs", np.asarray(np.arange(0.5,24.1,0.1), "d"))
hGenDplus = hGenDplus.Rebin(len(np.arange(0.5,24,0.1)), "hGenDplus", np.asarray(np.arange(0.5,24.1,0.1), "d"))

# Load cuts configuration
with open(infileCutset, 'r') as f:
    cuts_cfg = yaml.safe_load(f)

hRecDs = hGenDs.Clone("hRecDs")
hRecDplus = hGenDplus.Clone("hRecDplus")

for iPt in range(hRecDs.GetNbinsX()):
    hRecDs.SetBinContent(iPt+1, 0)
    hRecDs.SetBinError(iPt+1, 0)
    hRecDplus.SetBinContent(iPt+1, 0)
    hRecDplus.SetBinError(iPt+1, 0)

for iPt, (ptMin, ptMax) in enumerate(zip(ptMins,ptMaxs)):
    ModelHandl = ModelHandler()
    ModelHandl.load_model_handler(config['model_handlers'][iPt])
    cols_to_keep = ["fM", "fPt"]
    dfDs = apply_model_in_batches(ModelHandl, cols_to_keep, infileRecDs, f"fPt > {ptMin} and fPt < {ptMax}")
    dfDplus = apply_model_in_batches(ModelHandl, cols_to_keep, infileRecDplus, f"fPt > {ptMin} and fPt < {ptMax}")
    selToApply = ""
    for varName in cuts_cfg['cutvars']:
        if varName == 'InvMass' or varName == 'Pt':
            continue
        if selToApply != '':
            selToApply += ' & '
        selToApply += f"({cuts_cfg['cutvars'][varName]['min'][iPt]}<{cuts_cfg['cutvars'][varName]['name']}<{cuts_cfg['cutvars'][varName]['max'][iPt]})"

    for pt in dfDs.query(f"{selToApply}")['fPt']:
        hRecDs.Fill(pt)
    for pt in dfDplus.query(f"{selToApply}")['fPt']:
        hRecDplus.Fill(pt)

# Reweight generated and reconstructed with FONLL

hWeigthsDs = hFONLLDs.Clone("hWeigthsDs")
hWeigthsDs.Divide(hPYTHIADs)
hWeigthsDplus = hFONLLDplus.Clone("hWeigthsDplus")
hWeigthsDplus.Divide(hPYTHIADplus)

hGenDsReweighted = hGenDs.Clone("hGenDsReweighted")
hGenDsReweighted.Multiply(hWeigthsDs)
hRecDsReweighted = hRecDs.Clone("hRecDsReweighted")
hRecDsReweighted.Multiply(hWeigthsDs)
hGenDplusReweighted = hGenDplus.Clone("hGenDplusReweighted")
hGenDplusReweighted.Multiply(hWeigthsDplus)
hRecDplusReweighted = hRecDplus.Clone("hRecDplusReweighted")
hRecDplusReweighted.Multiply(hWeigthsDplus)

# Rebin histos for final efficiency calculation

hGenDsReb = hGenDs.Rebin(len(ptEdges)-1, "hGenDsReb", np.asarray(ptEdges, "d"))
hRecDsReb = hRecDs.Rebin(len(ptEdges)-1, "hRecDsReb", np.asarray(ptEdges, "d"))
hGenDsReweightedReb = hGenDsReweighted.Rebin(len(ptEdges)-1, "hGenDsReweightedReb", np.asarray(ptEdges, "d"))
hRecDsReweightedReb = hRecDsReweighted.Rebin(len(ptEdges)-1, "hRecDsReweightedReb", np.asarray(ptEdges, "d"))
hGenDplusReb = hGenDplus.Rebin(len(ptEdges)-1, "hGenDplusReb", np.asarray(ptEdges, "d"))
hRecDplusReb = hRecDplus.Rebin(len(ptEdges)-1, "hRecDplusReb", np.asarray(ptEdges, "d"))
hGenDplusReweightedReb = hGenDplusReweighted.Rebin(len(ptEdges)-1, "hGenDplusReweightedReb", np.asarray(ptEdges, "d"))
hRecDplusReweightedReb = hRecDplusReweighted.Rebin(len(ptEdges)-1, "hRecDplusReweightedReb", np.asarray(ptEdges, "d"))

# Evaluate efficiency and ratios

hEffDs = hRecDsReb.Clone("hEffDs")
hEffDs.Divide(hRecDsReb, hGenDsReb, 1, 1, "b")
hEffDsReweighted = hRecDsReweightedReb.Clone("hEffDsReweighted")
hEffDsReweighted.Divide(hRecDsReweightedReb, hGenDsReweightedReb, 1, 1, "b")
hEffDplus = hRecDplusReb.Clone("hEffDplus")
hEffDplus.Divide(hRecDplusReb, hGenDplusReb, 1, 1, "b")
hEffDplusReweighted = hRecDplusReweightedReb.Clone("hEffDplusReweighted")
hEffDplusReweighted.Divide(hRecDplusReweightedReb, hGenDplusReweightedReb, 1, 1, "b")

hRatio = hEffDs.Clone("hRatioEff")
hRatio.Divide(hEffDplus)
hRatioReweighted = hEffDsReweighted.Clone("hRatioEffReweighted")
hRatioReweighted.Divide(hEffDplusReweighted)

canvasEffDs = ROOT.TCanvas("canvasEffDs", "canvasEffDs", 800, 600)
canvasEffDs.DrawFrame(0, 1.e-5, 24, 1, ";#it{p}_{T} (GeV/#it{c});Efficiency")
hEffDs.SetMarkerColor(ROOT.kRed)
hEffDs.SetMarkerStyle(ROOT.kFullCircle)
hEffDs.SetLineColor(ROOT.kRed)
hEffDs.Draw("same")
hEffDsReweighted.SetMarkerColor(ROOT.kBlue)
hEffDsReweighted.SetMarkerStyle(ROOT.kFullCircle)
hEffDsReweighted.SetLineColor(ROOT.kBlue)
hEffDsReweighted.Draw("same")
canvasEffDs.SetLogy()
legEffDs = ROOT.TLegend(0.6, 0.3, 0.9, 0.4)
legEffDs.SetBorderSize(0)
legEffDs.SetFillStyle(0)
legEffDs.AddEntry(hEffDs, "Efficiency PYTHIA Mode 2", "p")
legEffDs.AddEntry(hEffDsReweighted, "Efficiency FONLL", "p")
legEffDs.Draw()

canvasEffDplus = ROOT.TCanvas("canvasEffDplus", "canvasEffDplus", 800, 600)
canvasEffDplus.DrawFrame(0, 1.e-5, 24, 1, ";#it{p}_{T} (GeV/#it{c});Efficiency")
hEffDplus.SetMarkerColor(ROOT.kRed)
hEffDplus.SetMarkerStyle(ROOT.kFullCircle)
hEffDplus.SetLineColor(ROOT.kRed)
hEffDplus.Draw("same")
hEffDplusReweighted.SetMarkerColor(ROOT.kBlue)
hEffDplusReweighted.SetMarkerStyle(ROOT.kFullCircle)
hEffDplusReweighted.SetLineColor(ROOT.kBlue)
hEffDplusReweighted.Draw("same")
canvasEffDplus.SetLogy()
legEffDplus = ROOT.TLegend(0.6, 0.3, 0.9, 0.4)
legEffDplus.SetBorderSize(0)
legEffDplus.SetFillStyle(0)
legEffDplus.AddEntry(hEffDplus, "Efficiency PYTHIA Mode 2", "p")
legEffDplus.AddEntry(hEffDplusReweighted, "Efficiency FONLL", "p")
legEffDplus.Draw()

canvasRatio = ROOT.TCanvas("canvasRatio", "canvasRatio", 800, 600)
canvasRatio.Divide(2)
canvasRatio.cd(1).DrawFrame(0, 0., 24, 1.5, ";#it{p}_{T} (GeV/#it{c});(Acc#times#varepsilon) (D_{s}^{+})/(Acc#times#varepsilon) (D^{+})")
hRatio.SetMarkerColor(ROOT.kRed)
hRatio.SetMarkerStyle(ROOT.kFullCircle)
hRatio.SetLineColor(ROOT.kRed)
hRatio.Draw("same")
hRatioReweighted.SetMarkerColor(ROOT.kBlue)
hRatioReweighted.SetMarkerStyle(ROOT.kFullCircle)
hRatioReweighted.SetLineColor(ROOT.kBlue)
hRatioReweighted.Draw("same")
legEffRatio = ROOT.TLegend(0.6, 0.3, 0.9, 0.4)
legEffRatio.SetBorderSize(0)
legEffRatio.SetFillStyle(0)
legEffRatio.AddEntry(hEffDplus, "PYTHIA Mode 2", "p")
legEffRatio.AddEntry(hEffDplusReweighted, "FONLL", "p")
legEffRatio.Draw()
canvasRatio.cd(2).DrawFrame(0, 0., 24, 1.5, ";#it{p}_{T} (GeV/#it{c});PYTHIA/FONLL [(Acc#times#varepsilon) (D_{s}^{+})/(Acc#times#varepsilon) (D^{+})]  ")
hRatioCentralOverReweighted = hRatio.Clone("hRatioCentralOverReweighted")
hRatioCentralOverReweighted.Divide(hRatioReweighted)
hRatioCentralOverReweighted.SetMarkerColor(ROOT.kRed)
hRatioCentralOverReweighted.SetMarkerStyle(ROOT.kFullCircle)
hRatioCentralOverReweighted.SetLineColor(ROOT.kRed)
hRatioCentralOverReweighted.Draw("same")

outfile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/PtShape/StudyPtShape.root", "recreate")
canvasEffDs.Write()
canvasEffDplus.Write()
canvasRatio.Write()
hEffDs.Write()
hEffDsReweighted.Write()
hEffDplus.Write()
hEffDplusReweighted.Write()
hRatio.Write()
hRatioReweighted.Write()
hFONLLDs.Write()
hFONLLDplus.Write()
hPYTHIADs.Write()
hPYTHIADplus.Write()
hGenDs.Write()
hRecDs.Write()
hGenDplus.Write()
hRecDplus.Write()
hGenDsReweighted.Write()
hRecDsReweighted.Write()
hGenDplusReweighted.Write()
hRecDplusReweighted.Write()
hWeigthsDs.Write()
hWeigthsDplus.Write()
hGenDsReb.Write()
hRecDsReb.Write()
hGenDsReweightedReb.Write()
hRecDsReweightedReb.Write()
hGenDplusReb.Write()
hRecDplusReb.Write()
hGenDplusReweightedReb.Write()
hRecDplusReweightedReb.Write()
outfile.Close()







