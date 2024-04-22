import ROOT 
import numpy as np

AnalysisResultName = '/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train198256/LHC22b1b_AnalysisResults.root'
ProjectionsName = '/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/MC/Projections_DsMC.root'
OutFileName = '/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/Tests/EffFromMcAsData/EfficiencyFromSparse.root'

ptMins = [0.5, 1., 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 12] 
ptMaxs = [1., 1.5, 2., 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 12, 24]
ptEdges = ptMins + [ptMaxs[-1]]

Projections = []
ProjectionsFile = ROOT.TFile.Open(ProjectionsName)
for i, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    hProjection = ProjectionsFile.Get(f'hMass_{ptMin*10:.0f}_{ptMax*10:.0f}')
    Projections.append(hProjection.GetEntries())

ProjectionsFile.Close()

genParticles = []
AnalysisResultFile = ROOT.TFile.Open(AnalysisResultName)
for i, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    hGenParticles = AnalysisResultFile.Get(f"hf-task-ds/MC/Ds/Prompt/hPtGen")                
    genParticles.append(hGenParticles.Integral(hGenParticles.FindBin(ptMin), hGenParticles.FindBin(ptMax)-1))

AnalysisResultFile.Close()

eff = [Projections[i]/genParticles[i] for i in range(len(ptMins))]

outFile = ROOT.TFile(OutFileName, 'recreate')
hEff = ROOT.TH1F('hEff', 'Efficiency', len(ptMins), np.asarray(ptEdges,'d'))
for iPt in range(len(ptMins)):
    hEff.SetBinContent(iPt+1, eff[iPt])
    hEff.SetBinError(iPt+1, np.sqrt(eff[iPt] * (1 - eff[iPt]) / genParticles[iPt]))
hEff.Write()

outFile.Close()