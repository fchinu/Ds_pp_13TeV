import pickle
import yaml
import ROOT
from itertools import product
from RawYieldSystematicsParallel import ProduceFigure

with open('/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/config_multi_trial.yml', 'rb') as f:
    config = yaml.safe_load(f)

with open("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/results/pt25.0_30.0.pkl", "rb") as f:
    d = pickle.load(f)

    #refFileName = config['reffilenames']['data']
    #refFile = ROOT.TFile.Open(refFileName)
    #hRawYieldsDsCentral = refFile.Get('hRawYields')
    #hRawYieldsDsCentral.SetDirectory(0)
    #hRawYieldsDplusCentral = refFile.Get('hRawYieldsSecondPeak')
    #hRawYieldsDplusCentral.SetDirectory(0)
    #hSigmaDsCentral = refFile.Get('hRawYieldsSigma')
    #hSigmaDsCentral.SetDirectory(0)
    #hSigmaDplusCentral = refFile.Get('hRawYieldsSigmaSecondPeak')
    #hSigmaDplusCentral.SetDirectory(0)
    #refFile.Close()
#
    #d['hRawYieldsDsCentral'] = hRawYieldsDsCentral
    #d['hRawYieldsDplusCentral'] = hRawYieldsDplusCentral
    #d['hSigmaDsCentral'] = hSigmaDsCentral
    #d['hSigmaDplusCentral'] = hSigmaDplusCentral

ProduceFigure(d, config['multitrial'], 2.5, 3.0)