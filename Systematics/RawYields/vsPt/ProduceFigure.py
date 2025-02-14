import pickle
import yaml
import ROOT
from itertools import product
from RawYieldSystematicsParallel import ProduceFigure
import argparse

parser = argparse.ArgumentParser(description='Produce figure')
parser.add_argument('ptmin', type=float, help='Minimum pT')
parser.add_argument('ptmax', type=float, help='Maximum pT')
parser.add_argument('configfile', type=str, help='multi trial config file')
args = parser.parse_args()

ptmin = float(args.ptmin)
ptmax = float(args.ptmax)


with open(args.configfile, 'rb') as f:
    config = yaml.safe_load(f)

with open(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/results/pt{ptmin*10:.1f}_{ptmax*10:.1f}.pkl", "rb") as f:
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

ProduceFigure(d, config['multitrial'], ptmin, ptmax)