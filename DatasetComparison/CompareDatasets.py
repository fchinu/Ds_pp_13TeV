import os
import sys
sys.path.append('/home/fchinu/DmesonAnalysis/utils')
import argparse
import pickle
import yaml
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

from DfUtils import LoadDfFromRootOrParquet
from ROOT import TFile

from hipe4ml import plot_utils
from hipe4ml.model_handler import ModelHandler, ModelHandlerNN
from hipe4ml.tree_handler import TreeHandler


def data_prep(inputCfg, iBin, PtBin, OutPutDirPt, PromptDf, FDDf, BkgDf): #pylint: disable=too-many-statements, too-many-branches
    '''
    function for data preparation
    '''
    nPrompt = len(PromptDf)
    nFD = len(FDDf)
    nBkg = len(BkgDf)
    if FDDf.empty:
        out = f'\n     Signal: {nPrompt}\n     Bkg: {nBkg}'
    else:
        out = f'\n     Prompt: {nPrompt}\n     FD: {nFD}\n     Bkg: {nBkg}'
    print(f'Number of available candidates in {PtBin[0]} < pT < {PtBin[1]} GeV/c:{out}')

    dataset_opt = inputCfg['data_prep']['dataset_opt']
    seed_split = inputCfg['data_prep']['seed_split']
    test_f = inputCfg['data_prep']['test_fraction']

    if dataset_opt == 'equal':
        if FDDf.empty:
            nCandToKeep = min([nPrompt, nBkg])
            out = 'signal'
            out2 = 'signal'
        else:
            nCandToKeep = min([nPrompt, nFD, nBkg, 10])
            out = 'prompt, FD'
            out2 = 'prompt'
        print((f'Keep same number of {out} and background (minimum) for training and '
               f'testing ({1 - test_f}-{test_f}): {nCandToKeep}'))
        print(f'Fraction of real data candidates used for ML: {nCandToKeep/nBkg:.5f}')

        if nPrompt > nCandToKeep:
            print((f'Remaining {out2} candidates ({nPrompt - nCandToKeep})'
                   'will be used for the efficiency together with test set'))
        if nFD > nCandToKeep:
            print((f'Remaining FD candidates ({nFD - nCandToKeep}) will be used for the '
                   'efficiency together with test set'))

        TotDf = pd.concat([BkgDf.iloc[:nCandToKeep], PromptDf.iloc[:nCandToKeep], FDDf.iloc[:nCandToKeep]], sort=True)
        if FDDf.empty:
            LabelsArray = np.array([0]*nCandToKeep + [1]*nCandToKeep)
        else:
            LabelsArray = np.array([0]*nCandToKeep + [1]*nCandToKeep + [2]*nCandToKeep)
        if test_f < 1:
            TrainSet, TestSet, yTrain, yTest = train_test_split(TotDf, LabelsArray, test_size=test_f,
                                                                random_state=seed_split)
        else:
            TrainSet = pd.DataFrame()
            TestSet = TotDf.copy()
            yTrain = pd.Series()
            yTest = LabelsArray.copy()

        TrainTestData = [TrainSet, yTrain, TestSet, yTest]
        PromptDfSelForEff = pd.concat([PromptDf.iloc[nCandToKeep:], TestSet[pd.Series(yTest).array == 1]], sort=False)
        if FDDf.empty:
            FDDfSelForEff = pd.DataFrame()
        else:
            FDDfSelForEff = pd.concat([FDDf.iloc[nCandToKeep:], TestSet[pd.Series(yTest).array == 2]], sort=False)
        del TotDf

    elif dataset_opt == 'max_signal':
        nCandBkg = round(inputCfg['data_prep']['bkg_mult'][iBin] * (nPrompt + nFD))
        out = 'signal' if FDDf.empty else 'prompt and FD'
        print((f'Keep all {out} and use {nCandBkg} bkg candidates for training and '
               f'testing ({1 - test_f}-{test_f})'))
        if nCandBkg >= nBkg:
            nCandBkg = nBkg
            print('\033[93mWARNING: using all bkg available, not good!\033[0m')
        print(f'Fraction of real data candidates used for ML: {nCandBkg/nBkg:.5f}')

        TotDf = pd.concat([BkgDf.iloc[:nCandBkg], PromptDf, FDDf], sort=True)
        if FDDf.empty:
            LabelsArray = np.array([0]*nCandBkg + [1]*nPrompt)
        else:
            LabelsArray = np.array([0]*nCandBkg + [1]*nPrompt + [2]*nFD)
        if test_f < 1:
            TrainSet, TestSet, yTrain, yTest = train_test_split(TotDf, LabelsArray, test_size=test_f,
                                                                random_state=seed_split)
        else:
            TrainSet = pd.DataFrame()
            TestSet = TotDf.copy()
            yTrain = pd.Series()
            yTest = LabelsArray.copy()

        TrainTestData = [TrainSet, yTrain, TestSet, yTest]
        PromptDfSelForEff = TestSet[pd.Series(yTest).array == 1]
        FDDfSelForEff = pd.DataFrame() if FDDf.empty else TestSet[pd.Series(yTest).array == 2]
        del TotDf

    else:
        print(f'\033[91mERROR: {dataset_opt} is not a valid option!\033[0m')
        sys.exit()

    # plots
    
    VarsToDraw = inputCfg['plots']['plotting_columns']
    if len(VarsToDraw)>6:
        chunks = [VarsToDraw[x:x+6] for x in range(0, len(VarsToDraw), 6)]
    else:
        chunks = [VarsToDraw]

    LegLabels = [inputCfg['output']['leg_labels']['Bkg'],
                    inputCfg['output']['leg_labels']['Prompt']]
    if inputCfg['output']['leg_labels']['FD'] is not None:
        LegLabels.append(inputCfg['output']['leg_labels']['FD'])
    OutputLabels = [inputCfg['output']['out_labels']['Bkg'],
                    inputCfg['output']['out_labels']['Prompt']]
    if inputCfg['output']['out_labels']['FD'] is not None:
        OutputLabels.append(inputCfg['output']['out_labels']['FD'])
    ListDf = [BkgDf, PromptDf] if FDDf.empty else [BkgDf, PromptDf, FDDf]
    #_____________________________________________
    for idx, chunk in enumerate(chunks):
        plot_utils.plot_distr(ListDf, chunk, 100, LegLabels, figsize=(12, 7),
                                alpha=0.3, log=True, grid=False, density=True)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(f'{OutPutDirPt}/DistributionsAll{idx}_pT_{PtBin[0]}_{PtBin[1]}.pdf')
        plt.close('all')
    #_____________________________________________
    CorrMatrixFig = plot_utils.plot_corr(ListDf, VarsToDraw, LegLabels)
    for Fig, Lab in zip(CorrMatrixFig, OutputLabels):
        plt.figure(Fig.number)
        plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
        Fig.savefig(f'{OutPutDirPt}/CorrMatrix{Lab}_pT_{PtBin[0]}_{PtBin[1]}.pdf')

    return TrainTestData, PromptDfSelForEff, FDDfSelForEff


def main(): #pylint: disable=too-many-statements
    # read config file
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileNameML.yml', help='config file name for ml')
    parser.add_argument("--train", help="perform only training and testing", action="store_true")
    parser.add_argument("--apply", help="perform only application", action="store_true")
    args = parser.parse_args()

    print('Loading analysis configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading analysis configuration: Done!')

    print('Loading and preparing data files: ...', end='\r')
    
    if not inputCfg['input']['merged']:
        DataFileName = inputCfg['input']['data']
        DataFile = TFile(DataFileName,'READ')
        DirIter = DataFile.GetListOfKeys()
        DataDirNames = [i.GetName() for i in DirIter if "parentFiles" not in i.GetName()]
        DataFileName = [DataFileName] * len(DataDirNames)
        DataTreeName = inputCfg['input']['treename']
        DataTreeName = [DataTreeName] * len(DataDirNames)
        DataFileName = DataFileName[:10]
        DataDirNames = DataDirNames[:10]
        DataTreeName = DataTreeName[:10]

        PromptFileName = inputCfg['input']['prompt']
        PromptFile = TFile(PromptFileName,'READ')
        DirIter = PromptFile.GetListOfKeys()
        PromptDirNames = [i.GetName() for i in DirIter if "parentFiles" not in i.GetName()]
        PromptFileName = [PromptFileName] * len(PromptDirNames)
        PromptTreeName = inputCfg['input']['treename']
        PromptTreeName = [PromptTreeName] * len(PromptDirNames)
        PromptFileName = PromptFileName[:10]
        PromptDirNames = PromptDirNames[:10]
        PromptTreeName = PromptTreeName[:10]

        if inputCfg['input']['FD'] is not None:
            FDFileName = inputCfg['input']['FD']
            FDFile = TFile(FDFileName,'READ')
            DirIter = FDFile.GetListOfKeys()
            FDDirNames = [i.GetName() for i in DirIter if "parentFiles" not in i.GetName()]
            FDFileName = [FDFileName] * len(FDDirNames)
            FDTreeName = inputCfg['input']['treename']
            FDTreeName = [FDTreeName] * len(FDDirNames)
            FDDf = LoadDfFromRootOrParquet(FDFileName, FDDirNames, FDTreeName)
            FDFileName = FDFileName[:10]
            FDDirNames = FDDirNames[:10]
            FDTreeName = FDTreeName[:10]
        else:
            FDDf = None

        DataDf = LoadDfFromRootOrParquet(DataFileName, DataDirNames, DataTreeName)
        PromptDf = LoadDfFromRootOrParquet(PromptFileName, PromptDirNames, PromptTreeName)
        #PromptDf = PromptDf.query('fFlagMcMatchRec == 4')

        if inputCfg['data_prep']['filt_bkg_mass']:
            DataDf = DataDf.query(inputCfg['data_prep']['filt_bkg_mass'])


        print(PromptDf['fPt'].max())
        print(PromptDf['fPt'].min())
        #print(PromptDf.query(f'{PtBin[0]} < fPt < {PtBin[1]}'))

        PtBins = [[a, b] for a, b in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max'])]

        print('Loading and preparing data files: Done!')
        for iBin, PtBin in enumerate(PtBins):
            print(f'\n\033[94mStarting ML analysis --- {PtBin[0]} < pT < {PtBin[1]} GeV/c\033[0m')

            OutPutDirPt = os.path.join(os.path.expanduser(inputCfg['output']['dir']), f'pt{PtBin[0]}_{PtBin[1]}')
            if os.path.isdir(OutPutDirPt):
                print((f'\033[93mWARNING: Output directory \'{OutPutDirPt}\' already exists,'
                    ' overwrites possibly ongoing!\033[0m'))
            else:
                os.makedirs(OutPutDirPt)

            # data preparation
            #_____________________________________________
            FDDfPt = pd.DataFrame() if FDDf is None else FDDf.query(f'{PtBin[0]} < fPt < {PtBin[1]} and (fFlagMcMatchRec == 0)')
            TrainTestData, PromptDfSelForEff, FDDfSelForEff = data_prep(inputCfg, iBin, PtBin, OutPutDirPt,
                                                                        PromptDf.query(f'{PtBin[0]} < fPt < {PtBin[1]} and (fFlagMcMatchRec == 0)'), FDDfPt,
                                                                        DataDf.query(f'{PtBin[0]} < fPt < {PtBin[1]}'))


    else:
        PromptHandler = TreeHandler(inputCfg['input']['prompt'], inputCfg['input']['treename'])
        FDHandler = None if inputCfg['input']['FD'] is None else TreeHandler(inputCfg['input']['FD'],
                                                                            inputCfg['input']['treename'])
        DataHandler = TreeHandler(inputCfg['input']['data'], inputCfg['input']['treename'])

        if inputCfg['data_prep']['filt_bkg_mass']:
            BkgHandler = DataHandler.get_subset(inputCfg['data_prep']['filt_bkg_mass'], frac=1.,
                                                rndm_state=inputCfg['data_prep']['seed_split'])
        else:
            BkgHandler = DataHandler

        PtBins = [[a, b] for a, b in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max'])]
        PromptHandler.slice_data_frame('fPt', PtBins, True)
        if FDHandler is not None:
            FDHandler.slice_data_frame('fPt', PtBins, True)
        DataHandler.slice_data_frame('fPt', PtBins, True)
        BkgHandler.slice_data_frame('fPt', PtBins, True)
        print('Loading and preparing data files: Done!')

        for iBin, PtBin in enumerate(PtBins):
            print(f'\n\033[94mStarting ML analysis --- {PtBin[0]} < pT < {PtBin[1]} GeV/c\033[0m')

            OutPutDirPt = os.path.join(os.path.expanduser(inputCfg['output']['dir']), f'pt{PtBin[0]}_{PtBin[1]}')
            if os.path.isdir(OutPutDirPt):
                print((f'\033[93mWARNING: Output directory \'{OutPutDirPt}\' already exists,'
                    ' overwrites possibly ongoing!\033[0m'))
            else:
                os.makedirs(OutPutDirPt)

            # data preparation
            #_____________________________________________
            FDDfPt = pd.DataFrame() if FDHandler is None else FDHandler.get_slice(iBin)
            TrainTestData, PromptDfSelForEff, FDDfSelForEff = data_prep(inputCfg, iBin, PtBin, OutPutDirPt,
                                                                        PromptHandler.get_slice(iBin), FDDfPt,
                                                                        BkgHandler.get_slice(iBin))

main()