import pandas as pd 
import numpy as np
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from DfUtils import read_parquet_in_batches

pTmins = [0, 1.5, 2, 3, 4, 5, 6, 8, 12] # list 
pTmaxs = [1.5, 2, 3, 4, 5, 6, 8, 12, 24] # list 

dfDspromptName = '/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train189890/LHC22b1b_PromptDs_Train.parquet'
dfDsFDName = '/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train189892/LHC22b1a_NonPromptDs_Train.parquet'
dfbkgName = '/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train189072/LHC22o_pass6_small.parquet'

for (pTmin, pTmax) in zip(pTmins, pTmaxs):
    dfDsprompt = read_parquet_in_batches(dfDspromptName, f'{pTmin} < fPt < {pTmax}')
    dfDsFD = read_parquet_in_batches(dfDsFDName, f'{pTmin} < fPt < {pTmax}')
    dfbkg = read_parquet_in_batches(dfbkgName, f'{pTmin} < fPt < {pTmax} and (1.7 < fM < 1.75 or  2.1 < fM < 2.15)')

    print(f'{pTmin}-{pTmax} & {len(dfDsprompt)} & {len(dfDsFD)} & {len(dfbkg)} \\\\')