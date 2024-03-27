import pandas as pd
import numpy as np
import sys
sys.path.append('/home/fchinu/Run3/ThesisUtils')
from DfUtils import read_parquet_in_batches


def GetNCandsKKPiAndPiKK(infile, pTbins):
    for idx, (ptMin, ptMax) in enumerate(zip(pTbins[:-1], pTbins[1:])):
        selections = f'fPt > {ptMin} and fPt < {ptMax}'
        df = read_parquet_in_batches(infile, selections)
        print(f'pT range: {ptMin} - {ptMax}')
        print("# PiKK candidates: ", len(df.query("fIsSelDsToPiKK>=3")))
        print("# KKPi candidates: ", len(df.query("fIsSelDsToKKPi>=3")))
        del df

if __name__ == "__main__":
   
    pTbins = [2,4,6,8,12]
    infile = '/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Train165702/Merged_LHC22b1b_Train165702_DsPrompt.parquet'
    print("Ds:")
    GetNCandsKKPiAndPiKK(infile, pTbins)
    print("")
    infile = '/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Train165702/Merged_LHC22b1b_Train165702_DplusPrompt.parquet'
    print("Dplus:")
    GetNCandsKKPiAndPiKK(infile, pTbins)
