from ROOT import TFile
import pandas as pd
import numpy as np
import sys
import uproot
sys.path.append('/home/fchinu/DmesonAnalysis/utils')
from pathlib import Path
from DfUtils import LoadDfFromRootOrParquet
from concurrent.futures import ThreadPoolExecutor
def ParquetFromRoot(input_file_path, output_file_path, treename, addNsigmaComb=False, preselection=None, n_parts=2, isMC=False):
    '''
    This function takes as input a root file and converts it to a parquet file.
    '''
    # Process the first half of the file
    InputFile = TFile(input_file_path,'READ')
    DirIter = InputFile.GetListOfKeys()
    DirNames = [i.GetName() for i in DirIter if "parentFiles" not in i.GetName()]
    FileName = [input_file_path] * len(DirNames)
    TreeName = treename
    TreeName = [TreeName] * len(DirNames)
    InputFile.Close()

    inputs = [f'{inFile}:{inDir}/{inTree}' for (inFile, inDir, inTree) in zip(FileName, DirNames, TreeName)]

    executor = ThreadPoolExecutor(3)
    iterator = uproot.iterate(inputs, library='pd', decompression_executor=executor,
                              interpretation_executor=executor)

    result = []
    for data in iterator:
        if preselection:
            data = data.query(preselection)
        result.append(data)

    result = pd.concat(result)
    return result

input_file_path = "/home/fchinu/Run3/Ds_pp_13TeV/TestO2/AO2D.root"
output_file_path = "/home/fchinu/Run3/Ds_pp_13TeV/TestO2/AO2D"   #No extension
treename = "O2hfcanddslite"
df = ParquetFromRoot(input_file_path, output_file_path, treename, n_parts=1, preselection = "")
df.to_parquet(output_file_path + ".parquet")