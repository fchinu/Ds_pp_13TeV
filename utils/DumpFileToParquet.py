from ROOT import TFile
import pandas as pd
import numpy as np
import sys
import uproot
sys.path.append('/home/fchinu/DmesonAnalysis/utils')
from pathlib import Path
from DfUtils import LoadDfFromRootOrParquet
from concurrent.futures import ThreadPoolExecutor

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def SplitTrainTest(df, train_size=0.8, random_state=42):
    trainset = df.sample(frac=train_size, random_state=random_state)
    testset = df.drop(trainset.index)
    return trainset, testset

def ParquetFromRoot(input_file_path, output_file_path, treename, addNsigmaComb=False, preselection=None, n_parts=2, isMC=False, train_size=1):
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
    if addNsigmaComb:
        result = AddNsigmaComb(result)

    if isMC:
        result = result.query('abs(fFlagMcMatchRec) == 4')
        PromptDs = result.query('fOriginMcRec == 1 and fFlagMcDecayChanRec == 1')
        NonPromptDs = result.query('fOriginMcRec == 2 and fFlagMcDecayChanRec == 1')
        PromptDplus = result.query('fOriginMcRec == 1 and fFlagMcDecayChanRec == 3')
        NonPromptDplus = result.query('fOriginMcRec == 2 and fFlagMcDecayChanRec == 3')
        del result
        if train_size < 1:
            # Ds
            trainsetPromptDs, testsetPromptDs = SplitTrainTest(PromptDs, train_size, random_state=42)
            trainsetPromptDs.to_parquet(output_file_path + "_DsPrompt_Train.parquet")
            testsetPromptDs.to_parquet(output_file_path + "_DsPrompt_Eff.parquet")
            trainsetNonPrompt, testsetNonPrompt = SplitTrainTest(NonPromptDs, train_size, random_state=42)
            trainsetNonPrompt.to_parquet(output_file_path + "_DsFD_Train.parquet")
            testsetNonPrompt.to_parquet(output_file_path + "_DsFD_Eff.parquet")
            del trainsetPromptDs, testsetPromptDs, trainsetNonPrompt, testsetNonPrompt, PromptDs, NonPromptDs
            # Dplus
            trainsetPromptDplus, testsetPromptDplus = SplitTrainTest(PromptDplus, train_size, random_state=42)
            trainsetPromptDplus.to_parquet(output_file_path + "_DplusPrompt_Train.parquet")
            testsetPromptDplus.to_parquet(output_file_path + "_DplusPrompt_Eff.parquet")
            trainsetNonPromptDplus, testsetNonPromptDplus = SplitTrainTest(NonPromptDplus, train_size, random_state=42)
            trainsetNonPromptDplus.to_parquet(output_file_path + "_DplusFD_Train.parquet")
            testsetNonPromptDplus.to_parquet(output_file_path + "_DplusFD_Eff.parquet")
            del trainsetPromptDplus, testsetPromptDplus, trainsetNonPromptDplus, testsetNonPromptDplus, PromptDplus, NonPromptDplus
        else:
            # Ds
            PromptDs.to_parquet(output_file_path + "_DsPrompt.parquet")
            NonPromptDs.to_parquet(output_file_path + "_DsFD.parquet")
            del PromptDs, NonPromptDs
            # Dplus
            PromptDplus.to_parquet(output_file_path + "_DplusPrompt.parquet")
            NonPromptDplus.to_parquet(output_file_path + "_DplusFD.parquet")
            del PromptDplus, NonPromptDplus
    else:
        result.to_parquet(output_file_path + ".parquet")
        del result


def GetNsigComb(row, particle, num):
    if not row[f'fNSigTpc{particle}{num}']:
        return row[f'fNSigTof{particle}{num}']
    if not row[f'fNSigTpc{particle}{num}']:
        return row[f'fNSigTof{particle}{num}']
    return np.sqrt((row[f'fNSigTpc{particle}{num}']**2 + row[f'fNSigTof{particle}{num}']**2)/2)

def AddNsigmaComb(df):
    df["fNSigCombPi0"] = df.apply(lambda row: GetNsigComb(row, "Pi", 0), axis=1)
    df["fNSigCombPi1"] = df.apply(lambda row: GetNsigComb(row, "Pi", 1), axis=1)
    df["fNSigCombPi2"] = df.apply(lambda row: GetNsigComb(row, "Pi", 2), axis=1)
    df["fNSigCombKa0"] = df.apply(lambda row: GetNsigComb(row, "Ka", 0), axis=1)
    df["fNSigCombKa1"] = df.apply(lambda row: GetNsigComb(row, "Ka", 1), axis=1)
    df["fNSigCombKa2"] = df.apply(lambda row: GetNsigComb(row, "Ka", 2), axis=1)
    return df

def AddNsigmaCombToParquet(input_file_path, output_file_path):
    df = pd.read_parquet(input_file_path)
    df = AddNsigmaComb(df)
    df.to_parquet(output_file_path + ".parquet")

def DivideIntoDsAndDplus(input_file_path, output_file_path):
    '''
    This function takes as input a pickle file and divides it into Ds and Dplus dataframes.
    '''
    df = pd.read_parquet(input_file_path)
    print("Processing Ds")
    DsDf = df.query('fFlagMcMatchRec == 4')
    DsDf.to_parquet(output_file_path + "_Ds.parquet")
    del DsDf

    print("Processing Dplus")
    DplusDf = df.query('fFlagMcMatchRec == 1')
    DplusDf.to_parquet(output_file_path + "_Dplus.parquet")
    del DplusDf

def GetDsDataframe(input_file_path, output_file_path):
    '''
    This function takes as input a pickle file and returns the Ds dataframe.
    '''
    df = pd.read_parquet(input_file_path)
    DsDf = df.query('fFlagMcMatchRec == 4')
    DsDf.to_parquet(output_file_path + "_Ds.parquet")
    del DsDf
    del df

def GetDplusDataframe(input_file_path, output_file_path):
    '''
    This function takes as input a pickle file and returns the Dplus dataframe.
    '''
    df = pd.read_parquet(input_file_path)
    DplusDf = df.query('fFlagMcMatchRec == 1')
    DplusDf.to_parquet(output_file_path + "_Dplus.parquet")
    del DplusDf
    del df

def MergeDataframes(input_file_paths, output_file_path):
    '''
    This function takes as input a list of parquet files and merges them into a single one.
    '''
    bin_edges = [1,2,3,4,5,6,7,8,10,12,16,24,36,50]

    PtBins = [[a, b] for a, b in zip(bin_edges[:-1], bin_edges[1:])]
    PtBinsLabels = [f"{PtBins[i][0]}_{PtBins[i][1]}" for i in range(len(PtBins))]


    for path in input_file_paths:
        print("Processing file " + path)
        print("Reading file " + path, end='\r')
        df = pd.read_parquet(path, filters = [("fPt", ">", bin_edges[0]), ("fPt", "<", bin_edges[-1])])
        # Create a new column 'fPt_bin' based on custom binning
        print("Cutting file " + path, end='\r')
        df['fPt_bin'] = pd.cut(df['fPt'], bins=bin_edges, labels=PtBinsLabels, right=False)
        # Write the dataframe to a parquet file
        print("Writing file " + path, end='\r')
        file_path = Path(output_file_path)
        if file_path.exists():
            df.to_parquet(file_path, engine='fastparquet', append=True, partition_cols=['fPt_bin'], compression='gzip')
        else:
            df.to_parquet(file_path, engine='fastparquet', partition_cols=['fPt_bin'], compression='gzip')
        del df
        print("Processing file " + path + " done!")

if __name__ == "__main__":
    # Define the input and output file paths
    input_file_path = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Train165702/Merged_LHC22b1b_Train165702_AO2D.root"
    output_file_path = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Train165702/Merged_LHC22b1b_Train165702"   #No extension
    treename = "O2hfcanddslite"
    selections = ''
    n_divisions = 1
    train_size = 0.8
    RecreateFile = True
    isMC = True



    ParquetFromRoot(input_file_path, output_file_path, treename, n_parts=n_divisions, preselection = selections, isMC=isMC, train_size=train_size)
    #AddNsigmaCombToParquet(output_file_path+".parquet",output_file_path+".parquet")
    #if isMC:
    #    for i in range (n_divisions):
    #        GetDsDataframe(output_file_path + f"_{i}.parquet", output_file_path + f"_{i}")
    #    MergeDataframes([output_file_path + f"_{i}_Ds.parquet" for i in range(n_divisions)], output_file_path + "_Ds.parquet")
    #    MergeDataframes([output_file_path + f"_{i}_Dplus.parquet" for i in range(n_divisions)], output_file_path + "_Dplus.parquet")
    #else:
    #    if RecreateFile:
    #        path = Path(output_file_path + ".parquet.gzip")
    #        if path.exists():
    #            os.remove(path)
    #    MergeDataframes([output_file_path + f"_{i}.parquet" for i in range(n_divisions)], output_file_path)
