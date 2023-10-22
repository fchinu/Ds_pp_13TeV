from ROOT import TFile
import pandas as pd
import sys
sys.path.append('/home/fchinu/DmesonAnalysis/utils')
from DfUtils import LoadDfFromRootOrParquet

# Define the input and output file paths
input_file_path = "/home/fchinu/Run3/Ds_pp_13TeV/ML/output/dataset/o2_Data_AO2D_Merged.root"
output_file_path = "/home/fchinu/Run3/Ds_pp_13TeV/ML/output/dataset/o2_Data_AO2D_Merged"   #No extension
treename = "O2hfcanddsfull"

def PickleFromRoot(input_file_path, output_file_path, treename):
    '''
    This function takes as input a root file and converts it to a pickle file.
    '''
    # Process the first half of the file
    InputFile = TFile(input_file_path,'READ')
    DirIter = InputFile.GetListOfKeys()
    DirNames = [i.GetName() for i in DirIter if "parentFiles" not in i.GetName()]
    FileName = [input_file_path] * len(DirNames)
    TreeName = treename
    TreeName = [TreeName] * len(DirNames)
    df = pd.DataFrame()
    for i in range(len(DirNames)):
        print("Processing file " + str(i) + " of " + str(len(DirNames)), end='\r')
        df_temp = LoadDfFromRootOrParquet(FileName[i], [DirNames[i]], TreeName[i])
        df = pd.concat([df, df_temp], ignore_index=True)
        del df_temp
    df.to_parquet(output_file_path + "_1.parquet")

    # Process the second half of the file
    start = len(DirNames)//2
    df = pd.DataFrame()
    for i in range(len(DirNames)//2):
        print("Processing file " + str(i) + " of " + str(len(DirNames)//2), end='\r')
        df_temp = LoadDfFromRootOrParquet(FileName[start + i], [DirNames[start + i]], TreeName[start + i])
        df = pd.concat([df, df_temp], ignore_index=True)
        del df_temp
    df.to_parquet(output_file_path + "_2.parquet")


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

def MergeDataframes(input_file_paths, output_file_path):
    '''
    This function takes as input two pickle files and merges them into a single one.
    '''
    Df = pd.concat([pd.read_parquet(path) for path in input_file_paths] , ignore_index=True)
    Df.to_parquet(output_file_path)
    del Df

# Merge the Ds and Dplus dataframes
#print("Processing Ds")
#DsDf = pd.concat([pd.read_parquet(output_file_path + "_Ds_1.parquet"), pd.read_parquet(output_file_path + "_Ds_2.parquet")] , ignore_index=True)
#DsDf.to_parquet(output_file_path + "_Ds.parquet")
#del DsDf
#print("Processing Dplus")
#DplusDf = pd.concat([pd.read_parquet(output_file_path + "_Dplus_1.parquet"), pd.read_parquet(output_file_path + "_Dplus_2.parquet")] , ignore_index=True)
#DplusDf.to_parquet(output_file_path + "_Dplus.parquet")
#del DplusDf

#PickleFromRoot(input_file_path, output_file_path, treename)
#DivideIntoDsAndDplus(output_file_path + "_1.parquet", output_file_path)
#DivideIntoDsAndDplus(output_file_path + "_2.parquet", output_file_path)
MergeDataframes([output_file_path + "_1.parquet", output_file_path + "_2.parquet"], output_file_path + ".parquet")