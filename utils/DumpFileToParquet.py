from ROOT import TFile
import pandas as pd
import sys
sys.path.append('/home/fchinu/DmesonAnalysis/utils')
from DfUtils import LoadDfFromRootOrParquet

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def ParquetFromRoot(input_file_path, output_file_path, treename, n_parts=2):
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
    for idx, (FileNamePart, DirNamesPart, TreeNamePart) in enumerate(zip(split_list(FileName, n_parts), split_list(DirNames, n_parts), split_list(TreeName, n_parts))):
        print("Processing part " + str(idx+1) + " of " + str(n_parts))
        df = pd.DataFrame()
        for i in range(len(DirNamesPart)):
            print("Processing file " + str(i) + " of " + str(len(DirNamesPart)), end='\r')
            df_temp = LoadDfFromRootOrParquet(FileNamePart[i], [DirNamesPart[i]], TreeNamePart[i])
            df = pd.concat([df, df_temp], ignore_index=True)
            del df_temp
        df.to_parquet(output_file_path + f"_{idx}.parquet")
        del df
        print("Processing part " + str(idx+1) + " of " + str(n_parts) + " done!")


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
    This function takes as input a list of pickle files and merges them into a single one.
    '''
    Df = pd.concat([pd.read_parquet(path) for path in input_file_paths] , ignore_index=True)
    Df.to_parquet(output_file_path)
    del Df

if __name__ == "__main__":
    # Define the input and output file paths
    input_file_path = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/o2_MC_big_Rebecca_Merged.root"
    output_file_path = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/o2_MC_big_Rebecca_Merged"   #No extension
    treename = "O2hfcanddsfull"
    n_divisions = 3
    isMC = False


    #ParquetFromRoot(input_file_path, output_file_path, treename, n_parts=n_divisions)
    if isMC:
        for i in range (n_divisions):
            DivideIntoDsAndDplus(output_file_path + f"_{i}.parquet", output_file_path + f"_{i}")
        MergeDataframes([output_file_path + f"_{i}_Ds.parquet" for i in range(n_divisions)], output_file_path + "_Ds.parquet")
        MergeDataframes([output_file_path + f"_{i}_Dplus.parquet" for i in range(n_divisions)], output_file_path + "_Dplus.parquet")
    else:
        MergeDataframes([output_file_path + f"_{i}.parquet" for i in range(n_divisions)], output_file_path + ".parquet")
