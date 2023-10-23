import pandas as pd

inputfiles = ["/home/fchinu/Run3/Ds_pp_13TeV/Datasets/o2_MC_big_Rebecca_Merged_Ds.parquet",
                "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/o2_Data_AO2D_Merged.parquet"]

query = "fCpa > 0.85 and fCpaXY > 0.85 and fDecayLengthXY < 0.04 and abs(fMaxNormalisedDeltaIP)<200"

for file in inputfiles:
    df = pd.read_parquet(file)
    dfSel = df.query(query)
    dfSel.to_parquet(file.replace(".parquet","_Skimmed.parquet"))
    del dfSel
    del df
