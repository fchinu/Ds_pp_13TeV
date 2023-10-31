import pandas as pd

inputfiles = ["/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Merged_LHC22o_Train130748_AO2D.parquet.gzip"]

OutputFileNames = ["/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/skimmed_data/Merged_LHC22o_Train130748_AO2D_Skimmed.parquet.gzip"]

query = "fCpa > 0.85 and fCpaXY > 0.85 and fDecayLengthXY < 1 and abs(fMaxNormalisedDeltaIP)<200 and fDeltaMassPhi < 0.015"

for (file, outfile) in zip(inputfiles, OutputFileNames):
    df = pd.read_parquet(file)
    dfSel = df.query(query)
    dfSel.to_parquet(outfile, compression="gzip")
    del dfSel
    del df
