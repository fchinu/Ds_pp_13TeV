# INFILE="/home/fchinu/Run3/Ds_pp_13TeV/utils/AODsForDownload.txt"
# AODLIST_OUTFILE="/home/fchinu/Run3/Ds_pp_13TeV/utils/fileList130748.txt"


# python3 /home/fchinu/Run3/Ds_pp_13TeV/utils/CreateListOfFilesFromGrid.py --input $INFILE --output $AODLIST_OUTFILE

# root -l -q /home/fchinu/Run3/Ds_pp_13TeV/utils/MergeFileOnGrid.C+

o2-aod-merger --input /home/fchinu/Run3/Ds_pp_13TeV/utils/inputs.txt --output /home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Merged_LHC22o_Train130748_AO2D.root
