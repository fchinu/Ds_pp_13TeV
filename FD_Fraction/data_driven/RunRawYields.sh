CutSetsDir="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/configs"
DataFile="/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train183601/LHC22o_pass6_AnalysisResults.root"
FitConfig="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/config_Ds_Fit_pp13TeV.yml"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

parallel -j10 python3 /home/fchinu/Run3/ThesisUtils/ProjectDataFromSparse.py ${DataFile} ${CutSetsDir}/cutset{1}.yml /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/RawYields/Projection_{1}.root ::: ${CutSets[@]}

parallel -j1 python3 /home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/GetRawYieldsDplusDsFlareFlyParallel.py ${FitConfig} /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/RawYields/Projection_{1}.root /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/RawYields/RawYields_{1}.root --batch ::: ${CutSets[@]} 

