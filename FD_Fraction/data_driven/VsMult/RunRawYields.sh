CutSetsDir="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/configs"
DataFile="/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train252059/LHC22o_AnalysisResults.root"
FitConfig="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/config_Ds_Fit_pp13TeV_multiplicity.yml"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

parallel -j10 python3 /home/fchinu/Run3/ThesisUtils/ProjectDataFromSparse.py ${DataFile} ${CutSetsDir}/cutset{1}.yml /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/Projections_RawYields/Projection_{1}.root ::: ${CutSets[@]}

parallel -j1 python3 /home/fchinu/Run3/ThesisUtils/GetRawYieldsDplusDsFlareFlyParallel.py /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/Projections_RawYields/Projection_{1}.root ${FitConfig} /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/Projections_RawYields/RawYields_{1}.root --batch ::: ${CutSets[@]} 

