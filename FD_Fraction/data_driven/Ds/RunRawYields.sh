CutSetsDir="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/configs"
DataFile="/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Train183601/LHC22o_pass6_AnalysisResults.root"
FitConfig="/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/config_Ds_Fit_pp13TeV.yml"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
do
    echo "Processing ${CutSets[$iCutSet]}"
    python3 /home/fchinu/Run3/ThesisUtils/ProjectDataFromSparse.py ${DataFile} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/RawYields/Projection_${CutSets[$iCutSet]}.root
    python3 /home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/GetRawYieldsDplusDs.py ${FitConfig} kpp13TeVPrompt /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/RawYields/Projection_${CutSets[$iCutSet]}.root /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/RawYields/RawYields_${CutSets[$iCutSet]}.root
done
