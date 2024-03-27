CutSetsDir="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/configs"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

configFile="/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/Config_Efficiency.yaml"

for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
do
    python3 /home/fchinu/Run3/ThesisUtils/EvaluateEfficiency.py -c ${configFile} -s ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml -o /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/Efficiency/Efficiency_${CutSets[$iCutSet]}.root
done