CutSetsDir="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/configs"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

configFile="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/ptSmearing/FD_Fraction/Ds/Efficiency/Config_Efficiency.yaml"

parallel -j10 python3 /home/fchinu/Run3/ThesisUtils/EvaluateEfficiency.py -c ${configFile} -s ${CutSetsDir}/cutset{1}.yml -o /home/fchinu/Run3/Ds_pp_13TeV/Systematics/Tracking/ptSmearing/FD_Fraction/Ds/Efficiency/Efficiency_{1}.root ::: ${CutSets[@]}