CutSetsDir="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/configs"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

configFile="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/Efficiency/Config_Efficiency.yaml"
configFile_newMC="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/Efficiency_newMC/Config_Efficiency.yaml"
configFile_newestMC="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/Efficiency_LHC24d3a/Config_Efficiency.yaml"

#parallel -j10 python3 /home/fchinu/Run3/ThesisUtils/EvaluateEfficiency.py -c ${configFile} -s ${CutSetsDir}/cutset{1}.yml -o /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/Efficiency/Efficiency_{1}.root ::: ${CutSets[@]}

#parallel -j10 python3 /home/fchinu/Run3/ThesisUtils/EvaluateEfficiency.py -c ${configFile_newMC} -s ${CutSetsDir}/cutset{1}.yml -o /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/Efficiency_newMC/Efficiency_{1}.root ::: ${CutSets[@]}

parallel -j10 python3 /home/fchinu/Run3/ThesisUtils/EvaluateEfficiency.py -c ${configFile_newestMC} -s ${CutSetsDir}/cutset{1}.yml -o /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/Efficiency_LHC24d3a/Efficiency_{1}.root ::: ${CutSets[@]}