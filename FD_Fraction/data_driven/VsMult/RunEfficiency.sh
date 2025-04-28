ConfigDir="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/mult_differential/2023/configs/efficiencies"
declare -a EffConfigs=()
for filename in ${ConfigDir}/*.yaml; do
    tmp_name="$(basename -- ${filename} .yaml)"
    tmp_name=${tmp_name:10}
    EffConfigs+=("${tmp_name}")
done

parallel -j10 python3 /home/fchinu/Run3/ThesisUtils/evaluate_efficiency_sparse.py ${ConfigDir}/config_eff{1}.yaml ::: ${EffConfigs[@]}
