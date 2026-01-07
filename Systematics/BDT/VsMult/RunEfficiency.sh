ConfigDir="/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/configs/efficiencies"
declare -a EffConfigs=()
for filename in ${ConfigDir}/*.yaml; do
    tmp_name="$(basename -- ${filename} .yaml)"
    tmp_name=${tmp_name:10}
    EffConfigs+=("${tmp_name}")
done

if [ ! -d "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/Efficiency" ]; then
    mkdir /home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/Efficiency
fi

parallel -j10 python3 /home/fchinu/Run3/Ds_Dp_ratio_PbPb/Efficiency/evaluate_efficiency_sparse.py ${ConfigDir}/config_eff{1}.yaml ::: ${EffConfigs[@]}
