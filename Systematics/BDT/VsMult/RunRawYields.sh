
# File with propmt enhanced projections (same bkg selections as for the cut variation)
ConfigFileDir="/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/configs/projections/data"
ConfigFileDirMc="/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/configs/projections/mc"
FitConfigDir="/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/configs/fits"
ProjectionsDir="/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/Projections_RawYields"

# Projections - Data
declare -a ProjectionsConfigs=()
for filename in ${ConfigFileDir}/*.yaml; do
    tmp_name="$(basename -- ${filename} .yaml)"
    tmp_name=${tmp_name:23}
    ProjectionsConfigs+=("${tmp_name}")
done

# Projections - MC
declare -a ProjectionsConfigsMc=()
for filename in ${ConfigFileDirMc}/*.yaml; do
    tmp_name="$(basename -- ${filename} .yaml)"
    tmp_name=${tmp_name:21}
    ProjectionsConfigsMc+=("${tmp_name}")
done

declare -a FitConfigs=()
for filename in ${FitConfigDir}/*.yaml; do
    tmp_name="$(basename -- ${filename} .yaml)"
    tmp_name=${tmp_name:11}
    FitConfigs+=("${tmp_name}")
done


nice -n 15 parallel -j2 python3 /home/fchinu/Run3/Ds_Dp_ratio_PbPb/Projections_RawYields/project_data_from_sparse.py ${ConfigFileDir}/config_projection_data_{1}.yaml ::: ${ProjectionsConfigs[@]}
declare -a ProjectionsDirToMerge=()
for dir in "$ProjectionsDir"/*; do
  if [[ -d $dir && $(basename "$dir") =~ [0-9]{2}$ ]]; then
      # ProjectionsDirToMerge+=("$(basename "$dir")")
      # echo "Merging ${ProjectionsDir}/projection_data_$(basename "$dir").root with ${ProjectionsDir}/$(basename "$dir")/*.root"
      hadd -f ${ProjectionsDir}/projection_data_$(basename "$dir").root ${ProjectionsDir}/$(basename "$dir")/*.root
  fi
done
nice -n 15 parallel -j2 python3 /home/fchinu/Run3/Ds_Dp_ratio_PbPb/Projections_RawYields/project_data_from_sparse.py ${ConfigFileDirMc}/config_projection_mc_{1}.yaml ::: ${ProjectionsConfigsMc[@]}
declare -a ProjectionsDirToMerge=()
for dir in "$ProjectionsDir"/MC/*; do
    if [[ -d $dir && $(basename "$dir") =~ [0-9]{2}$ ]]; then
      # ProjectionsDirToMerge+=("$(basename "$dir")")
      # echo "Merging ${ProjectionsDir}/projection_data_$(basename "$dir").root with ${ProjectionsDir}/$(basename "$dir")/*.root"
        hadd -f ${dir}/dplus.root ${dir}/dplus_nonprompt.root ${dir}/dplus_prompt.root
    fi
done
nice -n 15 parallel -j1 python3 /home/fchinu/Run3/Ds_Dp_ratio_PbPb/Projections_RawYields/get_raw_yields.py ${FitConfigDir}/config_fit_{1}.yaml ::: ${FitConfigs[@]} >& ${ProjectionsDir}/raw_yields.log

