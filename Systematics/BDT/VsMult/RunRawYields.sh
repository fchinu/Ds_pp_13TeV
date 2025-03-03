
# File with propmt enhanced projections (same bkg selections as for the cut variation)
ConfigFileDir="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/VsMult/configs/projections"
FitConfigDir="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/VsMult/configs/fits"
ProjectionsDir="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/VsMult/Projections_RawYields"
declare -a ProjectionsConfigs=()
for filename in ${ConfigFileDir}/*.yaml; do
    tmp_name="$(basename -- ${filename} .yaml)"
    tmp_name=${tmp_name:23}
    ProjectionsConfigs+=("${tmp_name}")
done

declare -a FitConfigs=()
for filename in ${FitConfigDir}/*.yaml; do
    tmp_name="$(basename -- ${filename} .yaml)"
    tmp_name=${tmp_name:11}
    FitConfigs+=("${tmp_name}")
done


nice -n 15 parallel -j2 python3 /home/fchinu/Run3/ThesisUtils/project_data_from_sparse.py ${ConfigFileDir}/config_projection_data_{1}.yaml ::: ${ProjectionsConfigs[@]}
declare -a ProjectionsDirToMerge=()
for dir in "$ProjectionsDir"/*; do
  if [[ -d $dir && $(basename "$dir") =~ [0-9]{2}$ ]]; then
      # ProjectionsDirToMerge+=("$(basename "$dir")")
      # echo "Merging ${ProjectionsDir}/projection_data_$(basename "$dir").root with ${ProjectionsDir}/$(basename "$dir")/*.root"
      hadd -f ${ProjectionsDir}/projection_data_$(basename "$dir").root ${ProjectionsDir}/$(basename "$dir")/*.root
  fi
done
nice -n 15 parallel -j1 python3 /home/fchinu/Run3/ThesisUtils/get_raw_yields_ds_dplus_flarefly.py ${FitConfigDir}/config_fit_{1}.yaml ::: ${FitConfigs[@]} >& ${ProjectionsDir}/raw_yields.log

