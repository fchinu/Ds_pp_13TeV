NPTBINS=13
PTMINS=(1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 8.0 12.0)
PTMAXS=(1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 8.0 12.0 24.0)
parallel -j $NPTBINS "python3 /home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/RawYieldSystematics.py /home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/config_multi_trial.yml -pmi {1} -pma {2}" ::: ${PTMINS[@]} :::+ ${PTMAXS[@]}
