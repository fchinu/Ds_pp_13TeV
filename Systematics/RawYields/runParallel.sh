#PTMINS=(4.0 4.5 5.0 5.5 6.0 8.0 12.0)
#PTMAXS=(4.5 5.0 5.5 6.0 8.0 12.0 24.0)
PTMINS=(4.0)
PTMAXS=(4.5)
parallel -j 1 "python3 /home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/RawYieldSystematicsParallel.py /home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/config_multi_trial_highpt.yml /home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/config_Ds_Fit_pp13TeV.yml -pmi {1} -pma {2}" ::: ${PTMINS[@]} :::+ ${PTMAXS[@]}
#parallel -j 5 "python3 /home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/ProduceFigure.py {1} {2} --highpt" ::: ${PTMINS[@]} :::+ ${PTMAXS[@]}

PTMINS=(0.5 1.0 1.5 2.0 2.5)
PTMAXS=(1.0 1.5 2.0 2.5 3.0)
#parallel -j 1 "python3 /home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/RawYieldSystematicsParallel.py /home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/config_multi_trial_lowpt.yml /home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/config_Ds_Fit_pp13TeV.yml -pmi {1} -pma {2}" ::: ${PTMINS[@]} :::+ ${PTMAXS[@]}
#parallel -j 5 "python3 /home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/ProduceFigure.py {1} {2}" ::: ${PTMINS[@]} :::+ ${PTMAXS[@]}