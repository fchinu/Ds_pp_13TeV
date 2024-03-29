#!/bin/bash
#steps to be performed
DoDataProjection=false
DoDataRawYields=true
DoEfficiency=false
DoRatios=false

#wheter you are projecting a tree or a sparse
ProjectTree=true

#PARAMETERS TO BE SET (use "" for parameters not needed)
################################################################################################
Particle="Ds" # # whether it is Dplus, D0, Ds, LctopK0s or LctopKpi analysis
Cent="kXic0pPb5TeVPrompt" # used also to asses prompt or non-prompt and system

cfgFileData="/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/config_Ds_pp_data_tree_13TeV_binary.yml"
cfgFileFit="/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/config_Ds_Fit_pp13TeV.yml"

#assuming cutsets config files starting with "cutset" and are .yml

CutSetsDir="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/configs"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}



OutDirRawyields="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/RawYields"
OutDirEfficiency="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/Efficiency"
OutDirCrossSec="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/CrossSection"
OutDirRatios="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/Ratios"
################################################################################################

if [ ! -d "${OutDirRawyields}" ]; then
  mkdir ${OutDirRawyields}
fi

if [ ! -d "${OutDirEfficiency}" ]; then
  mkdir ${OutDirEfficiency}
fi

if [ ! -d "${OutDirRatios}" ]; then
  mkdir ${OutDirRatios}
fi

if [ ! -d "${OutDirRaa}" ] && [ "${OutDirRaa}" != "" ]; then
  mkdir ${OutDirRaa}
fi

#project sparses or trees
ProjectScript="/home/fchinu/Run3/ThesisUtils/ProjectDataFromSparse.py"
SparseName="/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train183601/LHC22o_pass6_AnalysisResults.root"

if $DoDataProjection; then
    parallel -j 10 python3 ${ProjectScript} ${SparseName} ${CutSetsDir}/cutset{1}.yml ${OutDirRawyields}/Distr_${Particle}_data{1}.root ::: ${CutSets[@]}
fi

RawYieldsScript="/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/GetRawYieldsDplusDsFlareFly.py"
if $DoDataRawYields; then
  parallel -j 10 python3 ${RawYieldsScript} ${cfgFileFit} kpp13TeVPrompt ${OutDirRawyields}/Distr_${Particle}_data{1}.root ${OutDirRawyields}/RawYields${Meson}_data{1}.root --batch ::: ${CutSets[@]}
fi

configFileEff="/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/Config_Efficiency.yaml"
if $DoEfficiency; then
  for ((iCutSet=0; iCutSet<${arraylength}; iCutSet++));
  do
    python3 /home/fchinu/Run3/ThesisUtils/EvaluateEfficiency.py -c ${configFileEff} -s ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml -o ${OutDirEfficiency}/Efficiency_${CutSets[$iCutSet]}.root
  done
fi

RatioScript="/home/fchinu/Run3/Ds_pp_13TeV/Ratios/Ratio.py"
if $DoRatio; then
  for ((iCutSet=0; iCutSet<${arraylength}; iCutSet++));
  do
    python3 ${RatioScript} -r ${OutDirRawyields}/RawYields${Meson}_data${CutSets[$iCutSet]}.root -e ${OutDirEfficiency}/Efficiency_${CutSets[$iCutSet]}.root -o ${OutDirRatios}/Ratio_${CutSets[$iCutSet]}.root
  done
fi


#compute HFPtSpectrumRaa
#if $DoHFPtSpecRaa; then
#  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
#  do
#    echo $(tput setaf 4) Compute HFPtspectrumRaa $(tput sgr0)
#    echo '.x HFPtSpectrumRaa.C+("'${pprefFileName}'","'${OutDirCrossSec}'/HFPtSpectrum'${Particle}${CutSets[$iCutSet]}'.root","'${OutDirRaa}'/HFPtSpectrumRaa'$#{Particle}${CutSets[$iCutSet]}'.root",4,1,kNb,'${Cent}',k2018,k5dot023,1./3,3,6,false,1)' | root -l -b
#    echo '.q'
#  done
#fi
#
##compute corrected yield
#if $DoDmesonYield; then
#  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
#  do
#    echo $(tput setaf 4) Compute corrected yield $(tput sgr0)
#    echo '.x ComputeDmesonYield.C+(k'${Particle}','${Cent}',2,1,"'${pprefFileName}'","'${OutDirCrossSec}'/HFPtSpectrum'${Particle}${CutSets[$iCutSet]}'.root",#"","'${OutDirRaa}'/HFPtSpectrumRaa'${Particle}${CutSets[$iCutSet]}'.root","","'${OutDirCrossSec}'","'${CutSets[$iCutSet]}'",1,1./3,3,false,1)' | root -l -b
#    echo '.q'
#  done
#fi
#
#if $DoXic0Crosssec; then
#  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
#  do
#    echo $(tput setaf 4) Compute corrected yield $(tput sgr0)
#    python3 /home/fchinu/Xic0_pPb_5TeV/ML_CrossSection/Efficiency_and_Crosssection.py ${OutDirRawyields}/RawYields${Meson}_data${CutSets[$iCutSet]}.root $#{OutDirEfficiency}/Efficiency_${Particle}${Channel}${CutSets[$iCutSet]}.root ${OutDirCrossSec}/HFPtSpectrum${Particle}${CutSets[$iCutSet]}.root
#  done
#fi