#!/bin/bash
#steps to be performed
DoDataProjection=false
DoEfficiency=false
DoPromptFraction=false
DoDataRawYields=false
DoRatio=true

#wheter you are projecting a tree or a sparse
ProjectTree=true

#PARAMETERS TO BE SET (use "" for parameters not needed)
################################################################################################
Particle="Ds" # # whether it is Dplus, D0, Ds, LctopK0s or LctopKpi analysis
Cent="kXic0pPb5TeVPrompt" # used also to asses prompt or non-prompt and system

cfgFileData="/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/config_Ds_pp_data_tree_13TeV_binary.yml"
cfgFileFit="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/config_Ds_Fit_pp13TeV.yml"

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
OutDirEfficiency="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/Efficiency_LHC24d3a"
OutDirPromptFrac="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/Prompt_Fraction"
OutDirCrossSec="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/CrossSection"
OutDirRatios="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/Ratios_LHC24d3a"
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
SparseName="/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train191573/LHC22o_AnalysisResults.root"

if $DoDataProjection; then
  parallel -j 10 python3 ${ProjectScript} ${SparseName} ${CutSetsDir}/cutset{1}.yml ${OutDirRawyields}/Distr_${Particle}_data{1}.root ::: ${CutSets[@]}
fi

configFileEff="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/Config_Efficiency_LHC24d3a.yaml"
if $DoEfficiency; then
    parallel -j 10 python3 /home/fchinu/Run3/ThesisUtils/EvaluateEfficiency.py -c ${configFileEff} -s ${CutSetsDir}/cutset{1}.yml -o ${OutDirEfficiency}/Efficiency_{1}.root ::: ${CutSets[@]}
fi

PromptFractionScript="/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/GetPromptFrac.py"
DsCentralPromptFrac="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root"
DplusCentralPromptFrac="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root"
if $DoPromptFraction; then
  parallel -j 10 python3 ${PromptFractionScript} ${OutDirEfficiency}/Efficiency_{1}.root ${DsCentralPromptFrac} ${DplusCentralPromptFrac} ${OutDirPromptFrac}/PromptFracDs_{1}.root ${OutDirPromptFrac}/PromptFracDplus_{1}.root ::: ${CutSets[@]}
fi

RawYieldsScript="/home/fchinu/Run3/ThesisUtils/GetRawYieldsDplusDsFlareFlyParallel.py"
if $DoDataRawYields; then
  parallel -j 1 python3 ${RawYieldsScript} ${cfgFileFit} ${OutDirRawyields}/Distr_${Particle}_data{1}.root ${OutDirRawyields}/RawYields${Meson}_data{1}.root --batch ::: ${CutSets[@]}
fi

RatioScript="/home/fchinu/Run3/ThesisUtils/Ratio.py"
if $DoRatio; then
  for ((iCutSet=0; iCutSet<${arraylength}; iCutSet++));
  do
    python3 ${RatioScript} -r ${OutDirRawyields}/RawYields${Meson}_data${CutSets[$iCutSet]}.root -e ${OutDirEfficiency}/Efficiency_${CutSets[$iCutSet]}.root -fds ${OutDirPromptFrac}/PromptFracDs_${CutSets[$iCutSet]}.root -fdp ${OutDirPromptFrac}/PromptFracDplus_${CutSets[$iCutSet]}.root -o ${OutDirRatios}/Ratio_${CutSets[$iCutSet]}.root
  done
fi