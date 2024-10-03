CutSetsDir="/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/configs"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

infile="/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train253352/LHC24d3b_AnalysisResults.root"
weightsFile="/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/VsMult/Mult_weights.root"

parallel -j10 python3 /home/fchinu/Run3/ThesisUtils/EvaluateEfficiencySparse.py ${infile} -s ${CutSetsDir}/cutset{1}.yml -w ${weightsFile} -o /home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/Efficiency/Efficiency_{1}.root ::: ${CutSets[@]}