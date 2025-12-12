#!/bin/bash

# RunRawYields.sh
# Usage: bash RunRawYields.sh [TRIAL_DIR]

# Get trial directory (default to current directory if not provided)
TRIAL_DIR=${1:-.}

# Base paths - adjust these to your actual paths
PROJ_PATH="../../Projections_RawYields"

# Maximum number of retries
MAX_RETRIES=3

# Configuration paths (now relative to trial directory)
ConfigFitPromptEnhanced="${TRIAL_DIR}/config_fit_prompt_enhanced.yml"
FitConfigDir="${TRIAL_DIR}/configs/RawYields"
ProjectionsDir="${TRIAL_DIR}/Projections"
RawYieldsDir="${TRIAL_DIR}/RawYields"

# Create output directory if it doesn't exist
mkdir -p ${RawYieldsDir}

echo "=== Running Raw Yields Extraction ==="
echo "Trial directory: ${TRIAL_DIR}"
echo "Raw yields output: ${RawYieldsDir}"

# Build array of fit config files
declare -a FitConfigs=()
for filename in ${FitConfigDir}/*.yml; do
    if [[ -f "$filename" ]]; then
        tmp_name="$(basename -- ${filename} .yml)"
        tmp_name=${tmp_name:11}  # Adjust this number based on your naming convention
        FitConfigs+=("${tmp_name}")
    fi
done

echo "Found ${#FitConfigs[@]} fit configs"

# Run prompt enhanced fit
if [[ -f "${ConfigFitPromptEnhanced}" ]]; then
    echo "Running prompt enhanced fit..."

    for ((attempt=1; attempt<=MAX_RETRIES; attempt++)); do
        echo "Attempt $attempt of $MAX_RETRIES..."
        nice -n 15 python3 "${PROJ_PATH}/get_raw_yields.py" "${ConfigFitPromptEnhanced}" \
            >& "${RawYieldsDir}/raw_yields_prompt_enhanced.log"

        if [[ $? -eq 0 ]]; then
            echo "Prompt enhanced fit succeeded on attempt $attempt."
            break
        else
            echo "Prompt enhanced fit failed on attempt $attempt."
            if [[ $attempt -lt MAX_RETRIES ]]; then
                echo "Retrying..."
            else
                echo "All $MAX_RETRIES attempts failed. Check the log at ${RawYieldsDir}/raw_yields_prompt_enhanced.log"
                exit 1 
            fi
        fi
    done
else
    echo "Warning: ${ConfigFitPromptEnhanced} not found, skipping prompt enhanced fit"
fi

# Run parallel fits for all cutsets
if [[ ${#FitConfigs[@]} -gt 0 ]]; then
    echo "Running parallel fits for ${#FitConfigs[@]} configs..."
    nice -n 15 parallel --retries ${MAX_RETRIES} -j3 \
        --joblog "${RawYieldsDir}/parallel_joblog.txt" \
        --results "${RawYieldsDir}/logs" \
        python3 ${PROJ_PATH}/get_raw_yields.py ${FitConfigDir}/config_fit_{1}.yml \
        ::: "${FitConfigs[@]}"
    echo "Parallel fits log saved to ${RawYieldsDir}/raw_yields.log"
else
    echo "No fit configs found, skipping parallel fits"
fi

echo "=== Raw Yields Extraction Complete ==="