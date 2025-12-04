#!/usr/bin/env python3

import os
import copy
import yaml
import json
import argparse
import sys

def generate_configs(trial_dir):
    """Generate configuration files for each cutset."""
    
    # Base paths - adjust these to your actual file locations
    base_dir = os.path.abspath(os.path.dirname(__file__))
    
    # Input configuration files (adjust these paths as needed)
    proj_config_file = os.path.join(trial_dir, "config_projection.yml")
    proj_config_file_prompt_enhanced = os.path.join(trial_dir, "config_projection_prompt_enhanced.yml")
    fit_config_file = os.path.join(trial_dir, "config_fit.yml")
    fit_config_file_prompt_enhanced = os.path.join(trial_dir, "config_fit_prompt_enhanced.yml")
    eff_config_file = os.path.join(trial_dir, "config_efficiency.yml")
    cutvar_config_file = os.path.join(trial_dir, "config_cutvar.json")
    cutset_config_file = os.path.join(trial_dir, "cutset.yml")
    
    # Directories
    cutset_dir = os.path.join(trial_dir, "configs", "cutsets")
    output_dir = os.path.join(trial_dir, "configs")
    
    # Output directories for actual processing
    proj_out_dir = os.path.join(trial_dir, "Projections")
    proj_prompt_enhanced_out_dir = os.path.join(trial_dir, "Projections", "prompt_enhanced")
    fit_in_dir = os.path.join(trial_dir, "Projections")
    fit_out_dir = os.path.join(trial_dir, "RawYields")
    eff_out_dir = os.path.join(trial_dir, "Efficiencies")
    ds_frac_out_dir = os.path.join(trial_dir, "Fractions", "Ds")
    dplus_frac_out_dir = os.path.join(trial_dir, "Fractions", "Dplus")
    
    # Create output subdirectories
    os.makedirs(os.path.join(output_dir, "Projections"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "RawYields"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "Efficiencies"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "Fractions", "Ds"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "Fractions", "Dplus"), exist_ok=True)
    
    # Check if required files exist
    required_files = [proj_config_file, proj_config_file_prompt_enhanced, fit_config_file, fit_config_file_prompt_enhanced, eff_config_file, cutvar_config_file, cutset_config_file]
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"Error: Required file {file_path} not found!")
            return False
    
    if not os.path.exists(cutset_dir):
        print(f"Error: Cutset directory {cutset_dir} not found!")
        return False
    
    try:
        # Load configuration files
        with open(proj_config_file, "r", encoding="utf-8") as file:
            proj_config = yaml.safe_load(file)
        with open(proj_config_file_prompt_enhanced, "r", encoding="utf-8") as file:
            proj_config_prompt_enhanced = yaml.safe_load(file)
        with open(fit_config_file, "r", encoding="utf-8") as file:
            fit_config = yaml.safe_load(file)
        with open(fit_config_file_prompt_enhanced, "r", encoding="utf-8") as file:
            fit_config_prompt_enhanced = yaml.safe_load(file)
        with open(eff_config_file, "r", encoding="utf-8") as file:
            eff_config = yaml.safe_load(file)
        with open(eff_config_file, "r", encoding="utf-8") as file:
            eff_config = yaml.safe_load(file)
        with open(cutvar_config_file, "r", encoding="utf-8") as file:
            cutvar_config = json.load(file)
        with open(cutset_config_file, "r", encoding="utf-8") as file:
            cutset_config = yaml.safe_load(file)
    except Exception as e:
        print(f"Error loading configuration files: {e}")
        return False
    
    # Modify configurations for prompt enhanced
    proj_config_prompt_enhanced["inputs"]["cutset"] = os.path.join(trial_dir, "cutset_prompt_enhanced.yml")
    proj_config_prompt_enhanced["output"]["directory"] = proj_prompt_enhanced_out_dir
    fit_config_prompt_enhanced["inputs"]["data"] = os.path.join(fit_in_dir, "projection_prompt_enhanced.root")
    fit_config_prompt_enhanced["inputs"]["cutset"] = os.path.join(trial_dir, "cutset_prompt_enhanced.yml")
    fit_config_prompt_enhanced["output"]["directory"] = fit_out_dir
    fit_config_prompt_enhanced["output"]["save_all_fits"] = True
    fit_config_prompt_enhanced["output"]["suffix"] = "_prompt_enhanced"

    # Dump the modified prompt enhanced configurations
    with open(proj_config_file_prompt_enhanced, "w", encoding="utf-8") as file:
        yaml.dump(proj_config_prompt_enhanced, file, default_flow_style=False)
    with open(fit_config_file_prompt_enhanced, "w", encoding="utf-8") as file:
        yaml.dump(fit_config_prompt_enhanced, file, default_flow_style=False)

    # Process each cutset file
    cutset_files = [f for f in os.listdir(cutset_dir) if f.endswith('.yml')]
    if not cutset_files:
        print(f"Warning: No cutset files found in {cutset_dir}")
        return True

    for cutset_file in sorted(cutset_files):
        # Extract suffix from filename
        suffix = cutset_file.split(".yml")[0].split("cutset_ML_")[1]
        variation = suffix.split("_")[0]
        suffix = suffix.split(variation, 1)[1]
        cutset_file_path = os.path.join(cutset_dir, cutset_file)

        # Generate projection config
        proj_config_mod = copy.deepcopy(proj_config)
        proj_config_mod["inputs"]["cutset"] = cutset_file_path
        proj_config_mod["output"]["directory"] = os.path.join(proj_out_dir, suffix[1:])

        out_file_name = f"config_projection_data{suffix}.yml"
        out_file_path = os.path.join(output_dir, "Projections", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(proj_config_mod, file, default_flow_style=False)

        # Generate fit config
        fit_config_mod = copy.deepcopy(fit_config)
        data_proj_file = os.path.join(fit_in_dir, f"projection_data{suffix}.root")
        fit_config_mod["inputs"]["data"] = data_proj_file
        fit_config_mod["inputs"]["cutset"] = cutset_file_path
        fit_config_mod["fit_configs"]["signal"]["file_for_params_fix"] = os.path.join(fit_out_dir, "mass_fits_prompt_enhanced.root")
        fit_config_mod["output"]["directory"] = fit_out_dir
        fit_config_mod["output"]["save_all_fits"] = True
        fit_config_mod["output"]["suffix"] = suffix

        out_file_name = f"config_fit{suffix}.yml"
        out_file_path = os.path.join(output_dir, "RawYields", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(fit_config_mod, file, default_flow_style=False)

        # Generate efficiency config
        eff_config_mod = copy.deepcopy(eff_config)
        eff_config_mod["inputs"]["cutset"] = cutset_file_path
        eff_config_mod["output_dir"] = eff_out_dir
        eff_config_mod["suffix"] = suffix

        out_file_name = f"config_eff{suffix}.yml"
        out_file_path = os.path.join(output_dir, "Efficiencies", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(eff_config_mod, file, default_flow_style=False)

        print(f"Generated configs for cutset: {cutset_file}")
        
    # Generate cut variation configs
    cents = list(zip(cutset_config["cent"]["min"], cutset_config["cent"]["max"]))

    # Remove 0-100 if present in the cutset config
    index = cents.index((0, 100))
    cents.pop(index)

    # Minimum bias first
    config_ds_mb = copy.deepcopy(cutvar_config)
    config_ds_mb["rawyields"]["inputdir"] = fit_out_dir
    config_ds_mb["rawyields"]["histoname"] = "h_raw_yields_ds_0_100"
    config_ds_mb["efficiencies"]["inputdir"] = eff_out_dir
    config_ds_mb["efficiencies"]["histonames"] = {"prompt": "eff_DsPrompt_cent_0_100", "nonprompt": "eff_DsNonPrompt_cent_0_100"}
    config_ds_mb["output"]["directory"] = ds_frac_out_dir
    config_ds_mb["output"]["file"] = "CutVarDs_pp13TeV_MB.root"

    config_dplus_mb = copy.deepcopy(config_ds_mb)
    config_dplus_mb["rawyields"]["histoname"] = "h_raw_yields_dplus_0_100"
    config_dplus_mb["efficiencies"]["histonames"] = {"prompt": "eff_DplusPrompt_cent_0_100", "nonprompt": "eff_DplusNonPrompt_cent_0_100"}
    config_dplus_mb["output"]["directory"] = dplus_frac_out_dir
    config_dplus_mb["output"]["file"] = "CutVarDplus_pp13TeV_MB.root"

    out_file_path = os.path.join(output_dir, "Fractions", "Ds")
    with open(os.path.join(out_file_path, "config_ds_mb.json"), "w") as f:
        json.dump(config_ds_mb, f)
    out_file_path = os.path.join(output_dir, "Fractions", "Dplus")
    with open(os.path.join(out_file_path, "config_dplus_mb.json"), "w") as f:
        json.dump(config_dplus_mb, f)

    # Generate configurations for each centrality bin
    for cent_min, cent_max in cents:
        config_ds_cent = copy.deepcopy(config_ds_mb)
        config_ds_cent["rawyields"]["histoname"] = config_ds_mb["rawyields"]["histoname"].replace("0_100", f"{cent_min}_{cent_max}")
        config_ds_cent["efficiencies"]["histonames"] = {name: config_ds_mb["efficiencies"]["histonames"][name].replace("0_100", f"{cent_min}_{cent_max}") for name in config_ds_mb["efficiencies"]["histonames"]}
        config_ds_cent["output"]["file"] = config_ds_mb["output"]["file"].replace("MB", f"{cent_min}_{cent_max}")
        
        config_dplus_cent = copy.deepcopy(config_dplus_mb)
        config_dplus_cent["rawyields"]["histoname"] = config_dplus_mb["rawyields"]["histoname"].replace("0_100", f"{cent_min}_{cent_max}")
        config_dplus_cent["efficiencies"]["histonames"] = {name: config_dplus_mb["efficiencies"]["histonames"][name].replace("0_100", f"{cent_min}_{cent_max}") for name in config_dplus_mb["efficiencies"]["histonames"]}
        config_dplus_cent["output"]["file"] = config_dplus_mb["output"]["file"].replace("MB", f"{cent_min}_{cent_max}")

        # Save the configurations
        out_file_path = os.path.join(output_dir, "Fractions", "Ds")
        with open(os.path.join(out_file_path, f"config_ds_{cent_min}_{cent_max}.json"), "w") as f:
            json.dump(config_ds_cent, f)
        out_file_path = os.path.join(output_dir, "Fractions", "Dplus")
        with open(os.path.join(out_file_path, f"config_dplus_{cent_min}_{cent_max}.json"), "w") as f:
            json.dump(config_dplus_cent, f)
    
    print(f"Successfully generated configs for {len(cutset_files)} cutsets")
    return True

def main():
    parser = argparse.ArgumentParser(description="Generate configuration files for each cutset")
    parser.add_argument("--trial-dir", required=True, help="Trial directory path")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.trial_dir):
        print(f"Error: Trial directory {args.trial_dir} does not exist!")
        sys.exit(1)
    
    success = generate_configs(args.trial_dir)
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()