import os
import copy
import yaml

if __name__ == "__main__":
    proj_config_file = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/config_projection_data_FT0_2023.yaml"
    central_fit_file = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/data/VsMult/FT0M/2023_h1_fix_template_params_sigma/mass_fits.root" # to fix the sigma
    cutset_dir = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/VsMult/configs/cutsets"
    fit_config_file = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/config_fit_data_FT0_multiplicity_2023_fix_template_params_mc_sigma.yml"
    proj_out_dir = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/VsMult/Projections_RawYields"
    fit_in_dir = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/VsMult/Projections_RawYields" # different from proj_out_dir since we merge the projections
    fit_out_dir = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/VsMult/Projections_RawYields/doublecb" # different from proj_out_dir since we merge the projections
    eff_config_file = "/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/configs/config_efficiency_h1_FT0.yaml"
    eff_out_dir = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/VsMult/Efficiency"
    output_dir = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/VsMult/configs"

    with open(proj_config_file, "r", encoding="utf-8") as file:
        proj_config = yaml.safe_load(file)
    with open(fit_config_file, "r", encoding="utf-8") as file:
        fit_config = yaml.safe_load(file)
    with open(eff_config_file, "r", encoding="utf-8") as file:
        eff_config = yaml.safe_load(file)

    for i_cut, cutset_file in enumerate(sorted(os.listdir(cutset_dir))):
        suffix = cutset_file.split(".yml")[0].split("cutset_ML_")[1]
        cutset_file_path = os.path.join(cutset_dir, cutset_file)
        proj_config_mod = copy.deepcopy(proj_config)
        proj_config_mod["cut_set_file_name"] = cutset_file_path
        proj_config_mod["output"]["directory"] = proj_out_dir + f'/{suffix[1:]}'
        #proj_config_mod["output"]["file_names"] = f"projection_data{suffix}.root"
        out_file_name = f"config_projection_data{suffix}.yaml"
        out_file_path = os.path.join(output_dir, "projections", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(proj_config_mod, file, default_flow_style=False)

        fit_config_mod = copy.deepcopy(fit_config)
        data_proj_file = os.path.join(fit_in_dir, f"projection_data{suffix}.root")
        data_proj_file = os.path.join(fit_in_dir, f"projection_data{suffix}.root")
        fit_config_mod["inputs"]["data"]= data_proj_file
        fit_config_mod["inputs"]["cutset"] = cutset_file_path
        fit_config_mod["fit_configs"]["signal"]["fix_ds_sigma"] = [True]*6
        fit_config_mod["fit_configs"]["signal"]["fix_dplus_sigma"] = [True]*6
        fit_config_mod["fit_configs"]["signal"]["ds_sigma"] = central_fit_file
        fit_config_mod["fit_configs"]["signal"]["dplus_sigma"] = central_fit_file
        fit_config_mod["outputs"]["directory"] = fit_out_dir
        fit_config_mod["outputs"]["save_all_fits"] = True
        fit_config_mod["outputs"]["suffix"] = suffix
        out_file_name = f"config_fit{suffix}.yaml"
        out_file_path = os.path.join(output_dir, "fits", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(fit_config_mod, file, default_flow_style=False)
        
        eff_config_mod = copy.deepcopy(eff_config)
        eff_config_mod["inputs"]["cuts_file_name"] = cutset_file_path
        eff_config_mod["weights"]["apply"] = True
        eff_out_dir = os.path.join(eff_out_dir)
        eff_config_mod["output_dir"] = eff_out_dir
        eff_config_mod["suffix"] = suffix
        out_file_name = f"config_eff{suffix}.yaml"
        out_file_path = os.path.join(output_dir, "efficiencies", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(eff_config_mod, file, default_flow_style=False)
