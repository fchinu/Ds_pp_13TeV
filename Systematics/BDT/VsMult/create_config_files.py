import os
import copy
import yaml

if __name__ == "__main__":
    proj_config_file_data = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Projections_RawYields/data/doublecb/Projections/config_projection_data.yaml"
    proj_config_file_mc = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Projections_RawYields/MC/doublecb/Projections/w_bdt/config_projection_mc.yaml"
    central_fit_file = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Projections_RawYields/data/doublecb/RawYields/mass_fits.root" # to fix the sigma
    cutset_dir = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/configs/cutsets"
    fit_config_file = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Projections_RawYields/data/doublecb/RawYields/config_fit_data.yml"
    proj_out_dir_data = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/Projections_RawYields"
    proj_out_dir_mc = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/Projections_RawYields/MC"
    fit_in_dir = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/Projections_RawYields"
    fit_out_dir = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/Projections_RawYields/doublecb"
    eff_config_file = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Efficiency/efficiency/default/config_efficiency.yml"
    eff_out_dir = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/Efficiency"
    output_dir = "/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Systematics/BDT/VsMult/configs"

    with open(proj_config_file_data, "r", encoding="utf-8") as file:
        proj_config_data = yaml.safe_load(file)
    with open(proj_config_file_mc, "r", encoding="utf-8") as file:
        proj_config_mc = yaml.safe_load(file)
    with open(fit_config_file, "r", encoding="utf-8") as file:
        fit_config = yaml.safe_load(file)
    with open(eff_config_file, "r", encoding="utf-8") as file:
        eff_config = yaml.safe_load(file)
    
    #create output dirs
    if not os.path.exists(os.path.join(output_dir, "projections", "data")):
        os.makedirs(os.path.join(output_dir, "projections", "data"))
    if not os.path.exists(os.path.join(output_dir, "projections", "mc")):
        os.makedirs(os.path.join(output_dir, "projections", "mc"))
    if not os.path.exists(os.path.join(output_dir, "fits")):
        os.makedirs(os.path.join(output_dir, "fits"))
    if not os.path.exists(os.path.join(output_dir, "efficiencies")):
        os.makedirs(os.path.join(output_dir, "efficiencies"))

    for i_cut, cutset_file in enumerate(sorted(os.listdir(cutset_dir))):
        suffix = cutset_file.split(".yml")[0].split("cutset_ML_")[1]
        cutset_file_path = os.path.join(cutset_dir, cutset_file)

        # Projections - Data
        proj_config_mod = copy.deepcopy(proj_config_data)
        proj_config_mod["inputs"]["cutset"] = cutset_file_path
        proj_config_mod["output"]["directory"] = proj_out_dir_data + f'/{suffix[1:]}'
        out_file_name = f"config_projection_data{suffix}.yaml"
        out_file_path = os.path.join(output_dir, "projections", "data", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(proj_config_mod, file, default_flow_style=False)

        # Projections - MC
        proj_config_mod = copy.deepcopy(proj_config_mc)
        proj_config_mod["inputs"]["cutset"] = cutset_file_path
        proj_config_mod["output"]["directory"] = proj_out_dir_mc + f'/{suffix[1:]}'
        out_file_name = f"config_projection_mc{suffix}.yaml"
        out_file_path = os.path.join(output_dir, "projections", "mc", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(proj_config_mod, file, default_flow_style=False)

        # Fits
        fit_config_mod = copy.deepcopy(fit_config)
        data_proj_file = os.path.join(fit_in_dir, f"projection_data{suffix}.root")
        data_proj_file = os.path.join(fit_in_dir, f"projection_data{suffix}.root")
        fit_config_mod["inputs"]["data"]= data_proj_file
        fit_config_mod["inputs"]["cutset"] = cutset_file_path

        # fix all parameters to the central fit
        for func in fit_config_mod["fit_configs"]["signal"]["par_init_limit"]:
            for par in func:
                func[par]["fix_to_file"] = [True]*6
                func[par]["fix_to_config_value"] = [False]*6
                func[par]["fix_to_mb"] = [False]*6

        # Change the normalisation of the correlated backgrounds
        # Assuming here only one background component
        fit_config_mod["fit_configs"]["bkg"]["templ_norm"]["backgrounds"][0]["file_norm"] = os.path.join(proj_out_dir_mc, suffix[1:], f"dplus_bkg.root")
        fit_config_mod["fit_configs"]["bkg"]["templ_norm"]["signal"]["file_norm"] = os.path.join(proj_out_dir_mc, suffix[1:], f"dplus.root")

        fit_config_mod["fit_configs"]["signal"]["file_for_params_fix"] = central_fit_file
        fit_config_mod["output"]["directory"] = fit_out_dir
        fit_config_mod["output"]["save_all_fits"] = True
        fit_config_mod["output"]["suffix"] = suffix
        out_file_name = f"config_fit{suffix}.yaml"
        out_file_path = os.path.join(output_dir, "fits", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(fit_config_mod, file, default_flow_style=False)
        
        eff_config_mod = copy.deepcopy(eff_config)
        eff_config_mod["inputs"]["cutset"] = cutset_file_path
        eff_config_mod["weights"]["pt"]["apply"] = True
        eff_out_dir = os.path.join(eff_out_dir)
        eff_config_mod["output_dir"] = eff_out_dir
        eff_config_mod["suffix"] = suffix
        out_file_name = f"config_eff{suffix}.yaml"
        out_file_path = os.path.join(output_dir, "efficiencies", out_file_name)
        with open(out_file_path, "w", encoding="utf-8") as file:
            yaml.dump(eff_config_mod, file, default_flow_style=False)
