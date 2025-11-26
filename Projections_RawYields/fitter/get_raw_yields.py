import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position
from tqdm import tqdm
import uproot 
import yaml
from pypdf import PdfWriter
import tensorflow as tf
tf.config.threading.set_intra_op_parallelism_threads(20)
tf.config.threading.set_inter_op_parallelism_threads(20)
from fit_handler import FitHandler, BRInfo, CorrelatedBackground, CorrelatedBackgroundConfig, FitConfig


def merge_pdfs(cfg):
    """Merge individual PDF files into a single PDF."""
    if "pdf" not in cfg["output"]["formats"]:
        return

    # if output already exists, delete it
    if os.path.exists(os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits", "ds_mass_merged.pdf")):
        os.remove(os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits", "ds_mass_merged.pdf"))
    if os.path.exists(os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits", "ds_massres_merged.pdf")):
        os.remove(os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits", "ds_massres_merged.pdf"))

    output_dir = os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits")
    pdf_files = [f for f in os.listdir(output_dir) if f.endswith('.pdf') and f.startswith('ds_mass_')]
    pdf_files = sorted(pdf_files, key=lambda x: (
        int(x.split('_')[3]), 
        int(x.split('_')[4]) if 'cent' in x else -1
    ))

    pdf_files_residuals = [f for f in os.listdir(output_dir) if f.endswith('.pdf') and f.startswith('ds_massres_')]
    pdf_files_residuals = sorted(pdf_files_residuals, key=lambda x: (
        int(x.split('_')[3]),
        int(x.split('_')[4]) if 'cent' in x else -1
    ))

    merger = PdfWriter()
    for pdf in pdf_files:
        merger.append(os.path.join(output_dir, pdf))
    merger.write(os.path.join(output_dir, "ds_mass_merged.pdf"))
    merger.close()

    merger = PdfWriter()
    for pdf in pdf_files_residuals:
        merger.append(os.path.join(output_dir, pdf))
    merger.write(os.path.join(output_dir, "ds_massres_merged.pdf"))
    merger.close()

def run_fit(pt_info: Tuple[float, float], fit_config: FitConfig) -> Dict:
    """
    Run the fitting procedure based on the provided FitConfig.

    Parameters:
    fit_config (FitConfig): Configuration for the fit.

    Returns:
    Dict: Results of the fitting procedure.
    """
    pt_min, pt_max = pt_info
    fit_handler = FitHandler(fit_config)
    results = fit_handler.get_results()
    return (
        pt_min,
        pt_max,
        results
    )

def get_corr_bkg_config(cfg: Dict, i_pt: int) -> CorrelatedBackgroundConfig:
    """
    Create a CorrelatedBackgroundConfig object based on the provided configuration.

    Parameters:
    cfg (dict): Configuration dictionary.
    pt_min (float): Minimum pT value.
    pt_max (float): Maximum pT value.

    Returns:
    CorrelatedBackgroundConfig: Configured CorrelatedBackgroundConfig object.
    """
    bkg_cfg = cfg["fit_configs"]["bkg"]
    if bkg_cfg["use_bkg_templ"][i_pt] == False:
        return None
    
    bkg_norm = bkg_cfg["templ_norm"]

    return CorrelatedBackgroundConfig(
        fix_to_file=bkg_norm["fix_to_file"][i_pt],
        fix_to_mb=bkg_norm["fix_to_mb"][i_pt],
        fix_with_br=bkg_norm["fix_with_br"][i_pt],
        file_name_for_fix=bkg_norm["file_name_for_fix"],
        hist_name_for_fix=bkg_norm["hist_name_for_fix"],
        backgrounds=[
            CorrelatedBackground(
                name=bkg["name"],
                file_norm=bkg["file_norm"],
                norm_hist_name=bkg["norm_hist_name"],
                template_file=bkg["template_file"],
                template_hist_name=bkg["template_hist_name"],
                br=BRInfo(pdg=bkg["br"]["pdg"], simulations=bkg["br"]["simulations"])
            )
            for bkg in bkg_norm["backgrounds"]
        ],
        signal_norm_file=bkg_norm["signal"]["file_norm"],
        signal_hist_name=bkg_norm["signal"]["hist_name"],
        signal_br=BRInfo(pdg=bkg_norm["signal"]["br"]["pdg"], simulations=bkg_norm["signal"]["br"]["simulations"])
    )

def get_parameter_setup(cfg: List[Dict], i_pt: int, mb_result: Dict = None, result_sigma_fix: Dict = None) -> List[Dict]:
    """
    Get parameter setup based on the provided configuration.

    Parameters:
    cfg (list): List of configuration dictionaries.
    i_pt (int): Index of the pT bin.

    Returns:
    List[Dict]: Parameter setup dictionary.
    """
    param_cfg = cfg["fit_configs"]["signal"]["par_init_limit"]
    functions_setup = []
    for i_func, func_setup in enumerate(param_cfg):  # For each signal function
        param_setup = {}
        for par_name, par_values in func_setup.items():  # For each parameter in the function
            if par_values["fix_to_mb"][i_pt] and mb_result is not None:
                param_setup[par_name] = {
                    "init": mb_result[par_name][i_func][0],
                    "min": par_values["min"][i_pt],
                    "max": par_values["max"][i_pt],
                    "fix_to_config_value": True,
                    "fix_to_file": False,
                }
            elif par_name == "sigma" \
                and i_func == 1 \
                    and cfg["fit_configs"]["signal"]["fix_sigma_dplus_to_ds"][i_pt] \
                        and result_sigma_fix is not None:
                ratio_sigma_dplus_to_ds = -1.
                if isinstance(cfg["fit_configs"]["signal"]["ratio_sigma_dplus_to_ds"], list):
                    ratio_sigma_dplus_to_ds = cfg["fit_configs"]["signal"]["ratio_sigma_dplus_to_ds"][i_pt]
                elif isinstance(cfg["fit_configs"]["signal"]["ratio_sigma_dplus_to_ds"], (int, float)):
                    ratio_sigma_dplus_to_ds = cfg["fit_configs"]["signal"]["ratio_sigma_dplus_to_ds"]

                param_setup[par_name] = {
                    "init": result_sigma_fix[par_name][0][0] * ratio_sigma_dplus_to_ds,
                    "min": result_sigma_fix[par_name][0][0] * ratio_sigma_dplus_to_ds - 0.01,
                    "max": result_sigma_fix[par_name][0][0] * ratio_sigma_dplus_to_ds + 0.01,
                    "fix_to_config_value": True,
                    "fix_to_file": False
                }
            else:
                param_setup[par_name] = {
                    "init": par_values["init"][i_pt],
                    "min": par_values["min"][i_pt],
                    "max": par_values["max"][i_pt],
                    "fix_to_config_value": par_values["fix_to_config_value"][i_pt],
                    "fix_to_file": par_values["fix_to_file"][i_pt]
                }
        functions_setup.append(param_setup)
    return functions_setup

def get_config(
        cfg: Dict,
        pt_info: Tuple[int, float, float],
        cent_info: Tuple[List[float], List[float]],
        mb_results: Dict = None,
        result_sigma_fix: Dict = None
) -> FitConfig:
    """
    Create a FitConfig object based on the provided configuration.

    Parameters:
    cfg (dict): Configuration dictionary.
    pt_min (float): Minimum pT value.
    pt_max (float): Maximum pT value.

    Returns:
    FitConfig: Configured FitConfig object.
    """
    i_pt, pt_min, pt_max = pt_info
    cent_min, cent_max = cent_info

    ratio_sigma_dplus_to_ds = -1.
    if isinstance(cfg["fit_configs"]["signal"]["ratio_sigma_dplus_to_ds"], list):
        ratio_sigma_dplus_to_ds = cfg["fit_configs"]["signal"]["ratio_sigma_dplus_to_ds"][i_pt]
    elif isinstance(cfg["fit_configs"]["signal"]["ratio_sigma_dplus_to_ds"], (int, float)):
        ratio_sigma_dplus_to_ds = cfg["fit_configs"]["signal"]["ratio_sigma_dplus_to_ds"]

    return FitConfig(
        pt_min=pt_min,
        pt_max=pt_max,
        cent_min=cent_min,
        cent_max=cent_max,
        mass_range=[
            cfg["fit_configs"]["mass"]["mins"][i_pt],
            cfg["fit_configs"]["mass"]["maxs"][i_pt]
        ],
        signal_pdfs=cfg["fit_configs"]["signal"]["signal_funcs"][i_pt],
        bkg_pdfs=cfg["fit_configs"]["bkg"]["bkg_funcs"][i_pt],
        rebin=cfg["fit_configs"]["rebin"][i_pt],
        param_setup=get_parameter_setup(cfg, i_pt, mb_results, result_sigma_fix),
        data_path=cfg["inputs"]["data"],
        file_for_params_fix=cfg["fit_configs"]["signal"]["file_for_params_fix"],
        suffix_hist_for_params_fix=cfg["fit_configs"]["signal"]["suffix_hist_for_params_fix"],
        fix_dplus_sigma_to_ds=cfg["fit_configs"]["signal"]["fix_sigma_dplus_to_ds"][i_pt],
        ratio_sigma_dplus_to_ds=ratio_sigma_dplus_to_ds,
        correlated_bkg=get_corr_bkg_config(cfg, i_pt),
        draw_figures=cfg["output"]["save_all_fits"],
        draw_formats=cfg["output"]["formats"],
        output_dir=os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits")
    )

    # We fix parameters to mb sample if needed



def load_inputs(cfg, cut_set):
    """
    Load invariant event and collision histograms from a ROOT file.

    Parameters:
    cfg (dict): Configuration dictionary containing the path to the input data file.
    cut_set (dict): Dictionary containing the cuts applied to the data.

    Returns:
    tuple: A tuple containing:
        - h_ev (uproot.models.TH1): Event histogram.
        - h_collisions (uproot.models.TH1): Collisions histogram.
    """
    pt_mins = cut_set["pt"]["min"]
    pt_maxs = cut_set["pt"]["max"]
    cent_mins, cent_maxs = None, None
    if "cent" in cut_set:
        cent_mins = cut_set["cent"]["min"] 
        cent_maxs = cut_set["cent"]["max"] 

    with uproot.open(cfg["inputs"]["data"]) as f:
        h_ev = f["h_ev"]
        h_collisions = f["h_coll_rebinned"]

    return h_ev, h_collisions

def _mb_fit_is_needed(cfg, cut_set) -> bool:
    """Determines if a minimum-bias fit is needed based on the configuration."""
    if (0, 100) in list(zip(cut_set["cent"]["min"], cut_set["cent"]["max"])):
        return True
    for p_cfg in cfg["fit_configs"]["signal"]["par_init_limit"]:
        for _, settings in p_cfg.items():
            if settings["fix_to_mb"]:
                return True

    if cfg["fit_configs"]["correlated_bkg"] is not None:
        if cfg["fit_configs"]["correlated_bkg"]["fix_to_mb"]:
            return True
    return False

def fit(config_file_name):
    """
    Perform fitting procedure based on the provided configuration file.

    Parameters:
        config_file_name (str): Path to the YAML configuration file.
    """
    # load config
    with open(config_file_name, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    # load cut set
    with open(cfg["inputs"]["cutset"], "r", encoding="utf-8") as f:
        cut_set = yaml.safe_load(f)

    # zfit.run.set_cpus_explicit(
    #     intra=cfg['zfit_cpus']['intra'],
    #     inter=cfg['zfit_cpus']['inter']
    # )

    # load inputs
    h_ev, h_collisions = load_inputs(cfg, cut_set)
    n_ev = sum(h_ev.values())

    cent_mins, cent_maxs = None, None
    if "cent" in cut_set:
        cent_mins = cut_set["cent"]["min"]
        cent_maxs = cut_set["cent"]["max"]

    pt_mins = cut_set["pt"]["min"]
    pt_maxs = cut_set["pt"]["max"]

    # # create output handler
    # h_handler = HistHandler(
    #     cut_set["pt"]["min"], cut_set["pt"]["max"],
    #     cent_mins, cent_maxs
    # )
    # h_handler.set_n_ev(n_ev)
    # h_handler.set_h_collisions(h_collisions)

    # fit_configs = create_fit_configs(cfg, cut_set)

    # # recreate root file if needed
    # if "root" in cfg["output"]["formats"]:
    #     output_dir = os.path.join(
    #         os.path.expanduser(cfg["output"]["directory"]),
    #         "fits"
    #     )
    #     uproot.recreate(f"{output_dir}/fits_{cfg['output']['suffix']}.root")


    # if needed, run minimum-bias fit first
    mb_results = {}
    if _mb_fit_is_needed(cfg, cut_set):
        with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
            futures = []
            for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
                fit_config = get_config(cfg, (i_pt, pt_min, pt_max), (0, 100))
                futures.append(executor.submit(run_fit, (pt_min, pt_max), fit_config))

            for future in tqdm(as_completed(futures), total=len(futures)):
                pt_min, pt_max, results = future.result()
                mb_results[(pt_min, pt_max)] = results
    
    with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
        futures = []
        for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
            if cfg["fit_configs"]["signal"]["fix_sigma_dplus_to_ds"][i_pt]:
                fit_config = get_config(cfg, (i_pt, pt_min, pt_max), (0, 100), None, mb_results[(pt_min, pt_max)])
                futures.append(executor.submit(run_fit, (pt_min, pt_max), fit_config))

        for future in tqdm(as_completed(futures), total=len(futures)):
            pt_min, pt_max, results = future.result()
            mb_results[(pt_min, pt_max)] = results
            

    ## We get all the configurations for fits

    fit_configs = []
    for (pt_min, pt_max), mb_result in mb_results.items():
        for cent_min, cent_max in zip(cent_mins, cent_maxs):
            if cent_min == 0 and cent_max == 100:
                continue
            i_pt = pt_mins.index(pt_min)
            fit_config = get_config(cfg, (i_pt, pt_min, pt_max), (cent_min, cent_max), mb_result)
            fit_configs.append(fit_config)

    with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
        futures = []
        for fit_config in fit_configs:
            pt_min = fit_config.pt_min
            pt_max = fit_config.pt_max
            futures.append(executor.submit(
                run_fit, (pt_min, pt_max), fit_config
            ))

        for future in tqdm(as_completed(futures), total=len(futures)):
            pt_min, pt_max, results = future.result()
            # process results as needed

    with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
        futures = []
        for i_pt, ((pt_min, pt_max), mb_result) in enumerate(mb_results.items()):
            needs_further_fit = False
            if cfg["fit_configs"]["signal"]["fix_sigma_dplus_to_ds"][i_pt]:
                for i_func, func_par_setup in enumerate(cfg["fit_configs"]["signal"]["par_init_limit"]):
                    if "sigma" in func_par_setup and not func_par_setup["sigma"]["fix_to_mb"][i_pt]:
                        needs_further_fit = True
            if needs_further_fit:
                fit_config = get_config(cfg, (i_pt, pt_min, pt_max), (0, 100), None, mb_results[(pt_min, pt_max)])
                futures.append(executor.submit(run_fit, (pt_min, pt_max), fit_config))

            for future in tqdm(as_completed(futures), total=len(futures)):
                pt_min, pt_max, results = future.result()
                mb_results[(pt_min, pt_max)] = results

    merge_pdfs(cfg)

    #     mb_fit_config = get_config(
    #         cfg,
    #         (0, pt_mins[0], pt_maxs[0]),
    #         ( [cent_mins[mb_i_pt]], [cent_maxs[mb_i_pt]] )
    #     )
    #     print("Running minimum-bias fit first with config:", mb_fit_config)
    #     mb_fit_handler = FitHandler(mb_fit_config)
    #     mb_results = mb_fit_handler.get_results()
    #     print("Minimum-bias fit results:", mb_results)


    # # Perform fits in parallel
    # results = {}
    # with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
    #     futures = []
    #     for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
    #         fit_config = get_config(cfg, (i_pt, pt_min, pt_max), (cent_mins, cent_maxs))
    #         futures.append(executor.submit(run_fit, (i_pt, pt_min, pt_max), results, fit_config))

    #     for future in tqdm(as_completed(futures), total=len(futures)):
    #         print("results:", results)
    #         pt_min, pt_max, figs, axes, cents, results = future.result()
    #         # figs, axes, cents = results[(pt_min, pt_max)].draw_figures(
    #         #     style="ATLAS",
    #         #     show_extra_info=True,
    #         #     figsize=(8, 8), extra_info_loc=["lower left", "upper left"],
    #         #     legend_loc="upper right",
    #         #     axis_title=r"$M(\mathrm{KK\pi})$ GeV$/c^2$"
    #         # )
    #         for fig, cent in zip(figs, cents):
    #             fig.savefig(f"fit_{pt_min}_{pt_max}_{cent}.png")

    # # Perform fits sequentially
    # for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
    #     fit_config = get_config(cfg, (i_pt, pt_min, pt_max), (cent_mins, cent_maxs))
    #     print(fit_config)
    #     fit_handler = run_fit(fit_config)



    # if "root" in cfg["output"]["formats"]:
    #     # get all the partial root files in the output directory
    #     output_dir = os.path.join(
    #         os.path.expanduser(cfg["output"]["directory"]),
    #         "fits"
    #     )
    #     partial_files = [f for f in os.listdir(output_dir) if f.endswith('_partial.root') and f.startswith('fits_')]
    #     output_dir = os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits")
    #     out_root_file = f"{output_dir}/fits_{cfg['output']['suffix']}.root"
    #     print("hadd -f " + out_root_file + " " + " ".join([os.path.join(output_dir, f) for f in partial_files]))
    #     os.system("hadd -f " + out_root_file + " " + " ".join([os.path.join(output_dir, f) for f in partial_files]))
    #     os.system("rm " + " ".join([os.path.join(output_dir, f) for f in partial_files]))

    # # Save results
    # output_df = []
    # for result, fit_config in results:
    #     out_dict = result.result()
    #     out_dict.update(fit_config)
    #     output_df.append(out_dict)
    # for result, fit_config in results_fracs:
    #     out_dict = result.result()
    #     out_dict.update(fit_config)
    #     output_df.append(out_dict)
    # output_df = pd.DataFrame(output_df)
    
    # output_df.to_parquet(os.path.join(
    #     os.path.expanduser(cfg["output"]["directory"]),
    #     f'fit_results{cfg["output"]["suffix"]}.parquet'
    # ), index=False)

    # h_handler.set_histos(output_df)
    # h_handler.dump_to_root(os.path.join(
    #     os.path.expanduser(cfg["output"]["directory"]),
    #     f'mass_fits{cfg["output"]["suffix"]}.root'
    # ))

    # merge_pdfs(cfg)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('config_file_name', metavar='text', default='')
    args = parser.parse_args()



    fit(args.config_file_name)