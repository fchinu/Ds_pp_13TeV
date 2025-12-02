"""Script to perform systematic uncertainties estimation on raw yield extraction"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import dataclasses
import itertools
import multiprocessing
from typing import Dict, List, Tuple
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"  # FATAL only
os.environ["XLA_FLAGS"] = "--xla_gpu_cuda_data_dir=/nonexistent"
import sys
sys.path.append(os.path.abspath(os.path.join(__file__, '../../../../utils/fitter')))
from tqdm import tqdm
import uproot
import yaml
import contextlib
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import AnchoredText
import seaborn as sns
import tensorflow as tf
tf.config.threading.set_intra_op_parallelism_threads(80)
tf.config.threading.set_inter_op_parallelism_threads(80)
from fit_handler import (  # pylint: disable=import-error
    FitHandler,
    BRInfo,
    CorrelatedBackground,
    CorrelatedBackgroundConfig,
    FitConfig
)

def draw_multitrial(df_multitrial, cfg, is_mb):  # pylint: disable=too-many-locals, too-many-statements, too-many-branches  # noqa: 501
    """
    Produce a plot with the results of the multitrial procedure.

    Parameters:
    - df_multitrial (DataFrame): DataFrame containing the multitrial data.
    - cfg (dict): Configuration dictionary.
    - pt_min (float): Minimum pt value.
    - pt_max (float): Maximum pt value.
    - idx_assigned_syst (int): Index of the assigned systematic uncertainty.

    Returns:
    - None
    """
    multitrial_cfg = cfg["multitrial"]

    df_multitrial = df_multitrial.query(f"chi2 < {multitrial_cfg['quality_selections']['chi2']}")

    # Apply significance selection
    mask = df_multitrial['significance'].apply(
        lambda x:
            x[0][0] > multitrial_cfg['quality_selections']['significance'] and \
                x[1][0] > multitrial_cfg['quality_selections']['significance']
        )
    df_multitrial = df_multitrial[mask]

    if not is_mb:
        dump_results_to_root(df_multitrial, cfg)

    df_pt_cent_cfg = df_multitrial[
        ["pt_min_cfg", "pt_max_cfg", "cent_min_cfg", "cent_max_cfg"]
    ].drop_duplicates()
    for _, row in df_pt_cent_cfg.iterrows():
        pt_min = row["pt_min_cfg"]
        pt_max = row["pt_max_cfg"]
        cent_min = row["cent_min_cfg"]
        cent_max = row["cent_max_cfg"]

        df_pt_cent = df_multitrial.query(
            f"pt_min_cfg == {pt_min} and pt_max_cfg == {pt_max} and\
                cent_min_cfg == {cent_min} and cent_max_cfg == {cent_max}"
        )
        raw_yields_ds, raw_yields_ds_unc, raw_yields_dplus, raw_yields_dplus_unc = [], [], [], []
        raw_yields_ds_nsigma = [[] for _ in multitrial_cfg["bincounting_nsigma"]]
        raw_yields_ds_nsigma_unc = [[] for _ in multitrial_cfg["bincounting_nsigma"]]
        raw_yields_dplus_nsigma = [[] for _ in multitrial_cfg["bincounting_nsigma"]]
        raw_yields_dplus_nsigma_unc = [[] for _ in multitrial_cfg["bincounting_nsigma"]]
        sigma_ds, sigma_ds_unc, sigma_dplus, sigma_dplus_unc = [], [], [], []
        for _, row in df_pt_cent.iterrows():
            raw_yields_ds.append(row["raw_yields"][0][0])
            raw_yields_ds_unc.append(row["raw_yields"][0][1])
            raw_yields_dplus.append(row["raw_yields"][1][0])
            raw_yields_dplus_unc.append(row["raw_yields"][1][1])
            for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
                raw_yields_ds_nsigma[i_nsigma].append(
                    row[f"raw_yields_bincounting_{nsigma}"][0][0]
                )
                raw_yields_ds_nsigma_unc[i_nsigma].append(
                    row[f"raw_yields_bincounting_{nsigma}"][0][1]
                )
                raw_yields_dplus_nsigma[i_nsigma].append(
                    row[f"raw_yields_bincounting_{nsigma}"][1][0]
                )
                raw_yields_dplus_nsigma_unc[i_nsigma].append(
                    row[f"raw_yields_bincounting_{nsigma}"][1][1]
                )
            sigma_ds.append(row["sigma"][0][0])
            sigma_ds_unc.append(row["sigma"][0][1])
            sigma_dplus.append(row["sigma"][1][0])
            sigma_dplus_unc.append(row["sigma"][1][1])
        n_trials = len(df_pt_cent)
        n_bincounts = len(multitrial_cfg['bincounting_nsigma'])
        x_axis_range = n_trials * (n_bincounts + 1) + 1

        fig, axs = plt.subplots(2, 2, figsize=(20, 15))
        with uproot.open(cfg["reference_fits"]) as f:
            h_rawy_ds = f[f"h_raw_yields_ds_{cent_min}_{cent_max}"]
            h_sigma_ds = f["h_sigma_ds_0_100"]
            h_rawy_dplus = f[f"h_raw_yields_dplus_{cent_min}_{cent_max}"]
            h_sigma_dplus = f["h_sigma_dplus_0_100"]

        i_pt = np.digitize((pt_min + pt_max) / 2, h_rawy_ds.axis().edges()) - 1
        central_rawy_ds = h_rawy_ds.values()[i_pt]
        central_rawy_ds_unc = h_rawy_ds.errors()[i_pt]
        central_sigma_ds = h_sigma_ds.values()[i_pt]
        central_sigma_ds_unc = h_sigma_ds.errors()[i_pt]
        central_rawy_dplus = h_rawy_dplus.values()[i_pt]
        central_rawy_dplus_unc = h_rawy_dplus.errors()[i_pt]
        central_sigma_dplus = h_sigma_dplus.values()[i_pt]
        central_sigma_dplus_unc = h_sigma_dplus.errors()[i_pt]

        # Plot the results
        axs[0, 0].errorbar(
            x=range(1, n_trials + 1), y=raw_yields_ds,
            yerr=raw_yields_ds_unc, fmt='o', label=r'Fit ($\mathrm{D_s^+}$)',
            zorder=2
        )
        axs[0, 0].errorbar(
            x=range(1, n_trials + 1), y=raw_yields_dplus,
            yerr=raw_yields_dplus_unc, fmt='o', label=r'Fit ($\mathrm{D^+}$)',
            zorder=2
        )

        for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
            axs[0, 0].errorbar(
                x=list(range(n_trials * (i_nsigma + 1) + 1, n_trials * (i_nsigma + 2) + 1)),
                y=raw_yields_ds_nsigma[i_nsigma],
                yerr=raw_yields_ds_nsigma_unc[i_nsigma], fmt='o',
                label=fr'Bin counting {nsigma}$\sigma$',
                zorder=1
            )
            axs[0, 0].errorbar(
                x=list(range(n_trials * (i_nsigma + 1) + 1, n_trials * (i_nsigma + 2) + 1)),
                y=raw_yields_dplus_nsigma[i_nsigma],
                yerr=raw_yields_dplus_nsigma_unc[i_nsigma], fmt='o',
                label=fr'Bin counting {nsigma}$\sigma$',
                zorder=1
            )

        # Draw the central values
        axs[0, 0].axhline(y=central_rawy_ds, color='C0', linestyle='--')
        axs[0, 0].add_patch(
            Rectangle(
                (0, central_rawy_ds - central_rawy_ds_unc),
                x_axis_range, 2 * central_rawy_ds_unc,
                color='C0', alpha=0.3, zorder=0,
                label=r'Central value $\pm$ uncertainty ($\mathrm{D_s^+}$)'
            )
        )
        axs[0, 0].axhline(y=central_rawy_dplus, color='C1', linestyle='--')
        axs[0, 0].add_patch(
            Rectangle(
                (0, central_rawy_dplus - central_rawy_dplus_unc),
                x_axis_range, 2 * central_rawy_dplus_unc,
                color='C1', alpha=0.3, zorder=0,
                label=r'Central value $\pm$ uncertainty ($\mathrm{D^+}$)'
            )
        )

        axs[0, 0].set_xlim(0, x_axis_range)
        axs[0, 0].set_xlabel('Trial', fontsize=14)
        axs[0, 0].set_ylabel('Raw yield', fontsize=14)
        axs[0, 0].legend(fontsize=12)


        # Draw the ratio distribution
        ratio = np.array(raw_yields_ds) / np.array(raw_yields_dplus)
        ratio_nsigma = [np.array(raw_yields_ds_nsigma[i_nsigma]) / np.array(raw_yields_dplus_nsigma[i_nsigma])  # pylint: disable=line-too-long # noqa: 501
                        for i_nsigma in range(n_bincounts)]
        df_pt_cent["ratio"] = ratio
        axs[0, 1].hist(
            ratio, bins=30, alpha=0.7, label='Fit',
            histtype='stepfilled', ec="k", linewidth=2, zorder=2
        )

        for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
            axs[0, 1].hist(
                ratio_nsigma[i_nsigma],
                bins=30,
                alpha=0.3,
                label=fr'Bin Counting {nsigma}$\sigma$',
                histtype='stepfilled',
                ec="k",
                zorder=1
            )

        # Draw information
        info = 'Fit:\n'
        info += fr'$\mu =$ {np.mean(ratio):.3f}''\n'
        info += fr'$\sigma =$ {np.std(ratio):.3f}''\n'
        for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
            info += fr'{nsigma}$\sigma$ Bin counting:''\n'
            info += fr'$\mu =$ {np.mean(ratio_nsigma[i_nsigma]):.3f}''\n'
            info += fr'$\sigma =$ {np.std(ratio_nsigma[i_nsigma]):.3f}''\n'
        anchored_text_fit = AnchoredText(
            info,
            loc='upper left',
            frameon=False
        )

        # Draw the rms + shift from the central value
        central_ratio = central_rawy_ds / central_rawy_dplus
        axs[0, 1].axvline(x=central_ratio, color='r', linestyle='--')
        rms_shift = get_rms_shift_sum_quadrature(ratio, central_ratio, False)
        axs[0, 1].add_patch(
            Rectangle(
                (central_ratio - rms_shift, 0),
                2 * rms_shift, axs[0, 1].get_ylim()[1],
                color='r', alpha=0.3, zorder=0,
                label=r'$\mathrm{\sqrt{RMS^2 + \Delta^2}}$'f' = {rms_shift/central_ratio*100:.2f}%'
            )
        )
        axs[0, 1].add_artist(anchored_text_fit)

        # Draw the assigned systematic uncertainty
        axs[0, 1].set_xlabel('Ratio', fontsize=14)
        axs[0, 1].set_ylabel('Counts', fontsize=14)
        axs[0, 1].legend(fontsize=12, loc='upper right')

        x_min = min(
            ratio.min(),
            #*[df_pt_cent[f"rawy_bincounting_{nsigma}"].min()
            #    for nsigma in multitrial_cfg['bincounting_nsigma']],
            central_ratio - rms_shift
        )
        x_max = max(
            ratio.max(),
            #*[df_pt_cent[f"rawy_bincounting_{nsigma}"].max()
            #    for nsigma in multitrial_cfg['bincounting_nsigma']],
            central_ratio + rms_shift
        )
        axs[0, 1].set_xlim(x_min * 0.9, x_max * 1.1)

        # Draw the peak width
        axs[1, 0].errorbar(
            x=range(1, len(sigma_ds) + 1),
            y=sigma_ds,
            yerr=sigma_ds_unc, fmt='o',
            zorder=2
        )
        axs[1, 0].errorbar(
            x=range(1, len(sigma_dplus) + 1),
            y=sigma_dplus,
            yerr=sigma_dplus_unc, fmt='o',
            zorder=2
        )

        # Draw the central values
        axs[1, 0].axhline(y=central_sigma_ds, color='C0', linestyle='--')
        axs[1, 0].add_patch(
            Rectangle(
                (0, central_sigma_ds - central_sigma_ds_unc),
                x_axis_range, 2 * central_sigma_ds_unc,
                color='C0', alpha=0.3, zorder=0,
                label=r'Central value $\pm$ uncertainty ($\mathrm{D_s^+}$)'
            )
        )
        axs[1, 0].axhline(y=central_sigma_dplus, color='C1', linestyle='--')
        axs[1, 0].add_patch(
            Rectangle(
                (0, central_sigma_dplus - central_sigma_dplus_unc),
                x_axis_range, 2 * central_sigma_dplus_unc,
                color='C1', alpha=0.3, zorder=0,
                label=r'Central value $\pm$ uncertainty ($\mathrm{D^+}$)'
            )
        )

        axs[1, 0].set_xlim(0, x_axis_range)
        axs[1, 0].set_xlabel('Trial', fontsize=14)
        axs[1, 0].set_ylabel('Width ($GeV/c^2$)', fontsize=14)
        axs[1, 0].legend(fontsize=12)

        # Draw the chi2/ndf
        axs[1, 1].scatter(
            x=range(1, len(df_pt_cent["chi2"]) + 1),
            y=df_pt_cent["chi2"]
        )
        axs[1, 1].set_xlim(0, x_axis_range)
        axs[1, 1].set_xlabel('Trial', fontsize=14)
        axs[1, 1].set_ylabel(r'$\chi^2/$ndf', fontsize=14)

        plt.show()
        outfolder = os.path.join(cfg["output"]["dir"], "output", "syst")
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        fig.savefig(
            os.path.join(
                outfolder,
                f'fig_{pt_min*10:.0f}_{pt_max*10:.0f}_{cent_min}_{cent_max}.png'
            ), bbox_inches='tight'
        )

        bkg_funcs, sigmas = [], []
        for _, row in df_pt_cent.iterrows():
            bkg_funcs.append(row["bkg_pdfs_cfg"][0])
        df_pt_cent["bkg_pdfs_cfg"] = bkg_funcs
        variations = ["mins", "maxs", "rebins", "bkg_pdfs_cfg", "ratio_sigma_dplus_to_ds_cfg"]
        if not is_mb:
            for _, row in df_pt_cent.iterrows():
                sigmas.append(row["sigma_option"])
            df_pt_cent["sigma_option"] = sigmas
            variations.append("sigma_option")

        combinations = set(itertools.combinations(variations, 2))
        n_rows = 0
        if len(combinations) % 4 == 0:
            n_rows = len(combinations) // 4
        else:
            n_rows = len(combinations) // 4 + 1
        fig, axs = plt.subplots(n_rows, 4, figsize=(20, 5 * len(combinations)/4))
        for i_comb, combination in enumerate(combinations):
            sns.stripplot(
                data=df_pt_cent, x=combination[0], y=ratio, hue=combination[1],
                dodge=0.5, alpha=.5, legend=True, ax=axs[i_comb//4, i_comb%4], palette="tab10",
            )
            # sns.pointplot(
            #     data=df_pt_cent, x=combination[0], y=ratio, hue=combination[1],
            #     dodge=0.5, linestyle="none", errorbar=None, palette="tab10",
            #     marker="_", markersize=20, markeredgewidth=3, ax=axs[i_comb//4, i_comb%4]
            # )
        outfolder = os.path.join(cfg["output"]["dir"], "output", "ratio")
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        plt.savefig(
            os.path.join(
                outfolder,
                f'fig_ratio_{pt_min*10:.0f}_{pt_max*10:.0f}_{cent_min}_{cent_max}.png'
            ),
            bbox_inches='tight'
        )


def dump_results_to_root(df_multitrial, cfg):  # pylint: disable=too-many-locals
    """
    Dump the results to a ROOT file.

    Parameters:
        dfs (list of pandas.DataFrame): List of dataframes containing the data for each pt bin.
        cfg (dict): Configuration dictionary.
        cut_set (dict): Dictionary containing the cut sets.

    Returns:
        None
    """
    multitrial_cfg = cfg["multitrial"]

    df_multitrial = df_multitrial.query(f"chi2 < {multitrial_cfg['quality_selections']['chi2']}")

    # Apply significance selection
    mask = df_multitrial['significance'].apply(
        lambda x: x[0][0] > multitrial_cfg['quality_selections']['significance'] and \
                 x[1][0] > multitrial_cfg['quality_selections']['significance']
        )
    df_multitrial = df_multitrial[mask]

    cent_mins = list(df_multitrial["cent_min_cfg"].unique())
    cent_mins.sort()
    cent_maxs = list(df_multitrial["cent_max_cfg"].unique())
    cent_maxs.sort()

    pt_mins = list(df_multitrial["pt_min_cfg"].unique())
    pt_mins.sort()
    pt_maxs = list(df_multitrial["pt_max_cfg"].unique())
    pt_maxs.sort()

    pt_edges = np.array(pt_mins + [pt_maxs[-1]])

    histos_rms_shifts = {}
    histos_assigned_syst = {}

    df_cent_cfg = df_multitrial[
        ["cent_min_cfg", "cent_max_cfg"]
    ].drop_duplicates()
    for _, row in df_cent_cfg.iterrows():
        cent_min = row["cent_min_cfg"]
        cent_max = row["cent_max_cfg"]
        cent_str = f"{cent_min}_{cent_max}"

        df_cent = df_multitrial.query(
            f"cent_min_cfg == {cent_min} and cent_max_cfg == {cent_max}"
        )

        df_pt_cent_cfg = df_cent[
            ["pt_min_cfg", "pt_max_cfg"]
        ].drop_duplicates()

        histos_rms_shifts[cent_str] = [0.]*len(df_pt_cent_cfg)
        histos_assigned_syst[cent_str] = [0.]*len(df_pt_cent_cfg)
        for _, row in df_pt_cent_cfg.iterrows():
            pt_min = row["pt_min_cfg"]
            pt_max = row["pt_max_cfg"]
            pt_idx = pt_mins.index(pt_min)

            df_pt_cent = df_cent.query(
                f"pt_min_cfg == {pt_min} and pt_max_cfg == {pt_max}"
            )

            raw_yields_ds, raw_yields_dplus, = [], []
            for _, row in df_pt_cent.iterrows():
                raw_yields_ds.append(row["raw_yields"][0][0])
                raw_yields_dplus.append(row["raw_yields"][1][0])

            with uproot.open(cfg["reference_fits"]) as f:
                h_rawy_ds = f[f"h_raw_yields_ds_{cent_min}_{cent_max}"]
                h_rawy_dplus = f[f"h_raw_yields_dplus_{cent_min}_{cent_max}"]

            i_pt = np.digitize((pt_min + pt_max) / 2, h_rawy_ds.axis().edges()) - 1
            central_rawy_ds = h_rawy_ds.values()[i_pt]
            central_rawy_dplus = h_rawy_dplus.values()[i_pt]

            ratio = np.array(raw_yields_ds) / np.array(raw_yields_dplus)
            central_ratio = central_rawy_ds / central_rawy_dplus


            histos_rms_shifts[cent_str][pt_idx] = get_rms_shift_sum_quadrature(
                ratio,
                central_ratio,
                True
            )
            assigned_syst = cfg["assigned_syst"][pt_idx][cent_mins.index(cent_min)]
            histos_assigned_syst[cent_str][pt_idx] = assigned_syst


    if not os.path.exists(cfg["output"]["dir"]):
        os.makedirs(cfg["output"]["dir"])
    output_file_name = os.path.join(cfg["output"]["dir"], "raw_yields_systematic.root")
    with uproot.recreate(output_file_name) as f:
        for cent_min, cent_max in zip(cent_mins, cent_maxs):
            suffix = f"_{cent_min}_{cent_max}"
            f[f"rms_shifts_sum_quadrature{suffix}"] = (
                np.array(histos_rms_shifts[f"{cent_min}_{cent_max}"]),
                np.array(pt_edges)
            )
            f[f"assigned_syst{suffix}"] = (
                np.array(histos_assigned_syst[f"{cent_min}_{cent_max}"]),
                np.array(pt_edges)
            )


def get_rms_shift_sum_quadrature(ratio, central_ratio, rel=False):
    """
    Calculate the sum in quadrature of the RMS and shift from the central value for raw yields.

    Parameters:
        ratio (pandas.DataFrame): DataFrame containing the ratios.
        central_ratio (dict): value of the central ratio.
        rel (bool): If True, return the relative uncertainty.

    Returns:
        float: The sum in quadrature of the RMS and shift from the central value for raw yields.
    """
    if rel:
        return np.sqrt(
            np.std(ratio)**2 +
            (np.mean(ratio) - central_ratio)**2
        ) / central_ratio

    return np.sqrt(
        np.std(ratio)**2 +
        (np.mean(ratio) - central_ratio)**2
    )

def fitconfig_to_dict(cfg: FitConfig) -> dict:
    """Convert FitConfig dataclass to dictionary with modified keys."""
    base = dataclasses.asdict(cfg)
    return {f"{k}_cfg": v for k, v in base.items()}

def get_parameter_setup(
        param_cfg: List[Dict],
        trial: Dict,
    ) -> List[Dict]:
    """
    Get parameter setup based on the provided configuration.

    Parameters:
    cfg (list): List of configuration dictionaries.
    i_pt (int): Index of the pT bin.

    Returns:
    List[Dict]: Parameter setup dictionary.
    """
    i_pt = trial["pt_bins"]
    functions_setup = []
    for func_setup in param_cfg:  # For each signal function
        param_setup = {}
        for par_name, par_values in func_setup.items():  # For each parameter in the function
            param_setup[par_name] = {
                "init": par_values["init"][i_pt],
                "min": par_values["min"][i_pt],
                "max": par_values["max"][i_pt],
                "fix_to_config_value": par_values["fix_to_config_value"][i_pt],
                "fix_to_file": par_values["fix_to_file"][i_pt]
            }
        functions_setup.append(param_setup)
    return functions_setup

def get_corr_bkg_config(cfg: Dict) -> CorrelatedBackgroundConfig:
    """
    Create a CorrelatedBackgroundConfig object based on the provided configuration.

    Parameters:
    cfg (dict): Configuration dictionary.
    pt_min (float): Minimum pT value.
    pt_max (float): Maximum pT value.

    Returns:
    CorrelatedBackgroundConfig: Configured CorrelatedBackgroundConfig object.
    """
    bkg_norm = cfg["multitrial"]["backgrounds"]
    sgn_norm = cfg["multitrial"]["signal"]

    return CorrelatedBackgroundConfig(
        fix_to_file=False,
        fix_to_mb=False,
        fix_with_br=True,
        file_name_for_fix=None,
        hist_name_for_fix=None,
        backgrounds=[
            CorrelatedBackground(
                name=bkg["name"],
                file_norm=bkg["file_norm"],
                norm_hist_name=bkg["norm_hist_name"],
                template_file=bkg["template_file"],
                template_hist_name=bkg["template_hist_name"],
                br=BRInfo(pdg=bkg["br"]["pdg"], simulations=bkg["br"]["simulations"])
            )
            for bkg in bkg_norm
        ],
        signal_norm_file=sgn_norm["file_norm"],
        signal_hist_name=sgn_norm["hist_name"],
        signal_br=BRInfo(
            pdg=sgn_norm["br"]["pdg"],
            simulations=sgn_norm["br"]["simulations"]
        )
    )

def get_sigma_ds(result_sigma_fix: Dict, cfg_trial: Dict) -> float:  # pylint: disable=too-many-return-statements
    """
    Get the sigma value for Ds from the reference result depending on the configuration.

    Parameters:
    result_sigma_fix (Dict): Reference result dictionary.
    cfg_trial (Dict): Trial configuration dictionary.

    Returns:
    float: Sigma value for Ds.
    """
    if cfg_trial["sigma_option"] == "free":
        return None
    if cfg_trial["sigma_option"] == "fixed":
        return result_sigma_fix["sigma"][0][0]
    if cfg_trial["sigma_option"] == "fixed_plus_unc":
        return result_sigma_fix["sigma"][0][0] + result_sigma_fix["sigma"][0][1]
    if cfg_trial["sigma_option"] == "fixed_minus_unc":
        return result_sigma_fix["sigma"][0][0] - result_sigma_fix["sigma"][0][1]
    if "fixed_plus_" in cfg_trial["sigma_option"] and "_perc" in cfg_trial["sigma_option"]:
        perc = float(cfg_trial["sigma_option"].split("_")[2]) / 100.0
        return result_sigma_fix["sigma"][0][0] * (1 + perc)
    if "fixed_minus_" in cfg_trial["sigma_option"] and "_perc" in cfg_trial["sigma_option"]:
        perc = float(cfg_trial["sigma_option"].split("_")[2]) / 100.0
        return result_sigma_fix["sigma"][0][0] * (1 - perc)
    return None

def get_config(
        cfg: Dict,
        cfg_trial: Dict,
        pt_info: Tuple[int, float, float],
        ref_result: Dict = None
) -> FitConfig:
    """
    Create a FitConfig object based on the provided configuration.

    Parameters:
    cfg (dict): Configuration dictionary.
    cfg_trial (dict): Trial configuration dictionary.
    pt_min (float): Minimum pT value.
    pt_max (float): Maximum pT value.
    cent_min (float): Minimum centrality value.
    cent_max (float): Maximum centrality value.

    Returns:
    FitConfig: Configured FitConfig object.
    """
    i_pt, pt_mins, pt_maxs = pt_info
    cent_min, cent_max = cfg_trial["cent_min_cfg"], cfg_trial["cent_max_cfg"]
    par_setup = cfg_trial["par_init_limit"]

    sigma = None
    if ref_result is not None:
        if not (cent_min == 0 and cent_max == 100):
            # Fix ds sigma to the MB result if needed
            sigma = get_sigma_ds(ref_result, cfg_trial)
            if sigma is not None:
                # Sigma not free, set it for Ds
                par_setup[0]["sigma"] = {
                    "init": sigma,
                    "min": sigma - 0.01,
                    "max": sigma + 0.01,
                    "fix_to_config_value": True,
                    "fix_to_file": False
                }
                # Set also ratio for D+
                ratio_sigma_dplus_to_ds = cfg_trial["fix_ratio_ds_dplus_width"]
                par_setup[1]["sigma"] = {
                    "init": sigma * ratio_sigma_dplus_to_ds,
                    "min": sigma * ratio_sigma_dplus_to_ds - 0.01,
                    "max": sigma * ratio_sigma_dplus_to_ds + 0.01,
                    "fix_to_config_value": True,
                    "fix_to_file": False
                }
            else:
                ratio_sigma_dplus_to_ds = cfg_trial["fix_ratio_ds_dplus_width"]
                par_setup[1]["sigma"] = {
                    "init": ref_result["sigma"][0][0] * ratio_sigma_dplus_to_ds,
                    "min": ref_result["sigma"][0][0] * ratio_sigma_dplus_to_ds - 0.01,
                    "max": ref_result["sigma"][0][0] * ratio_sigma_dplus_to_ds + 0.01,
                    "fix_to_config_value": True,
                    "fix_to_file": False
                }
        else:
            # For MB we only want to fix D+ sigma to Ds if needed
            ratio_sigma_dplus_to_ds = cfg_trial["fix_ratio_ds_dplus_width"]
            par_setup[1]["sigma"] = {
                "init": ref_result["sigma"][0][0] * ratio_sigma_dplus_to_ds,
                "min": ref_result["sigma"][0][0] * ratio_sigma_dplus_to_ds - 0.01,
                "max": ref_result["sigma"][0][0] * ratio_sigma_dplus_to_ds + 0.01,
                "fix_to_config_value": True,
                "fix_to_file": False
            }

    return FitConfig(
        pt_min=pt_mins[i_pt],
        pt_max=pt_maxs[i_pt],
        cent_min=cent_min,
        cent_max=cent_max,
        mass_range=[
            cfg_trial["mins"],
            cfg_trial["maxs"]
        ],
        signal_pdfs=cfg_trial["sgn_funcs"],
        bkg_pdfs=cfg_trial["bkg_funcs"],
        rebin=cfg_trial["rebins"],
        param_setup=par_setup,
        data_path=cfg["inputs"]["data"],
        file_for_params_fix=cfg["multitrial"]["file_for_params_fix"],
        suffix_hist_for_params_fix=cfg["multitrial"]["suffix_hist_for_params_fix"],
        fix_dplus_sigma_to_ds=True,
        ratio_sigma_dplus_to_ds=cfg_trial["fix_ratio_ds_dplus_width"],
        correlated_bkg=get_corr_bkg_config(cfg),
        draw_figures=cfg["output"]["save_all_fits"],
        draw_formats=cfg["output"]["formats"],
        output_dir=os.path.join(
            os.path.expanduser(cfg["output"]["dir"]), cfg["output"]["dir_fits"]
        ),
        fig_suffix=str(cfg_trial["idx"]),
        nsigma_bincounting=cfg["multitrial"]["bincounting_nsigma"]
    )


def get_trials(  # pylint: disable=too-many-locals
        cfg: Dict,
        is_mb: bool,
        cent_mins: List[float] = None,
        cent_maxs: List[float] = None
    ) -> List[Dict]:
    """Generate all trial configurations for MB analysis."""
    base_cfg_keys = ["pt_bins", "mins", "maxs", "rebins", "bkg_funcs", "fix_ratio_ds_dplus_width"]
    if not is_mb:
        base_cfg_keys.append("sigma_option")
    trial_params = {key: cfg[key] for key in base_cfg_keys}

    # First we create the base trials without the signal setup
    keys, values = zip(*trial_params.items())
    base_trials = [dict(zip(keys, v)) for v in itertools.product(*values)]

    # We add the signal setup
    trials = []

    sgn_setups = cfg["sgn_funcs"]
    par_setups = cfg["par_init_limit"]

    for t in base_trials:
        for sgn_funcs, par_block in zip(sgn_setups, par_setups):
            if not is_mb:
                for cent_min, cent_max in zip(cent_mins, cent_maxs):
                    trial = t.copy()
                    trial["sgn_funcs"] = sgn_funcs
                    trial["par_init_limit"] = get_parameter_setup(par_block, trial)
                    trial["cent_min_cfg"] = cent_min
                    trial["cent_max_cfg"] = cent_max
                    trials.append(trial)
            else:
                trial = t.copy()
                trial["sgn_funcs"] = sgn_funcs
                trial["par_init_limit"] = get_parameter_setup(par_block, trial)
                trial["cent_min_cfg"] = 0
                trial["cent_max_cfg"] = 100
                trial["sigma_option"] = None
                trials.append(trial)


    # Add index to each trial for pT bin
    counters = [0] * len(cfg["pt_bins"])
    for trial in trials:
        trial["idx"] = counters[cfg["pt_bins"].index(trial["pt_bins"])]
        counters[cfg["pt_bins"].index(trial["pt_bins"])] += 1

    return trials

def get_matching_mb_result(mb_results: List[Dict], cent_trial: Dict) -> Dict:
    """
    Find the matching MB result for a given centrality trial.

    Parameters:
    mb_results (List[Dict]): List of MB results.
    cent_trial (Dict): Centrality trial configuration.

    Returns:
    Dict: Matching MB result.
    """
    base_cfg_keys = ["pt_bins", "mins", "maxs", "rebins", "bkg_funcs", "fix_ratio_ds_dplus_width"]
    for mb_result in mb_results:
        for key in base_cfg_keys:
            if mb_result["trial"][key] != cent_trial[key]:
                break
        else:
            return mb_result["results"][1]

    return None

def run_fit(cfg, fit_config: FitConfig) -> Dict:
    """
    Run the fitting procedure based on the provided FitConfig.

    Parameters:
    fit_config (FitConfig): Configuration for the fit.

    Returns:
    Dict: Results of the fitting procedure.
    """
    with open(os.path.join(cfg["output"]["dir"], "fit_log.txt"), "a") as f, contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
        fit_handler = FitHandler(fit_config)
        results = fit_handler.get_results()
        return (
            fit_config,
            results
        )

def multi_trial(config_file_name: str, draw=False):  # pylint: disable=too-many-locals, too-many-statements, too-many-branches
    """
    Perform multiple trials based on the given configuration file.

    Parameters:
        config_file_name (str): The path to the configuration file.
    Returns:
        None
    """
    # load config
    with open(config_file_name, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    # load cut set
    with open(cfg["inputs"]["cutset"], "r", encoding="utf-8") as f:
        cut_set = yaml.safe_load(f)

    pt_mins = cut_set["pt"]["min"]
    pt_maxs = cut_set["pt"]["max"]
    cent_mins = cut_set["cent"]["min"]
    cent_maxs = cut_set["cent"]["max"]

    index = [*zip(cent_mins, cent_maxs)].index((0, 100))
    cent_mins.pop(index)
    cent_maxs.pop(index)

    if not os.path.exists(cfg["output"]["dir"]):
        os.makedirs(cfg["output"]["dir"])

    if cfg["output"]["save_all_fits"]:
        fit_dir = os.path.join(cfg["output"]["dir"], cfg["output"]["dir_fits"])
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)

    log_file = os.path.join(cfg["output"]["dir"], "fit_log.txt")
    if os.path.exists(log_file):
        os.remove(log_file)

    multitrial_cfg = cfg["multitrial"]
    if not draw:
        mb_results = []
        cent_results = []
        # define all MB trials
        mb_trials = get_trials(multitrial_cfg, is_mb=True)
        cent_trials = get_trials(
            multitrial_cfg,
            is_mb=False,
            cent_mins=cent_mins,
            cent_maxs=cent_maxs
        )
        for ipt in cfg["multitrial"]["pt_bins"]:
            print(f"Starting MB fits for pT bin {ipt}...")
            with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
                future_to_trial_first_stage = {}  # Map futures to trials
                future_to_trial_second_stage = {}  # Map futures to trials
                for trial in mb_trials:
                    if trial["pt_bins"] != ipt:
                        continue
                    fit_config = get_config(cfg, trial, (trial["pt_bins"], pt_mins, pt_maxs))
                    future = executor.submit(run_fit, cfg, fit_config)
                    future_to_trial_first_stage[future] = trial

                for future in tqdm(as_completed(future_to_trial_first_stage), total=len(future_to_trial_first_stage),
                        desc="1st stage MB fits"):
                    trial = future_to_trial_first_stage[future]
                    try:
                        fit_cfg, results = future.result()
                    except Exception as e:
                        print(f"Error in trial {trial}: {e}")
                        continue
                    fit_config = get_config(
                        cfg,
                        trial,
                        (trial["pt_bins"], pt_mins, pt_maxs),
                        results
                    )

                    feature2 = executor.submit(run_fit, cfg, fit_config)
                    future_to_trial_second_stage[feature2] = trial

                print("Finalising MB fits...")
                # Collect second-stage MB fits
                for future in tqdm(as_completed(future_to_trial_second_stage), total=len(future_to_trial_second_stage),
                        desc="2nd stage MB fits"):
                    trial = future_to_trial_second_stage[future]
                    try:
                        fit_cfg, results = future.result()
                    except Exception as e:
                        print(f"Error in trial {trial}: {e}")
                        continue
                    mb_results.append({
                        "trial": trial,
                        "results": (fit_cfg, results)
                    })

                print("MB fits completed.")

            # Now we fit the centrality bins
            # Spawn a new process for each centrality bin
            # This avoids issues with multiprocessing getting stuck
            print(f"Starting centrality fits for pT bin {ipt}...")
            mp_ctx = multiprocessing.get_context('spawn')
            with ProcessPoolExecutor(max_workers=cfg["max_workers"], mp_context=mp_ctx) as executor:
                future_to_trial = {}
                for trial in cent_trials:
                    if trial["pt_bins"] != ipt:
                        continue
                    fit_config = get_config(
                        cfg,
                        trial,
                        (trial["pt_bins"], pt_mins, pt_maxs),
                        ref_result=get_matching_mb_result(mb_results, trial)
                    )
                    future = executor.submit(run_fit, cfg, fit_config)
                    future_to_trial[future] = trial

                for future in tqdm(as_completed(future_to_trial), total=len(future_to_trial),
                        desc="1st stage centrality fits"):
                    trial = future_to_trial[future]
                    try:
                        fit_cfg, results = future.result()
                    except Exception as e:
                        print(f"Error in trial {trial}: {e}")
                        continue
                    if trial["sigma_option"] == "free":
                        fit_config = get_config(
                            cfg,
                            trial,
                            (trial["pt_bins"], pt_mins, pt_maxs),
                            ref_result=results
                        )
                        feature2 = executor.submit(run_fit, cfg, fit_config)
                        future_to_trial_second_stage[feature2] = trial

                    else:
                        cent_results.append({
                            "trial": trial,
                            "results": (fit_cfg, results)
                        })
                
                print("Finalising centrality fits...")
                # Collect second-stage centrality fits
                for future in tqdm(as_completed(future_to_trial_second_stage), total=len(future_to_trial_second_stage),
                        desc="2nd stage centrality fits"):
                    trial = future_to_trial_second_stage[future]
                    try:
                        fit_cfg, results = future.result()
                    except Exception as e:
                        print(f"Error in trial {trial}: {e}")
                        continue
                    cent_results.append({
                        "trial": trial,
                        "results": (fit_cfg, results)
                    })
                print("Centrality fits completed.")

        # Save results to parquet
        mb_df = pd.DataFrame([{
            **trial["trial"],
            **fitconfig_to_dict(trial["results"][0]),
            **trial["results"][1]
        } for trial in mb_results])

        cent_df = pd.DataFrame([{
            **trial["trial"],
            **fitconfig_to_dict(trial["results"][0]),
            **trial["results"][1]
        } for trial in cent_results])

        mb_df.to_parquet(os.path.join(
            cfg["output"]["dir"],
            f"mb_results{cfg['output']['suffix']}.parquet"
        ))
        cent_df.to_parquet(os.path.join(
            cfg["output"]["dir"],
            f"cent_results{cfg['output']['suffix']}.parquet"
        ))


    else:
        if cfg["files_to_draw"]["mb"] is not None:
            mb_df, cent_df = [], []
            for file in cfg["files_to_draw"]["mb"]:
                mb_df.append(pd.read_parquet(file))
            mb_df = pd.concat(mb_df)
            for file in cfg["files_to_draw"]["cent"]:
                cent_df.append(pd.read_parquet(file))
            cent_df = pd.concat(cent_df)
        else:
            mb_df = pd.read_parquet(os.path.join(
                cfg["output"]["dir"],
                f"mb_results{cfg['output']['suffix']}.parquet"
            ))
            cent_df = pd.read_parquet(os.path.join(
                cfg["output"]["dir"],
                f"cent_results{cfg['output']['suffix']}.parquet"
            ))

    draw_multitrial(cent_df, cfg, is_mb=False)
    draw_multitrial(mb_df, cfg, is_mb=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fit Ds')
    parser.add_argument('configFile', type=str, help='Path to the configuration file')
    parser.add_argument('--draw', action='store_true', help='Only draw the results')
    args = parser.parse_args()

    multi_trial(args.configFile, args.draw)
