"""
This script performs a multitrial procedure to get the raw yields systematic uncertainty.

Usage:
    python get_ry_systematic.py configFile
"""
# pylint: disable=too-many-lines
import re
import argparse
import itertools
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import yaml
import uproot
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
from flarefly.utils import Logger
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Rectangle
import seaborn as sns
import pandas as pd
#zfit.run.set_n_cpu(1)


def draw_multitrial(df_multitrial, cfg): #, pt_min, pt_max, idx_assigned_syst):  # pylint: disable=too-many-locals, too-many-statements # noqa: 501
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
                yerr=raw_yields_ds_nsigma_unc, fmt='o',
                label=fr'Bin counting {nsigma}$\sigma$',
                zorder=1
            )
            axs[0, 0].errorbar(
                x=list(range(n_trials * (i_nsigma + 1) + 1, n_trials * (i_nsigma + 2) + 1)),
                y=raw_yields_dplus_nsigma[i_nsigma],
                yerr=raw_yields_dplus_nsigma_unc, fmt='o',
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
        rms_shift = get_rms_shift_sum_quadrature(ratio, central_ratio,i_pt)
        axs[0, 1].add_patch(
            Rectangle(
                (central_ratio - rms_shift, 0),
                2 * rms_shift, axs[0, 1].get_ylim()[1],
                color='r', alpha=0.3, zorder=0,
                label=r'$\mathrm{\sqrt{RMS^2 + \Delta^2}}$'
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
        fig.savefig(
            os.path.join(
                cfg["output"]["dir"], "output", "syst",
                f'fig_{pt_min*10:.0f}_{pt_max*10:.0f}_{cent_min}_{cent_max}.png'
            ), bbox_inches='tight'
        )

        bkg_funcs = []
        for _, row in df_pt_cent.iterrows():
            bkg_funcs.append(row["bkg_funcs_cfg"][0])
        df_pt_cent["bkg_funcs_cfg"] = bkg_funcs
        variations = ["mins_cfg", "maxs_cfg", "rebins_cfg", "bkg_funcs_cfg"]
        combinations = set(itertools.combinations(variations, 2))
        n_rows = 0
        if len(combinations) % 2 == 0:
            n_rows = len(combinations) // 2
        else:
            n_rows = len(combinations) // 2 + 1
        fig, axs = plt.subplots(n_rows, 2, figsize=(20, 5 * len(combinations)/2))
        for i_comb, combination in enumerate(combinations):
            sns.stripplot(
                data=df_pt_cent, x=combination[0], y=ratio, hue=combination[1],
                dodge=0.5, alpha=.5, legend=False, ax=axs[i_comb//2, i_comb%2], palette="tab10",
            )
            sns.pointplot(
                data=df_pt_cent, x=combination[0], y=ratio, hue=combination[1],
                dodge=0.5, linestyle="none", errorbar=None, palette="tab10",
                marker="_", markersize=20, markeredgewidth=3, ax=axs[i_comb//2, i_comb%2]
            )
        plt.savefig(
            os.path.join(
                cfg["output"]["dir"], "output", "ratio",
                f'fig_ratio_{pt_min*10:.0f}_{pt_max*10:.0f}_{cent_min}_{cent_max}.png'
            ),
            bbox_inches='tight'
        )

def get_input_data(cfg, pt_mins, pt_maxs, bdt_cut_mins, bdt_cut_maxs):  # pylint: disable=too-many-locals # noqa: 501
    """
    Retrieve input data based on the given configuration and selection criteria.

    Parameters:
    - cfg (dict): Configuration dictionary.
    - pt_mins (list): List of minimum pt values.
    - pt_maxs (list): List of maximum pt values.
    - bdt_cut_mins (list): List of minimum BDT cut values.
    - bdt_cut_maxs (list): List of maximum BDT cut values.

    Returns:
    - tuple: A tuple containing:
        - df (pandas.DataFrame): DataFrame containing the selected data.
        - df_mc_sig (pandas.DataFrame): DataFrame containing the selected signal MC data.
        - df_mc_prd_bkg (pandas.DataFrame): DataFrame containing the selected
            partially reco decays background MC data.
    """
    selection_string = ""
    for ipt, (pt_min, pt_max, bdt_cut_min, bdt_cut_max) in enumerate(zip(pt_mins, pt_maxs, bdt_cut_mins, bdt_cut_maxs)):  # pylint: disable=line-too-long # noqa: 501
        if ipt == 0:
            selection_string += f"({pt_min} < fPt < {pt_max} and\
                {bdt_cut_min} < ML_output < {bdt_cut_max})"
        else:
            selection_string += f" or ({pt_min} < fPt < {pt_max} and\
                {bdt_cut_min} < ML_output < {bdt_cut_max})"

    # load data
    df = pd.DataFrame()
    for file in cfg["inputs"]["data"]:
        df = pd.concat([df, pd.read_parquet(file)])
    df.query(selection_string, inplace=True)

    # load mc
    df_mc = pd.DataFrame()
    for file in cfg["inputs"]["mc"]:
        df_mc = pd.concat([df_mc, pd.read_parquet(file)])
    df_mc.query(selection_string, inplace=True)

    df_mc_sig = df_mc.query("fFlagMcMatchRec == -1 or fFlagMcMatchRec == 1")
    if any(cfg["multitrial"]["use_bkg_templ"]):
        df_mc_prd_bkg = df_mc.query("fFlagMcMatchRec == 8")  # prd = partly reco decays
    else:
        df_mc_prd_bkg = None

    return df, df_mc_sig, df_mc_prd_bkg


def get_fixed_parameter(par, unc, string):
    """
    Return the value of a fixed parameter based on the given string.

    Parameters:
    - par (float): The value of the parameter.
    - unc (float): The uncertainty of the parameter.
    - string (str): The string indicating the type of fixed parameter.

    Returns:
    - float: The value of the fixed parameter.

    Raises:
    - FATAL: If the string does not match any of the patterns.

    """
    if string == "fixed":
        return par
    if string == "fixed_plus_unc":
        return par + unc
    if string == "fixed_minus_unc":
        return par - unc

    # Check if the string matches the pattern fixed_plus_XX_perc
    pattern = r'^fixed_plus_\d+_perc$'
    if re.match(pattern, string):
        # Extract the number
        number = int(string.split('_')[-2])
        return par + par * number / 100

    # Check if the string matches the pattern fixed_minus_XX_perc
    pattern = r'^fixed_minus_\d+_perc$'
    if re.match(pattern, string):
        # Extract the number
        number = int(string.split('_')[-2])
        return par - par * number / 100

    Logger(
        f"The string {string} does not match any of the available patterns", "FATAL"
    )
    return None


def get_mean_sigma_for_data(mean_with_unc, sigma_with_unc, trial):
    """
    Calculate the mean and sigma values for data based on the option chosen for the trial.

    Parameters:
    - mean_with_unc (tuple): A tuple containing the mean value and its uncertainty.
    - sigma_with_unc (tuple): A tuple containing the sigma value and its uncertainty.
    - trial (dict): A dictionary containing the trial parameters (mins, maxs, sgn_funcs,
        bkg_funcs, sigma, mean, use_bkg_templ, bincounting_nsigma).

    Returns:
    - mean_for_data (float or None): The mean value for data,
        or None if trial["mean"] == "free".
    - sigma_for_data (float or None): The sigma value for data,
        or None if trial["sigma"] == "free".
    """
    if trial["sigma"] == "free":
        sigma_for_data = None
    else:
        sigma_for_data = get_fixed_parameter(sigma_with_unc[0], sigma_with_unc[1], trial["sigma"])
    if trial["mean"] == "free":
        mean_for_data = None
    else:
        mean_for_data = get_fixed_parameter(mean_with_unc[0], mean_with_unc[1], trial["mean"])
    return mean_for_data, sigma_for_data


def get_mc_sigma_width_with_unc(trial, data_hdl_mc_sig, fitter_suffix):
    """
    Calculate the mean and sigma with uncertainties for the Monte Carlo (MC) signal.

    Parameters:
    - trial (dict): A dictionary containing the trial parameters (mins, maxs, sgn_funcs,
        bkg_funcs, sigma, mean, use_bkg_templ, bincounting_nsigma).
    - data_hdl_mc_sig (flarefly.DataHandler): The data handler for the MC signal.
    - fitter_suffix (str): A suffix to be added to the fitter name.

    Returns:
    - mean_with_unc: The mean of the MC signal with its uncertainty.
    - sigma_with_unc: The sigma of the MC signal with its uncertainty.
    - None: If the fit did not converge.
    """
    fitter_mc_pt = F2MassFitter(
        data_hdl_mc_sig,
        trial["sgn_funcs"],
        ["nobkg"],
        name=f"mc_{fitter_suffix}",
        label_signal_pdf=[r"$\mathrm{B}^{0}$ signal"],
        verbosity=0
    )
    fitter_mc_pt.set_particle_mass(0, pdg_id=511, fix=False)
    fitter_mc_pt.set_signal_initpar(0, "sigma", 0.01, fix=False, limits=[0.001, 0.1])

    result = fitter_mc_pt.mass_zfit()
    if result.converged:
        mean_with_unc = fitter_mc_pt.get_mass(0)
        sigma_with_unc = fitter_mc_pt.get_sigma(0)

        return mean_with_unc, sigma_with_unc

    # If the fit did not converge, return None
    Logger("The fit for the MC signal did not converge, skipping this pt interval", "WARNING")
    return None, None


def build_fitter(
        trial, data_hdl, data_hdl_prd_bkg,
        mean_with_unc, sigma_with_unc, fitter_suffix
    ):  # pylint: disable=too-many-arguments, too-many-positional-arguments  # noqa: 121, 125
    """
    Build a flarefly mass fitter for fitting the data candidate distribution.

    Parameters:
    - trial (dict): A dictionary containing the trial parameters (mins, maxs, sgn_funcs, bkg_funcs,
        bkg_funcs, sigma, mean, use_bkg_templ, bincounting_nsigma).
    - data_hdl (flarefly.DataHandler): The data handler for the data.
    - data_hdl_prd_bkg (flarefly.DataHandler): The data handler for the MC partly reco decays.
    - mean_with_unc (tuple): the mean value of the data with uncertainty.
    - sigma_with_unc (tuple): the sigma value of the data with uncertainty.
    - fitter_suffix (str): A suffix to be added to the fitter name.

    Returns:
    - fitter (flarefly.F2MassFitter): The mass fitter object.
    """
    bkg_funcs = trial["bkg_funcs"].copy()
    if not isinstance(bkg_funcs, list):
        bkg_funcs = [bkg_funcs]
    sgn_funcs = trial["sgn_funcs"]
    if not isinstance(sgn_funcs, list):
        sgn_funcs = [sgn_funcs]
    label_bkg_pdf = ["Comb. bkg"]
    if trial["use_bkg_templ"]:
        label_bkg_pdf.insert(0, "Partly reco decays")
        bkg_funcs.insert(0, "kde_grid")

    mean_for_data, sigma_for_data = get_mean_sigma_for_data(mean_with_unc, sigma_with_unc, trial)

    fitter = F2MassFitter(
        data_hdl,
        sgn_funcs,
        bkg_funcs,
        name=f"b0_{fitter_suffix}",
        label_signal_pdf=[r"$\mathrm{B}^{0}$ signal"],
        label_bkg_pdf=label_bkg_pdf,
        verbosity=0
    )
    if trial["use_bkg_templ"]:
        fitter.set_background_kde(0, data_hdl_prd_bkg)
        fitter.set_background_initpar(0, "frac", 0.01, limits=[0., 1.])

    if mean_for_data is not None:  # if mean is fixed, set the particle mass
        fitter.set_particle_mass(0, mass=mean_for_data, fix=True)
    else:  # mean is free
        fitter.set_particle_mass(0, pdg_id=511, fix=False)

    if sigma_for_data is not None:  # if sigma is fixed, set the sigma
        fitter.set_signal_initpar(0, "sigma", sigma_for_data, fix=True)
    else:  # sigma is free
        fitter.set_signal_initpar(0, "sigma", 0.04, fix=False, limits=[0.001, 0.1])

    if trial["bkg_funcs"] == ["chebpol2"]:
        fitter.set_background_initpar(1, "c0", 1.)
        fitter.set_background_initpar(1, "c1", -0.3)
        fitter.set_background_initpar(1, "c2", 0.02)

    return fitter


def fit(fitter, cfg, i_trial, suffix):  # pylint: disable=too-many-locals
    """
    Perform a fit using the given fitter object and configuration.

    Parameters:
    - fitter (flarefly.F2MassFitter): The mass fitter object.
    - cfg (dict): Configuration dictionary.
    - i_trial (int): The trial number.
    - suffix (str): The suffix for the output file.

    Returns:
    - output_dict (dict): A dictionary containing the fit results.
    """
    result = fitter.mass_zfit()
    if result.converged:
        rawy, rawy_unc = fitter.get_raw_yield(0)
        if cfg["multitrial"]["bincounting_nsigma"]:  # if there is at least one nsigma
            rawy_bincounting, rawy_bincounting_unc = zip(
                *[fitter.get_raw_yield_bincounting(0, nsigma=nsigma)
                    for nsigma in cfg["multitrial"]["bincounting_nsigma"]]
            )
        else:
            rawy_bincounting, rawy_bincounting_unc = None, None
        significance, significance_unc = fitter.get_significance(0)
        soverb, soverb_unc = fitter.get_signal_over_background(0)
        mean, mean_unc = fitter.get_mass(0)
        sigma, sigma_unc = fitter.get_sigma(0)
        chi2_ndf = fitter.get_chi2_ndf()
        if cfg["save_all_fits"]:
            fig, _ = fitter.plot_mass_fit(
                style="ATLAS",
                figsize=(8, 8),
                axis_title=r"$M(\mathrm{D^-\pi^+})$ (GeV/$c^2$)"
            )
            if not os.path.exists(os.path.join(cfg["output_dir"], cfg["output_dir_fits"])):
                os.makedirs(os.path.join(cfg["output_dir"], cfg["output_dir_fits"]))
            fig.savefig(
                os.path.join(cfg["output_dir"], cfg["output_dir_fits"], f"mass_fit_{suffix}.pdf")
            )
    else:
        rawy, rawy_unc = None, None
        rawy_bincounting, rawy_bincounting_unc = [None] * len(cfg["multitrial"]["bincounting_nsigma"]), [None] * len(cfg["multitrial"]["bincounting_nsigma"])  # pylint: disable=line-too-long # noqa: 501
        significance, significance_unc = None, None
        soverb, soverb_unc = None, None
        mean, mean_unc = None, None
        sigma, sigma_unc = None, None
        chi2_ndf = None
        Logger(
            f"The fit for the trial {i_trial} did not converge, skipping this trial",
            "WARNING"
        )

    output_dict = {
        "rawy": rawy, "rawy_unc": rawy_unc, "significance": significance,
        "significance_unc": significance_unc, "soverb": soverb, "soverb_unc": soverb_unc,
        "mean": mean, "mean_unc": mean_unc, "sigma": sigma, "sigma_unc": sigma_unc,
        "chi2_ndf": chi2_ndf
    }
    for i_nsigma, nsigma in enumerate(cfg["multitrial"]["bincounting_nsigma"]):
        output_dict[f"rawy_bincounting_{nsigma}"] = rawy_bincounting[i_nsigma]
        output_dict[f"rawy_bincounting_{nsigma}_unc"] = rawy_bincounting_unc[i_nsigma]

    return output_dict


def get_rms_shift_sum_quadrature(ratio, central_ratio, rel=False):
    """
    Calculate the sum in quadrature of the RMS and shift from the central value for raw yields.

    Parameters:
        df (pandas.DataFrame): DataFrame containing the raw yields.
        cfg (dict): Configuration dictionary.
        i_pt (int): Index of pt.
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


def initialise_pars(fitter, fit_config, params):  # pylint: disable=too-many-branches
    """
    Initialise the parameters of the fitter based on the given configuration.

    Parameters:
    - fitter (flarefly.F2MassFitter): The mass fitter object.
    - fit_config (dict): The configuration dictionary.
    - params (dict): A dictionary containing the parameters for the fit.

    Returns:
    - None
    """
    if params is not None:
        signal_params = params["signal"]
        bkg_params = params["bkg"]
        for i_func, _ in enumerate(signal_params):
            fitter.set_signal_initpar(i_func, "frac", 0.1, limits=[0.0002, 1.])
            #for param, value in signal_params[i_func].items():
            #    fitter.set_signal_initpar(i_func, param, value)
            if fit_config['fix_sigma_to_mb']:
                if fit_config['sigma'] == "fixed_plus_unc":
                    fitter.set_signal_initpar(
                        i_func, "sigma",
                        signal_params[i_func]['sigma'] + signal_params[i_func]['sigma_unc'],
                        fix=True
                    )
                elif fit_config['sigma'] == "fixed_minus_unc":
                    fitter.set_signal_initpar(
                        i_func, "sigma",
                        signal_params[i_func]['sigma'] - signal_params[i_func]['sigma_unc'],
                        fix=True
                    )
                else:
                    fitter.set_signal_initpar(
                        i_func, "sigma", signal_params[i_func]['sigma'], fix=True
                    )
        #for i_func, _ in enumerate(bkg_params):
        #    for param, value in bkg_params[i_func].items():
        #        fitter.set_background_initpar(i_func, param, value)
        #        #fitter.set_signal_initpar(i_func, "frac", 0.1, limits=[0.0002, 1.])
        if fit_config['fix_corr_bkg_to_mb']:
            fitter.fix_bkg_frac_to_signal_pdf(
                0, 1, bkg_params[0]['frac'] / signal_params[1]['frac']
            )
    else:
        fitter.set_particle_mass(0, pdg_id=431)
        fitter.set_particle_mass(1, pdg_id=411)
        for i_func, _ in enumerate(fitter.get_name_signal_pdf()):
            fitter.set_signal_initpar(i_func, "sigma", 0.01, limits=[0.001, 0.05])
        for i_func, bkg_func in enumerate(fitter.get_name_background_pdf()):
            if bkg_func == "expo":
                fitter.set_background_initpar(i_func, "lam", -2)
            elif bkg_func == "chebpol2":
                fitter.set_background_initpar(i_func, "c0", 0.6)
                fitter.set_background_initpar(i_func, "c1", -0.2)
                fitter.set_background_initpar(i_func, "c2", 0.01)
            elif bkg_func == "chebpol3":
                fitter.set_background_initpar(i_func, "c0", 0.4)
                fitter.set_background_initpar(i_func, "c1", -0.2)
                fitter.set_background_initpar(i_func, "c2", -0.01)
                fitter.set_background_initpar(i_func, "c3", 0.01)

def dump_results_to_root(dfs, cfg, cut_set):  # pylint: disable=too-many-locals
    """
    Dump the results to a ROOT file.

    Parameters:
        dfs (list of pandas.DataFrame): List of dataframes containing the data for each pt bin.
        cfg (dict): Configuration dictionary.
        cut_set (dict): Dictionary containing the cut sets.

    Returns:
        None
    """
    pt_mins = cut_set["pt"]["mins"]
    pt_maxs = cut_set["pt"]["maxs"]
    pt_edges = np.asarray(pt_mins + [pt_maxs[-1]], "d")
    pt_bins = cfg["multitrial"]["pt_bins"]
    if pt_bins is None:
        pt_bins = list(range(len(pt_mins)))

    rms_shifts = []
    assigned_syst = []

    cols_to_save = [
        "rawy", "rawy_unc", "significance", "significance_unc",
        "soverb", "soverb_unc", "mean", "mean_unc", "sigma", "sigma_unc", "chi2_ndf"
    ]

    idx_assigned_syst = 0
    for i_pt in range(len(pt_mins)):
        if i_pt not in pt_bins:
            rms_shifts.append(0)
            assigned_syst.append(0)
            continue
        rms_shifts.append(get_rms_shift_sum_quadrature(dfs[i_pt], cfg, i_pt, rel=True))
        assigned_syst.append(cfg["assigned_syst"][idx_assigned_syst])
        idx_assigned_syst += 1

    if not os.path.exists(cfg["output_dir"]):
        os.makedirs(cfg["output_dir"])
    output_file_name = os.path.join(cfg["output_dir"], "raw_yields_systematic.root")
    with uproot.recreate(output_file_name) as f:
        for pt_bin in pt_bins:
            suffix = f"_{pt_mins[pt_bin] * 10:.0f}_{pt_maxs[pt_bin] * 10:.0f}"
            f[f"df{suffix}"] = dfs[pt_bin][cols_to_save]

        f["rms_shifts_sum_quadrature"] = (np.array(rms_shifts), pt_edges)
        f["assigned_syst"] = (np.array(assigned_syst), pt_edges)



def do_fit(fit_config, cfg, params=None):  # pylint: disable=too-many-locals, too-many-branches, too-many-statements # noqa: E501
    """Fit the invariant mass spectrum for a given configuration."""
    pt_min = fit_config["pt_min"]
    pt_max = fit_config["pt_max"]

    cent_min, cent_max = None, None
    if "cent_min" in fit_config and "cent_max" in fit_config:
        cent_min = fit_config["cent_min"]
        cent_max = fit_config["cent_max"]

        data_hdl = DataHandler(
            data=cfg["inputs"]["data"],
            histoname=f'h_mass_{pt_min*10:.0f}_{pt_max*10:.0f}_cent_{cent_min:.0f}_{cent_max:.0f}',
            limits=[fit_config["mins"], fit_config["maxs"]], rebin=fit_config["rebins"]
        )
    else:
        data_hdl = DataHandler(
            data=cfg["inputs"]["data"],
            histoname=f'h_mass_{pt_min*10:.0f}_{pt_max*10:.0f}',
            limits=[fit_config["mins"], fit_config["maxs"]], rebin=fit_config["rebins"]
        )

    bkg_funcs = fit_config["bkg_funcs"].copy()
    label_signal_pdfs = [r"$\mathrm{D_{s}^{+}}$ signal", r"$\mathrm{D^{+}}$ signal"]
    label_bkg_pdfs = ["Combinatorial background"]
    if fit_config["use_bkg_templ"]:
        bkg_funcs.insert(0, "hist")
        label_bkg_pdfs.insert(
            0,
            r"$\mathrm{D^{+}}\rightarrow K^{-}\pi^{+}\pi^{+}$"
            "\ncorrelated background"
        )

        data_corr_bkg = DataHandler(
            data=cfg["inputs"]["template"]["file"],
            histoname=cfg["inputs"]["template"]["hist_name"].format(
                pt_min=f"{pt_min*10:.0f}", pt_max=f"{pt_max*10:.0f}",
                cent_min=cent_min, cent_max=cent_max
            ),
            limits=[fit_config["mins"], fit_config["maxs"]], rebin=fit_config["rebins"]
        )
        fitter = F2MassFitter(
            data_hdl, name_signal_pdf=fit_config["sgn_funcs"],
            name_background_pdf=bkg_funcs,
            name=f"ds_pt{pt_min*10:.0f}_{pt_max*10:.0f}", chi2_loss=True,
            label_signal_pdf=label_signal_pdfs,
            label_bkg_pdf=label_bkg_pdfs,
            verbosity=7, tol=1.e-1
        )
        fitter.set_background_template(0, data_corr_bkg)
        fitter.set_background_initpar(0, "frac", 0.01, limits=[0., 1.])

        fitter.set_signal_initpar(0, "sigma", 0.01, limits=[0.001, 0.05])
        fitter.set_signal_initpar(1, "sigma", 0.01, limits=[0.001, 0.05])
        fitter.set_signal_initpar(0, "frac", 0.1, limits=[0.0002, 1.])
        fitter.set_signal_initpar(1, "frac", 0.1, limits=[0.0002, 1.])

    else:
        fitter = F2MassFitter(
            data_hdl, name_signal_pdf=fit_config["sgn_funcs"],
            name_background_pdf=fit_config["bkg_funcs"],
            name=f"ds_pt{pt_min*10:.0f}_{pt_max*10:.0f}", chi2_loss=True,
            label_signal_pdf=label_signal_pdfs,
            label_bkg_pdf=label_bkg_pdfs,
            verbosity=1, tol=1.e-1
        )


    fitter.set_background_initpar(1, "c0", 0.6)
    fitter.set_background_initpar(1, "c1", -0.2)
    fitter.set_background_initpar(1, "c2", 0.01)
    # signals initialisation
    fitter.set_particle_mass(0, pdg_id=431)
    fitter.set_particle_mass(1, pdg_id=411)
    initialise_pars(fitter, fit_config, params)

    n_signal = len(fit_config["sgn_funcs"])
    fit_result = fitter.mass_zfit()
    #if fit_result.converged:
    if cfg["output"]["save_all_fits"]:
        output_dir = os.path.join(
            os.path.expanduser(cfg["output"]["dir"]),
            cfg["output"]["dir_fits"]
        )
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        loc = ["lower left", "upper left"]
        ax_title = r"$M(\mathrm{KK\pi})$ GeV$/c^2$"
        fig, _ = fitter.plot_mass_fit(
            style="ATLAS",
            show_extra_info = (bkg_funcs != ["nobkg"]),
            figsize=(8, 8), extra_info_loc=loc,
            extra_info_massrange=[1.8, 2.2],
            axis_title=ax_title
        )
        figres = fitter.plot_raw_residuals(
            figsize=(8, 8), style="ATLAS",
            extra_info_loc=loc, axis_title=ax_title
        )
        for frmt in cfg["output"]["formats"]:
            if cent_min is not None and cent_max is not None:
                suffix = f"_{pt_min * 10:.0f}_{pt_max * 10:.0f}_cent_{cent_min:.0f}_{cent_max:.0f}_"  # pylint: disable=line-too-long # noqa: E501
            else:
                suffix = f"_{pt_min * 10:.0f}_{pt_max * 10:.0f}_"
            suffix += cfg["output"]["suffix"]
            fig.savefig(f"{output_dir}/ds_mass_pt{suffix}.{frmt}")
            figres.savefig(f"{output_dir}/ds_massres_pt{suffix}.{frmt}")
            if frmt == "root":
                fitter.dump_to_root(
                    f"{output_dir}/fits_{cfg['output']['suffix']}.{frmt}",
                    option="update", suffix=suffix, num=5000
                )
        plt.close(fig)
        plt.close(figres)

        fracs = fitter._F2MassFitter__get_all_fracs() # pylint: disable=protected-access
        corr_bkg_frac = fracs[1][0]
        corr_bkg_frac_err = fracs[4][0]
        corr_bkg_over_dplus_signal = fracs[1][0] / fracs[0][1]
        corr_bkg_over_dplus_signal_err = np.sqrt(
            (fracs[4][0] / fracs[1][0])**2 + (fracs[3][1] / fracs[0][1])**2
        ) * corr_bkg_over_dplus_signal
        out_dict = {
            "raw_yields": [fitter.get_raw_yield(i) for i in range(n_signal)],
            "mean": [fitter.get_mass(i) for i in range(n_signal)],
            "chi2": float(fitter.get_chi2_ndf()),
            "significance": [fitter.get_significance(i, min=1.8, max=2.2) for i in range(n_signal)],
            "signal": [fitter.get_signal(i, min=1.8, max=2.) for i in range(n_signal)],
            "background": [fitter.get_background(i, min=1.8, max=2.) for i in range(n_signal)],
            "bkg_frac": (corr_bkg_frac, corr_bkg_frac_err),
            "fracs": fracs,
            "converged": fit_result.converged, 
            "corr_bkg_over_dplus_signal": (corr_bkg_over_dplus_signal, corr_bkg_over_dplus_signal_err)  # pylint: disable=line-too-long # noqa: E501
        }

        if fit_config["sgn_funcs"][0] == "doublegaus":
            out_dict["sigma1"] = [
                fitter.get_signal_parameter(i, "sigma1") for i in range(n_signal)
                ]
            out_dict["sigma2"] = [
                fitter.get_signal_parameter(i, "sigma2") for i in range(n_signal)
                ]
            out_dict["frac1"] = [
                fitter.get_signal_parameter(i, "frac1") for i in range(n_signal)
                ]
        elif fit_config["sgn_funcs"][0] == "gaussian":
            out_dict["sigma"] = [
                fitter.get_sigma(i) for i in range(n_signal)]
        elif fit_config["sgn_funcs"][0] == "doublecb":
            out_dict["sigma"] = [
                fitter.get_signal_parameter(i, "sigma") for i in range(n_signal)
                ]
            out_dict["alphar"] = [
                fitter.get_signal_parameter(i, "alphar") for i in range(n_signal)
                ]
            out_dict["alphal"] = [
                fitter.get_signal_parameter(i, "alphal") for i in range(n_signal)
                ]
            out_dict["nl"] = [
                fitter.get_signal_parameter(i, "nl") for i in range(n_signal)
                ]
            out_dict["nr"] = [
                fitter.get_signal_parameter(i, "nr") for i in range(n_signal)
                ]
        elif fit_config["sgn_funcs"][0] == "doublecbsymm":
            out_dict["sigma"] = [
                fitter.get_signal_parameter(i, "sigma") for i in range(n_signal)
                ]
            out_dict["alpha"] = [
                fitter.get_signal_parameter(i, "alpha") for i in range(n_signal)
                ]
            out_dict["n"] = [
                fitter.get_signal_parameter(i, "n") for i in range(n_signal)
                ]
        elif fit_config["sgn_funcs"][0] == "genergausexptailsymm":
            out_dict["sigma"] = [
                fitter.get_signal_parameter(i, "sigma") for i in range(n_signal)
                ]
            out_dict["alpha"] = [
                fitter.get_signal_parameter(i, "alpha") for i in range(n_signal)
                ]
    else:
        out_dict = {
            "raw_yields": [None] * n_signal,
            "sigma": [None] * n_signal,
            "mean": [None] * n_signal,
            "chi2": None,
            "significance": [None] * n_signal,
            "signal": [None] * n_signal,
            "background": [None] * n_signal,
            "fracs": None,
            "converged": False
        }
    return out_dict, fitter


def get_matching_trial(trial_to_match, trials, cent_min, cent_max, pt_min, pt_max):  # pylint: disable=too-many-arguments, too-many-positional-arguments
    """
    Find a trial in the list of trials that matches the given trial_to_match 
    based on the specified centrality range (cent_min to cent_max).

    Args:
        trial_to_match (dict): The trial dictionary to match against.
        trials (list of dict): A list of trial dictionaries to search through.
        cent_min (float): The minimum centrality value for matching.
        cent_max (float): The maximum centrality value for matching.
        pt_min (float): The minimum pt value for matching.
        pt_max (float): The maximum

    Returns:
        dict: The matching trial dictionary from the list of trials.

    Raises:
        ValueError: If no matching trial is found within the specified centrality range.
    """
    for trial in trials:
        if trial["cent_min"] != cent_min or trial["cent_max"] != cent_max\
                or trial["pt_min"] != pt_min or trial["pt_max"] != pt_max:
            continue
        match = True
        for key in trial:
            if key not in ["cent_min", "cent_max"]:
                if trial[key] != trial_to_match[key]:
                    match = False
                    break
        if match:
            return trial
    raise ValueError("No matching trial found. Trial to match: ", trial_to_match)


def run_pt_trial(mb_trial, trials, cfg, pt_min, pt_max, cent_mins, cent_maxs):  # pylint: disable=too-many-locals, too-many-positional-arguments, too-many-arguments
    """
    Perform a trial run by fitting the MB sample and then fitting the different centrality bins.

    Parameters:
    - mb_trial (dict): The minimum bias trial data.
    - trials (list): A list of trial data for different centrality bins.
    - cfg (dict): Configuration parameters for the fitting process.
    - pt_min (float): The minimum pt value for the trial.
    - pt_max (float): The maximum pt value for the trial.
    - cent_mins (list): A list of minimum centrality values for the bins.
    - cent_maxs (list): A list of maximum centrality values for the bins.

    Returns:
    - tuple: A tuple containing:
        - mb_result (dict): The result of the fit for the MB sample.
        - results (list): A list of results for the fits of the trials 
          for different centrality bins.
    """
    # we firstly fit the MB sample
    mb_result, mb_fitter = do_fit(mb_trial, cfg)
    mb_params = {}
    mb_params['signal'] = mb_fitter.get_signal_pars()
    mb_params['bkg'] = mb_fitter.get_bkg_pars()
    mb_params['signal'][0]['sigma_unc'] = mb_fitter.get_sigma(0)[1]
    mb_params['signal'][1]['sigma_unc'] = mb_fitter.get_sigma(1)[1]
    results = []
    for cent_min, cent_max in zip(cent_mins, cent_maxs):
        trial_cent = get_matching_trial(mb_trial, trials, cent_min, cent_max, pt_min, pt_max)
        result, _ = do_fit(trial_cent, cfg, mb_params)
        trial_cent_cfg = trial_cent.copy()
        for key in trial_cent:
            trial_cent_cfg[key + "_cfg"] = trial_cent_cfg.pop(key)
        result.update(trial_cent_cfg)
        results.append(result)
    return mb_result, results

def multi_trial(config_file_name: str, draw=False):  # pylint: disable=too-many-locals, too-many-statements
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

    multitrial_cfg = cfg["multitrial"]

    # define all the trials
    trials = list(itertools.product(*(multitrial_cfg[var] for var in MULTITRIAL_PARAMS)))
    trials = itertools.product(trials, [*zip(cent_mins, cent_maxs)], [*zip(pt_mins, pt_maxs)])
    trials = [(*trial, *cent, *pt) for trial, cent, pt in trials]
    trials = [dict(
        zip(MULTITRIAL_PARAMS + ["cent_min", "cent_max", "pt_min", "pt_max"], trial)
    ) for trial in trials]

    mb_trials = list(itertools.product(*(multitrial_cfg[var] for var in MULTITRIAL_PARAMS)))
    mb_trials = itertools.product(mb_trials, [(0, 100)], [*zip(pt_mins, pt_maxs)])
    mb_trials = [(*trial, *cent, *pt) for trial, cent, pt in mb_trials]
    mb_trials = [dict(
        zip(MULTITRIAL_PARAMS + ["cent_min", "cent_max", "pt_min", "pt_max"], trial)
    ) for trial in mb_trials]

    if not draw:
        futures = []
        with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
            for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
                if multitrial_cfg["pt_bins"] is not None and i_pt not in multitrial_cfg["pt_bins"]:
                    continue

                for mb_trial in mb_trials:
                    if mb_trial["pt_min"] != pt_min or mb_trial["pt_max"] != pt_max:
                        continue
                    futures.append(executor.submit(
                        run_pt_trial, mb_trial, trials, cfg,
                        pt_min,pt_max, cent_mins, cent_maxs
                    ))

        df_mb = []
        df_cent = []
        for future in futures:
            mb_result, cent_result = future.result()
            df_mb.append(mb_result)
            for result in cent_result:
                df_cent.append(result)
        df_mb = pd.DataFrame(df_mb)
        df_cent = pd.DataFrame(df_cent)
        df_mb.to_parquet(os.path.join(cfg["output"]["dir"], "mb_results.parquet"))
        df_cent.to_parquet(os.path.join(cfg["output"]["dir"], "cent_results.parquet"))
    else:
        df_mb = pd.read_parquet(os.path.join(cfg["output"]["dir"], "mb_results.parquet"))
        df_cent = pd.read_parquet(os.path.join(cfg["output"]["dir"], "cent_results.parquet"))

    draw_multitrial(df_cent, cfg)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fit Ds')
    parser.add_argument('configFile', type=str, help='Path to the configuration file')
    parser.add_argument('--draw', action='store_true', help='Only draw the results')
    args = parser.parse_args()

    # "bincounting_nsigma" removed so that we do not fit multiple times for each nsigma
    MULTITRIAL_PARAMS = [
        "mins", "maxs", "rebins", "sgn_funcs", "bkg_funcs", "sigma",
        "fix_sigma_to_mb", "fix_corr_bkg_to_mb", "mean", "use_bkg_templ"
    ]
    multi_trial(args.configFile, args.draw)
