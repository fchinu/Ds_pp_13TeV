"""
This script performs a multitrial procedure to get the raw yields systematic uncertainty.

Usage:
    python get_ry_systematic.py configFile
"""

import re
import argparse
import itertools
import os
import numpy as np
import yaml
import uproot
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
from flarefly.utils import Logger
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Rectangle
import pandas as pd


def draw_multitrial(df_multitrial, cfg, pt_min, pt_max, idx_assigned_syst):  # pylint: disable=too-many-locals, too-many-statements # noqa: 501
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
    n_trials = len(df_multitrial)
    n_bincounts = len(multitrial_cfg['bincounting_nsigma'])
    x_axis_range = n_trials * (n_bincounts + 1) + 1

    fig, axs = plt.subplots(2, 2, figsize=(20, 15))
    with uproot.open(cfg["reference_fits"]) as f:
        h_rawy = f["h_rawyields"]
        h_sigma = f["h_sigmas"]

    i_pt = np.digitize((pt_min + pt_max) / 2, h_rawy.axis().edges()) - 1
    central_rawy = h_rawy.values()[i_pt]
    central_rawy_unc = h_rawy.errors()[i_pt]
    central_sigma = h_sigma.values()[i_pt]
    central_sigma_unc = h_sigma.errors()[i_pt]

    # Plot the results
    axs[0, 0].errorbar(
        x=range(1, n_trials + 1), y=df_multitrial["rawy"],
        yerr=df_multitrial["rawy_unc"], fmt='o', label='Fit',
        zorder=2
    )

    for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
        axs[0, 0].errorbar(
            x=list(range(n_trials * (i_nsigma + 1) + 1, n_trials * (i_nsigma + 2) + 1)),
            y=df_multitrial[f"rawy_bincounting_{nsigma}"],
            yerr=df_multitrial[f"rawy_bincounting_{nsigma}_unc"], fmt='o',
            label=fr'Bin counting {nsigma}$\sigma$',
            zorder=1
        )

    # Draw the central values
    axs[0, 0].axhline(y=central_rawy, color='r', linestyle='--')
    axs[0, 0].add_patch(
        Rectangle(
            (0, central_rawy - central_rawy_unc),
            x_axis_range, 2 * central_rawy_unc,
            color='r', alpha=0.3, zorder=0,
            label=r'Central value $\pm$ uncertainty'
        )
    )

    axs[0, 0].set_xlim(0, x_axis_range)
    axs[0, 0].set_xlabel('Trial', fontsize=14)
    axs[0, 0].set_ylabel('Raw yield', fontsize=14)
    axs[0, 0].legend(fontsize=12)

    # Draw the raw yields distribution
    axs[0, 1].hist(
        df_multitrial["rawy"], bins=30, alpha=0.7, label='Fit',
        histtype='stepfilled', ec="k", linewidth=2, zorder=2
    )

    for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
        axs[0, 1].hist(
            df_multitrial[f"rawy_bincounting_{nsigma}"],
            bins=30,
            alpha=0.3,
            label=fr'Bin Counting {nsigma}$\sigma$',
            histtype='stepfilled',
            ec="k",
            zorder=1
        )

    # Draw information
    info = 'Fit:\n'
    info += fr'$\mu =$ {np.mean(df_multitrial["rawy"]):.3f}''\n'
    info += fr'$\sigma =$ {np.std(df_multitrial["rawy"]):.3f}''\n'
    for i_nsigma, nsigma in enumerate(multitrial_cfg['bincounting_nsigma']):
        info += fr'{nsigma}$\sigma$ Bin counting:''\n'
        info += fr'$\mu =$ {np.mean(df_multitrial[f"rawy_bincounting_{nsigma}"]):.3f}''\n'
        info += fr'$\sigma =$ {np.std(df_multitrial[f"rawy_bincounting_{nsigma}"]):.3f}''\n'
    anchored_text_fit = AnchoredText(
        info,
        loc='upper left',
        frameon=False
    )

    # Draw the rms + shift from the central value
    rms_shift = get_rms_shift_sum_quadrature(df_multitrial, cfg, i_pt)
    axs[0, 1].axvline(x=central_rawy, color='r', linestyle='--')
    axs[0, 1].add_patch(
        Rectangle(
            (central_rawy - rms_shift, 0),
            2 * rms_shift, axs[0, 1].get_ylim()[1],
            color='r', alpha=0.3, zorder=0,
            label=r'$\mathrm{\sqrt{RMS^2 + \Delta^2}}$'
        )
    )
    axs[0, 1].add_artist(anchored_text_fit)

    # Draw the assigned systematic uncertainty
    axs[0, 1].add_patch(
        Rectangle(
            (central_rawy - cfg["assigned_syst"][idx_assigned_syst] * central_rawy, 0),
            2 * cfg["assigned_syst"][idx_assigned_syst] * central_rawy, axs[0, 1].get_ylim()[1],
            color='limegreen', alpha=0.3, zorder=0,
            label='Assigned syst.'
        )
    )

    axs[0, 1].set_xlabel('Raw yields', fontsize=14)
    axs[0, 1].set_ylabel('Counts', fontsize=14)
    axs[0, 1].legend(fontsize=12, loc='upper right')

    x_min = min(
        df_multitrial["rawy"].min(),
        *[df_multitrial[f"rawy_bincounting_{nsigma}"].min()
            for nsigma in multitrial_cfg['bincounting_nsigma']],
        central_rawy - rms_shift
    )
    x_max = max(
        df_multitrial["rawy"].max(),
        *[df_multitrial[f"rawy_bincounting_{nsigma}"].max()
            for nsigma in multitrial_cfg['bincounting_nsigma']],
        central_rawy + rms_shift
    )
    axs[0, 1].set_xlim(x_min * 0.9, x_max * 1.1)

    # Draw the peak width
    axs[1, 0].errorbar(
        x=range(1, len(df_multitrial["sigma"]) + 1),
        y=df_multitrial["sigma"],
        yerr=df_multitrial["sigma_unc"], fmt='o',
        zorder=2
    )

    # Draw the central values
    axs[1, 0].axhline(y=central_sigma, color='r', linestyle='--')
    axs[1, 0].add_patch(
        Rectangle(
            (0, central_sigma - central_sigma_unc),
            x_axis_range, 2 * central_sigma_unc,
            color='r', alpha=0.3, zorder=0,
            label=r'Central value $\pm$ uncertainty'
        )
    )

    axs[1, 0].set_xlim(0, x_axis_range)
    axs[1, 0].set_xlabel('Trial', fontsize=14)
    axs[1, 0].set_ylabel('Width ($GeV/c^2$)', fontsize=14)
    axs[1, 0].legend(fontsize=12)

    # Draw the chi2/ndf
    axs[1, 1].scatter(
        x=range(1, len(df_multitrial["chi2_ndf"]) + 1),
        y=df_multitrial["chi2_ndf"]
    )
    axs[1, 1].set_xlim(0, x_axis_range)
    axs[1, 1].set_xlabel('Trial', fontsize=14)
    axs[1, 1].set_ylabel(r'$\chi^2/$ndf', fontsize=14)

    plt.show()
    fig.savefig(
        os.path.join(cfg["output_dir"], f'fig_{pt_min*10:.0f}_{pt_max*10:.0f}.png'),
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


def build_data_handlers(trial, cfg, pt_min, pt_max, cent_min, cent_max):
    """
    Build data handlers for the given trial.

    Parameters:
    - trial (dict): A dictionary containing the trial parameters (mins, maxs, sgn_funcs, bkg_funcs,
        sigma, mean, use_bkg_templ, bincounting_nsigma).
    - df_pt (pandas.DataFrame): The data frame for the data.
    - df_mc_sig_pt (pandas.DataFrame): The data frame for the MC signal.
    - df_mc_prd_bkg_pt (pandas.DataFrame): The data frame for the MC partly reco decays.
        Can be None if use_bkg_templ is False.

    Returns:
    - data_hdl (flarefly.DataHandler): The data handler for the data.
    - data_hdl_mc (flarefly.DataHandler): The data handler for the MC signal.
    - data_hdl_bkg (flarefly.DataHandler): The data handler for the MC partly reco decays,
        or None if use_bkg_templ is False.
    """
    # data
    data_hdl = DataHandler(
        cfg["inputs"]["data"], var_name="fM",
        histoname=f"h_mass_{pt_min*10:.0f}_{pt_max*10:.0f}_{cent_min:.0f}_{cent_max:.0f}",
        limits=[trial["mins"], trial["maxs"]]
    )

    # mc signal
    data_hdl_mc = DataHandler(df_mc_sig_pt, var_name="fM",
                              limits=[trial["mins"], trial["maxs"]])

    # mc partly reco decays
    if trial["use_bkg_templ"]:
        data_hdl_bkg = DataHandler(
            df_mc_prd_bkg_pt, var_name="fM",
            limits=[trial["mins"], trial["maxs"]]
        )
    else:
        data_hdl_bkg = None

    return data_hdl, data_hdl_mc, data_hdl_bkg


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


def get_rms_shift_sum_quadrature(df, cfg, i_pt, rel=False):
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
    with uproot.open(cfg["reference_fits"]) as f:
        h_rawy = f["h_rawyields"]
    central_rawy = h_rawy.values()[i_pt]

    if rel:
        return np.sqrt(
            np.std(df["rawy"])**2 +
            (np.mean(df["rawy"]) - central_rawy)**2
        ) / central_rawy

    return np.sqrt(
        np.std(df["rawy"])**2 +
        (np.mean(df["rawy"]) - central_rawy)**2
    )


def dump_results_to_root(dfs, cfg, cut_set):
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


def multi_trial(config_file_name: str):  # pylint: disable=too-many-locals, too-many-statements
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

    multitrial_cfg = cfg["multitrial"]

    # define all the trials
    trials = list(itertools.product(*(multitrial_cfg[var] for var in MULTITRIAL_PARAMS)))
    trials = [dict(zip(MULTITRIAL_PARAMS, trial)) for trial in trials]

    rms_shifts = []  # sum in quadrature of the rms and shifts of the raw yield distributions
    dfs = []

    idx_assigned_syst = 0
    for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        if multitrial_cfg["pt_bins"] is not None and i_pt not in multitrial_cfg["pt_bins"]:
            rms_shifts.append(0)
            dfs.append(None)
            continue

        if any(multitrial_cfg["use_bkg_templ"]):
            df_mc_prd_bkg_pt = df_mc_prd_bkg.query(f"{pt_min} < fPt < {pt_max}")
        else:
            df_mc_prd_bkg_pt = None

        trial_rows = []
        for i_trial, trial in enumerate(trials):
            print(trial)
            data_hdl, data_hdl_mc_sig, data_hdl_prd_bkg = build_data_handlers(
                trial, cfg, pt_min, pt_max
            )

            suffix = f"{pt_min*10:.0f}_{pt_max*10:.0f}_{i_trial}"
            if "free" not in trial["mean"] or "free" not in trial["sigma"]:
                mean_with_unc, sigma_with_unc = get_mc_sigma_width_with_unc(
                    trial, data_hdl_mc_sig, suffix)
                # if the fit did not converge, skip the trial
                if mean_with_unc is None and sigma_with_unc is None:
                    continue
            else:
                mean_with_unc, sigma_with_unc = None, None

            fitter = build_fitter(
                trial, data_hdl, data_hdl_prd_bkg, mean_with_unc, sigma_with_unc, suffix
            )
            trial_dict = fit(fitter, cfg, i_trial, suffix)
            trial_renamed = trial.copy()
            trial_renamed['sigma_type'] = trial_renamed.pop('sigma')
            trial_renamed['mean_type'] = trial_renamed.pop('mean')
            trial_dict.update(trial_renamed)
            trial_rows.append(trial_dict)

        # save the results as pandas dataframe
        df_trials = pd.DataFrame(trial_rows)
        if not os.path.exists(cfg["output_dir"]):
            os.makedirs(cfg["output_dir"])
        df_trials.to_parquet(
            os.path.join(cfg["output_dir"], f"raw_yields_{pt_min*10:.0f}_{pt_max*10:.0f}.parquet")
        )
        dfs.append(df_trials)

        draw_multitrial(df_trials, cfg, pt_min, pt_max, idx_assigned_syst)
        idx_assigned_syst += 1
    dump_results_to_root(dfs, cfg, cut_set)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fit Ds')
    parser.add_argument('configFile', type=str, help='Path to the configuration file')
    args = parser.parse_args()

    # "bincounting_nsigma" removed so that we do not fit multiple times for each nsigma
    MULTITRIAL_PARAMS = ["mins", "maxs", "sgn_funcs", "bkg_funcs", "sigma", "mean", "use_bkg_templ"]
    multi_trial(args.configFile)