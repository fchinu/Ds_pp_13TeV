"""
Script for fitting D+, D0 and Ds+ invariant-mass spectra.

run: python get_raw_yields_ds_dplus_flarefly.py config.yaml
"""

import argparse
import itertools
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.offsetbox import AnchoredText  # noqa: E402
import dataclasses  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import uproot  # noqa: E402
import yaml  # noqa: E402
import ROOT  # noqa: E402
import tensorflow as tf
import zfit
from flarefly.data_handler import DataHandler  # noqa: E402
from flarefly.fitter import F2MassFitter  # noqa: E402
from flarefly.utils import Logger  # noqa: E402
# pylint: disable=no-member


import dataclasses
import numpy as np
import ROOT

@dataclasses.dataclass
class BinsHelper:
    """
    Helper class for binning.

    Parameters:
    - mins (list or array-like): A list or array of minimum values for each bin.
    - maxs (list or array-like): A list or array of maximum values for each bin.

    Attributes:
    - mins (list or array-like): Stores the minimum values for each bin.
    - maxs (list or array-like): Stores the maximum values for each bin.
    - bins (list of tuples): A list of tuples where each tuple contains
        the minimum and maximum values for a bin.
    - edges (numpy.ndarray): An array of bin edges.
    - n_bins (int): The number of bins.
    """
    mins: list
    maxs: list
    
    def __post_init__(self):
        self.bins = [*zip(self.mins, self.maxs)]
        self.edges = np.asarray(self.mins + [self.maxs[-1]], 'd')
        self.n_bins = len(self.mins)


class HistHandler:  # pylint: disable=too-many-instance-attributes
    """
    Class designed to handle the creation, manipulation, and storage
    of histograms for various observables.

    Parameters:
    - pt_mins (list or array-like): Minimum values for pT bins.
    - pt_maxs (list or array-like): Maximum values for pT bins.
    - cent_mins (list or array-like, optional): Minimum values for centrality bins.
    - cent_maxs (list or array-like, optional): Maximum values for centrality bins.
    """

    def __init__(self, pt_mins, pt_maxs, cent_mins=None, cent_maxs=None):
        self._pt_info = BinsHelper(pt_mins, pt_maxs)
        self._cent_info = BinsHelper(cent_mins, cent_maxs) if cent_mins is not None else None
        self._pt_axis_title = '#it{p}_{T} (GeV/#it{c})'
        self._n_ev = None
        self._histos = {}
        
        self._observable_config = self._get_observable_config()
        self._build_histos()

    def _get_observable_config(self):
        """Get configuration for all observables and their axis titles."""
        return {
            # Common observables (have both _ds and _dplus versions)
            "common": {
                "raw_yields": "Raw yields",
                "sigma": "Width (GeV/#it{c}^{2})",
                "sigma1": "Width_{1} (GeV/#it{c}^{2})",
                "sigma2": "Width_{2} (GeV/#it{c}^{2})",
                "frac1": "Gaussian fraction",
                "mean": "Mean (GeV/#it{c}^{2})",
                "raw_yield_over_ev": "Raw yields / N_{ev}",
                "significance": "Significance (3#sigma)",
                "significance_over_sqrt_ev": "Significance / #sqrt{N_{ev}}",
                "s_over_b": "S/B (3#sigma)",
                "signal": "Signal (3#sigma)",
                "background": "Background (3#sigma)",
                "alphal": "#alpha_{l}",
                "alphar": "#alpha_{r}",
                "nl": "n_{l}",
                "nr": "n_{r}",
                "alpha": "#alpha",
                "n": "n"
            },
            # Non-common observables (single version only)
            "not_common": {
                "chi2": "#chi^{2}/#it{ndf}",
                "sigma_ratio_second_first_peak": "Width second peak / width first peak",
                "corr_bkg_frac_over_dplus_0": "Corr. bkg / D^{+} signal",
                "corr_bkg_frac_0": "Corr. bkg fraction"
            }
        }

    def _create_histo_name(self, obs, cent_min=None, cent_max=None):
        """Generate histogram name based on observable and centrality."""
        base_name = f'h_{obs}'
        if cent_min is not None and cent_max is not None:
            base_name += f'_{cent_min:.0f}_{cent_max:.0f}'
        return base_name

    def _create_histogram(self, obs, y_title):
        """Create histograms for a given observable."""
        histos = []
        
        if self._cent_info is None:
            hist = ROOT.TH1D(
                self._create_histo_name(obs),
                f';{self._pt_axis_title};{y_title}',
                self._pt_info.n_bins, self._pt_info.edges
            )
            hist.SetDirectory(0)
            histos.append(hist)
        else:
            for cent_min, cent_max in zip(self._cent_info.mins, self._cent_info.maxs):
                hist = ROOT.TH1D(
                    self._create_histo_name(obs, cent_min, cent_max),
                    f';{self._pt_axis_title};{y_title}',
                    self._pt_info.n_bins, self._pt_info.edges
                )
                hist.SetDirectory(0)
                histos.append(hist)
        
        return histos

    def _build_histos(self):
        """Build histograms for all observables."""
        # Create common observables (both _ds and _dplus versions)
        for obs, y_title in self._observable_config["common"].items():
            self._histos[f"{obs}_ds"] = self._create_histogram(f"{obs}_ds", y_title)
            self._histos[f"{obs}_dplus"] = self._create_histogram(f"{obs}_dplus", y_title)
        
        # Create non-common observables
        for obs, y_title in self._observable_config["not_common"].items():
            self._histos[obs] = self._create_histogram(obs, y_title)

    def set_n_ev(self, n_ev):
        """Set the number of events."""
        self._n_ev = n_ev

    def _get_centrality_index(self, row):
        """Get centrality index from row data."""
        if "cent_min" in row and "cent_max" in row:
            cent_tuple = (row["cent_min"], row["cent_max"])
            return self._cent_info.bins.index(cent_tuple)
        return 0

    def _set_histogram_values(self, hist, i_pt, value, error=None):
        """Set histogram bin content and error."""
        hist.SetBinContent(i_pt + 1, value)
        if error is not None:
            hist.SetBinError(i_pt + 1, error)

    def _get_available_observables(self, row):
        """Determine which observables are available in the row data."""
        observables = ["raw_yields", "mean", "significance", "signal", "background"]
        
        # Check sigma variants
        if "sigma" in row:
            observables.append("sigma")
        else:
            observables.extend(["sigma1", "sigma2", "frac1"])
        
        # Check Crystal Ball parameters
        if "alphal" in row:
            observables.extend(["alphal", "alphar", "nl", "nr"])
        if "alpha" in row:
            observables.extend(["alpha", "n"])
        
        return observables

    def _process_common_observables(self, row, i_pt, i_cent, observables):
        """Process common observables (those with _ds and _dplus versions)."""
        for obs in observables:
            if obs not in row:
                continue
                
            # Set ds values
            self._set_histogram_values(
                self._histos[f"{obs}_ds"][i_cent], i_pt,
                row[obs][0][0], row[obs][0][1]
            )
            
            # Set dplus values if available
            if len(row[obs]) > 1:
                self._set_histogram_values(
                    self._histos[f"{obs}_dplus"][i_cent], i_pt,
                    row[obs][1][0], row[obs][1][1]
                )

    def _process_derived_observables(self, row, i_pt, i_cent):
        """Process derived observables that require calculation."""
        # Raw yield over events
        for suffix in ["_ds", "_dplus"]:
            idx = 0 if suffix == "_ds" else 1
            if len(row["raw_yields"]) > idx:
                value = row["raw_yields"][idx][0] / self._n_ev
                error = row["raw_yields"][idx][1] / self._n_ev
                self._set_histogram_values(
                    self._histos[f"raw_yield_over_ev{suffix}"][i_cent],
                    i_pt, value, error
                )

        # Significance over sqrt(events)
        for suffix in ["_ds", "_dplus"]:
            idx = 0 if suffix == "_ds" else 1
            if len(row["significance"]) > idx:
                value = row["significance"][idx][0] / np.sqrt(self._n_ev)
                error = row["significance"][idx][1] / np.sqrt(self._n_ev)
                self._set_histogram_values(
                    self._histos[f"significance_over_sqrt_ev{suffix}"][i_cent],
                    i_pt, value, error
                )

        # Signal over background
        for suffix in ["_ds", "_dplus"]:
            idx = 0 if suffix == "_ds" else 1
            if (len(row["signal"]) > idx and len(row["background"]) > idx and 
                row["background"][idx][0] != 0):
                value = row["signal"][idx][0] / row["background"][idx][0]
                error = row["signal"][idx][1] / row["background"][idx][0]
                self._set_histogram_values(
                    self._histos[f"s_over_b{suffix}"][i_cent],
                    i_pt, value, error
                )

        # Sigma ratio (second peak / first peak)
        if ("sigma" in row and len(row["sigma"]) > 1 and 
            row["sigma"][0][0] != 0):
            value = row["sigma"][1][0] / row["sigma"][0][0]
            self._set_histogram_values(
                self._histos["sigma_ratio_second_first_peak"][i_cent],
                i_pt, value
            )

    def _process_not_common_observables(self, row, i_pt, i_cent):
        """Process non-common observables."""
        # Chi2
        if "chi2" in row:
            value = row["chi2"] if not isinstance(row["chi2"], (tuple, list)) else row["chi2"][0]
            error = None if not isinstance(row["chi2"], (tuple, list)) else row["chi2"][1]
            self._set_histogram_values(self._histos["chi2"][i_cent], i_pt, value, error)

        # Correlated background observables
        corr_bkg_cols = [col for col in row.keys() if "corr_bkg_frac" in col and "_cfg" not in col]
        for obs in corr_bkg_cols:
            if obs in self._histos:
                value = row[obs] if not isinstance(row[obs], (tuple, list)) else row[obs][0]
                error = None if not isinstance(row[obs], (tuple, list)) else row[obs][1]
                self._set_histogram_values(self._histos[obs][i_cent], i_pt, value, error)

    def set_histos(self, df):
        """Set histogram bin contents and errors based on the provided DataFrame."""
        for _, row in df.iterrows():
            i_pt = row["i_pt"]
            i_cent = self._get_centrality_index(row)
            
            # Get available observables for this row
            available_observables = self._get_available_observables(row)
            
            # Process different types of observables
            self._process_common_observables(row, i_pt, i_cent, available_observables)
            self._process_derived_observables(row, i_pt, i_cent)
            self._process_not_common_observables(row, i_pt, i_cent)

    def dump_to_root(self, output_file):
        """Dump histograms to a ROOT file."""
        with ROOT.TFile(output_file, "RECREATE") as outfile:
            for key in self._histos:
                outfile.mkdir(key)
            for histos in self._histos.values():
                for hist in histos:
                    hist.Write()

    @property
    def obs_common(self):
        """Get list of common observable names for backward compatibility."""
        return list(self._observable_config["common"].keys())
    
    @property
    def axes_titles_common(self):
        """Get list of common observable axis titles for backward compatibility."""
        return list(self._observable_config["common"].values())
    
    @property
    def obs_not_common(self):
        """Get list of non-common observable names for backward compatibility."""
        return list(self._observable_config["not_common"].keys())
    
    @property
    def axes_titles_not_common(self):
        """Get list of non-common observable axis titles for backward compatibility."""
        return list(self._observable_config["not_common"].values())


# pylint: disable=too-many-arguments
def add_info_on_canvas(axs, loc, pt_min, pt_max, cent_min=None, cent_max=None, fitter=None):
    """
    Helper method to add text on flarefly mass fit plot

    Parameters
    ----------
    - axs: matplotlib.figure.Axis
        Axis instance of the mass fit figure

    - loc: str
        Location of the info on the figure

    - pt_min: float
        Minimum pT value in the pT range

    - pt_max: float
        Maximum pT value in the pT range
    
    - cent_min: float
        Minimum centrality value (default: None)

    - cent_max: float
        Maximum centrality value (default: None)

    - fitter: F2MassFitter
        Fitter instance allowing to access chi2 and ndf if wanted
    """
    xspace = " "
    text = xspace
    if fitter is not None:
        chi2 = fitter.get_chi2()
        ndf = fitter.get_ndf()
        text += fr"$\chi^2 / \mathrm{{ndf}} =${chi2:.2f} / {ndf} $\simeq$ {chi2/ndf:.2f}""\n"

    text += "\n\n"
    text += xspace + fr"{pt_min:.0f} < $p_{{\mathrm{{T}}}}$ < {pt_max:.0f} GeV/$c$, $|y|$ < 0.5""\n"
    if cent_min is not None and cent_max is not None:
        text += xspace + fr"{cent_min:.0f}% < Cent. < {cent_max:.0f}%" + "\n"

    anchored_text = AnchoredText(text, loc=loc, frameon=False)
    axs.add_artist(anchored_text)


def load_inputs(cfg, cut_set):
    """
    Load invariant mass histograms and event histogram from a ROOT file.

    Parameters:
    cfg (dict): Configuration dictionary containing the path to the input data file.
    cut_set (dict): Dictionary containing the cuts applied to the data.

    Returns:
    tuple: A tuple containing:
        - h_mass (list): List of invariant mass histograms.
        - h_ev (uproot.models.TH1): Event histogram.
    """
    pt_mins = cut_set["pt"]["min"]
    pt_maxs = cut_set["pt"]["max"]
    cent_mins, cent_maxs = None, None
    if "cent" in cut_set:
        cent_mins = cut_set["cent"]["min"] 
        cent_maxs = cut_set["cent"]["max"] 

    # load inv-mass histos
    h_mass = []
    with uproot.open(cfg["inputs"]["data"]) as f:
        h_ev = f["h_ev"]

        for pt_min, pt_max in zip(pt_mins, pt_maxs):
            if cent_mins is not None and cent_maxs is not None:
                for cent_min, cent_max in zip(cent_mins, cent_maxs):
                    suffix = f"{pt_min*10:.0f}_{pt_max*10:.0f}_cent_{cent_min:.0f}_{cent_max:.0f}"
                    h_mass.append(f[f'h_mass_{suffix}'])
            else:
                suffix = f"{pt_min*10:.0f}_{pt_max*10:.0f}"
                h_mass.append(f[f'h_mass_{suffix}'])

    return h_mass, h_ev


def create_fit_configs(cfg, cut_set):
    """
    Generate a list of fit configurations based on the provided configuration and cut sets.

    Args:
        cfg (dict): Configuration dictionary containing fit settings.
        cut_set (dict): Dictionary containing cut sets for pt and optionally cent.

    Returns:
        - fit_configs (list): A list of dictionaries (one per bin of pt-centrality),
            each representing a fit configuration with keys:
            - "cent_min" (optional): Minimum centrality cut.
            - "cent_max" (optional): Maximum centrality cut.
            - "pt_min": Minimum pt cut.
            - "pt_max": Maximum pt cut.
            - "mass_min": Minimum mass cut.
            - "mass_max": Maximum mass cut.
            - "rebin": Rebin factor.
            - "signal_func": Signal function.
            - "fix_ds_sigma": Fix Ds sigma.
            - "fix_dplus_sigma": Fix Dplus sigma.
            - "fix_sigma_dplus_to_ds": Fix sigma ratio of Dplus to Ds.
            - "bkg_func": Background function.
            - "use_bkg_templ": Use background template.
    """
    fit_configs = [(idx, *cfgs) for idx, cfgs in enumerate(zip(
        cut_set["pt"]["min"],
        cut_set["pt"]["max"],
        cfg["fit_configs"]["mass"]["mins"],
        cfg["fit_configs"]["mass"]["maxs"],
        cfg["fit_configs"]["rebin"],
        cfg["fit_configs"]["signal"]["signal_funcs"],
        cfg["fit_configs"]["signal"]["fix_ds_sigma"],
        cfg["fit_configs"]["signal"]["fix_dplus_sigma"],
        cfg["fit_configs"]["signal"]["fix_sigma_dplus_to_ds"],
        cfg["fit_configs"]["bkg"]["bkg_funcs"],
        cfg["fit_configs"]["bkg"]["use_bkg_templ"],
        *[sgn_func[par]["init"] for sgn_func in cfg["fit_configs"]["signal"]["par_init_limit"] for par in sgn_func],
        *[sgn_func[par]["min"] for sgn_func in cfg["fit_configs"]["signal"]["par_init_limit"] for par in sgn_func],
        *[sgn_func[par]["max"] for sgn_func in cfg["fit_configs"]["signal"]["par_init_limit"] for par in sgn_func],
        *[sgn_func[par]["fix"] for sgn_func in cfg["fit_configs"]["signal"]["par_init_limit"] for par in sgn_func]
    ))]
    if "cent" in cut_set:
        fit_configs = itertools.product(
            [*zip(cut_set["cent"]["min"], cut_set["cent"]["max"])],
            fit_configs
        )
        fit_configs = [(*cent, *config) for cent, config in fit_configs]
        config_keys = [
            "cent_min", "cent_max", "i_pt", "pt_min", "pt_max", "mass_min", "mass_max",
            "rebin", "signal_func", "fix_ds_sigma", "fix_dplus_sigma", "fix_sigma_dplus_to_ds",
            "bkg_func", "use_bkg_templ"
        ]
    else:
        config_keys = [
            "i_pt", "pt_min", "pt_max", "mass_min", "mass_max",
            "rebin", "signal_func", "fix_ds_sigma", "fix_dplus_sigma", "fix_sigma_dplus_to_ds",
            "bkg_func", "use_bkg_templ"
        ]
    config_keys += [
        *[f"{par}_init_{i_func}" for i_func, sgn_func in enumerate(cfg["fit_configs"]["signal"]["par_init_limit"]) for par in sgn_func],
        *[f"{par}_min_{i_func}" for i_func, sgn_func in enumerate(cfg["fit_configs"]["signal"]["par_init_limit"]) for par in sgn_func],
        *[f"{par}_max_{i_func}" for i_func, sgn_func in enumerate(cfg["fit_configs"]["signal"]["par_init_limit"]) for par in sgn_func],
        *[f"{par}_fix_{i_func}" for i_func, sgn_func in enumerate(cfg["fit_configs"]["signal"]["par_init_limit"]) for par in sgn_func]
    ]
    fit_configs = [dict(zip(config_keys, config)) for config in fit_configs]
    return fit_configs

def get_sigma_from_cfg(cfg, fit_config, particle_name):
    """
    Retrieve the sigma value for a given particle from the configuration.

    Parameters:
    - cfg (dict): Configuration dictionary containing input data and sigma settings.
    - fit_config (dict): Dictionary containing fit configuration parameters
        such as mass limits and rebinning.
    - particle_name (str): Name of the particle, either "ds" or "dplus".

    Returns:
    - float: The sigma value for the specified particle.

    Raises:
    - ValueError: If the particle_name is not "ds" or "dplus".
    """
    if particle_name not in ('ds', 'dplus'):
        raise ValueError(f"Invalid particle name: {particle_name}")
    
    idx_signal = 0 if particle_name == "ds" else 1
    cfg_sigma_fix = cfg["fit_configs"]["signal"][f"{particle_name}_sigma"]
    signal_func = fit_config["signal_func"][idx_signal]
    
    # If cfg_sigma_fix is False, fit MC to get sigma
    if not cfg_sigma_fix:
        return _fit_mc_for_sigma(cfg, fit_config, particle_name, idx_signal, signal_func)
    
    # If cfg_sigma_fix is a string (ROOT file path), read from file
    if isinstance(cfg_sigma_fix, str):
        return _read_sigma_from_root(cfg_sigma_fix, fit_config, particle_name, signal_func)
    
    # If cfg_sigma_fix is a list, return the appropriate pt bin value
    if isinstance(cfg_sigma_fix, list):
        return cfg_sigma_fix[fit_config["i_pt"]]
    
    # Otherwise, return the fixed value for all pt bins
    return cfg_sigma_fix


def _fit_mc_for_sigma(cfg, fit_config, particle_name, idx_signal, signal_func):
    """Fit MC data to extract sigma parameters."""
    pt_suffix = f"{fit_config['pt_min'] * 10:.0f}_{fit_config['pt_max'] * 10:.0f}"
    suffix = pt_suffix
    
    if "cent_min" in fit_config and "cent_max" in fit_config:
        suffix += f"_cent_{fit_config['cent_min']}_{fit_config['cent_max']}"
    
    data_hdl_mc = DataHandler(
        data=cfg["inputs"]["mc"][particle_name],
        histoname=f"h_mass_{suffix}",
        limits=[fit_config["mass_min"], fit_config["mass_max"]], 
        rebin=fit_config["rebin"]
    )
    
    fitter_mc = F2MassFitter(
        data_hdl_mc, name_signal_pdf=signal_func,
        name_background_pdf="nobkg",
        name=f"{particle_name}_pt{pt_suffix}_for_sigma", 
        chi2_loss=True, verbosity=1, tol=1.e-1
    )
    
    pdg_id = 431 if particle_name == "ds" else 411
    fitter_mc.set_particle_mass(0, pdg_id=pdg_id)
    
    _configure_signal_parameters(fitter_mc, signal_func)
    fitter_mc.mass_zfit()
    
    return _extract_fit_results(fitter_mc, signal_func)


def _configure_signal_parameters(fitter_mc, signal_func):
    """Configure initial parameters for different signal functions."""
    param_configs = {
        "doublegaus": [
            ("sigma1", 0.01, [0.001, 0.050]),
            ("sigma2", 0.02, [0.001, 0.050]),
            ("frac1", 0.02, [0., 1.])
        ],
        "gaussian": [
            ("sigma", 0.01, [0.001, 0.050])
        ],
        "doublecb": [
            ("sigma", 0.01, [0.001, 0.050]),
            ("alphal", 1.5, [1., 3.]),
            ("alphar", 1.5, [1., 3.]),
            ("nl", 50, [30., 100.]),
            ("nr", 50, [30., 100.])
        ],
        "doublecbsymm": [
            ("sigma", 0.01, [0.001, 0.050]),
            ("alpha", 1.5, [1., 3.]),
            ("n", 50, [30., 100.])
        ]
    }
    
    if signal_func not in param_configs:
        Logger("Signal function not supported", "FATAL")
    
    for param_name, init_val, limits in param_configs[signal_func]:
        fitter_mc.set_signal_initpar(0, param_name, init_val, limits=limits)


def _extract_fit_results(fitter_mc, signal_func):
    """Extract fit results based on signal function type."""
    if signal_func == "doublegaus":
        return (
            fitter_mc.get_signal_parameter(0, "sigma1")[0],
            fitter_mc.get_signal_parameter(0, "sigma2")[0],
            fitter_mc.get_signal_parameter(0, "frac1")[0]
        )
    elif signal_func == "gaussian":
        return fitter_mc.get_sigma(0)[0]
    elif signal_func == "doublecb":
        return (
            fitter_mc.get_signal_parameter(0, "sigma")[0],
            fitter_mc.get_signal_parameter(0, "alphar")[0],
            fitter_mc.get_signal_parameter(0, "alphal")[0],
            fitter_mc.get_signal_parameter(0, "nl")[0],
            fitter_mc.get_signal_parameter(0, "nr")[0]
        )
    elif signal_func == "doublecbsymm":
        return (
            fitter_mc.get_signal_parameter(0, "sigma")[0],
            fitter_mc.get_signal_parameter(0, "alpha")[0],
            fitter_mc.get_signal_parameter(0, "n")[0]
        )


def _read_sigma_from_root(file_path, fit_config, particle_name, signal_func):
    """Read sigma values from ROOT file."""
    with uproot.open(file_path) as f:
        has_centrality = "cent_min" in fit_config and "cent_max" in fit_config
        cent_suffix = f'_{fit_config["cent_min"]}_{fit_config["cent_max"]}' if has_centrality else ''
        
        histogram_names = {
            "doublegaus": [f'h_sigma1_{particle_name}{cent_suffix}', 
                          f'h_sigma2_{particle_name}{cent_suffix}', 
                          f'h_frac1_{particle_name}{cent_suffix}'],
            "gaussian": [f'h_sigma_{particle_name}{cent_suffix}'],
            "doublecb": [f'h_sigma_{particle_name}{cent_suffix}', 
                        f'h_alphar_{particle_name}{cent_suffix}',
                        f'h_alphal_{particle_name}{cent_suffix}',
                        f'h_nl_{particle_name}{cent_suffix}',
                        f'h_nr_{particle_name}{cent_suffix}'],
            "doublecbsymm": [f'h_sigma_{particle_name}{cent_suffix}',
                             f'h_alpha_{particle_name}{cent_suffix}',
                             f'h_n_{particle_name}{cent_suffix}']
        }
        
        hist_names = histogram_names[signal_func]
        values = [f[name].values()[fit_config["i_pt"]] for name in hist_names]
        
        return values[0] if len(values) == 1 else tuple(values)


def get_sigma_dplus_to_ds(cfg, fit_config, fitter):
    """
    Retrieve the sigma of the D+ peak using that of the Ds peak.

    Parameters:
    - cfg (dict): Configuration dictionary containing the ratio of sigma D+ to Ds.
    - fit_config (dict): Configuration dictionary for the fit.
    - fitter (flarefly.F2MassFitter): Fitter object used to perform the mass fit
        and extract sigma values.

    Returns:
    - float: The sigma value for the D+ peak.
    """
    # get sigma of ds peak
    fitter.mass_zfit()
    sigma = fitter.get_sigma(0)[0]
    cfg_sigma_ratio = cfg["fit_configs"]["signal"]["ratio_sigma_dplus_to_ds"]

    if not cfg_sigma_ratio:
        return sigma

    if isinstance(cfg_sigma_ratio, list):
        return sigma * cfg_sigma_ratio[fit_config["i_pt"]]

    if isinstance(cfg_sigma_ratio, str):
        with uproot.open(cfg_sigma_ratio) as f:
            h_sigma_ratio = f["h_sigma_ratio"]
            return sigma * h_sigma_ratio.values()[fit_config["i_pt"]]

    return sigma * cfg_sigma_ratio


def initialise_signal(fitter, fit_config, idx):
    """
    Initialise the signal parameters for the fitter based on the provided configuration.

    Parameters:
    - fitter (object): The fitter object that will be used to set the signal parameters.
    - fit_config (dict): Configuration for the signal function and parameters.
    - idx (int): Index of the signal component (0 for Ds, 1 for D+).
    """

    def set_param(param):
        fitter.set_signal_initpar(
            idx,
            param,
            fit_config.get(f"{param}_init_{idx}", defaults[param]["init"]),
            limits=[
                fit_config.get(f"{param}_min_{idx}", defaults[param]["min"]),
                fit_config.get(f"{param}_max_{idx}", defaults[param]["max"]),
            ],
            fix=fit_config.get(f"{param}_fix_{idx}", False)
        )

    func_type = fit_config["signal_func"][idx]

    # Always set 'mu' if present
    if f"mu_init_{idx}" in fit_config:
        set_param("mu")

    # Define parameter defaults per function type
    PARAMS_BY_FUNC = {
        "gaussian": ["sigma"],
        "doublegaus": ["sigma1", "sigma2", "frac1"],
        "doublecb": ["sigma", "alphar", "alphal", "nl", "nr"],
        "doublecbsymm": ["sigma", "alpha", "n"],
        "genergausexptailsymm": ["sigma", "alpha"]
    }

    defaults = {
        "mu":      {"init": 1.9, "min": 0., "max": 10.},
        "sigma":   {"init": 0.01, "min": 0.001, "max": 0.03},
        "sigma1":  {"init": 0.01, "min": 0.001, "max": 0.03},
        "sigma2":  {"init": 0.01, "min": 0.001, "max": 0.03},
        "frac1":   {"init": 0.01, "min": 0.0, "max": 1.0},
        "alpha":   {"init": 5.0, "min": 0.5, "max": 10.0},
        "alphar":  {"init": 0.5, "min": 0.0, "max": 10.0},
        "alphal":  {"init": 0.5, "min": 0.0, "max": 10.0},
        "n":       {"init": 10.0, "min": 5.0, "max": 100.0},
        "nl":      {"init": 1.0, "min": 0.0, "max": 10.0},
        "nr":      {"init": 1.0, "min": 0.0, "max": 10.0},
    }

    for param in PARAMS_BY_FUNC.get(func_type, []):
        set_param(param)

    # Always set the global signal frac
    fitter.set_signal_initpar(idx, "frac", 0.1, limits=[0.0002, 1.0])


def fix_signal_parameters(fitter, fit_config, idx):
    particle = "ds" if idx == 0 else "dplus"
    if fit_config["signal_func"][idx] == "gaussian":
        sigma_ds = get_sigma_from_cfg(cfg, fit_config, particle)
        fitter.set_signal_initpar(idx, "sigma", sigma_ds, fix=True)
    elif fit_config["signal_func"][idx] == "doublegaus":
        sigma_ds = get_sigma_from_cfg(cfg, fit_config, particle)
        fitter.set_signal_initpar(idx, "sigma1", sigma_ds[0], fix=True)
        fitter.set_signal_initpar(idx, "sigma2", sigma_ds[1], fix=True)
        fitter.set_signal_initpar(idx, "frac1", sigma_ds[2], fix=True)
    elif fit_config["signal_func"][idx] == "doublecb":
        sigma_ds = get_sigma_from_cfg(cfg, fit_config, particle)
        fitter.set_signal_initpar(idx, "sigma", sigma_ds, fix=True)
        fitter.set_signal_initpar(idx, "alphal", 1.5, limits=[1., 3.])
        fitter.set_signal_initpar(idx, "alphar", 1.5, limits=[1., 3.])
        fitter.set_signal_initpar(idx, "nl", 50, limits=[30., 100.])
        fitter.set_signal_initpar(idx, "nr", 50, limits=[30., 100.])

def do_fit(fit_config, cfg):  # pylint: disable=too-many-locals, too-many-branches, too-many-statements # noqa: E501
    """Fit the invariant mass spectrum for a given configuration."""
    pt_min = fit_config["pt_min"]
    pt_max = fit_config["pt_max"]
    pt_cent_suffix = f"{pt_min*10:.0f}_{pt_max*10:.0f}"

    cent_min, cent_max = None, None
    if "cent_min" in fit_config and "cent_max" in fit_config:
        cent_min = fit_config["cent_min"]
        cent_max = fit_config["cent_max"]
        pt_cent_suffix += f"_cent_{cent_min:.0f}_{cent_max:.0f}"

    # Create data handler based on centrality
    if cent_min is not None and cent_max is not None:
        histoname = f'h_mass_{pt_min*10:.0f}_{pt_max*10:.0f}_cent_{cent_min:.0f}_{cent_max:.0f}'
    else:
        histoname = f'h_mass_{pt_min*10:.0f}_{pt_max*10:.0f}'
    
    data_hdl = DataHandler(
        data=cfg["inputs"]["data"],
        histoname=histoname,
        limits=[fit_config["mass_min"], fit_config["mass_max"]], 
        rebin=fit_config["rebin"]
    )

    bkg_funcs = fit_config["bkg_func"]
    label_signal_pdfs = [r"$\mathrm{D_{s}^{+}}$ signal", r"$\mathrm{D^{+}}$ signal"]
    label_bkg_pdfs = ["Combinatorial background"]
    
    # Handle background templates if needed
    if fit_config["use_bkg_templ"]:
        fitter = _setup_fitter_with_templates(
            cfg, fit_config, data_hdl, bkg_funcs, label_signal_pdfs, 
            label_bkg_pdfs, pt_cent_suffix, pt_min, pt_max, cent_min, cent_max
        )
    else:
        fitter = F2MassFitter(
            data_hdl, name_signal_pdf=fit_config["signal_func"],
            name_background_pdf=fit_config["bkg_func"],
            name=f"ds_pt_{pt_cent_suffix}", chi2_loss=True,
            label_signal_pdf=label_signal_pdfs,
            label_bkg_pdf=label_bkg_pdfs,
            verbosity=1, tol=1.e-1
        )

    # signals initialisation
    fitter.set_particle_mass(0, pdg_id=431)
    initialise_signal(fitter, fit_config, 0)

    if len(fit_config["signal_func"]) > 1:
        fitter.set_particle_mass(1, pdg_id=411)
        initialise_signal(fitter, fit_config, 1)

    # bkg initialisation
    _initialize_background_functions(fitter, bkg_funcs)

    # Handle fixed sigma configurations
    _apply_fixed_sigma_configs(fitter, cfg, fit_config)

    n_signal = len(fit_config["signal_func"])
    fit_result = fitter.mass_zfit()
    
    if cfg["output"]["save_all_fits"]:
        _save_fit_outputs(cfg, fitter, fit_config, bkg_funcs, pt_min, pt_max, cent_min, cent_max)
        
        fracs = fitter._F2MassFitter__get_all_fracs() # pylint: disable=protected-access
        corr_bkg_frac_dict, corr_bkg_frac_over_dplus_dict = _extract_background_fractions(fracs)
        
        out_dict = {
            "raw_yields": [fitter.get_raw_yield(i) for i in range(n_signal)],
            "mean": [fitter.get_mass(i) for i in range(n_signal)],
            "chi2": float(fitter.get_chi2_ndf()),
            "significance": [fitter.get_significance(i, min=1.8, max=2.2) for i in range(n_signal)],
            "signal": [fitter.get_signal(i, min=1.8, max=2.) for i in range(n_signal)],
            "background": [fitter.get_background(i, min=1.8, max=2.) for i in range(n_signal)],
            **corr_bkg_frac_dict,
            "fracs": fracs,
            "converged": fit_result.converged, 
            **corr_bkg_frac_over_dplus_dict
        }

        # Add signal function specific parameters
        _add_signal_parameters(out_dict, fitter, fit_config, n_signal)
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
    return out_dict


def _setup_fitter_with_templates(cfg, fit_config, data_hdl, bkg_funcs, label_signal_pdfs, 
                                 label_bkg_pdfs, pt_cent_suffix, pt_min, pt_max, cent_min, cent_max):
    """Setup fitter with background templates."""
    bkg_cfg = cfg["fit_configs"]["bkg"]["templ_norm"]["backgrounds"]
    data_corr_bkgs = []
    
    for i_bkg, bkg in enumerate(bkg_cfg):
        bkg_funcs.insert(i_bkg, "hist")
        label_bkg_pdfs.insert(i_bkg, bkg["name"] + "\ncorrelated background")
        data_corr_bkgs.append(DataHandler(
            data=bkg["template_file"],
            histoname=bkg["template_hist_name"].format(
                pt_min=f"{pt_min*10:.0f}", pt_max=f"{pt_max*10:.0f}",
                cent_min=cent_min, cent_max=cent_max
            ),
            limits=[fit_config["mass_min"], fit_config["mass_max"]], 
            rebin=fit_config["rebin"]
        ))

    fitter = F2MassFitter(
        data_hdl, name_signal_pdf=fit_config["signal_func"],
        name_background_pdf=bkg_funcs,
        name=f"ds_pt_{pt_cent_suffix}", chi2_loss=True,
        label_signal_pdf=label_signal_pdfs,
        label_bkg_pdf=label_bkg_pdfs,
        verbosity=1, tol=1.e-1
    )

    templ_norm_cfg = cfg["fit_configs"]["bkg"]["templ_norm"]
    if templ_norm_cfg["fix_with_br"][fit_config['i_pt']]:
        data_signal_reference = DataHandler(
            data=templ_norm_cfg["signal"]["file_norm"],
            histoname=templ_norm_cfg["signal"]["hist_name"].format(
                pt_min=f"{pt_min*10:.0f}", pt_max=f"{pt_max*10:.0f}",
                cent_min=cent_min, cent_max=cent_max
            ),
            limits=[fit_config["mass_min"], fit_config["mass_max"]], 
            rebin=fit_config["rebin"]
        )

    for i_func, data_corr_bkg in enumerate(data_corr_bkgs):
        fitter.set_background_template(i_func, data_corr_bkg)
        fitter.set_background_initpar(i_func, "frac", 0.01, limits=[0., 1.])
        
        if templ_norm_cfg["fix_with_br"][fit_config['i_pt']]:
            data_corr_bkg_reference = DataHandler(
                data=bkg_cfg[i_func]["file_norm"],
                histoname=bkg_cfg[i_func]["norm_hist_name"].format(
                    pt_min=f"{pt_min*10:.0f}", pt_max=f"{pt_max*10:.0f}",
                    cent_min=cent_min, cent_max=cent_max
                ),
                limits=[fit_config["mass_min"], fit_config["mass_max"]], 
                rebin=fit_config["rebin"]
            )
            norm_ratio = (data_corr_bkg_reference.get_norm() * 
                         bkg_cfg[i_func]["br"]["pdg"] / bkg_cfg[i_func]["br"]["simulations"]) / \
                        (data_signal_reference.get_norm() * 
                         templ_norm_cfg["signal"]["br"]["pdg"] / templ_norm_cfg["signal"]["br"]["simulations"])
            fitter.fix_bkg_frac_to_signal_pdf(i_func, 1, norm_ratio)
        elif f"corr_bkg_frac_over_dplus_{i_func}_cfg" in fit_config:
            fitter.fix_bkg_frac_to_signal_pdf(
                i_func, 1, fit_config[f"corr_bkg_frac_over_dplus_{i_func}_cfg"]
            )
    
    return fitter


def _initialize_background_functions(fitter, bkg_funcs):
    """Initialize background function parameters."""
    for i_func, bkg_func in enumerate(bkg_funcs):
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


def _apply_fixed_sigma_configs(fitter, cfg, fit_config):
    """Apply fixed sigma configurations for both Ds and D+ signals."""
    if fit_config["fix_ds_sigma"]:
        _fix_signal_sigma(fitter, cfg, fit_config, 0, "ds")
    
    if fit_config["fix_dplus_sigma"]:
        _fix_signal_sigma(fitter, cfg, fit_config, 1, "dplus")
    
    if fit_config["fix_sigma_dplus_to_ds"]:
        sigma_dplus_ratio = get_sigma_dplus_to_ds(cfg, fit_config, fitter)
        fitter.set_signal_initpar(1, "sigma", sigma_dplus_ratio, 
                                 limits=[sigma_dplus_ratio-0.001, sigma_dplus_ratio+0.001], fix=True)


def _fix_signal_sigma(fitter, cfg, fit_config, signal_idx, particle_type):
    """Fix signal parameters for a specific signal and particle type."""
    func_type = fit_config["signal_func"][signal_idx]
    sigma_values = get_sigma_from_cfg(cfg, fit_config, particle_type)
    
    if func_type == "gaussian":
        fitter.set_signal_initpar(signal_idx, "sigma", sigma_values, fix=True)
    elif func_type == "doublegaus":
        fitter.set_signal_initpar(signal_idx, "sigma1", sigma_values[0], fix=True)
        fitter.set_signal_initpar(signal_idx, "sigma2", sigma_values[1], fix=True)
        fitter.set_signal_initpar(signal_idx, "frac1", sigma_values[2], fix=True)
    elif func_type == "doublecbsymm":
        fitter.set_signal_initpar(signal_idx, "sigma", sigma_values[0], fix=True)
        fitter.set_signal_initpar(signal_idx, "alpha", sigma_values[1], fix=True)
        fitter.set_signal_initpar(signal_idx, "n", sigma_values[2], fix=True)


def _extract_background_fractions(fracs):
    """Extract background fraction dictionaries from fit fractions."""
    corr_bkg_frac_dict = {}
    corr_bkg_frac_over_dplus_dict = {}
    
    if len(fracs[1]) > 0:
        for i_corr_bkg, (corr_bkg_frac, corr_bkg_err) in enumerate(zip(fracs[1][:-1], fracs[4][:-1])):
            corr_bkg_frac_dict[f"corr_bkg_frac_{i_corr_bkg}"] = [corr_bkg_frac, corr_bkg_err]
            ratio = corr_bkg_frac / fracs[0][1]
            ratio_err = ratio * np.sqrt((corr_bkg_err / corr_bkg_frac)**2 + (fracs[3][1] / fracs[0][1])**2)
            corr_bkg_frac_over_dplus_dict[f"corr_bkg_frac_over_dplus_{i_corr_bkg}"] = [ratio, ratio_err]
    else:
        corr_bkg_frac_dict["corr_bkg_frac_0"] = [0, 0]
        corr_bkg_frac_over_dplus_dict["corr_bkg_frac_over_dplus_0"] = [0, 0]
    
    return corr_bkg_frac_dict, corr_bkg_frac_over_dplus_dict


def _add_signal_parameters(out_dict, fitter, fit_config, n_signal):
    """Add signal function specific parameters to output dictionary."""
    func_type = fit_config["signal_func"][0]
    
    if func_type == "doublegaus":  #TODO: generalise in case different signal functions are used
        out_dict["sigma1"] = [fitter.get_signal_parameter(i, "sigma1") for i in range(n_signal)]
        out_dict["sigma2"] = [fitter.get_signal_parameter(i, "sigma2") for i in range(n_signal)]
        out_dict["frac1"] = [fitter.get_signal_parameter(i, "frac1") for i in range(n_signal)]
    elif func_type == "gaussian":
        out_dict["sigma"] = [fitter.get_sigma(i) for i in range(n_signal)]
    elif func_type == "doublecb":
        out_dict["sigma"] = [fitter.get_signal_parameter(i, "sigma") for i in range(n_signal)]
        out_dict["alphar"] = [fitter.get_signal_parameter(i, "alphar") for i in range(n_signal)]
        out_dict["alphal"] = [fitter.get_signal_parameter(i, "alphal") for i in range(n_signal)]
        out_dict["nl"] = [fitter.get_signal_parameter(i, "nl") for i in range(n_signal)]
        out_dict["nr"] = [fitter.get_signal_parameter(i, "nr") for i in range(n_signal)]
    elif func_type == "doublecbsymm":
        out_dict["sigma"] = [fitter.get_signal_parameter(i, "sigma") for i in range(n_signal)]
        out_dict["alpha"] = [fitter.get_signal_parameter(i, "alpha") for i in range(n_signal)]
        out_dict["n"] = [fitter.get_signal_parameter(i, "n") for i in range(n_signal)]
    elif func_type == "genergausexptailsymm":
        out_dict["sigma"] = [fitter.get_signal_parameter(i, "sigma") for i in range(n_signal)]
        out_dict["alpha"] = [fitter.get_signal_parameter(i, "alpha") for i in range(n_signal)]


def _save_fit_outputs(cfg, fitter, fit_config, bkg_funcs, pt_min, pt_max, cent_min, cent_max):
    """Save fit plots and ROOT outputs."""
    output_dir = os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    loc = ["lower left", "upper left"]
    ax_title = r"$M(\mathrm{KK\pi})$ GeV$/c^2$"
    fig, ax = fitter.plot_mass_fit(
        style="ATLAS",
        show_extra_info=(bkg_funcs != ["nobkg"]),
        figsize=(8, 8), extra_info_loc=loc,
        legend_loc="upper right",
        axis_title=ax_title
    )
    add_info_on_canvas(ax, "center right", pt_min, pt_max, cent_min, cent_max)
    figres = fitter.plot_raw_residuals(
        figsize=(8, 8), style="ATLAS",
        extra_info_loc=loc, axis_title=ax_title
    )
    
    for frmt in cfg["output"]["formats"]:
        if cent_min is not None and cent_max is not None:
            suffix = f"_{pt_min * 10:.0f}_{pt_max * 10:.0f}_cent_{cent_min:.0f}_{cent_max:.0f}_"
        else:
            suffix = f"_{pt_min * 10:.0f}_{pt_max * 10:.0f}_"
        suffix += cfg["output"]["suffix"]
        
        if frmt == "root":
            fitter.dump_to_root(
                f"{output_dir}/fits_{cfg['output']['suffix']}_{suffix}_partial.root", 
                option="recreate", suffix=suffix, num=5000
            )
        else:
            fig.savefig(f"{output_dir}/ds_mass_pt{suffix}.{frmt}")
            figres.savefig(f"{output_dir}/ds_massres_pt{suffix}.{frmt}")
            plt.close(fig)
            plt.close(figres)


def _apply_signal_constraints(fit_config, mb_result):
    """Apply signal parameter constraints from MB fit results."""
    for signal_idx in [0, 1]:
        func_type = fit_config["signal_func"][signal_idx]
        
        if func_type == 'gaussian':
            _fix_parameter(fit_config, "sigma", mb_result["sigma"][signal_idx][0], signal_idx)
            
        elif func_type == 'doublegaus':
            _fix_parameter(fit_config, "sigma1", mb_result["sigma1"][signal_idx][0], signal_idx)
            _fix_parameter(fit_config, "sigma2", mb_result["sigma2"][signal_idx][0], signal_idx)
            _fix_parameter(fit_config, "frac1", mb_result["frac1"][signal_idx][0], signal_idx)
            
        elif func_type == 'doublecbsymm':
            _fix_parameter(fit_config, "sigma", mb_result["sigma"][signal_idx][0], signal_idx)
            _fix_parameter(fit_config, "alpha", mb_result["alpha"][signal_idx][0], signal_idx)
            _fix_parameter(fit_config, "n", mb_result["n"][signal_idx][0], signal_idx)


def _fix_parameter(fit_config, param_name, value, signal_idx, tolerance=0.001):
    """Fix a parameter to a specific value with small tolerance."""
    suffix = f"_{signal_idx}"
    fit_config[f"{param_name}_fix{suffix}"] = True
    fit_config[f"{param_name}_init{suffix}"] = value
    fit_config[f"{param_name}_min{suffix}"] = value - tolerance
    fit_config[f"{param_name}_max{suffix}"] = value + tolerance


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

    zfit.run.set_cpus_explicit(
        intra=cfg['zfit_cpus']['intra'],
        inter=cfg['zfit_cpus']['inter']
    )

    # load inputs
    _, h_ev = load_inputs(cfg, cut_set)
    n_ev = sum(h_ev.values())

    cent_mins, cent_maxs = None, None
    if "cent" in cut_set:
        cent_mins = cut_set["cent"]["min"]
        cent_maxs = cut_set["cent"]["max"]

    # create output handler
    h_handler = HistHandler(
        cut_set["pt"]["min"], cut_set["pt"]["max"],
        cent_mins, cent_maxs
    )
    h_handler.set_n_ev(n_ev)

    fit_configs = create_fit_configs(cfg, cut_set)

    # recreate root file if needed
    if "root" in cfg["output"]["formats"]:
        output_dir = os.path.join(
            os.path.expanduser(cfg["output"]["directory"]),
            "fits"
        )
        uproot.recreate(f"{output_dir}/fits_{cfg['output']['suffix']}.root")

    bkg_cfg = cfg["fit_configs"]["bkg"]
    if (any(bkg_cfg["use_bkg_templ"]) and any(bkg_cfg["templ_norm"]["fix_to_mb"])) or\
            any(cfg["fit_configs"]["signal"]["fix_sigma_to_mb"]):
        
        # Execute MB fits first
        results_fracs = []
        with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
            for fit_config in fit_configs:
                if "cent" in cut_set and fit_config["cent_min"] == 0 and fit_config["cent_max"] == 100:
                    results_fracs.append((executor.submit(do_fit, fit_config, cfg), fit_config))
        
        # Apply MB constraints to other fits
        for result, fit_cfg_result in results_fracs:
            res = result.result()
            sig_fracs, bkg_fracs, _, _, _, _ = res["fracs"]
            
            for fit_config in fit_configs:
                if (fit_config["pt_min"] == fit_cfg_result["pt_min"] and
                    fit_config["pt_max"] == fit_cfg_result["pt_max"] and
                    not (fit_config["cent_max"] == 100 and fit_config["cent_min"] == 0)):
                    
                    # Apply background constraints
                    if bkg_cfg["templ_norm"]["fix_to_mb"][fit_config["i_pt"]]:
                        for i_corr_bkg, corr_bkg_frac in enumerate(bkg_fracs[:-1]):
                            fit_config[f"corr_bkg_frac_over_dplus_{i_corr_bkg}_cfg"] = corr_bkg_frac / sig_fracs[1]
                    
                    # Apply signal constraints
                    if cfg["fit_configs"]["signal"]["fix_sigma_to_mb"][fit_config["i_pt"]]:
                        _apply_signal_constraints(fit_config, res)

    # Handle external background fraction constraints
    if bkg_cfg["templ_norm"]["fix_to_file_name"] is not None:
        with uproot.open(bkg_cfg["templ_norm"]["fix_to_file_name"]) as f:
            h_frac1 = f[bkg_cfg["templ_norm"]["hist_name"].format(
                cent_min=0, cent_max=100 # TODO: generalise
            )]
            for fit_config in fit_configs:
                fit_config["corr_bkg_frac_over_dplus_0"] = h_frac1.values()[fit_config["i_pt"]]

    # Execute all fits
    results = []
    with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
        for fit_config in fit_configs:
            results.append((executor.submit(do_fit, fit_config, cfg), fit_config))

    if "root" in cfg["output"]["formats"]:
        # get all the partial root files in the output directory
        output_dir = os.path.join(
            os.path.expanduser(cfg["output"]["directory"]),
            "fits"
        )
        partial_files = [f for f in os.listdir(output_dir) if f.endswith('_partial.root') and f.startswith('fits_')]
        output_dir = os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits")
        out_root_file = f"{output_dir}/fits_{cfg['output']['suffix']}.root"
        print("hadd -f " + out_root_file + " " + " ".join([os.path.join(output_dir, f) for f in partial_files]))
        os.system("hadd -f " + out_root_file + " " + " ".join([os.path.join(output_dir, f) for f in partial_files]))
        os.system("rm " + " ".join([os.path.join(output_dir, f) for f in partial_files]))

    # Save results
    output_df = []
    for result, fit_config in results:
        out_dict = result.result()
        out_dict.update(fit_config)
        output_df.append(out_dict)
    output_df = pd.DataFrame(output_df)
    
    output_df.to_parquet(os.path.join(
        os.path.expanduser(cfg["output"]["directory"]),
        f'fit_results{cfg["output"]["suffix"]}.parquet'
    ), index=False)

    h_handler.set_histos(output_df)
    h_handler.dump_to_root(os.path.join(
        os.path.expanduser(cfg["output"]["directory"]),
        f'mass_fits{cfg["output"]["suffix"]}.root'
    ))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('config_file_name', metavar='text', default='')
    args = parser.parse_args()

    fit(args.config_file_name)
