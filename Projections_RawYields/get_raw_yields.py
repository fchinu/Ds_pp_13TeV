"""Script to perform fitting of invariant mass distributions for Ds and D+ mesons"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import dataclasses
from typing import Dict, List, Tuple
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""  # pylint: disable=wrong-import-position
import sys
sys.path.append(os.path.abspath(os.path.join(__file__, '../../utils/fitter')))
from tqdm import tqdm
import uproot
import yaml
import pandas as pd
import numpy as np
import ROOT
from pypdf import PdfWriter
import tensorflow as tf
tf.config.threading.set_intra_op_parallelism_threads(20)
tf.config.threading.set_inter_op_parallelism_threads(20)
from fit_handler import (
    FitHandler,
    BRInfo,
    CorrelatedBackground,
    CorrelatedBackgroundConfig,
    FitConfig
)

# pylint: disable=no-member  # (ROOT dynamic members)

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
        self._h_collisions = None

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

    def set_h_collisions(self, h_collisions):
        """Set the histogram for collisions."""
        self._h_collisions = h_collisions

    def _get_centrality_index(self, row):
        """Get centrality index from row data."""
        if "cent_min_cfg" in row and "cent_max_cfg" in row:
            cent_tuple = (row["cent_min_cfg"], row["cent_max_cfg"])
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
            i_pt = self._pt_info.mins.index(row["pt_min_cfg"])
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
        with uproot.update(output_file) as f:
            f["h_coll_rebinned"] = self._h_collisions

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


def fitconfig_to_dict(cfg: FitConfig) -> dict:
    """Convert FitConfig dataclass to dictionary with modified keys."""
    base = dataclasses.asdict(cfg)
    return {f"{k}_cfg": v for k, v in base.items()}

def merge_pdfs(cfg):
    """Merge individual PDF files into a single PDF."""
    if "pdf" not in cfg["output"]["formats"]:
        return

    # if output already exists, delete it
    fits_out_path = os.path.join(
        os.path.expanduser(cfg["output"]["directory"]),
        "fits",
        "ds_mass_merged.pdf"
    )
    if os.path.exists(fits_out_path):
        os.remove(fits_out_path)

    residuals_out_path = os.path.join(
        os.path.expanduser(cfg["output"]["directory"]),
        "fits",
        "ds_massres_merged.pdf"
    )
    if os.path.exists(residuals_out_path):
        os.remove(residuals_out_path)

    output_dir = os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits")
    files = os.listdir(output_dir)
    pdf_files = [f for f in files if f.endswith('.pdf') and f.startswith('ds_mass_')]
    pdf_files = sorted(pdf_files, key=lambda x: (
        int(x.split('_')[3]),
        int(x.split('_')[4]) if 'cent' in x else -1
    ))

    pdf_files_residuals = [f for f in files if f.endswith('.pdf') and f.startswith('ds_massres_')]
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

def run_fit(fit_config: FitConfig) -> Dict:
    """
    Run the fitting procedure based on the provided FitConfig.

    Parameters:
    fit_config (FitConfig): Configuration for the fit.

    Returns:
    Dict: Results of the fitting procedure.
    """
    fit_handler = FitHandler(fit_config)
    results = fit_handler.get_results()
    return (
        fit_config,
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
    if not bkg_cfg["use_bkg_templ"][i_pt]:
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
        signal_br=BRInfo(
            pdg=bkg_norm["signal"]["br"]["pdg"],
            simulations=bkg_norm["signal"]["br"]["simulations"]
        )
    )

def get_parameter_setup(
        cfg: List[Dict],
        i_pt: int,
        mb_result: Dict = None,
        result_sigma_fix: Dict = None
    ) -> List[Dict]:
    """
    Get parameter setup based on the provided configuration.

    Parameters:
    cfg (list): List of configuration dictionaries.
    i_pt (int): Index of the pT bin.

    Returns:
    List[Dict]: Parameter setup dictionary.
    """
    sgn_cfg = cfg["fit_configs"]["signal"]
    param_cfg = sgn_cfg["par_init_limit"]
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
                    and sgn_cfg["fix_sigma_dplus_to_ds"][i_pt] \
                        and result_sigma_fix is not None:
                ratio_sigma_dplus_to_ds = -1.
                if isinstance(sgn_cfg["ratio_sigma_dplus_to_ds"], list):
                    ratio_sigma_dplus_to_ds = sgn_cfg["ratio_sigma_dplus_to_ds"][i_pt]
                elif isinstance(sgn_cfg["ratio_sigma_dplus_to_ds"], (int, float)):
                    ratio_sigma_dplus_to_ds = sgn_cfg["ratio_sigma_dplus_to_ds"]

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


def load_ev_collisions(cfg):
    """
    Load invariant event and collision histograms from a ROOT file.

    Parameters:
    cfg (dict): Configuration dictionary containing the path to the input data file.

    Returns:
    tuple: A tuple containing:
        - h_ev (uproot.models.TH1): Event histogram.
        - h_collisions (uproot.models.TH1): Collisions histogram.
    """

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

def fit(config_file_name):  # pylint: disable=too-many-locals, too-many-statements, too-many-branches
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

    # Create fit directory if it doesn't exist
    os.makedirs(os.path.join(os.path.expanduser(cfg["output"]["directory"]), "fits"), exist_ok=True)

    # zfit.run.set_cpus_explicit(
    #     intra=cfg['zfit_cpus']['intra'],
    #     inter=cfg['zfit_cpus']['inter']
    # )

    # load inputs
    h_ev, h_collisions = load_ev_collisions(cfg)
    n_ev = sum(h_ev.values())

    cent_mins, cent_maxs = None, None
    if "cent" in cut_set:
        cent_mins = cut_set["cent"]["min"]
        cent_maxs = cut_set["cent"]["max"]

    pt_mins = cut_set["pt"]["min"]
    pt_maxs = cut_set["pt"]["max"]

    # create output handler
    h_handler = HistHandler(
        pt_mins, pt_maxs,
        cent_mins, cent_maxs
    )
    h_handler.set_n_ev(n_ev)
    h_handler.set_h_collisions(h_collisions)

    # recreate root file if needed
    if "root" in cfg["output"]["formats"]:
        output_dir = os.path.join(
            os.path.expanduser(cfg["output"]["directory"]),
            "fits"
        )
        uproot.recreate(f"{output_dir}/fits_{cfg['output']['suffix']}.root")


    # if needed, run minimum-bias fit first
    mb_results = {}
    if _mb_fit_is_needed(cfg, cut_set):
        with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
            futures = []
            for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
                fit_config = get_config(cfg, (i_pt, pt_min, pt_max), (0, 100))
                futures.append(executor.submit(run_fit, fit_config))

            for future in tqdm(as_completed(futures), total=len(futures)):
                fit_cfg, results = future.result()
                pt_min = fit_cfg.pt_min
                pt_max = fit_cfg.pt_max
                cent_min = fit_cfg.cent_min
                cent_max = fit_cfg.cent_max
                mb_results[(pt_min, pt_max), (cent_min, cent_max)] = (fit_cfg, results)

    with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
        futures = []
        for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
            if cfg["fit_configs"]["signal"]["fix_sigma_dplus_to_ds"][i_pt]:
                fit_config = get_config(
                    cfg,
                    (i_pt, pt_min, pt_max),
                    (0, 100),
                    None,
                    mb_results[(pt_min, pt_max), (0, 100)][1]
                )
                futures.append(executor.submit(run_fit, fit_config))

        for future in tqdm(as_completed(futures), total=len(futures)):
            fit_cfg, results = future.result()
            pt_min = fit_cfg.pt_min
            pt_max = fit_cfg.pt_max
            cent_min = fit_cfg.cent_min
            cent_max = fit_cfg.cent_max
            mb_results[(pt_min, pt_max), (cent_min, cent_max)] = (fit_cfg, results)

    ## We get all the configurations for fits

    fit_configs = []
    for ((pt_min, pt_max), (_, _)), (fit_cfg, mb_result) in mb_results.items():
        for cent_min, cent_max in zip(cent_mins, cent_maxs):
            if cent_min == 0 and cent_max == 100:
                continue
            i_pt = pt_mins.index(pt_min)
            fit_config = get_config(cfg, (i_pt, pt_min, pt_max), (cent_min, cent_max), mb_result)
            fit_configs.append(fit_config)

    results = {}
    with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
        futures = []
        for fit_config in fit_configs:
            pt_min = fit_config.pt_min
            pt_max = fit_config.pt_max
            futures.append(executor.submit(
                run_fit, fit_config
            ))

        for future in tqdm(as_completed(futures), total=len(futures)):
            fit_cfg, res = future.result()
            pt_min = fit_cfg.pt_min
            pt_max = fit_cfg.pt_max
            cent_min = fit_cfg.cent_min
            cent_max = fit_cfg.cent_max
            results[(pt_min, pt_max), (cent_min, cent_max)] = (fit_cfg, res)

    with ProcessPoolExecutor(max_workers=cfg["max_workers"]) as executor:
        futures = []
        for ((pt_min, pt_max), (cent_min, cent_max)), (fit_cfg, mb_result) in results.items():
            i_pt = pt_mins.index(pt_min)
            needs_further_fit = False
            if cfg["fit_configs"]["signal"]["fix_sigma_dplus_to_ds"][i_pt]:
                for func_par_setup in cfg["fit_configs"]["signal"]["par_init_limit"]:
                    if "sigma" in func_par_setup and not func_par_setup["sigma"]["fix_to_mb"][i_pt]:
                        needs_further_fit = True
            if needs_further_fit:
                fit_config = get_config(
                    cfg,
                    (i_pt, pt_min, pt_max),
                    (cent_min, cent_max),
                    None,
                    results[(pt_min, pt_max), (cent_min, cent_max)][1]
                )
                futures.append(executor.submit(run_fit, fit_config))

            for future in tqdm(as_completed(futures), total=len(futures)):
                fit_cfg, res = future.result()
                pt_min = fit_cfg.pt_min
                pt_max = fit_cfg.pt_max
                cent_min = fit_cfg.cent_min
                cent_max = fit_cfg.cent_max
                results[(pt_min, pt_max), (cent_min, cent_max)] = (fit_cfg, res)

    merge_pdfs(cfg)

    # Save results
    output_df = []
    for (fit_config, result) in results.values():
        out_dict = result.copy()
        out_dict.update(fitconfig_to_dict(fit_config))
        output_df.append(out_dict)
    for (fit_config, result) in mb_results.values():
        out_dict = result.copy()
        out_dict.update(fitconfig_to_dict(fit_config))
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
