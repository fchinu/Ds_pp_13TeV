import os
# Use only one thread per worker
for _threads_var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_threads_var, "1")

from pathlib import Path
from functools import partial
import multiprocessing as mp
import ROOT
import numpy as np
import yaml
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.stats import norm, spearmanr
from tqdm import tqdm


ROOT.TH1.AddDirectory(False)

infile = Path(__file__).parent / "ds_over_dplus_ratio_vmult_ptint.root"
infile_oo = Path(__file__).parent / "ds_over_dp_ratio_ptint_oo.root"
# rhos = {
#     "rawyield": {
#         "pp": 0.3,
#         "pbpb": 0.3,
#         "pp_pbpb": 0.2,
#     },
#     "np_frac": {
#         "pp": 0.8,
#         "pbpb": 0.8,
#         "pp_pbpb": 0.5,
#     },
#     "bdt_eff": {
#         "pp": 0.8,
#         "pbpb": 0.8,
#         "pp_pbpb": 0.7,
#     },
#     "pt_shape": {
#         "pp": 0.8,
#         "pbpb": 0.8,
#         "pp_pbpb": 0.,
#     },
#     "extrapolation": {
#         "pp": 0.6,
#         "pbpb": 0.6,
#         "pp_pbpb": 0.3,
#     },
# }



# rhos = {
#     "rawyield": {
#         "pp": 0.3,
#         "pbpb": 0.3,
#         "oo": 0.3,
#         "pp_pbpb": 0.2,
#         "pp_pbpb_oo": 0.2,
#     },
#     "np_frac": {
#         "pp": 1.,
#         "pbpb": 1.,
#         "oo": 1.,
#         "pp_pbpb": 1.,
#         "pp_pbpb_oo": 1.,
#     },
#     "bdt_eff": {
#         "pp": 1.0,
#         "pbpb": 1.0,
#         "oo": 1.0,
#         "pp_pbpb": 0.4,
#         "pp_pbpb_oo": 0.4,
#     },
#     "pt_shape": {
#         "pp": 1.0,
#         "pbpb": 1.0,
#         "oo": 1.0,
#         "pp_pbpb": 0.,
#         "pp_pbpb_oo": 0.,
#     },
#     "extrapolation": {
#         "pp": 1.,
#         "pbpb": 1.,
#         "oo": 1.,
#         "pp_pbpb": 0.3,
#         "pp_pbpb_oo": 0.3,
#     },
# }

rhos = {
    "rawyield": {
        "pp": 0.3,
        "pbpb": 0.3,
        "oo": 0.3,
        "pp_pbpb": 0.2,
        "pp_pbpb_oo": 0.2,
    },
    "np_frac": {
        "pp": 0.8,
        "pbpb": 0.8,
        "oo": 0.8,
        "pp_pbpb": 0.8,
        "pp_pbpb_oo": 0.8,
    },
    "bdt_eff": {
        "pp": 0.8,
        "pbpb": 0.8,
        "oo": 0.8,
        "pp_pbpb": 0.4,
        "pp_pbpb_oo": 0.4,
    },
    "pt_shape": {
        "pp": 0.8,
        "pbpb": 0.8,
        "oo": 0.8,
        "pp_pbpb": 0.,
        "pp_pbpb_oo": 0.,
    },
    "extrapolation": {
        "pp": 0.6,
        "pbpb": 0.6,
        "oo": 0.6,
        "pp_pbpb": 0.3,
        "pp_pbpb_oo": 0.3,
    },
}

# rhos = {
#     "rawyield": {
#         "pp": 0,
#         "pbpb": 0,
#         "pp_pbpb": 0,
#     },
#     "np_frac": {
#         "pp": 0,
#         "pbpb": 0,
#         "pp_pbpb": 0
#     },
#     "bdt_eff": {
#         "pp": 0,
#         "pbpb": 0,
#         "pp_pbpb": 0,
#     },
#     "pt_shape": {
#         "pp": 0,
#         "pbpb": 0,
#         "pp_pbpb": 0
#     },
#     "extrapolation": {
#         "pp": 0,
#         "pbpb": 0,
#         "pp_pbpb": 0,
#     },
# }

rhos_confidence_intervals = {
    "rawyield": {
        "pp": [0., 0.5],
        "pbpb": [0.1, 0.5],
        "pp_pbpb": [0.1, 0.4],
    },
    "np_frac": {
        "pp": [0.6, 1.0],
        "pbpb": [0.6, 1.0],
        "pp_pbpb": [0., 1.0],
    },
    "bdt_eff": {
        "pp": [0.6, 1.0],
        "pbpb": [0.6, 1.0],
        "pp_pbpb": [0.4, 1.0],
    },
    "pt_shape": {
        "pp": [0.6, 1.0],
        "pbpb": [0.6, 1.0],
        "pp_pbpb": [0., 1.0],
    },
    "extrapolation": {
        "pp": [0.4, 1.0],
        "pbpb": [0.4, 1.0],
        "pp_pbpb": [0.1, 1.0],
    },
}

def get_systematics_database():
    syst_database_pp = Path(__file__).parent.parent.parent / "Systematics/database_systematics_pp.yml"
    with open(syst_database_pp, 'r') as f:
        syst_pp = yaml.safe_load(f)
    syst_database_pbpb = Path(__file__).parent.parent.parent / "Systematics/database_systematics_pbpb.yml"
    with open(syst_database_pbpb, 'r') as f:
        syst_pbpb = yaml.safe_load(f)
    syst_database_oo = Path(__file__).parent.parent.parent / "Systematics/database_systematics_oo.yml"
    with open(syst_database_oo, 'r') as f:
        syst_oo = yaml.safe_load(f)
    
    return syst_pp, syst_pbpb, syst_oo

def get_data():
    """
    Get the ratio and uncertainties from the ROOT file.

    Returns:
        tuple: A tuple containing two tuples, one for pp and one for PbPb. Each tuple contains three lists: x values, y values, and uncertainties.
    """
    with ROOT.TFile.Open(str(infile.resolve())) as f:
        dir = f.Get("single_ratios")
        dir.cd()
        g_stat_pp = dir.Get("graph_pp_stat")
        g_syst_pp = dir.Get("graph_pp_syst_nobr")
        g_stat_pbpb = dir.Get("graph_pbpb_stat")
        g_syst_pbpb = dir.Get("graph_pbpb_syst_nobr")

    x_pp, y_pp = [], []
    stat_pp, syst_pp = [], []
    for i in range(g_stat_pp.GetN()):
        x_pp.append(g_stat_pp.GetPointX(i))
        y_pp.append(g_stat_pp.GetPointY(i))
        stat_pp.append(g_stat_pp.GetErrorY(i))
        syst_pp.append(g_syst_pp.GetErrorY(i))

    x_pbpb, y_pbpb = [], []
    stat_pbpb, syst_pbpb = [], []
    for i in range(g_stat_pbpb.GetN()):
        x_pbpb.append(g_stat_pbpb.GetPointX(i))
        y_pbpb.append(g_stat_pbpb.GetPointY(i))
        stat_pbpb.append(g_stat_pbpb.GetErrorY(i))
        syst_pbpb.append(g_syst_pbpb.GetErrorY(i))

    return (x_pp, y_pp, stat_pp, syst_pp), (x_pbpb, y_pbpb, stat_pbpb, syst_pbpb)

def get_data_oo():
    """
    Get the ratio and uncertainties from the ROOT file.

    Returns:
        tuple: A tuple containing two tuples, one for pp and one for PbPb. Each tuple contains three lists: x values, y values, and uncertainties.
    """
    with ROOT.TFile.Open(str(infile_oo.resolve())) as f:
        g_stat_oo = f.Get("gstat_ds_over_dp_ptint")
        g_syst_oo = f.Get("gsyst_tot_nobr_ds_over_dp_ptint")

    x_oo, y_oo = [], []
    stat_oo, syst_oo = [], []
    for i in range(g_stat_oo.GetN()):
        x_oo.append(g_stat_oo.GetPointX(i))
        y_oo.append(g_stat_oo.GetPointY(i))
        stat_oo.append(g_stat_oo.GetErrorY(i))
        syst_oo.append(g_syst_oo.GetErrorY(i))

    return (x_oo, y_oo, stat_oo, syst_oo)

def get_systematics(syst_database, x_data, y_data):
    """
    Match the systematics from the database to the data points based on x values.

    Args:
        syst_database (dict): The systematics database.
        x_data (list): The x values of the data points.
        y_data (list): The y values of the data points.

    Returns:
        list: A list of systematic uncertainties corresponding to the data points.
    """
    syst = {
        "rawyield": [],
        "np_frac": [],
        "bdt_eff": [],
        "pt_shape": [],
        "extrapolation": [],
    }
    cents = [k.replace("cent_", "") for k in syst_database["systematics"].keys() if "cent" in k]
    cents = [(int(k.split("_")[0]), int(k.split("_")[1])) for k in cents]
    cents = sorted(cents, key=lambda x: (x[0], x[1]))
    if len(cents) != len(x_data):
        raise ValueError("The number of centrality bins in the systematics database does not match the number of data points.")
    for i, (cent_low, cent_high) in enumerate(cents):
        cent_key = f"cent_{cent_low}_{cent_high}"
        for syst_name in syst.keys():
            # We take the systematic uncertainty from the lowest pt bin
            if not isinstance(syst_database["systematics"][cent_key][syst_name], list):
                syst[syst_name].append(np.mean(syst_database["systematics"][cent_key][syst_name]) * y_data[i])
            elif isinstance(syst_database["systematics"][cent_key][syst_name][0], list):
                # TODO: We take the mean for now
                syst[syst_name].append(np.mean(syst_database["systematics"][cent_key][syst_name][0]) * y_data[i])
            else:
                syst[syst_name].append(syst_database["systematics"][cent_key][syst_name][0] * y_data[i])
    return syst

def build_covariance_matrix(stat_pp, syst_pp, stat_pbpb, syst_pbpb, rhos):
    cov_pp = np.zeros((len(stat_pp), len(stat_pp)))
    # add the statistical and systematic uncertainty as a fully uncorrelated component
    for i in range(len(stat_pp)):
        cov_pp[i][i] = stat_pp[i]**2
        for syst_name in syst_pp.keys():
            cov_pp[i][i] += syst_pp[syst_name][i]**2
    # add the systematic uncertainty as a fully correlated component
    for i in range(len(stat_pp)):
        for j in range(len(stat_pp)):
            if i != j:
                for syst_name, rho_pp in rhos.items():
                    cov_pp[i][j] += rho_pp["pp"] * syst_pp[syst_name][i] * syst_pp[syst_name][j]

    cov_pbpb = np.zeros((len(stat_pbpb), len(stat_pbpb)))
    # add the statistical and systematic uncertainty as a fully uncorrelated component
    for i in range(len(stat_pbpb)):
        cov_pbpb[i][i] = stat_pbpb[i]**2
        for syst_name in syst_pbpb.keys():
            cov_pbpb[i][i] += syst_pbpb[syst_name][i]**2
    # add the systematic uncertainty as a fully correlated component
    for i in range(len(stat_pbpb)):
        for j in range(len(stat_pbpb)):
            if i != j:
                for syst_name, rho_pbpb in rhos.items():
                    cov_pbpb[i][j] += rho_pbpb["pbpb"] * syst_pbpb[syst_name][i] * syst_pbpb[syst_name][j]

    cov_total = np.zeros((len(stat_pp) + len(stat_pbpb), len(stat_pp) + len(stat_pbpb)))
    cov_total[:len(stat_pp), :len(stat_pp)] = cov_pp
    cov_total[len(stat_pp):, len(stat_pp):] = cov_pbpb

    # add the correlation between pp and PbPb
    for i in range(len(stat_pp)):
        for j in range(len(stat_pbpb)):
            for syst_name, rho_pp_pbpb in rhos.items():
                cov_total[i][len(stat_pp) + j] += rho_pp_pbpb["pp_pbpb"] * syst_pp[syst_name][i] * syst_pbpb[syst_name][j]
            cov_total[len(stat_pp) + j][i] = cov_total[i][len(stat_pp) + j] # symmetric matrix

    return cov_total

def build_covariance_matrix_with_oo(stat_pp, syst_pp, stat_pbpb, syst_pbpb, stat_oo, syst_oo, rhos):
    cov_pp = np.zeros((len(stat_pp), len(stat_pp)))
    # add the statistical and systematic uncertainty as a fully uncorrelated component
    for i in range(len(stat_pp)):
        cov_pp[i][i] = stat_pp[i]**2
        for syst_name in syst_pp.keys():
            cov_pp[i][i] += syst_pp[syst_name][i]**2
    # add the systematic uncertainty as a fully correlated component
    for i in range(len(stat_pp)):
        for j in range(len(stat_pp)):
            if i != j:
                for syst_name, rho_pp in rhos.items():
                    cov_pp[i][j] += rho_pp["pp"] * syst_pp[syst_name][i] * syst_pp[syst_name][j]

    cov_pbpb = np.zeros((len(stat_pbpb), len(stat_pbpb)))
    # add the statistical and systematic uncertainty as a fully uncorrelated component
    for i in range(len(stat_pbpb)):
        cov_pbpb[i][i] = stat_pbpb[i]**2
        for syst_name in syst_pbpb.keys():
            cov_pbpb[i][i] += syst_pbpb[syst_name][i]**2
    # add the systematic uncertainty as a fully correlated component
    for i in range(len(stat_pbpb)):
        for j in range(len(stat_pbpb)):
            if i != j:
                for syst_name, rho_pbpb in rhos.items():
                    cov_pbpb[i][j] += rho_pbpb["pbpb"] * syst_pbpb[syst_name][i] * syst_pbpb[syst_name][j]

    cov_oo = np.zeros((len(stat_oo), len(stat_oo)))
    # add the statistical and systematic uncertainty as a fully uncorrelated component
    for i in range(len(stat_oo)):
        cov_oo[i][i] = stat_oo[i]**2
        for syst_name in syst_oo.keys():
            cov_oo[i][i] += syst_oo[syst_name][i]**2
    # add the systematic uncertainty as a fully correlated component
    for i in range(len(stat_oo)):
        for j in range(len(stat_oo)):
            if i != j:
                for syst_name, rho_oo in rhos.items():
                    cov_oo[i][j] += rho_oo["oo"] * syst_oo[syst_name][i] * syst_oo[syst_name][j]

    cov_total = np.zeros((len(stat_pp) + len(stat_pbpb) + len(stat_oo), len(stat_pp) + len(stat_pbpb) + len(stat_oo)))
    cov_total[:len(stat_pp), :len(stat_pp)] = cov_pp
    cov_total[len(stat_pp):len(stat_pp) + len(stat_pbpb), len(stat_pp):len(stat_pp) + len(stat_pbpb)] = cov_pbpb
    cov_total[len(stat_pp) + len(stat_pbpb):, len(stat_pp) + len(stat_pbpb):] = cov_oo

    # add the correlation between pp, PbPb, OO. The three cross blocks have different shapes
    # (n_pp x n_pbpb, n_pp x n_oo and n_pbpb x n_oo), so each one needs its own loop bounds
    off_pbpb = len(stat_pp)
    off_oo = len(stat_pp) + len(stat_pbpb)

    # pp-PbPb
    for i in range(len(stat_pp)):
        for j in range(len(stat_pbpb)):
            for syst_name, rho_pp_pbpb_oo in rhos.items():
                cov_total[i][off_pbpb + j] += rho_pp_pbpb_oo["pp_pbpb_oo"] * syst_pp[syst_name][i] * syst_pbpb[syst_name][j]
            cov_total[off_pbpb + j][i] = cov_total[i][off_pbpb + j] # symmetric matrix

    # pp-OO
    for i in range(len(stat_pp)):
        for j in range(len(stat_oo)):
            for syst_name, rho_pp_pbpb_oo in rhos.items():
                cov_total[i][off_oo + j] += rho_pp_pbpb_oo["pp_pbpb_oo"] * syst_pp[syst_name][i] * syst_oo[syst_name][j]
            cov_total[off_oo + j][i] = cov_total[i][off_oo + j] # symmetric matrix

    # PbPb-OO
    for i in range(len(stat_pbpb)):
        for j in range(len(stat_oo)):
            for syst_name, rho_pp_pbpb_oo in rhos.items():
                cov_total[off_pbpb + i][off_oo + j] += rho_pp_pbpb_oo["pp_pbpb_oo"] * syst_pbpb[syst_name][i] * syst_oo[syst_name][j]
            cov_total[off_oo + j][off_pbpb + i] = cov_total[off_pbpb + i][off_oo + j] # symmetric matrix

    return cov_total

def pol0_weights(cov):
    """
    Get the linear weights of the pol0 fit.

    Args:
        cov (np.ndarray): The covariance matrix of the data points.

    Returns:
        tuple: The weight vector w and the uncertainty on the fitted constant.
    """
    cov_inv = np.linalg.inv(cov)
    A = np.ones((cov.shape[0], 1))
    param_cov = np.linalg.inv(A.T @ cov_inv @ A)
    w = np.asarray(param_cov @ (A.T @ cov_inv)).flatten()

    return w, np.sqrt(param_cov[0][0])


def sample_multivariate(mean, cov, n_toys, rng):
    """
    Draw n_toys pseudo-measurements from a multivariate Gaussian with the given covariance.

    Uses a Cholesky decomposition when possible, falling back to an eigenvalue
    decomposition with clipping if the covariance is not positive definite (which can
    happen for inconsistent choices of the rho's).

    Args:
        mean (np.ndarray): The central values around which to fluctuate.
        cov (np.ndarray): The covariance matrix.
        n_toys (int): Number of toy experiments.
        rng (np.random.Generator): The random number generator.

    Returns:
        np.ndarray: Array of shape (len(mean), n_toys) with the pseudo-measurements stored
            as columns, so that the pol0 weights can be applied as w @ y_toys.
    """
    try:
        chol = scipy.linalg.cholesky(cov, lower=True)
    except scipy.linalg.LinAlgError:
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        n_negative = np.sum(eigenvalues < 0)
        if n_negative:
            print(f"  WARNING: covariance matrix is not positive definite "
                  f"({n_negative} negative eigenvalue(s), min = {eigenvalues.min():.2e}). "
                  f"Clipping to zero: check the consistency of the rho's.")
        chol = eigenvectors @ np.diag(np.sqrt(np.clip(eigenvalues, 0., None)))

    return np.asarray(mean)[:, None] + chol @ rng.standard_normal((n_toys, len(mean))).T


def fit_pol0(y, cov):
    """Fit a constant (pol0) to the data using the covariance matrix."""
    # Invert the covariance matrix (V^-1)
    cov_inv = np.linalg.inv(cov)

    # build the data vector
    data = np.array(y)
    # Define the Design Matrix 'A' for a pol0 model (y = c)
    # This is a column vector of ones with the same length as your data
    A = np.ones((len(data), 1))

    # Calculate the covariance of the parameters: (A^T * V^-1 * A)^-1
    # For pol0, this is a 1x1 matrix containing the variance of the constant 'c'
    param_cov = np.linalg.inv(A.T @ cov_inv @ A)

    # Calculate the best fit parameter values: Param_cov * (A^T * V^-1 * Y)
    # The @ symbol is the matrix multiplication operator in Python
    theta_hat = param_cov @ (A.T @ cov_inv @ data)

    # Extract the value and error
    c_val = theta_hat[0]          # The fitted constant value
    c_err = np.sqrt(param_cov[0][0]) # The uncertainty on the fit

    return c_val, c_err


def get_significance_from_fit_to_all_data(pp_data, pbpb_data, fit_name):
    """
    Fit a constant to the combined pp and PbPb data and calculate the significance of the fit result.

    Args:
        pp_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for pp data.
        pbpb_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for PbPb data.
    """
    x_pp, y_pp, stat_pp, syst_pp, syst_tot_pp = pp_data
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb, syst_tot_pbpb = pbpb_data

    cov_total = build_covariance_matrix(stat_pp, syst_pp, stat_pbpb, syst_pbpb, rhos)
    cov_inv = np.linalg.inv(cov_total)

    data = np.array(y_pp + y_pbpb)
    c_val, c_err = fit_pol0(data, cov_total)


    print(f"Fit Result (pol0): {c_val:.4f} +/- {c_err:.4f}")

    # Calculate the residuals vector (Data - Model)
    # Since the model is just the constant c_val
    residuals = data - c_val

    # Chi2 = Residuals^T * V^-1 * Residuals
    chi2 = residuals.T @ cov_inv @ residuals
    ndf = len(data) - 1 # Number of points - Number of parameters

    print(f"Chi2 / NDF: {chi2:.2f} / {ndf} = {chi2/ndf:.2f}")
    # p-value calculation
    p_value = ROOT.TMath.Prob(chi2, ndf)
    print(f"p-value: {p_value:.4f}")

    # get significance of the fit result from the inverse of the normal distribution (one-tailed test)
    significance = norm.ppf(1 - p_value)
    print(f"Significance: {significance:.2f} sigma")

    fig, ax = plt.subplots(figsize=(8, 6))
    # Plotting pp and PbPb separately for clarity
    plt.errorbar(x_pp, y_pp, yerr=stat_pp, fmt='o', label='pp Data (Stat)')
    plt.errorbar(x_pbpb, y_pbpb, yerr=stat_pbpb, fmt='s', label='PbPb Data (Stat)')
    # Draw the systematic uncertainties as shaded boxes
    # Define a function to draw systematic boxes to keep code clean
    def draw_syst_boxes(ax, x, y, y_err, x_width, color, label):
        for i in range(len(x)):
            # Create a rectangle: (lower_left_x, lower_left_y), width, height
            rect = plt.Rectangle(
                (x[i] - x[i]*x_width/2, y[i] - y_err[i]), 
                x[i]*x_width, 
                2 * y_err[i],
                linewidth=0, 
                facecolor=color, 
                alpha=0.3,
                label=label if i == 0 else "" # Only add label once for the legend
            )
            ax.add_patch(rect)
    box_width = 0.2
    draw_syst_boxes(ax, x_pp, y_pp, syst_tot_pp, box_width, 'blue', 'pp Syst')
    draw_syst_boxes(ax, x_pbpb, y_pbpb, syst_tot_pbpb, box_width, 'orange', 'PbPb Syst')

    # Plot the fit band
    # We plot a horizontal line at c_val
    total_x = x_pp + x_pbpb
    plt.axhline(c_val, color='r', linestyle='--', label=f'Combined Fit: {c_val:.2f}')
    plt.fill_between([min(total_x), max(total_x)], c_val - c_err, c_val + c_err, color='r', alpha=0.3)

    # add text with the fit result and chi2/ndf
    plt.text(0.05, 0.95, f'Fit Result: {c_val:.4f} +/- {c_err:.4f}\nChi2/NDF: {chi2:.2f}/{ndf} = {chi2/ndf:.2f}\np-value: {p_value:.4f}\nSignificance: {significance:.2f} sigma', transform=plt.gca().transAxes, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    plt.xlabel("Multiplicity")
    plt.ylabel("Ds/D+ Ratio")
    # set log scale for x axis
    plt.xscale('log')
    # set y limits
    plt.ylim(0, 0.9)
    plt.legend()
    plt.savefig((Path(__file__).parent / f"{fit_name}.png").resolve())

def get_significance_two_fits(pp_data, pbpb_data, fit_name):
    x_pp, y_pp, stat_pp, syst_pp, syst_tot_pp = pp_data
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb, syst_tot_pbpb = pbpb_data

    cov_total = build_covariance_matrix(stat_pp, syst_pp, stat_pbpb, syst_pbpb, rhos)
    cov_pp = cov_total[:len(x_pp), :len(x_pp)]
    cov_pbpb = cov_total[len(x_pp):, len(x_pp):]
    cov_pp_inv = np.linalg.inv(cov_pp)
    cov_pbpb_inv = np.linalg.inv(cov_pbpb)

    c_val_pp, c_err_pp = fit_pol0(y_pp, cov_pp)
    print(f"Fit Result pp (pol0): {c_val_pp:.4f} +/- {c_err_pp:.4f}")
    c_val_pbpb, c_err_pbpb = fit_pol0(y_pbpb, cov_pbpb)
    print(f"Fit Result PbPb (pol0): {c_val_pbpb:.4f} +/- {c_err_pbpb:.4f}")

    # Chi2 = Residuals^T * V^-1 * Residuals, with ndf = n_points - 1 (one free parameter)
    residuals_pp = np.array(y_pp) - c_val_pp
    chi2_pp = residuals_pp.T @ cov_pp_inv @ residuals_pp
    ndf_pp = len(y_pp) - 1
    print(f"Chi2 / NDF pp: {chi2_pp:.2f} / {ndf_pp} = {chi2_pp/ndf_pp:.2f}, "
          f"p-value: {ROOT.TMath.Prob(chi2_pp, ndf_pp):.4f}")

    residuals_pbpb = np.array(y_pbpb) - c_val_pbpb
    chi2_pbpb = residuals_pbpb.T @ cov_pbpb_inv @ residuals_pbpb
    ndf_pbpb = len(y_pbpb) - 1
    print(f"Chi2 / NDF PbPb: {chi2_pbpb:.2f} / {ndf_pbpb} = {chi2_pbpb/ndf_pbpb:.2f}, "
          f"p-value: {ROOT.TMath.Prob(chi2_pbpb, ndf_pbpb):.4f}")

    significance_uncorr = abs(c_val_pp - c_val_pbpb) / np.sqrt(c_err_pp**2 + c_err_pbpb**2)
    print(f"Significance (uncorr.): {significance_uncorr:.2f} sigma")

    fig, ax = plt.subplots(figsize=(8, 6))
    # Plotting pp and PbPb separately for clarity
    plt.errorbar(x_pp, y_pp, yerr=stat_pp, fmt='o', label='pp Data (Stat)')
    plt.errorbar(x_pbpb, y_pbpb, yerr=stat_pbpb, fmt='s', label='PbPb Data (Stat)')
    # Draw the systematic uncertainties as shaded boxes
    # Define a function to draw systematic boxes to keep code clean
    def draw_syst_boxes(ax, x, y, y_err, x_width, color, label):
        for i in range(len(x)):
            # Create a rectangle: (lower_left_x, lower_left_y), width, height
            rect = plt.Rectangle(
                (x[i] - x[i]*x_width/2, y[i] - y_err[i]), 
                x[i]*x_width, 
                2 * y_err[i],
                linewidth=0, 
                facecolor=color, 
                alpha=0.3,
                label=label if i == 0 else "" # Only add label once for the legend
            )
            ax.add_patch(rect)
    box_width = 0.2
    draw_syst_boxes(ax, x_pp, y_pp, syst_tot_pp, box_width, 'blue', 'pp Syst')
    draw_syst_boxes(ax, x_pbpb, y_pbpb, syst_tot_pbpb, box_width, 'orange', 'PbPb Syst')

    # Plot the fit band
    # We plot a horizontal line at c_val_pp and c_val_pbpb
    plt.axhline(c_val_pp, color='r', linestyle='--', label=f'pp fit: {c_val_pp:.2f}')
    plt.fill_between([min(x_pp), max(x_pp)], c_val_pp - c_err_pp, c_val_pp + c_err_pp, color='r', alpha=0.3)
    plt.axhline(c_val_pbpb, color='g', linestyle='--', label=f'PbPb fit: {c_val_pbpb:.2f}')
    plt.fill_between([min(x_pbpb), max(x_pbpb)], c_val_pbpb - c_err_pbpb, c_val_pbpb + c_err_pbpb, color='g', alpha=0.3)

    # add text with the fit result and chi2/ndf
    plt.text(0.05, 0.95, f'pp fit result: {c_val_pp:.4f} +/- {c_err_pp:.4f} (Chi2/NDF: {chi2_pp:.2f}/{ndf_pp} = {chi2_pp/ndf_pp:.2f})\nPbPb fit result: {c_val_pbpb:.4f} +/- {c_err_pbpb:.4f} (Chi2/NDF: {chi2_pbpb:.2f}/{ndf_pbpb} = {chi2_pbpb/ndf_pbpb:.2f})\nSignificance (uncorr.): {significance_uncorr:.2f} sigma', transform=plt.gca().transAxes, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    plt.xlabel("Multiplicity")
    plt.ylabel("Ds/D+ Ratio")
    # set log scale for x axis
    plt.xscale('log')
    # set y limits
    plt.ylim(0, 0.9)
    plt.legend()
    plt.savefig((Path(__file__).parent / f"{fit_name}.png").resolve())

def toy_mc_rho_confidence_intervals(pp_data, pbpb_data, n_toys=1000):
    """
    Perform a toy Monte Carlo study to get the significance assuming pbpb points to have same values as pp points.

    Args:
        pp_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for pp data.
        pbpb_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for PbPb data.
        n_toys (int): Number of toy experiments to perform.
    """

    x_pp, y_pp, stat_pp, syst_pp, syst_tot_pp = pp_data
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb, syst_tot_pbpb = pbpb_data

    significance_distribution = []
    for i in range(n_toys):
        # Sample the rhos from the confidence intervals
        rhos_sampled = {}
        for syst_name, syst_rhos in rhos_confidence_intervals.items():
            rhos_sampled[syst_name] = {}
            for key, (low, high) in syst_rhos.items():
                rhos_sampled[syst_name][key] = np.random.uniform(low, high)

        cov_total = build_covariance_matrix(stat_pp, syst_pp, stat_pbpb, syst_pbpb, rhos_sampled)
        cov_pp = cov_total[:len(x_pp), :len(x_pp)]
        cov_pbpb = cov_total[len(x_pp):, len(x_pp):]

        c_val_pp, c_err_pp = fit_pol0(y_pp, cov_pp)
        c_val_pbpb, c_err_pbpb = fit_pol0(y_pbpb, cov_pbpb)

        significance_uncorr = abs(c_val_pp - c_val_pbpb) / np.sqrt(c_err_pp**2 + c_err_pbpb**2)
        significance_distribution.append(significance_uncorr)

    # Draw the histogram of the significance distribution
    plt.figure(figsize=(8, 6))
    plt.hist(significance_distribution, bins=30, density=True, alpha=0.7, color='blue')
    plt.xlabel("Significance")
    plt.ylabel("Probability Density")
    plt.title(fr"Significance Distribution from {n_toys} Toy MC Experiments (with different $\rho$)")
    plt.grid()
    plt.savefig((Path(__file__).parent / f"significance_distribution_{n_toys}.png").resolve())
    plt.close()


def toy_mc_mean_pbpb_pp_fixed(pp_data, pbpb_data, label, n_toys=100000, seed=42):
    """
    Test H0 (true PbPb ratio = true pp ratio) with the pp reference assumed exact.

    Args:
        pp_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for pp data.
        pbpb_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for PbPb data.
        label (str): Label for the plot files.
        n_toys (int): Number of toy experiments to perform.
        seed (int): Seed of the random number generator, for reproducibility.

    Returns:
        tuple: toys, observed value, H0 reference, p-value and significance.
    """
    rng = np.random.default_rng(seed)
    x_pp, y_pp, stat_pp, syst_pp, _ = pp_data
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb, _ = pbpb_data
    n_pp, n_pbpb = len(x_pp), len(x_pbpb)

    cov_total = build_covariance_matrix(stat_pp, syst_pp, stat_pbpb, syst_pbpb, rhos)
    cov_pp = cov_total[:n_pp, :n_pp]
    cov_pbpb = cov_total[n_pp:, n_pp:]

    # These are always the same for all fits
    w_pp, err_pp = pol0_weights(cov_pp)
    w_pbpb, err_pbpb = pol0_weights(cov_pbpb)

    c_val_pp = w_pp @ np.array(y_pp)
    c_val_pbpb = w_pbpb @ np.array(y_pbpb)

    y_pbpb_toys = sample_multivariate(np.full(n_pbpb, c_val_pp),
                                      cov_pbpb, n_toys, rng)
    c_pbpb_toys = w_pbpb @ y_pbpb_toys

    n_extreme = np.sum(np.abs(c_pbpb_toys - c_val_pp) >= abs(c_val_pbpb - c_val_pp))
    p_value = n_extreme / len(c_pbpb_toys)
    significance = norm.ppf(1 - p_value/2)  # two-tailed test
    print(f"[pp fixed] p-value = {p_value:.2e} "
          f"-> significance = {significance:.2f} sigma")

    fig, ax = plt.subplots(figsize=(8, 6))
    counts, _, _ = ax.hist(c_pbpb_toys, bins=100, density=True, log=True, alpha=0.7, label=f'Pol0 fit to {n_toys} PbPb toys (using pp reference)')
    ax.axvline(c_val_pbpb, color='k', linestyle='--', lw=2, label=f'Observed: {c_val_pbpb:.4f}')
    ax.axvline(c_val_pp, color='g', linestyle=':', lw=2, label=f'pp: {c_val_pp:.4f}')
    ax.set_xlabel("pol0 coefficient")
    ax.set_ylabel("Probability Density")
    ax.set_title("pp fixed")
    # Headroom on the log scale so the legend and the p-value box do not cover the toys
    positive_counts = counts[counts > 0]
    ax.set_ylim(positive_counts.min() / 2., positive_counts.max() * 10)
    ax.text(0.03, 0.97,
            f'p-value = {p_value:.2e}\nsignificance = {significance:.2f}$\\sigma$',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig((Path(__file__).parent / f"{label}_{n_toys}.png").resolve())
    plt.close(fig)


def toy_mc_mean_pbpb_both_fluctuate(pp_data, pbpb_data, label, n_toys=100000, seed=42):
    """
    Test H0 (true PbPb ratio = true pp ratio) letting both systems fluctuate.

    Args:
        pp_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for pp data.
        pbpb_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for PbPb data.
        label (str): Label for the plot files.
        n_toys (int): Number of toy experiments to perform.
        seed (int): Seed of the random number generator, for reproducibility.

    Returns:
        tuple: toys, observed value, H0 reference, p-value and significance.
    """
    rng = np.random.default_rng(seed)
    x_pp, y_pp, stat_pp, syst_pp, _ = pp_data
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb, _ = pbpb_data
    n_pp, n_pbpb = len(x_pp), len(x_pbpb)

    cov_total = build_covariance_matrix(stat_pp, syst_pp, stat_pbpb, syst_pbpb, rhos)
    cov_pp = cov_total[:n_pp, :n_pp]
    cov_pbpb = cov_total[n_pp:, n_pp:]

    # These are always the same for all fits
    w_pp, err_pp = pol0_weights(cov_pp)
    w_pbpb, err_pbpb = pol0_weights(cov_pbpb)

    c_val_pp = w_pp @ np.array(y_pp)
    c_val_pbpb = w_pbpb @ np.array(y_pbpb)

    delta_obs = c_val_pbpb - c_val_pp
    print(f"Observed difference PbPb - pp: {delta_obs:+.4f}")

    # The common true value is irrelevant for delta
    y_toys = sample_multivariate(np.full(n_pp + n_pbpb, c_val_pp),
                                 cov_total, n_toys, rng)
    delta_toys = w_pbpb @ y_toys[n_pp:, :] - w_pp @ y_toys[:n_pp, :]

    n_extreme = np.sum(np.abs(delta_toys) >= abs(delta_obs))
    p_value = n_extreme / len(delta_toys)
    significance = norm.ppf(1 - p_value/2)  # two-tailed test
    print(f"[both fluctuating] p-value = {p_value:.2e} "
          f"-> significance = {significance:.2f} sigma")


    fig, ax = plt.subplots(figsize=(8, 6))
    counts, _, _ = ax.hist(delta_toys, bins=100, density=True, log=True, alpha=0.7, label=f'Pol0 fits to {n_toys} pp and PbPb toys (using pp reference)')
    ax.axvline(delta_obs, color='k', linestyle='--', lw=2, label=f'Observed: {delta_obs:.4f}')
    ax.axvline(0., color='g', linestyle=':', lw=2, label='$H_0$: 0')
    ax.set_xlabel("difference in pol0 coefficients")
    ax.set_ylabel("Probability Density")
    ax.set_title("pp also fluctuating")
    # Headroom on the log scale so the legend and the p-value box do not cover the toys
    positive_counts = counts[counts > 0]
    ax.set_ylim(positive_counts.min() / 2., positive_counts.max() * 10)
    ax.text(0.03, 0.97,
            f'p-value = {p_value:.2e}\nsignificance = {significance:.2f}$\\sigma$',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig((Path(__file__).parent / f"{label}_{n_toys}.png").resolve())
    plt.close(fig)

    return significance


def _rho_significance_worker(task, y_pp, y_pbpb, n_toys, chunk_size=1_000_000):
    """
    Significance of the pp/PbPb difference for one sampled set of rho's.

    Args:
        task (tuple): The covariance matrix and the seed for this iteration.
        y_pp (list): The pp central values.
        y_pbpb (list): The PbPb central values.
        n_toys (int): Number of toy experiments.
        chunk_size (int): Number of toys drawn per chunk.

    Returns:
        float: the significance in sigma.
    """
    cov_total, seed = task

    rng = np.random.default_rng(seed)
    n_pp, n_pbpb = len(y_pp), len(y_pbpb)
    w_pp, _ = pol0_weights(cov_total[:n_pp, :n_pp])
    w_pbpb, _ = pol0_weights(cov_total[n_pp:, n_pp:])

    c_val_pp = w_pp @ np.array(y_pp)
    delta_obs = w_pbpb @ np.array(y_pbpb) - c_val_pp

    # The common true value is irrelevant for delta
    mean = np.full(n_pp + n_pbpb, c_val_pp)
    n_extreme = 0
    for start in range(0, n_toys, chunk_size):
        y_toys = sample_multivariate(mean, cov_total, min(chunk_size, n_toys - start), rng)
        delta_toys = w_pbpb @ y_toys[n_pp:, :] - w_pp @ y_toys[:n_pp, :]
        n_extreme += np.sum(np.abs(delta_toys) >= abs(delta_obs))

    # With no extreme toy the p-value is only bounded from above: quote the bound
    p_value = max(n_extreme / n_toys, 1. / n_toys)

    return norm.ppf(1 - p_value/2)  # two-tailed test


def toy_mc_mean_pbpb_both_fluctuate_rho_confidence_intervals(
        pp_data, pbpb_data, label, n_rho_samples=10000, n_toys=50000000, seed=42,
        n_workers=None):
    """
    Spread of the pp/PbPb significance when the rho's are varied in their allowed ranges.

    For each of the n_rho_samples the rho's are drawn uniformly in their confidence
    intervals, the full covariance is rebuilt with those rho's and the significance of the
    pp/PbPb difference is recomputed. The rho samples are independent, so they are
    distributed over a pool of worker processes.

    Args:
        pp_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for pp data.
        pbpb_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for PbPb data.
        label (str): Label for the plot files.
        n_rho_samples (int): Number of sampled sets of rho's.
        n_toys (int): Toys per rho sample. This caps the significance that can be resolved,
            since a p-value below 1/n_toys cannot be measured: 1e6 toys saturate at about
            4.9 sigma, 1e7 at about 5.3 sigma.
        seed (int): Seed of the random number generator, for reproducibility.
        n_workers (int or None): Number of worker processes.
    """
    _, y_pp, stat_pp, syst_pp, _ = pp_data
    _, y_pbpb, stat_pbpb, syst_pbpb, _ = pbpb_data

    # Sample all the rho sets up front, from a single seeded generator, so that the result
    # is reproducible and independent of how the work is spread over the workers
    rng = np.random.default_rng(seed)
    child_seeds = rng.bit_generator.seed_seq.spawn(n_rho_samples)
    tasks, n_rejected = [], 0
    for child_seed in child_seeds:
        while True:
            rhos_sampled = {}
            for syst_name, syst_rhos in rhos_confidence_intervals.items():
                rhos_sampled[syst_name] = {}
                for key, (low, high) in syst_rhos.items():
                    rhos_sampled[syst_name][key] = rng.uniform(low, high)
            cov_total = build_covariance_matrix(stat_pp, syst_pp, stat_pbpb, syst_pbpb,
                                                rhos_sampled)
            try:
                scipy.linalg.cholesky(cov_total, lower=True)
                break
            except scipy.linalg.LinAlgError:
                n_rejected += 1
        tasks.append((cov_total, child_seed))

    n_drawn = n_rho_samples + n_rejected
    print(f"[rho scan] {n_rejected}/{n_drawn} ({100.*n_rejected/n_drawn:.1f}%) of the drawn "
          f"rho sets were rejected as not positive definite")

    worker = partial(_rho_significance_worker, y_pp=y_pp, y_pbpb=y_pbpb, n_toys=n_toys)

    if n_workers is None:
        n_workers = 100
    n_workers = max(1, min(n_workers, n_rho_samples))

    if n_workers == 1:
        significance_distribution = [worker(task) for task in tqdm(tasks, desc="rho scan")]
    else:
        with mp.get_context("fork").Pool(n_workers) as pool:
            chunksize = max(1, n_rho_samples // (8 * n_workers))
            significance_distribution = list(tqdm(
                pool.imap(worker, tasks, chunksize=chunksize),
                total=n_rho_samples, desc="rho scan", unit="rho"))

    significance_distribution = np.array(significance_distribution)
    low, median, high = np.percentile(significance_distribution, [2.5, 50., 97.5])
    print(f"[rho scan, {n_rho_samples} samples] significance = {median:.2f} sigma (median), "
          f"95% range [{low:.2f}, {high:.2f}], "
          f"min = {significance_distribution.min():.2f}, "
          f"max = {significance_distribution.max():.2f}")

    # Draw the histogram of the significance distribution
    plt.figure(figsize=(8, 6))
    plt.hist(significance_distribution, bins=30, density=True, alpha=0.7, color='blue')
    plt.axvline(median, color='k', linestyle='--', lw=2, label=f'median: {median:.2f}$\\sigma$')
    plt.axvline(low, color='k', linestyle=':', lw=1.5, label=f'95% range: [{low:.2f}, {high:.2f}]$\\sigma$')
    plt.axvline(high, color='k', linestyle=':', lw=1.5)
    plt.xlabel("Significance")
    plt.ylabel("Probability Density")
    plt.title(fr"Significance Distribution from {n_rho_samples} Toy MC Experiments (with different $\rho$)")
    plt.legend()
    plt.grid()
    plt.savefig((Path(__file__).parent / f"significance_distribution_two_fits_both_fluctuate_{label}_{n_rho_samples}.png").resolve())
    plt.close()

def kendalltau(y):
    """Get Kenadll's tau from an ordered set of data"""
    y = np.asarray(y)
    is_single_set = y.ndim == 1
    if is_single_set:
        y = y[:, None]

    n = y.shape[0]
    num = np.zeros(y.shape[1], dtype=np.int32)
    den = np.zeros(y.shape[1], dtype=np.int32)
    for i in range(n):
        for j in range(i + 1, n):
            concordant = y[j] > y[i]
            discordant = y[j] < y[i]
            num += concordant
            num -= discordant
            den += concordant
            den += discordant

    tau = num / den

    return tau[0] if is_single_set else tau

def spearman_rho(y):
    """Get Spearman's rho of an ordered set of data against its position"""
    y = np.asarray(y)
    if y.ndim == 1:
        return spearmanr(np.arange(y.shape[0]), y).statistic

    n = y.shape[0]
    # rank of each point = how many points of the same set are below it
    ranks = np.zeros(y.shape, dtype=np.int16)
    for i in range(n):
        for j in range(i + 1, n):
            below = y[j] < y[i]
            ranks[i] += below
            ranks[j] += ~below

    d = ranks - np.arange(n, dtype=np.int16)[:, None]
    d2 = np.einsum('ij,ij->j', d, d, dtype=np.int64)

    return 1. - 6. * d2 / (n * (n * n - 1))

def order_data_by_x(data):
    """Order the data by x values for Kendall's tau calculation."""
    x, y, stat, syst, _ = data
    order = np.argsort(x)
    x_ordered = np.array(x)[order]
    y_ordered = np.array(y)[order]
    stat_ordered = np.array(stat)[order]
    syst_ordered = {k: np.array(v)[order] for k, v in syst.items()}
    return x_ordered, y_ordered, stat_ordered, syst_ordered

def toy_mc_kendall_single_system(data, label, n_toys=100000, seed=42):
    """
    Perform a toy Monte Carlo study to get the significance assuming pbpb points to have same values as pp points.

    Args:
        data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties.
        label (str): Label for the plot files.
        n_toys (int): Number of toy experiments to perform.
        seed (int): Seed of the random number generator, for reproducibility.
    """
    rng = np.random.default_rng(seed)
    x, y, stat, syst, _ = data
    # we order the points in increasing multiplicity for Kendall's tau
    x, y, stat, syst = order_data_by_x(data)

    n = len(x)

    cov = build_covariance_matrix(stat, syst, stat, syst, rhos)[:n, :n]

    # These are always the same for all fits
    w, err = pol0_weights(cov)

    c_val = w @ np.array(y)

    tau_distribution = []
    y_toys = sample_multivariate(np.full(n , c_val),
                                cov, n_toys, rng)
    tau_distribution = kendalltau(y_toys)

    observed_tau = kendalltau(y)
    print(f"Observed Kendall's tau: {observed_tau:.4f}")
    n_extreme = np.sum(abs(np.array(tau_distribution)) >= abs(observed_tau))
    p_value = n_extreme / len(tau_distribution)
    significance = norm.ppf(1 - p_value/2)  # two-tailed test
    print(f"[both fluctuating] p-value = {p_value:.2e} "
            f"-> significance = {significance:.2f} sigma")

    # Draw the histogram of the Kendall's tau distribution
    fig, ax = plt.subplots(figsize=(8, 6))
    counts, _, _ = ax.hist(tau_distribution, bins=30, density=True, alpha=0.7, color='blue')
    ax.axvline(observed_tau, color='k', linestyle='--', lw=2, label=f'Observed: {observed_tau:.4f}')
    ax.axvline(0., color='g', linestyle=':', lw=2, label='$H_0$: 0')
    # Headroom on the log scale so the legend and the p-value box do not cover the toys
    positive_counts = counts[counts > 0]
    ax.set_ylim(positive_counts.min() / 2., positive_counts.max() * 1.2)
    ax.text(0.03, 0.97,
            f'p-value = {p_value:.2e}\nsignificance = {significance:.2f}$\\sigma$',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)

    plt.xlabel("Kendall's tau")
    plt.ylabel("Probability Density")
    plt.title(fr"Kendall's tau Distribution from {n_toys} Toy MC Experiments")
    plt.savefig((Path(__file__).parent / f"kendall_tau_distribution_{label}_{n_toys}.png").resolve())
    plt.close()

def toy_mc_kendall_tau(pp_data, pbpb_data, label, n_toys=100000, seed=42):
    """
    Perform a toy Monte Carlo study to get the significance assuming pbpb points to have same values as pp points.

    Args:
        pp_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for pp data.
        pbpb_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for PbPb data.
        label (str): Label for the plot files.
        n_toys (int): Number of toy experiments to perform.
        seed (int): Seed of the random number generator, for reproducibility.
    """
    rng = np.random.default_rng(seed)
    x_pp, y_pp, stat_pp, syst_pp, _ = pp_data
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb, _ = pbpb_data
    # we order the points in increasing multiplicity for Kendall's tau
    x_pp, y_pp, stat_pp, syst_pp = order_data_by_x(pp_data)
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb = order_data_by_x(pbpb_data)

    n_pp = len(x_pp)
    n_pbpb = len(x_pbpb)

    cov_total = build_covariance_matrix(stat_pp, syst_pp, stat_pbpb, syst_pbpb, rhos)
    cov_pp = cov_total[:n_pp, :n_pp]

    # These are always the same for all fits
    w_pp, err_pp = pol0_weights(cov_pp)

    c_val_pp = w_pp @ np.array(y_pp)

    tau_distribution = []
    y_toys = sample_multivariate(np.full(n_pp + n_pbpb, c_val_pp),
                                cov_total, n_toys, rng)
    tau_distribution = kendalltau(y_toys)

    observed_tau = kendalltau(np.concatenate([y_pp, y_pbpb]))
    print(f"Observed Kendall's tau: {observed_tau:.4f}")
    n_extreme = np.sum(abs(np.array(tau_distribution)) >= abs(observed_tau))
    p_value = n_extreme / len(tau_distribution)
    significance = norm.ppf(1 - p_value/2)  # two-tailed test
    print(f"[both fluctuating] p-value = {p_value:.2e} "
            f"-> significance = {significance:.2f} sigma")

    # Draw the histogram of the Kendall's tau distribution
    fig, ax = plt.subplots(figsize=(8, 6))
    counts, _, _ = ax.hist(tau_distribution, bins=30, density=True, alpha=0.7, color='blue')
    ax.axvline(observed_tau, color='k', linestyle='--', lw=2, label=f'Observed: {observed_tau:.4f}')
    ax.axvline(0., color='g', linestyle=':', lw=2, label='$H_0$: 0')
    # Headroom on the log scale so the legend and the p-value box do not cover the toys
    positive_counts = counts[counts > 0]
    ax.set_ylim(positive_counts.min() / 2., positive_counts.max() * 1.2)
    ax.text(0.03, 0.97,
            f'p-value = {p_value:.2e}\nsignificance = {significance:.2f}$\\sigma$',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)

    plt.xlabel("Kendall's tau")
    plt.ylabel("Probability Density")
    plt.title(fr"Kendall's tau Distribution from {n_toys} Toy MC Experiments")
    plt.savefig((Path(__file__).parent / f"kendall_tau_distribution_{label}_{n_toys}.png").resolve())
    plt.close()

def toy_mc_kendall_tau_oo(pp_data, pbpb_data, oo_data, label, n_toys=100000, seed=42):
    """
    Perform a toy Monte Carlo study to get the significance assuming pbpb points to have same values as pp points.

    Args:
        pp_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for pp data.
        pbpb_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for PbPb data.
        oo_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for oo data.
        label (str): Label for the plot files.
        n_toys (int): Number of toy experiments to perform.
        seed (int): Seed of the random number generator, for reproducibility.
    """
    rng = np.random.default_rng(seed)
    x_pp, y_pp, stat_pp, syst_pp, _ = pp_data
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb, _ = pbpb_data
    x_oo, y_oo, stat_oo, syst_oo, _ = oo_data
    # we order the points in increasing multiplicity for Kendall's tau
    x_pp, y_pp, stat_pp, syst_pp = order_data_by_x(pp_data)
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb = order_data_by_x(pbpb_data)
    x_oo, y_oo, stat_oo, syst_oo = order_data_by_x(oo_data)

    n_pp = len(x_pp)
    n_pbpb = len(x_pbpb)
    n_oo = len(x_oo)

    cov_total = build_covariance_matrix_with_oo(stat_pp, syst_pp, stat_pbpb, syst_pbpb, stat_oo, syst_oo, rhos)
    cov_pp = cov_total[:n_pp, :n_pp]

    # These are always the same for all fits
    w_pp, err_pp = pol0_weights(cov_pp)

    c_val_pp = w_pp @ np.array(y_pp)

    order_tot = np.argsort(np.concatenate([x_pp, x_pbpb, x_oo]))
    x_tot = np.concatenate([x_pp, x_pbpb, x_oo])[order_tot]
    y_tot = np.concatenate([y_pp, y_pbpb, y_oo])[order_tot]
    stat_tot = np.concatenate([stat_pp, stat_pbpb, stat_oo])[order_tot]
    # order the systematic uncertainties in the same way
    syst_tot = {}
    for k in syst_pp.keys():
        syst_tot[k] = np.concatenate([syst_pp[k], syst_pbpb[k], syst_oo[k]])[order_tot]

    cov_sorted = cov_total[np.ix_(order_tot, order_tot)]

    y_toys = sample_multivariate(np.full(n_pp + n_pbpb + n_oo, c_val_pp),
                                cov_sorted, n_toys, rng)
    tau_distribution = kendalltau(y_toys)
    observed_tau = kendalltau(y_tot)
    print(f"Observed Kendall's tau: {observed_tau:.4f}")
    n_extreme = np.sum(abs(np.array(tau_distribution)) >= abs(observed_tau))
    p_value = n_extreme / len(tau_distribution)
    significance = norm.ppf(1 - p_value/2)  # two-tailed test
    print(f"[both fluctuating] p-value = {p_value:.2e} "
            f"-> significance = {significance:.2f} sigma")

    # Draw the histogram of the Kendall's tau distribution
    fig, ax = plt.subplots(figsize=(8, 6))
    counts, _, _ = ax.hist(tau_distribution, bins=30, density=True, log=True, alpha=0.7)
    ax.axvline(observed_tau, color='k', linestyle='--', lw=2, label=f'Observed: {observed_tau:.4f}')
    ax.axvline(0., color='g', linestyle=':', lw=2, label='$H_0$: 0')
    # Headroom on the log scale so the legend and the p-value box do not cover the toys
    positive_counts = counts[counts > 0]
    ax.set_ylim(positive_counts.min() / 2., positive_counts.max() * 1.2)
    ax.text(0.03, 0.97,
            f'p-value = {p_value:.2e}\nsignificance = {significance:.2f}$\\sigma$',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)

    plt.xlabel("Kendall's tau")
    plt.ylabel("Probability Density")
    plt.title(fr"Kendall's tau Distribution from {n_toys} Toy MC Experiments")
    plt.savefig((Path(__file__).parent / f"kendall_tau_distribution_{label}_{n_toys}.png").resolve())
    plt.close()

def toy_mc_spearman_tau_oo(pp_data, pbpb_data, oo_data, label, n_toys=100000, seed=42):
    """
    Perform a toy Monte Carlo study to get the significance assuming pbpb points to have same values as pp points.

    Args:
        pp_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for pp data.
        pbpb_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for PbPb data.
        oo_data (tuple): A tuple containing x values, y values, statistical uncertainties, systematic uncertainties, and total systematic uncertainties for oo data.
        label (str): Label for the plot files.
        n_toys (int): Number of toy experiments to perform.
        seed (int): Seed of the random number generator, for reproducibility.
    """
    rng = np.random.default_rng(seed)
    x_pp, y_pp, stat_pp, syst_pp, _ = pp_data
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb, _ = pbpb_data
    x_oo, y_oo, stat_oo, syst_oo, _ = oo_data
    # we order the points in increasing multiplicity for Spearman's rho
    x_pp, y_pp, stat_pp, syst_pp = order_data_by_x(pp_data)
    x_pbpb, y_pbpb, stat_pbpb, syst_pbpb = order_data_by_x(pbpb_data)
    x_oo, y_oo, stat_oo, syst_oo = order_data_by_x(oo_data)

    n_pp = len(x_pp)
    n_pbpb = len(x_pbpb)
    n_oo = len(x_oo)

    cov_total = build_covariance_matrix_with_oo(stat_pp, syst_pp, stat_pbpb, syst_pbpb, stat_oo, syst_oo, rhos)
    cov_pp = cov_total[:n_pp, :n_pp]

    # These are always the same for all fits
    w_pp, err_pp = pol0_weights(cov_pp)

    c_val_pp = w_pp @ np.array(y_pp)

    order_tot = np.argsort(np.concatenate([x_pp, x_pbpb, x_oo]))
    x_tot = np.concatenate([x_pp, x_pbpb, x_oo])[order_tot]
    y_tot = np.concatenate([y_pp, y_pbpb, y_oo])[order_tot]
    stat_tot = np.concatenate([stat_pp, stat_pbpb, stat_oo])[order_tot]
    # order the systematic uncertainties in the same way
    syst_tot = {}
    for k in syst_pp.keys():
        syst_tot[k] = np.concatenate([syst_pp[k], syst_pbpb[k], syst_oo[k]])[order_tot]

    cov_sorted = cov_total[np.ix_(order_tot, order_tot)]

    y_toys = sample_multivariate(np.full(n_pp + n_pbpb + n_oo, c_val_pp),
                                cov_sorted, n_toys, rng)
    rho_distribution = spearman_rho(y_toys)
    observed_rho = spearman_rho(y_tot)
    print(f"Observed Spearman's rho: {observed_rho:.4f}")
    n_extreme = np.sum(abs(np.array(rho_distribution)) >= abs(observed_rho))
    p_value = n_extreme / len(rho_distribution)
    significance = norm.ppf(1 - p_value/2)  # two-tailed test
    print(f"[both fluctuating] p-value = {p_value:.2e} "
            f"-> significance = {significance:.2f} sigma")

    # Draw the histogram of the Spearman's rho distribution
    fig, ax = plt.subplots(figsize=(8, 6))
    counts, _, _ = ax.hist(rho_distribution, bins=30, density=True, log=True, alpha=0.7)
    ax.axvline(observed_rho, color='k', linestyle='--', lw=2, label=f'Observed: {observed_rho:.4f}')
    ax.axvline(0., color='g', linestyle=':', lw=2, label='$H_0$: 0')
    # Headroom on the log scale so the legend and the p-value box do not cover the toys
    positive_counts = counts[counts > 0]
    ax.set_ylim(positive_counts.min() / 2., positive_counts.max() * 1.2)
    ax.text(0.03, 0.97,
            f'p-value = {p_value:.2e}\nsignificance = {significance:.2f}$\\sigma$',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)

    plt.xlabel("Spearman's rho")
    plt.ylabel("Probability Density")
    plt.title(fr"Spearman's rho Distribution from {n_toys} Toy MC Experiments")
    plt.savefig((Path(__file__).parent / f"spearman_rho_distribution_{label}_{n_toys}.png").resolve())
    plt.close()

if __name__ == "__main__":
    (x_pp, y_pp, stat_pp, syst_tot_pp), (x_pbpb, y_pbpb, stat_pbpb, syst_tot_pbpb) = get_data()
    (x_oo, y_oo, stat_oo, syst_tot_oo) = get_data_oo()
    database_syst_pp, database_syst_pbpb, database_syst_oo = get_systematics_database()
    syst_pp = get_systematics(database_syst_pp, x_pp, y_pp)
    syst_pbpb = get_systematics(database_syst_pbpb, x_pbpb, y_pbpb)
    syst_oo = get_systematics(database_syst_oo, x_oo, y_oo)

    pp_data = (x_pp, y_pp, stat_pp, syst_pp, syst_tot_pp)
    pbpb_data = (x_pbpb, y_pbpb, stat_pbpb, syst_pbpb, syst_tot_pbpb)
    oo_data = (x_oo, y_oo, stat_oo, syst_oo, syst_tot_oo)
    syst_pbpb_central = {k: v[:len(x_pbpb)-3] for k, v in syst_pbpb.items()}
    pbpb_data_central = (x_pbpb[:len(x_pbpb)-3], y_pbpb[:len(x_pbpb)-3], stat_pbpb[:len(x_pbpb)-3], syst_pbpb_central, syst_tot_pbpb[:len(x_pbpb)-3])

    syst_pbpb_peripheral = {k: v[-3:] for k, v in syst_pbpb.items()}
    pbpb_data_peripheral = (x_pbpb[-3:], y_pbpb[-3:], stat_pbpb[-3:], syst_pbpb_peripheral, syst_tot_pbpb[-3:])

    print("Fitting a constant to the combined pp and PbPb data...")
    get_significance_from_fit_to_all_data(pp_data, pbpb_data, "combined_fit")
    print()

    print("Fitting a constant to the combined pp and most central PbPb data...")
    get_significance_from_fit_to_all_data(pp_data, pbpb_data_central, "combined_fit_central")
    print()

    print("Fitting a constant to pp and another to PbPb data...")
    get_significance_two_fits(pp_data, pbpb_data, "two_fits")
    print()

    print("Fitting a constant to pp and another to central PbPb data...")
    get_significance_two_fits(pp_data, pbpb_data_central, "two_fits_central")
    print()

    print("Performing toy Monte Carlo study to get significance distribution given the different correlations...")
    toy_mc_rho_confidence_intervals(pp_data, pbpb_data, n_toys=10000)    
    print()

    print("Performing toy Monte Carlo study under H0 (true PbPb ratio = true pp ratio), pp reference fixed...")
    toy_mc_mean_pbpb_pp_fixed(pp_data, pbpb_data, "pp_fixed", n_toys=10000000)
    print()

    print("Performing toy Monte Carlo study under H0 (true PbPb ratio = true pp ratio), both systems fluctuating...")
    toy_mc_mean_pbpb_both_fluctuate(pp_data, pbpb_data, "pp_also_fluctuating", n_toys=10000000)
    print()

    print("Performing toy Monte Carlo study under H0 only with central pbpb data (true PbPb ratio = true pp ratio), pp reference fixed...")
    toy_mc_mean_pbpb_pp_fixed(pp_data, pbpb_data_central, "pp_fixed_central", n_toys=10000000)
    print()

    print("Performing toy Monte Carlo study under H0 only with central pbpb data (true PbPb ratio = true pp ratio), both systems fluctuating...")
    toy_mc_mean_pbpb_both_fluctuate(pp_data, pbpb_data_central, "pp_also_fluctuating_central", n_toys=10000000)
    print()

    # print("Performing toy Monte Carlo study under H0 only with central pbpb data (true PbPb ratio = true pp ratio), both systems fluctuating...")
    # toy_mc_mean_pbpb_both_fluctuate_rho_confidence_intervals(pp_data, pbpb_data_central, "pp_also_fluctuating_central", n_rho_samples=10000)
    # print()

    print("Performing toy Monte Carlo for Kendall's tau")
    toy_mc_kendall_tau(pp_data, pbpb_data, "kendalls", n_toys=100000000)
    print()

    print("Performing toy Monte Carlo for Kendall's tau")
    toy_mc_kendall_tau_oo(pp_data, pbpb_data, oo_data, "kendalls_oo", n_toys=100000000)
    print()

    print("Performing toy Monte Carlo for Kendall's tau (pp only)")
    toy_mc_kendall_single_system(pp_data, "kendalls_pp_only", n_toys=100000000)
    print()

    print("Performing toy Monte Carlo for Kendall's tau (PbPb only)")
    toy_mc_kendall_single_system(pbpb_data, "kendalls_pbpb_only", n_toys=100000000)
    print()

    print("Performing toy Monte Carlo for Kendall's tau (OO only)")
    toy_mc_kendall_single_system(oo_data, "kendalls_oo_only", n_toys=100000000)
    print()

    print("Performing toy Monte Carlo for Spearman's coefficient")
    toy_mc_spearman_tau_oo(pp_data, pbpb_data, oo_data, "spearman_oo", n_toys=100000000)
    print()


