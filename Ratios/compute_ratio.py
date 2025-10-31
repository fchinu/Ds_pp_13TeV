import ROOT
import numpy as np
import argparse
import yaml


def tsallis(pt, pars):
    """
    Tsallis function for pT extrapolation

    Parameters:
    - pt (list): 1-D list for pT (variable)
    - pars (list): list of parameters

    Returns:
    - double: power law function evaluated for a given pt

    """
    c_par = (pars[1] - 1) * (pars[1] - 2) / (pars[1] * pars[2] * (pars[1] * pars[2] + 2 * (pars[1] - 2) * pars[3]))
    mt_par = np.sqrt(pars[3]*pars[3] + pt[0]*pt[0]) # pars[3] is the mass of the particle

    return pars[0] * c_par * pt[0] * np.power((1 + (mt_par - pars[3]) / pars[1] / pars[2]), -pars[1])


# definition of shared parameters for Ds, D+ combined Tsallis fit
pars_tsallis_ds = np.array(
    [
        0, # integral Ds
        1, # n
        2, # T
        3 # m(Ds)
    ],
    dtype=np.int32
)  # exp amplitude in B histo and exp common parameter


pars_tsallis_dplus = np.array(
    [
        4,  # integral D+
        1,  # n
        2,  # T
        5  # m(D+)
    ],
    dtype=np.int32,
)

class GlobalChi2(object):
    """
    Class for combined Chi2 in case of simultaneous Ds, D+ Tsallis fit
    """
    def __init__(self, f1, f2):
        self._f1 = f1
        self._f2 = f2

    def __call__(self, par):
        # parameter vector is first background (in common 1 and 2) and then is
        # signal (only in 2)

        # the zero-copy way to get a numpy array from a double *
        par_arr = np.frombuffer(par, dtype=np.float64, count=6)

        p1 = par_arr[pars_tsallis_ds]
        p2 = par_arr[pars_tsallis_dplus]

        return self._f1(p1) + self._f2(p2)


def get_raw_yields(config, cent_bin):
    """
    Retrieve raw yield histograms from a ROOT file based on centrality ranges.

    Parameters:
    - config (dict): Configuration dictionary containing input file paths.
    - cent_bin (int): centrality bin

    Returns:
    - tuple: A tuple containing two ROOT histograms:
        - h_raw_yields_ds: Histogram for raw yields of ds.
        - h_raw_yields_dplus: Histogram for raw yields of dplus.
    """
    with ROOT.TFile.Open(config['inputs']['files']['raw_yields']) as in_file_raw_yields:
        h_raw_yields_ds = in_file_raw_yields.Get(
            config['inputs']['histo_names']['raw_yields_ds'][cent_bin])
        h_raw_yields_ds.SetDirectory(0)
        h_raw_yields_dplus = in_file_raw_yields.Get(
            config['inputs']['histo_names']['raw_yields_dp'][cent_bin])
        h_raw_yields_dplus.SetDirectory(0)

    return h_raw_yields_ds, h_raw_yields_dplus


def get_efficiencies(config, cent_bin):
    """
    Retrieve efficiency histograms for prompt Ds and D+ from a ROOT file.

    Parameters:
    - config (dict): Configuration dictionary containing input file paths.
    - cent_bin (int): centrality bin

    Returns:
    - tuple: A tuple containing two ROOT histograms:
        - h_eff_prompt_ds: Efficiency histogram for prompt Ds.
        - h_eff_prompt_dplus: Efficiency histogram for prompt D+.
    """
    with ROOT.TFile.Open(config['inputs']['files']['eff']) as eff_file:
        h_eff_prompt_ds = eff_file.Get(
            config['inputs']['histo_names']['eff_ds'][cent_bin])
        h_eff_prompt_ds.SetDirectory(0)
        h_eff_prompt_dplus = eff_file.Get(
            config['inputs']['histo_names']['eff_dp'][cent_bin])
        h_eff_prompt_dplus.SetDirectory(0)

    return h_eff_prompt_ds, h_eff_prompt_dplus


def get_prompt_fractions(config, cent_bin):
    """
    Retrieves the prompt fractions histograms for Ds and D+ from the specified ROOT files.

    Args:
    - config (dict): Configuration dictionary containing input file paths.
    - cent_bin (int): centrality bin

    Returns:
    - tuple: A tuple containing two ROOT histograms:
        - h_prompt_frac_ds: Prompt fraction for Ds.
        - h_prompt_frac_dplus: Prompt fraction for D+.
    """
    if config['inputs']['files']['frac_ds'][cent_bin] is None or config['inputs']['files']['frac_dp'][cent_bin] is None:
        return None, None

    with ROOT.TFile.Open(config['inputs']['files']['frac_ds'][cent_bin]) as in_file_prompt_frac_ds:
        h_prompt_frac_ds = in_file_prompt_frac_ds.Get("hRawFracPrompt")
        h_prompt_frac_ds.SetDirectory(0)

    with ROOT.TFile.Open(config['inputs']['files']['frac_dp'][cent_bin]) as in_file_prompt_frac_dplus:
        h_prompt_frac_dplus = in_file_prompt_frac_dplus.Get("hRawFracPrompt")
        h_prompt_frac_dplus.SetDirectory(0)

    return h_prompt_frac_ds, h_prompt_frac_dplus


def get_collisions(config):
    """
    Retrieve the collisions histogram from a ROOT file.

    Parameters:
    - config (dict): Configuration dictionary containing input file paths.
    - cent_bin (int): centrality bin

    Returns:
    - float: Number of collisions for the specified centrality bin.
    """
    with ROOT.TFile.Open(config['inputs']['files']['raw_yields']) as in_file_raw_yields:
        h_collisions = in_file_raw_yields.Get('h_coll_rebinned')
        h_collisions.SetDirectory(0)

    return h_collisions


def get_ratio_vs_pt(
    config, h_raw_yields_ds, h_raw_yields_dplus, h_eff_prompt_ds,
    h_eff_prompt_dplus, h_prompt_frac_ds, h_prompt_frac_dplus
    ):
    """
    Calculate the ratio of corrected yields of D_s^+ to D^+.

    Parameters:
    - config (dict): Configuration dictionary.
    - h_raw_yields_ds (ROOT.TH1): Histogram of raw yields for D_s^+.
    - h_raw_yields_dplus (ROOT.TH1): Histogram of raw yields for D^+.
    - h_eff_prompt_ds (ROOT.TH1): Histogram of prompt efficiency for D_s^+.
    - h_eff_prompt_dplus (ROOT.TH1): Histogram of prompt efficiency for D^+.
    - h_prompt_frac_ds (ROOT.TH1): Histogram of prompt fraction for D_s^+.
    - h_prompt_frac_dplus (ROOT.TH1): Histogram of prompt fraction for D^+.

    Returns:
        ROOT.TH1: Histogram of the ratio of corrected yields of D_s^+ to D^+.
    """
    h_corrected_yields_ds = h_raw_yields_ds.Clone("h_corrected_yields_ds")
    h_corrected_yields_ds.Divide(h_eff_prompt_ds)
    h_corrected_yields_ds.Scale(1 / config['br']['ds_to_phipi_to_kkpi'], "width")
    if h_prompt_frac_ds is not None:
        h_corrected_yields_ds.Multiply(h_prompt_frac_ds)
    else:
        print("WARNING: computing ratios without feeddown correction")

    h_corrected_yields_dplus = h_raw_yields_dplus.Clone("h_corrected_yields_dplus")
    h_corrected_yields_dplus.Divide(h_eff_prompt_dplus)
    h_corrected_yields_dplus.Scale(1 / config['br']['dplus_to_phipi_to_kkpi'], "width")
    if h_prompt_frac_dplus is not None:
        h_corrected_yields_dplus.Multiply(h_prompt_frac_dplus)
    else:
        print("WARNING: computing ratios without feeddown correction")

    # pt extrapolation
    pt_min = 0.
    pt_max_ds = 24.
    pt_max_dplus = 12. # poor precision at high pT
    f_corrected_yields_ds = [
        ROOT.TF1("f_corrected_yields_ds_tsallis", tsallis, pt_min, pt_max_ds, 4),
        ROOT.TF1("f_corrected_yields_ds_tsallissim", tsallis, pt_min, pt_max_ds, 4)
    ]
    f_corrected_yields_dplus = [
        ROOT.TF1("f_corrected_yields_dplus_tsallis", tsallis, pt_min, pt_max_dplus, 4),
        ROOT.TF1("f_corrected_yields_dplus_tsallissim", tsallis, pt_min, pt_max_dplus, 4)
    ]
    # indpendent fits with Tsallis
    f_corrected_yields_ds[0].SetParameters(h_corrected_yields_ds.Integral(), 10, 0.5, 1.968)
    f_corrected_yields_dplus[0].SetParameters(h_corrected_yields_dplus.Integral(), 10, 0.5, 1.870)
    f_corrected_yields_ds[0].FixParameter(3, 1.968)
    f_corrected_yields_dplus[0].FixParameter(3, 1.870)
    f_corrected_yields_ds[0].SetParLimits(1, 1.5, 15.)
    f_corrected_yields_dplus[0].SetParLimits(1, 1.5, 15.)
    fitres_ds = h_corrected_yields_ds.Fit(f_corrected_yields_ds[0], "IMEQ0S")
    fitres_dplus = h_corrected_yields_dplus.Fit(f_corrected_yields_dplus[0], "IMEQ0S")
    # We need to retrieve the covariance matrices to compute the uncertainties on the integrals
    cov_mat_ds = fitres_ds.GetCovarianceMatrix()
    cov_mat_dplus = fitres_dplus.GetCovarianceMatrix()

    # let's also do a simultaneous fit with Tsallis
    wf_corrected_yields_ds = ROOT.Math.WrappedMultiTF1(f_corrected_yields_ds[1], 1)
    wf_corrected_yields_dplus = ROOT.Math.WrappedMultiTF1(f_corrected_yields_dplus[1], 1)

    opt = ROOT.Fit.DataOptions()
    range_ds = ROOT.Fit.DataRange()
    range_ds.SetRange(pt_min, pt_max_ds)
    data_ds = ROOT.Fit.BinData(opt, range_ds)
    ROOT.Fit.FillData(data_ds, h_corrected_yields_ds)
    range_dplus = ROOT.Fit.DataRange()
    range_dplus.SetRange(pt_min, pt_max_ds)
    data_dplus = ROOT.Fit.BinData(opt, range_dplus)
    ROOT.Fit.FillData(data_dplus, h_corrected_yields_dplus)

    chi2_ds = ROOT.Fit.Chi2Function(data_ds, wf_corrected_yields_ds)
    chi2_dplus = ROOT.Fit.Chi2Function(data_dplus, wf_corrected_yields_dplus)
    global_chi2 = GlobalChi2(chi2_ds, chi2_dplus)
    fitter_tsallis = ROOT.Fit.Fitter()

    pars_all = np.array(
        [h_corrected_yields_ds.Integral(), 10, 0.5, 1.968,
         h_corrected_yields_dplus.Integral(), 1.870]
    )
    fitter_tsallis.Config().SetParamsSettings(6, pars_all)
    # fix the masses
    fitter_tsallis.Config().ParSettings(3).Fix()
    fitter_tsallis.Config().ParSettings(5).Fix()
    # set limits
    fitter_tsallis.Config().ParSettings(1).SetLimits(1.5, 15.)
    fitter_tsallis.Config().MinimizerOptions().SetPrintLevel(0)
    fitter_tsallis.Config().SetMinimizer("Minuit2", "Migrad")

    # we can't pass the Python object global_chi2 directly to FitFCN.
    # It needs to be wrapped in a ROOT::Math::Functor.
    global_chi2_functor = ROOT.Math.Functor(global_chi2, 6)

    # fit FCN function
    # (specify optionally data size and flag to indicate that is a chi2 fit)
    fitter_tsallis.FitFCN(global_chi2_functor, 0, data_ds.Size() + data_dplus.Size(), 1)
    result_tsallis = fitter_tsallis.Result()
    result_tsallis.Print(ROOT.std.cout)
    f_corrected_yields_ds[1].SetFitResult(result_tsallis, pars_tsallis_ds)
    f_corrected_yields_dplus[1].SetFitResult(result_tsallis, pars_tsallis_dplus)

    # We need to retrieve the covariance matrices to compute the uncertainties on the integrals
    cov_mat_sim = ROOT.TMatrixDSym(6)
    result_tsallis.GetCovarianceMatrix(cov_mat_sim)

    h_ratio = h_corrected_yields_ds.Clone("h_ratio")
    h_ratio.Divide(h_corrected_yields_dplus)
    h_ratio.SetTitle(";#it{p}_{T} (GeV/#it{c});D_{s}^{+}/D^{+} Ratio")

    h_ratio_ptint = [
        ROOT.TH1D("h_ratio_ptint_tsallis", ";;D_{s}^{+}/D^{+} Ratio", 1, pt_min, pt_max_ds),
        ROOT.TH1D("h_ratio_ptint_tsallissim", ";;D_{s}^{+}/D^{+} Ratio", 1, pt_min, pt_max_ds)
    ]

    for ifunc in range(2):
        func_covmat_ds = cov_mat_ds if ifunc == 0 else cov_mat_sim
        func_fitres_ds = fitres_ds if ifunc == 0 else result_tsallis
        func_covmat_dplus = cov_mat_dplus if ifunc == 0 else cov_mat_sim
        func_fitres_dplus = fitres_dplus if ifunc == 0 else result_tsallis
        
        # No need to multiply by bin width (TF1 integral already does that)
        pt_min_meas = h_corrected_yields_ds.GetBinLowEdge(1)
        yield_ds = f_corrected_yields_ds[ifunc].Integral(0., pt_min_meas)
        yield_dplus = f_corrected_yields_dplus[ifunc].Integral(0., pt_min_meas)
        unc_yield_ds = f_corrected_yields_ds[ifunc].IntegralError(
            0., pt_min_meas, 
            func_fitres_ds.GetParams(), 
            func_covmat_ds.GetMatrixArray())**2
        
        unc_yield_dplus = f_corrected_yields_dplus[ifunc].IntegralError(
            0., pt_min_meas, 
            func_fitres_dplus.GetParams(), 
            func_covmat_dplus.GetMatrixArray())**2
        for ipt in range(1, h_corrected_yields_ds.GetNbinsX()+1):
            yield_ds += h_corrected_yields_ds.GetBinContent(ipt)
            yield_dplus += h_corrected_yields_dplus.GetBinContent(ipt)
            unc_yield_ds += h_corrected_yields_ds.GetBinError(ipt)**2
            unc_yield_dplus += h_corrected_yields_dplus.GetBinError(ipt)**2
        unc_yield_ds = np.sqrt(unc_yield_ds)
        unc_yield_dplus = np.sqrt(unc_yield_dplus)
        ratio = yield_ds/yield_dplus
        unc_ratio = np.sqrt((unc_yield_ds/yield_ds)**2+(unc_yield_dplus/yield_dplus)**2) * ratio
        h_ratio_ptint[ifunc].SetBinContent(1, ratio)
        h_ratio_ptint[ifunc].SetBinError(1, unc_ratio)

    n_points = 99
    g_ratio_tsallis = [ROOT.TGraph(n_points), ROOT.TGraph(n_points)]
    g_ratio_tsallis[0].SetName("g_ratio_tsallis")
    g_ratio_tsallis[1].SetName("g_ratio_tsallissim")
    ipt=0
    for pt in np.arange(pt_min+0.01, pt_max_ds, (pt_max_ds-pt_min-0.01)/n_points):
        for ifunc in range(2):
            g_ratio_tsallis[ifunc].SetPoint(
                ipt, pt,
                f_corrected_yields_ds[ifunc].Eval(pt)/f_corrected_yields_dplus[ifunc].Eval(pt))
        ipt += 1

    return h_ratio, h_corrected_yields_ds, h_corrected_yields_dplus, \
        f_corrected_yields_ds, f_corrected_yields_dplus, g_ratio_tsallis, h_ratio_ptint


def get_ratios_vs_cent(h_ratios_vs_pt, h_ratios_pt_int, cent_mins, cent_maxs):
    """
    Generate histograms of ratios versus centrality from input histograms of ratios versus pT.

    Parameters:
    - h_ratios_vs_pt (list of ROOT.TH1): List of histograms representing ratios
        versus pT for different centrality bins.
    - h_ratios_pt_int (list of list of ROOT.TH1): List of list of histograms representing pT-integrated ratios 
        for different centrality bins.
    - cent_mins (list of float): List of minimum centrality values for each bin.
    - cent_maxs (list of float): List of maximum centrality values for each bin.

    Returns:
    - histos (list of ROOT.TH1D): List of histograms representing ratios versus centrality.
    """
    h_ratios_vs_pt = h_ratios_vs_pt.copy()
    cent_tuple = list(zip(cent_mins, cent_maxs))
    if (0, 100) in cent_tuple:
        idx_zero_hundred = cent_tuple.index((0, 100))
        cent_mins.pop(idx_zero_hundred)
        cent_maxs.pop(idx_zero_hundred)
        h_ratios_vs_pt.pop(idx_zero_hundred)
        h_ratios_pt_int.pop(idx_zero_hundred)
    cent_edges = np.asarray(cent_mins + [cent_maxs[-1]], "d")
    histos = []
    for i in range(h_ratios_vs_pt[0].GetNbinsX()):
        suffix = f"{h_ratios_vs_pt[0].GetBinLowEdge(i+1)*10:.0f}_{h_ratios_vs_pt[0].GetBinLowEdge(i+2)*10:.0f}"
        histos.append(ROOT.TH1D(
            f"h_ratio_cent_{suffix}",
            ";Centrality percentile;D_{s}^{+}/D^{+} Ratio",
            len(cent_edges)-1, cent_edges
        ))
        histos[-1].SetDirectory(0)
        for i_cent in range(len(cent_mins)):
            histos[-1].SetBinContent(i_cent+1, h_ratios_vs_pt[i_cent].GetBinContent(i+1))
            histos[-1].SetBinError(i_cent+1, h_ratios_vs_pt[i_cent].GetBinError(i+1))

    for ifunc, _ in enumerate(h_ratios_pt_int[0]):
        func_str = h_ratios_pt_int[0][ifunc].GetName().split(sep="_")[-1]
        suffix = f"pt_int_{func_str}"
        histos.append(ROOT.TH1D(
            f"h_ratio_cent_{suffix}",
            ";Centrality percentile;D_{s}^{+}/D^{+} Ratio",
            len(cent_edges)-1, cent_edges
        ))
        for i_cent in range(len(cent_mins)):
            histos[-1].SetBinContent(i_cent+1, h_ratios_pt_int[i_cent][ifunc].GetBinContent(1))
            histos[-1].SetBinError(i_cent+1, h_ratios_pt_int[i_cent][ifunc].GetBinError(1))

    return histos


def get_ratios_vs_dndeta(h_ratios_vs_cent, cent_mins, cent_maxs, mult_values, mult_unc_low, mult_unc_high):
    cent_tuple = list(zip(cent_mins, cent_maxs))
    if (0, 100) in cent_tuple:
        idx_zero_hundred = cent_tuple.index((0, 100))
        cent_mins.pop(idx_zero_hundred)
        cent_maxs.pop(idx_zero_hundred)
        mult_values.pop(idx_zero_hundred)
        mult_unc_low.pop(idx_zero_hundred)
        mult_unc_high.pop(idx_zero_hundred)
        h_ratios_vs_cent.pop(idx_zero_hundred)
    graphs = []
    for histo in h_ratios_vs_cent:
        y, y_err = [], []
        for i in range(histo.GetNbinsX()):
            y.append(histo.GetBinContent(i+1))
            y_err.append(histo.GetBinError(i+1))
        x = np.asarray(mult_values, "d")
        y = np.asarray(y, "d")
        y_err = np.asarray(y_err, "d")
        x_err_low = np.asarray(mult_unc_low, "d")
        x_err_high = np.asarray(mult_unc_high, "d")
        suffix = histo.GetName().split("h_ratio_cent")[-1]
        graphs.append(ROOT.TGraphAsymmErrors(len(x), x, y, x_err_low, x_err_high, y_err, y_err))
        graphs[-1].SetName(f"g_ratio_dndeta{suffix}")

    return graphs


def evaluate_ratio(config_file_name):
    """
    Evaluate the Ds+/D+ ratio and save the histograms to a ROOT file.

    Parameters:
        config_file_name (str): Path to the configuration YAML file.
    """
    with open(config_file_name, 'r', encoding="utf8") as file:
        config = yaml.safe_load(file)

    h_ratios_vs_pt, h_corry_ds_vs_pt, h_corry_dplus_vs_pt, f_corry_ds_vs_pt, f_corry_dplus_vs_pt, \
        g_ratio_tsallis_vs_pt, h_ratio_ptint = ([] for _ in range(7))
    cent_mins = config["centrality"]["mins"]
    cent_maxs = config["centrality"]["maxs"]
    with ROOT.TFile(config['output_file'], "recreate") as output_file:
        for icent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
            h_rawy_ds, h_rawy_dplus = get_raw_yields(config, icent)
            h_eff_ds_prompt, h_eff_dplus_prompt = get_efficiencies(config, icent)
            h_fprompt_ds, h_fprompt_dplus = get_prompt_fractions(config, icent)
            ratio_out = get_ratio_vs_pt(
                config, h_rawy_ds, h_rawy_dplus, h_eff_ds_prompt,
                h_eff_dplus_prompt, h_fprompt_ds, h_fprompt_dplus
            )
            h_ratios_vs_pt.append(ratio_out[0])
            h_corry_ds_vs_pt.append(ratio_out[1])
            h_corry_dplus_vs_pt.append(ratio_out[2])
            f_corry_ds_vs_pt.append(ratio_out[3])
            f_corry_dplus_vs_pt.append(ratio_out[4])
            g_ratio_tsallis_vs_pt.append(ratio_out[5])
            h_ratio_ptint.append(ratio_out[6])
            if cent_min is not None and cent_max is not None:
                output_file.mkdir(f"centrality_{cent_min}_{cent_max}")
                output_file.cd(f"centrality_{cent_min}_{cent_max}")

            h_rawy_ds.Write()
            h_rawy_dplus.Write()
            h_eff_ds_prompt.Write()
            h_eff_dplus_prompt.Write()
            if h_fprompt_ds is not None:
                h_fprompt_ds.Write()
            if h_fprompt_dplus is not None:
                h_fprompt_dplus.Write()
            h_ratios_vs_pt[icent].Write()
            h_corry_ds_vs_pt[icent].Write()
            h_corry_dplus_vs_pt[icent].Write()
            for func in f_corry_ds_vs_pt[icent]:
                func.Write()
            for func in f_corry_dplus_vs_pt[icent]:
                func.Write()
            for graph in g_ratio_tsallis_vs_pt[icent]:
                graph.Write()
            for histo in h_ratio_ptint[icent]:
                histo.Write()

        if cent_mins is not None and cent_maxs is not None:
            h_ratios_vs_cent = get_ratios_vs_cent(h_ratios_vs_pt, h_ratio_ptint, cent_mins, cent_maxs)
            output_file.cd()
            for h in h_ratios_vs_cent:
                h.Write()
            mult_values = config["mult"]["values"]
            mult_unc_low = config["mult"]["unc_low"]
            mult_unc_high = config["mult"]["unc_high"]
            g_ratios_vs_dndeta = get_ratios_vs_dndeta(h_ratios_vs_cent, cent_mins, cent_maxs,
                                                      mult_values, mult_unc_low, mult_unc_high)
            for g in g_ratios_vs_dndeta:
                g.Write()

        h_collisions = get_collisions(config)
        output_file.cd()
        h_collisions.Write()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate Ds+/D+ ratio')
    parser.add_argument('config_file', metavar='text', help='Path to the config file')
    args = parser.parse_args()

    evaluate_ratio(args.config_file)
