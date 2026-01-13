import yaml
import ROOT
import argparse
import numpy as np

ROOT.TH1.AddDirectory(False)

def get_data(particle, is_prompt, cent_min=None, cent_max=None):
    """
    Get data histogram from unified data ROOT file

    Parameters:
    - particle (str): "Ds" or "D0"
    - is_prompt (bool): True for prompt, False for non-prompt
    - cent_min (int): minimum centrality
    - cent_max (int): maximum centrality

    Returns:
    - TH1: data histogram
    """
    is_pp = cent_min is None and cent_max is None
    if is_pp:
        name = "data/Data_pp5TeV_Run2.root"
    else:
        name = "data/Data_PbPb5TeV_Run2.root"
    with ROOT.TFile.Open(name) as f:
        dir_particle = f.Get(particle)
        dir_type = dir_particle.Get("Prompt" if is_prompt else "NonPrompt")
        if is_pp:
            h_stat = dir_type.Get("h_stat")
            g_syst = dir_type.Get("g_syst")
            h_combined = dir_type.Get("h_combined")
        else:
            dir_cent = dir_type.Get(f"{cent_min}{cent_max}")
            h_stat = dir_cent.Get("h_stat")
            g_syst = dir_cent.Get("g_syst")
            h_combined = dir_cent.Get("h_combined")

    h_stat_plussyst = h_stat.Clone("h_stat_plussyst")
    h_stat_minussyst = h_stat.Clone("h_stat_minussyst")
    for ipt in range(1, h_stat.GetNbinsX()+1):
        cent = h_stat.GetBinContent(ipt)
        sys_high = g_syst.GetErrorYhigh(ipt-1)
        sys_low = g_syst.GetErrorYlow(ipt-1)
        h_stat_plussyst.SetBinContent(ipt, cent + sys_high)
        h_stat_minussyst.SetBinContent(ipt, cent - sys_low)
        
    return h_stat, g_syst, h_combined, h_stat_plussyst, h_stat_minussyst


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

def power_law(pt, pars):
    """
    Power law function for pT extrapolation

    Parameters:
    - pt (list): 1-D list for pT (variable)
    - pars (list): list of parameters

    Returns:
    - double: power law function evaluated for a given pt

    """
    return pars[0] * pt[0] / ((1 + (pt[0] / pars[1])**pars[3])**pars[2])

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

def fit_histogram(config, hist, initial_pars):
    """
    Fit function to histogram

    Parameters:
    - config (dict): configuration dictionary
    - hist (TH1): histogram to fit
    - initial_pars (list): initial parameters for the fit

    Returns:
    - TF1: fitted Tsallis function
    """
    if config["fit_function"] == "tsallis":
        func = ROOT.TF1("func", tsallis, 0, 16, len(initial_pars))

        # func.SetParLimits(0, 1.e-10, 1.e20)
        func.SetParLimits(1, 1.5, 15.)
        func.SetParLimits(2, 1.e-10, 15.)
        func.FixParameter(3, initial_pars[3])  # Fix mass parameter
    elif config["fit_function"] == "powerlaw":
        func = ROOT.TF1("func", power_law, 0, 16, len(initial_pars))

        func.SetParLimits(0, .1, 1000)
        func.SetParLimits(1, 0.5, 10.)
        func.SetParLimits(2, 0.5, 10.)
        func.SetParLimits(3, 0.5, 20.0)


    for i, par in enumerate(initial_pars):
        func.SetParameter(i, par)
    
    fit_res = hist.Fit(func, "IME0S")
    
    return func, fit_res

def simultaneous_fit(h_corrected_yields_ds, h_corrected_yields_dplus, fitter_tsallis):
    """
    Simultaneous fit of Ds and D+ with Tsallis

    Parameters:
    - config (dict): configuration dictionary
    - hist (TH1): histogram to fit
    - initial_pars (list): initial parameters for the fit

    Returns:
    - TF1: fitted Tsallis functions
    """
    pt_min = 0.
    pt_max_ds = h_corrected_yields_ds.GetXaxis().GetBinUpEdge(h_corrected_yields_ds.GetNbinsX())
    pt_max_dplus = h_corrected_yields_dplus.GetXaxis().GetBinUpEdge(h_corrected_yields_dplus.GetNbinsX())
    f_corrected_yields_ds = ROOT.TF1("f_corrected_yields_ds_tsallissim", tsallis, pt_min, pt_max_ds, 4)
    f_corrected_yields_dplus = ROOT.TF1("f_corrected_yields_dplus_tsallissim", tsallis, pt_min, pt_max_dplus, 4)

    wf_corrected_yields_ds = ROOT.Math.WrappedMultiTF1(f_corrected_yields_ds, 1)
    wf_corrected_yields_dplus = ROOT.Math.WrappedMultiTF1(f_corrected_yields_dplus, 1)

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
    f_corrected_yields_ds.SetFitResult(result_tsallis, pars_tsallis_ds)
    f_corrected_yields_dplus.SetFitResult(result_tsallis, pars_tsallis_dplus)

    return f_corrected_yields_ds, f_corrected_yields_dplus, result_tsallis


def get_confidence_intervals():
    h = ROOT.TH1D("h_confidence_intervals", "", 1000, 0, 16)
    virt_fitter = ROOT.TVirtualFitter.GetFitter()
    virt_fitter.GetConfidenceIntervals(h, 0.683)
    return h

def get_yield_1_2(tsallis_func, fit_res):
    """
    Get yields in 1-2 GeV/c from fitted Tsallis function

    Parameters:
    - tsallis_func (TF1): fitted Tsallis function
    - fit_res (TFitResultPtr): fit result

    Returns:
    - double: yield in 1-2 GeV/c
    - double: uncertainty on yield in 1-2 GeV/c
    """
    yield_1_2 = tsallis_func.Integral(1.0, 2.0)
    unc_yield_1_2 = tsallis_func.IntegralError(1.0, 2.0, fit_res.GetParams(), fit_res.GetCovarianceMatrix().GetMatrixArray())
    h_yield_1_2 = ROOT.TH1D("h_yield_1_2", "Yield 1-2 GeV/c", 1, 1, 2)
    h_yield_1_2.SetBinContent(1, yield_1_2)
    h_yield_1_2.SetBinError(1, unc_yield_1_2)

    return h_yield_1_2

def get_yields_from_ci(config, h_errors):
    h_errors_upper = h_errors.Clone("h_errors_upper")
    h_errors_lower = h_errors.Clone("h_errors_lower")
    for i in range(1, h_errors.GetNbinsX() + 1):
        h_errors_upper.SetBinContent(i, h_errors.GetBinContent(i) + h_errors.GetBinError(i))
        h_errors_lower.SetBinContent(i, h_errors.GetBinContent(i) - h_errors.GetBinError(i))

    pt_mins = config["pt"]["min"]
    pt_maxs = config["pt"]["max"]

    g_yield_1_2_from_ci = ROOT.TGraphAsymmErrors(len(pt_mins))

    for i_pt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        bin_center = (pt_min + pt_max) / 2.0

        y = h_errors.Integral(h_errors.FindBin(pt_min), h_errors.FindBin(pt_max), "width") / (pt_max - pt_min)
        unc_y_upper = (h_errors_upper.Integral(h_errors_upper.FindBin(pt_min), h_errors_upper.FindBin(pt_max), "width") / (pt_max - pt_min)) - y
        unc_y_lower = y - (h_errors_lower.Integral(h_errors_lower.FindBin(pt_min), h_errors_lower.FindBin(pt_max), "width") / (pt_max - pt_min))

        g_yield_1_2_from_ci.SetPoint(i_pt, bin_center, y)
        g_yield_1_2_from_ci.SetPointError(i_pt, (pt_max - pt_min) / 2.0, (pt_max - pt_min) / 2.0, unc_y_lower, unc_y_upper)

    return g_yield_1_2_from_ci

def get_efficiencies_ds(eff_infile, cent_min, cent_max):
    """
    Get Ds prompt and non-prompt efficiencies from ROOT file

    Parameters:
    - eff_infile (str): input ROOT file with efficiencies
    - cent_min (int): minimum centrality
    - cent_max (int): maximum centrality

    Returns:
    - TH1: Ds prompt efficiency
    - TH1: Ds non-prompt efficiency
    """
    with ROOT.TFile.Open(eff_infile) as f:
        h_eff_ds_prompt = f.Get(f"eff_DsPrompt_cent_{cent_min}_{cent_max}")
        h_eff_ds_nonprompt = f.Get(f"eff_DsNonPrompt_cent_{cent_min}_{cent_max}")

        h_eff_ds_prompt.SetDirectory(0)
        h_eff_ds_nonprompt.SetDirectory(0)

    return h_eff_ds_prompt, h_eff_ds_nonprompt

def get_efficiencies_dplus(eff_infile, cent_min, cent_max):
    """
    Get D+ prompt and non-prompt efficiencies from ROOT file

    Parameters:
    - eff_infile (str): input ROOT file with efficiencies
    - cent_min (int): minimum centrality
    - cent_max (int): maximum centrality

    Returns:
    - TH1: Ds prompt efficiency
    - TH1: Ds non-prompt efficiency
    """
    with ROOT.TFile.Open(eff_infile) as f:
        h_eff_ds_prompt = f.Get(f"eff_DplusPrompt_cent_{cent_min}_{cent_max}")
        h_eff_ds_nonprompt = f.Get(f"eff_DplusNonPrompt_cent_{cent_min}_{cent_max}")

        h_eff_ds_prompt.SetDirectory(0)
        h_eff_ds_nonprompt.SetDirectory(0)

    return h_eff_ds_prompt, h_eff_ds_nonprompt

def get_fraction_corrected(g_yield_prompt, g_yield_nonprompt, g_raa_prompt, g_raa_nonprompt):
    """
    Get corrected prompt fraction

    Parameters:
    - g_yield_prompt (TGraphAsymmErrors): prompt yield in 1-2 GeV/c
    - g_yield_nonprompt (TGraphAsymmErrors): non-prompt yield in 1-2 GeV/c

    Returns:
    - TGraphAsymmErrors: fraction of Ds prompt over total Ds
    """

    g_frac = ROOT.TGraphAsymmErrors(g_yield_prompt.GetN())
    g_frac.SetName("g_frac_corrected")
    g_frac.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    g_frac.GetYaxis().SetTitle("Prompt fraction (corrected)")
    for i_pt in range(g_frac.GetN()):
        bin_center = g_yield_prompt.GetX()[i_pt]
        n_prompt = g_yield_prompt.GetY()[i_pt]
        n_nonprompt = g_yield_nonprompt.GetY()[i_pt]
        if g_raa_prompt is not None and g_raa_nonprompt is not None:
            n_prompt = n_prompt * g_raa_prompt.Eval(bin_center)
            n_nonprompt = n_nonprompt * g_raa_nonprompt.Eval(bin_center)

        frac_prompt = n_prompt / (n_prompt + n_nonprompt)

        # Uncertainty propagation
        err_prompt = g_yield_prompt.GetErrorY(i_pt)
        err_nonprompt = g_yield_nonprompt.GetErrorY(i_pt)

        frac_err = err_prompt**2 * n_nonprompt**2 / n_prompt**2
        frac_err += err_nonprompt**2
        frac_err = np.sqrt(frac_err) * frac_prompt / (n_prompt + n_nonprompt)

        g_frac.SetPoint(i_pt, bin_center, frac_prompt)
        g_frac.SetPointError(i_pt, g_yield_prompt.GetErrorX(i_pt), g_yield_prompt.GetErrorX(i_pt), frac_err, frac_err)

    return g_frac

def get_fraction_raw(h_eff_prompt, h_eff_nonprompt, g_yield_prompt, g_yield_nonprompt, g_raa_prompt, g_raa_nonprompt):
    """
    Get raw prompt fraction

    Parameters:
    - h_eff_prompt (TH1): prompt efficiency
    - h_eff_nonprompt (TH1): non-prompt efficiency
    - g_yield_prompt (TGraphAsymmErrors): prompt yield in 1-2 GeV/c
    - g_yield_nonprompt (TGraphAsymmErrors): non-prompt yield in 1-2 GeV/c

    Returns:
    - TGraphAsymmErrors: fraction of Ds prompt over total Ds
    """

    g_frac = ROOT.TGraphAsymmErrors(g_yield_prompt.GetN())
    g_frac.SetName("g_frac_corrected")
    g_frac.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    g_frac.GetYaxis().SetTitle("Prompt fraction (raw)")
    for i_pt in range(g_frac.GetN()):
        bin_center = g_yield_prompt.GetX()[i_pt]
        n_prompt = g_yield_prompt.GetY()[i_pt]
        n_nonprompt = g_yield_nonprompt.GetY()[i_pt]
        if g_raa_prompt is not None and g_raa_nonprompt is not None:
            n_prompt = n_prompt * g_raa_prompt.Eval(bin_center)
            n_nonprompt = n_nonprompt * g_raa_nonprompt.Eval(bin_center)

        eff_prompt = h_eff_prompt.GetBinContent(h_eff_prompt.FindBin(bin_center))
        eff_nonprompt = h_eff_nonprompt.GetBinContent(h_eff_nonprompt.FindBin(bin_center))
        raw_yield_prompt = n_prompt * eff_prompt
        raw_yield_nonprompt = n_nonprompt * eff_nonprompt

        frac_prompt = raw_yield_prompt / (raw_yield_prompt + raw_yield_nonprompt)

        # Uncertainty propagation
        err_prompt = g_yield_prompt.GetErrorY(i_pt)
        err_nonprompt = g_yield_nonprompt.GetErrorY(i_pt)
        err_eff_prompt = h_eff_prompt.GetBinError(h_eff_prompt.FindBin(bin_center))
        err_eff_nonprompt = h_eff_nonprompt.GetBinError(h_eff_nonprompt.FindBin(bin_center))

        eff_prompt_term = (eff_nonprompt*n_nonprompt / eff_prompt)**2
        n_prompt_term = (eff_nonprompt*n_nonprompt / n_prompt)**2
        eff_nonprompt_term = n_nonprompt**2
        n_nonprompt_term = eff_nonprompt**2

        frac_err = err_eff_prompt**2 * eff_prompt_term
        frac_err += err_prompt**2 * n_prompt_term
        frac_err += err_eff_nonprompt**2 * eff_nonprompt_term
        frac_err += err_nonprompt**2 * n_nonprompt_term
        frac_err = np.sqrt(frac_err) * frac_prompt / (raw_yield_prompt + raw_yield_nonprompt)

        g_frac.SetPoint(i_pt, bin_center, frac_prompt)
        g_frac.SetPointError(i_pt, g_yield_prompt.GetErrorX(i_pt), g_yield_prompt.GetErrorX(i_pt), frac_err, frac_err)

    return g_frac


def get_spectrum(config, particle, is_prompt, cent_min=None, cent_max=None):
    """
    Get pT spectrum from ROOT file

    Parameters:
    - particle (str): "Ds" or "D0"
    - is_prompt (bool): True for prompt, False for non-prompt
    - cent_min (int): minimum centrality
    - cent_max (int): maximum centrality

    Returns:
    - TH1: pT spectrum histogram
    """
    with open(config["inputs"]["cutset"], "r") as f:
        cutset = yaml.safe_load(f)

    h_stat, g_syst, h_combined, h_stat_plussyst, h_stat_minussys = get_data(particle, is_prompt, cent_min, cent_max)
    mass = 1.968 if particle == "Ds" else 1.865
    if config["fit_function"] == "tsallis":
        initial_pars = [h_stat.Integral(), 5, 2, mass]  # Example initial parameters
    elif config["fit_function"] == "powerlaw":
        initial_pars = [h_stat.Integral()/2, 2.3, 1.70, 3.0]  # Example initial parameters
    func, fit_res = fit_histogram(config, h_stat, initial_pars)
    h_errors = get_confidence_intervals()
    if cent_min is not None and cent_max is not None:
        name = f"{config["fit_function"]}_{particle}_{'Prompt' if is_prompt else 'NonPrompt'}_{cent_min}_{cent_max}"
    else:
        name = f"{config["fit_function"]}_{particle}_{'Prompt' if is_prompt else 'NonPrompt'}"
    h_errors.SetTitle(name)

    g_yield_prompt_from_ci = get_yields_from_ci(cutset, h_errors)

    with ROOT.TFile.Open(config["output"], "UPDATE") as output_file:
        h_stat.Write(f"{particle}{'Prompt' if is_prompt else 'NonPrompt'}_Yield_Stat")
        g_syst.Write(f"{particle}{'Prompt' if is_prompt else 'NonPrompt'}_Yield_Syst")
        h_combined.Write(f"{particle}{'Prompt' if is_prompt else 'NonPrompt'}_Yield_Combined")
        func.Write(f"{config["fit_function"]}_{particle}{'Prompt' if is_prompt else 'NonPrompt'}")
        h_errors.Write(f"{config["fit_function"]}_{particle}{'Prompt' if is_prompt else 'NonPrompt'}_ConfidenceIntervals")
        g_yield_prompt_from_ci.Write(f"{particle}{'Prompt' if is_prompt else 'NonPrompt'}_Yields_FromCI")

    return g_yield_prompt_from_ci

def get_spectrum_simfit(is_prompt, cent_min=None, cent_max=None):
    """
    Get pT spectrum from ROOT file

    Parameters:
    - is_prompt (bool): True for prompt, False for non-prompt
    - cent_min (int): minimum centrality
    - cent_max (int): maximum centrality

    Returns:
    - TH1: pT spectrum histogram
    """
    with open(config["inputs"]["cutset"], "r") as f:
        cutset = yaml.safe_load(f)

    h_stat_ds, g_syst_ds, h_combined_ds, h_stat_plussys_ds, h_stat_minussys_ds = get_data("Ds", is_prompt, cent_min, cent_max)
    h_stat_dplus, g_syst_dplus, h_combined_dplus, h_stat_plussys_dplus, h_stat_minussys_dplus = get_data("Dplus", is_prompt, cent_min, cent_max)

    fitter_tsallis = ROOT.Fit.Fitter()

    func_plussyst_ds, func_plussyst_dplus, result_plussyst_tsallis = simultaneous_fit(h_stat_plussys_ds, h_stat_plussys_dplus, fitter_tsallis)
    cov_mat_plussyst = ROOT.TMatrixDSym(6)
    result_plussyst_tsallis.GetCovarianceMatrix(cov_mat_plussyst)
    params_plussyst = result_plussyst_tsallis.GetParams()

    func_minussyst_ds, func_minussyst_dplus, result_minussyst_tsallis = simultaneous_fit(h_stat_minussys_ds, h_stat_minussys_dplus, fitter_tsallis)
    cov_mat_minussyst = ROOT.TMatrixDSym(6)
    result_minussyst_tsallis.GetCovarianceMatrix(cov_mat_minussyst)
    params_minussyst = result_minussyst_tsallis.GetParams()

    func_ds, func_dplus, result_tsallis = simultaneous_fit(h_stat_ds, h_stat_dplus, fitter_tsallis)
    cov_mat = ROOT.TMatrixDSym(6)
    result_tsallis.GetCovarianceMatrix(cov_mat)
    params = result_tsallis.GetParams()

    pt_mins = config["pt"]["min"]
    pt_maxs = config["pt"]["max"]
    g_yield_prompt_from_ci_ds = ROOT.TGraphAsymmErrors(len(pt_mins))
    g_yield_prompt_from_ci_dplus = ROOT.TGraphAsymmErrors(len(pt_mins))

    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        pt_cent = (pt_min + pt_max) / 2
        delta_pt = (pt_max - pt_min)
        pt_unc = delta_pt / 2
        int_ds = func_ds.Integral(pt_min, pt_max) / delta_pt
        int_dplus = func_dplus.Integral(pt_min, pt_max) / delta_pt
        g_yield_prompt_from_ci_ds.SetPoint(ipt, pt_cent, int_ds)
        g_yield_prompt_from_ci_dplus.SetPoint(ipt, pt_cent, int_dplus)
        unc_ds = func_ds.IntegralError(pt_min, pt_max, params, cov_mat.GetMatrixArray()) / delta_pt
        unc_dplus = func_dplus.IntegralError(pt_min, pt_max, params, cov_mat.GetMatrixArray()) / delta_pt
        unc_syshigh_ds = func_plussyst_ds.Integral(pt_min, pt_max) / delta_pt - int_ds
        unc_syslow_ds = int_ds - func_minussyst_ds.Integral(pt_min, pt_max) / delta_pt
        unc_syshigh_dplus = func_plussyst_dplus.Integral(pt_min, pt_max) / delta_pt - int_dplus
        unc_syslow_dplus = int_dplus - func_minussyst_dplus.Integral(pt_min, pt_max) / delta_pt
        unc_totlow_ds = np.sqrt(unc_ds**2 + unc_syslow_ds**2)
        unc_tothigh_ds = np.sqrt(unc_ds**2 + unc_syshigh_ds**2)
        unc_totlow_dplus = np.sqrt(unc_dplus**2 + unc_syslow_dplus**2)
        unc_tothigh_dplus = np.sqrt(unc_dplus**2 + unc_syshigh_dplus**2)
        g_yield_prompt_from_ci_ds.SetPointError(ipt, pt_unc, pt_unc, unc_totlow_ds, unc_tothigh_ds)
        g_yield_prompt_from_ci_dplus.SetPointError(ipt, pt_unc, pt_unc, unc_totlow_dplus, unc_tothigh_dplus)

    with ROOT.TFile.Open(config["output"], "UPDATE") as output_file:
        h_stat_ds.Write(f"Ds{'Prompt' if is_prompt else 'NonPrompt'}_Yield_Stat")
        g_syst_ds.Write(f"Ds{'Prompt' if is_prompt else 'NonPrompt'}_Yield_Syst")
        h_combined_ds.Write(f"Ds{'Prompt' if is_prompt else 'NonPrompt'}_Yield_Combined")
        func_ds.Write(f"Tsallissim_Ds{'Prompt' if is_prompt else 'NonPrompt'}")
        g_yield_prompt_from_ci_ds.Write(f"Ds{'Prompt' if is_prompt else 'NonPrompt'}_Yields_FromCI")
        h_stat_dplus.Write(f"Dplus{'Prompt' if is_prompt else 'NonPrompt'}_Yield_Stat")
        g_syst_dplus.Write(f"Dplus{'Prompt' if is_prompt else 'NonPrompt'}_Yield_Syst")
        h_combined_dplus.Write(f"Dplus{'Prompt' if is_prompt else 'NonPrompt'}_Yield_Combined")
        func_dplus.Write(f"Tsallissim_Dplus{'Prompt' if is_prompt else 'NonPrompt'}")
        g_yield_prompt_from_ci_dplus.Write(f"Dplus{'Prompt' if is_prompt else 'NonPrompt'}_Yields_FromCI")

    func_ds.SetLineColor(ROOT.kRed-2)
    func_dplus.SetLineColor(ROOT.kAzure-3)
    h_combined_ds.SetMarkerStyle(ROOT.kFullCircle)
    h_combined_ds.SetLineColor(ROOT.kRed)
    h_combined_ds.SetMarkerColor(ROOT.kRed)
    h_combined_dplus.SetMarkerStyle(ROOT.kFullCircle)
    h_combined_dplus.SetLineColor(ROOT.kAzure-5)
    h_combined_dplus.SetMarkerColor(ROOT.kAzure-5)
    c = ROOT.TCanvas("c", "c", 800, 600)
    c.SetLogy()
    c.DrawFrame(0, 1.e-4, 24, 100, " ;#it{p}_{T} (GeV/#it{c});d#it{#sigma}/d#it{p}_{T} #mub(GeV/#it{c})^{-1}")

    leg = ROOT.TLegend(0.55, 0.65, 0.85, 0.85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_combined_ds, f"{'Prompt' if is_prompt else 'NonPrompt'} D_{{s}}^{{+}} data", "pl")
    leg.AddEntry(func_ds, f"{'Prompt' if is_prompt else 'NonPrompt'} D_{{s}}^{{+}} fit", "l")
    leg.AddEntry(h_combined_dplus, f"{'Prompt' if is_prompt else 'NonPrompt'} D^{{+}} data", "pl")
    leg.AddEntry(func_dplus, f"{'Prompt' if is_prompt else 'NonPrompt'} D^{{+}} fit", "l")
    leg.Draw()

    func_ds.Draw("same")
    h_combined_ds.Draw("E same")
    func_dplus.Draw("same")
    h_combined_dplus.Draw("E same")

    c.SaveAs(config["output"].replace(".root", f"_{'Prompt' if is_prompt else 'NonPrompt'}.pdf"))

    return g_yield_prompt_from_ci_ds, g_yield_prompt_from_ci_dplus 

def fc_from_pb_pb_data(config):
    # Load cutset from config
    with open(config["inputs"]["cutset"], "r") as f:
        cutset = yaml.safe_load(f)

    cent_min = config["cent"]["min"]
    cent_max = config["cent"]["max"]

    # Create output ROOT file
    with ROOT.TFile.Open(config["output"], "RECREATE") as output_file:
        pass

    print("Fitting Ds prompt spectrum...")
    g_yield_prompt_from_ci = get_spectrum(config, "Ds", True, cent_min, cent_max)
    print("Fitting Ds non-prompt spectrum...")
    g_yield_nonprompt_from_ci = get_spectrum(config, "Ds", False, cent_min, cent_max)

    # h_yield_prompt_1_2 = get_yield_1_2(tsallis_prompt, fit_res_prompt)
    # h_yield_nonprompt_1_2 = get_yield_1_2(tsallis_nonprompt, fit_res_nonprompt)
    # h_yield_prompt_1_2.SetName("h_yield_prompt_1_2")
    # h_yield_nonprompt_1_2.SetName("h_yield_nonprompt_1_2")

    h_eff_ds_prompt, h_eff_ds_nonprompt = get_efficiencies_ds(config["inputs"]["efficiency"], cent_min, cent_max)

    g_frac_corr_ds = get_fraction_corrected(g_yield_prompt_from_ci, g_yield_nonprompt_from_ci)
    g_frac_raw_ds = get_fraction_raw(h_eff_ds_prompt, h_eff_ds_nonprompt, g_yield_prompt_from_ci, g_yield_nonprompt_from_ci)
    
    print("Fitting D0 prompt spectrum...")
    g_yield_prompt_from_ci = get_spectrum(config, "D0", True, cent_min, cent_max)
    print("Fitting D0 non-prompt spectrum...")
    g_yield_nonprompt_from_ci = get_spectrum(config, "D0", False, cent_min, cent_max)

    # h_yield_prompt_1_2 = get_yield_1_2(tsallis_prompt, fit_res_prompt)
    # h_yield_nonprompt_1_2 = get_yield_1_2(tsallis_nonprompt, fit_res_nonprompt)
    # h_yield_prompt_1_2.SetName("h_yield_prompt_1_2")
    # h_yield_nonprompt_1_2.SetName("h_yield_nonprompt_1_2")

    h_eff_dp_prompt, h_eff_dp_nonprompt = get_efficiencies_dplus(config["inputs"]["efficiency"], cent_min, cent_max)

    g_frac_corr_dp = get_fraction_corrected(g_yield_prompt_from_ci, g_yield_nonprompt_from_ci)
    g_frac_raw_dp = get_fraction_raw(h_eff_dp_prompt, h_eff_dp_nonprompt, g_yield_prompt_from_ci, g_yield_nonprompt_from_ci)

    output_file = ROOT.TFile.Open(config["output"], "UPDATE")
    output_file.cd()
    g_frac_corr_ds.Write("DsPrompt_CorrFraction")
    g_frac_raw_ds.Write("DsPrompt_RawFraction")
    g_frac_corr_dp.Write("D0Prompt_CorrFraction")
    g_frac_raw_dp.Write("D0Prompt_RawFraction")
    output_file.Close()

def fc_from_pp_data(config):
    # Load cutset from config
    with open(config["inputs"]["cutset"], "r") as f:
        cutset = yaml.safe_load(f)

    cent_min = config["cent"]["min"]
    cent_max = config["cent"]["max"]

    # Create output ROOT file
    with ROOT.TFile.Open(config["output"], "RECREATE") as output_file:
        pass

    if config["fit_function"] != "tsallissim":
        print("Fitting Ds prompt fraction...")
        g_yield_prompt_ds_from_ci = get_spectrum(config, "Ds", True)
        print("Fitting Ds non-prompt fraction...")
        g_yield_nonprompt_ds_from_ci = get_spectrum(config, "Ds", False)

        # h_yield_prompt_1_2 = get_yield_1_2(tsallis_prompt, fit_res_prompt)
        # h_yield_nonprompt_1_2 = get_yield_1_2(tsallis_nonprompt, fit_res_nonprompt)
        # h_yield_prompt_1_2.SetName("h_yield_prompt_1_2")
        # h_yield_nonprompt_1_2.SetName("h_yield_nonprompt_1_2")

        print("Fitting D+ prompt spectrum...")
        g_yield_prompt_dplus_from_ci = get_spectrum(config, "Dplus", True)
        print("Fitting D+ non-prompt spectrum...")
        g_yield_nonprompt_dplus_from_ci = get_spectrum(config, "Dplus", False)

        # h_yield_prompt_1_2 = get_yield_1_2(tsallis_prompt, fit_res_prompt)
        # h_yield_nonprompt_1_2 = get_yield_1_2(tsallis_nonprompt, fit_res_nonprompt)
        # h_yield_prompt_1_2.SetName("h_yield_prompt_1_2")
        # h_yield_nonprompt_1_2.SetName("h_yield_nonprompt_1_2")
    else:
        print("Fitting prompt spectra...")
        g_yield_prompt_ds_from_ci, g_yield_prompt_dplus_from_ci = get_spectrum_simfit(True)
        print("Fitting nonprompt spectra...")
        g_yield_nonprompt_ds_from_ci, g_yield_nonprompt_dplus_from_ci = get_spectrum_simfit(False)

    h_eff_ds_prompt, h_eff_ds_nonprompt = get_efficiencies_ds(config["inputs"]["efficiency"], cent_min, cent_max)
    h_eff_dp_prompt, h_eff_dp_nonprompt = get_efficiencies_dplus(config["inputs"]["efficiency"], cent_min, cent_max)

    g_raa_dp_p, g_raa_ds_p, g_raa_dp_np, g_raa_ds_np = None, None, None, None
    if config["raa"]["use"]:
        model = config["raa"]["model"]
        infile_raa = ROOT.TFile.Open(config["raa"]["file"])
        g_raa_dp_p = infile_raa.Get(f"gr_d_raa_{model}_{cent_min}_{cent_max}")
        g_raa_ds_p = infile_raa.Get(f"gr_ds_raa_{model}_{cent_min}_{cent_max}")
        g_raa_dp_np = infile_raa.Get(f"gr_npd_raa_{model}_{cent_min}_{cent_max}")
        g_raa_ds_np = infile_raa.Get(f"gr_npds_raa_{model}_{cent_min}_{cent_max}")
        infile_raa.Close()

    g_frac_corr_ds = get_fraction_corrected(g_yield_prompt_ds_from_ci, g_yield_nonprompt_ds_from_ci, g_raa_ds_p, g_raa_ds_np)
    g_frac_raw_ds = get_fraction_raw(h_eff_ds_prompt, h_eff_ds_nonprompt, g_yield_prompt_ds_from_ci, g_yield_nonprompt_ds_from_ci, g_raa_ds_p, g_raa_ds_np)

    g_frac_corr_dp = get_fraction_corrected(g_yield_prompt_dplus_from_ci, g_yield_nonprompt_dplus_from_ci, g_raa_dp_p, g_raa_dp_np)
    g_frac_raw_dp = get_fraction_raw(h_eff_dp_prompt, h_eff_dp_nonprompt, g_yield_prompt_dplus_from_ci, g_yield_nonprompt_dplus_from_ci, g_raa_dp_p, g_raa_dp_np)

    output_file = ROOT.TFile.Open(config["output"], "UPDATE")
    output_file.cd()
    g_frac_corr_ds.Write("DsPrompt_CorrFraction")
    g_frac_raw_ds.Write("DsPrompt_RawFraction")
    g_frac_corr_dp.Write("DplusPrompt_CorrFraction")
    g_frac_raw_dp.Write("DplusPrompt_RawFraction")
    output_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Ds prompt and non-prompt yields with Tsallis fit")
    parser.add_argument("config_file", type=str, help="Input configuration file")
    args = parser.parse_args()

    # Load configuration
    with open(args.config_file, "r") as f:
        config = yaml.safe_load(f)

    if config["data"] == "pp":
        fc_from_pp_data(config)
    else:
        fc_from_pb_pb_data(config)
