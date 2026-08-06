import argparse
import sys
import yaml
import numpy as np
import ROOT


def tsallis(pt, pars):
    """
    Tsallis function for pT extrapolation

    Parameters:
    - pt (list): 1-D list for pT (variable)
    - pars (list): list of parameters

    Returns:
    - double: Tsallis function evaluated for a given pt

    """
    c_par = (pars[1] - 1) * (pars[1] - 2) / (pars[1] * pars[2] * (pars[1] * pars[2] + 2 * (pars[1] - 2) * pars[3]))
    mt_par = np.sqrt(pars[3]*pars[3] + pt[0]*pt[0]) # pars[3] is the mass of the particle

    return pars[0] * c_par * pt[0] * np.power((1 + (mt_par - pars[3]) / pars[1] / pars[2]), -pars[1])


def power_law(pt, pars):
    """
    power law function for pT extrapolation

    Parameters:
    - pt (list): 1-D list for pT (variable)
    - pars (list): list of parameters

    Returns:
    - double: power law function evaluated for a given pt

    """

    return pars[0] * pt[0] / np.power((1 + np.power(pt[0] / pars[1], pars[3])), pars[2])

def bw_integrand(r, pars):
    """
    Blast-wave integral function for pT extrapolation

    Parameters:
    - r (list): radius
    - pars (list): list of parameters

    Returns:
    - double: BW integrand function evaluated for a given r

    """
    mt = pars[0]
    pt = pars[1]
    beta_max = pars[2]
    temp_1 = 1./pars[3]
    n = pars[4]

    beta = beta_max * np.power(r[0], n)
    if beta > 0.9999999999999999:
        beta = 0.9999999999999999
    rho = np.arctanh(beta)
    arg_I0 = pt * np.sinh(rho) * temp_1
    if arg_I0 > 700.:
        arg_I0 = 700.
    arg_K1 = mt * np.cosh(rho) * temp_1
 
    return r[0] * mt * ROOT.TMath.BesselI0(arg_I0) * ROOT.TMath.BesselK1(arg_K1)

# we need to have this global
func_bw_integrand = ROOT.TF1("func_bw_integrand", bw_integrand, 0., 1., 5)

def bw_func(pt, pars):
    """
    Blast-wave function for pT extrapolation

    Parameters:
    - r (list): radius
    - pars (list): list of parameters

    Returns:
    - double: BW integrand function evaluated for a given r

    """

    norm = pars[0]
    mass = pars[1]
    mt = np.sqrt(pt[0]**2 + mass**2)
    beta_max = pars[2]
    t = pars[3]
    n = pars[4]

    func_bw_integrand.SetParameters(mt, pt[0], beta_max, t, n)

    return norm * pt[0] * func_bw_integrand.Integral(0, 1) 

# definition of shared parameters for Ds, D+ combined Tsallis and BW fits
pars_tsallis_ds = np.array(
    [
        0, # integral Ds
        1, # n
        2, # T
        3  # m(Ds)
    ],
    dtype=np.int32
)

pars_tsallis_dplus = np.array(
    [
        4, # integral D+
        1, # n
        2, # T
        5  # m(D+)
    ],
    dtype=np.int32,
)

SYST_SOURCES = ["rawyield", "np_frac", "bdt_eff", "pt_shape"]

pars_bw_ds = np.array(
    [
        0, # integral Ds
        1, # m(Ds)
        2, # beta_max
        3, # T
        4  # n
    ],
    dtype=np.int32
)

pars_bw_dplus = np.array(
    [
        5, # integral D+
        6, # m(D+)
        2, # beta_max
        3, # T
        4  # n
    ],
    dtype=np.int32
)

class GlobalChi2Tsallis(object):
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

class GlobalChi2BW(object):
    """
    Class for combined Chi2 in case of simultaneous Ds, D+ Blast-wave fit
    """
    def __init__(self, f1, f2):
        self._f1 = f1
        self._f2 = f2

    def __call__(self, par):
        # parameter vector is first background (in common 1 and 2) and then is
        # signal (only in 2)

        # the zero-copy way to get a numpy array from a double *
        par_arr = np.frombuffer(par, dtype=np.float64, count=7)

        p1 = par_arr[pars_bw_ds]
        p2 = par_arr[pars_bw_dplus]

        return self._f1(p1) + self._f2(p2)

def compute_extrap_factor(yield_tot, yield_vis, unc_yield_tot, unc_yield_vis):
    """
    """
    extrap_factor = yield_tot/yield_vis
    rho = np.sqrt(yield_vis/yield_tot)
    if rho > 1:
        rho = 1./rho

    unc_extrap_factor = np.sqrt(unc_yield_tot**2/yield_vis**2 + unc_yield_vis**2 * yield_tot**2/yield_vis**4 - 2*rho*unc_yield_tot*unc_yield_vis*yield_tot/yield_vis**3)

    return extrap_factor, unc_extrap_factor


colors = [
    ROOT.kRed+2,
    ROOT.kRed+1,
    ROOT.kOrange+10,
    ROOT.kOrange-3,
    ROOT.kSpring-5,
    ROOT.kGreen+2,
    ROOT.kCyan+2,
    ROOT.kAzure+4,
    ROOT.kBlue+2
]

colors_models = [
    ROOT.kRed+1,
    ROOT.kAzure+4,
    ROOT.kGreen+2
]

def set_obj_style(obj, color, marker=ROOT.kFullCircle, fillstyle=-1, alpha=1.1):
    """
    Set style
    """
    if isinstance(obj, ROOT.TH1):
        obj.SetDirectory(0)
    obj.SetLineWidth(2)
    obj.SetLineColor(color)
    if isinstance(obj, ROOT.TF1):
        return
    obj.SetMarkerColor(color)
    obj.SetMarkerStyle(marker)
    if fillstyle > -1:
        obj.SetFillStyle(fillstyle)
        obj.SetFillColorAlpha(color, alpha)


def extapolate(config_file_name):
    """
    Main function for extrapolation
    """

    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetPadTopMargin(0.035)
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetTitleSize(0.045, "xy")
    ROOT.gStyle.SetLabelSize(0.045, "xy")
    ROOT.gStyle.SetTitleOffset(1.3, "x")
    ROOT.gStyle.SetOptStat(0)

    with open(config_file_name, 'r', encoding="utf8") as file:
        config = yaml.safe_load(file)

    with open(config["extrap"]["sys_unc_db"], 'r', encoding="utf8") as file:
        db_sys = yaml.safe_load(file)

    hist_corry_ds, hist_corry_dp = {}, {}
    func_tsallis_corry_ds, func_powlaw_corry_ds, func_tsallissim_corry_ds, func_bwsim_corry_ds = ({} for _ in range(4))
    func_tsallis_corry_dp, func_powlaw_corry_dp, func_tsallissim_corry_dp, func_bwsim_corry_dp = ({} for _ in range(4))
    graph_extrapfactor_tsallis_ds, graph_extrapfactor_tsallis_dp, graph_extrapfactor_tsallis_ds_over_dp = (ROOT.TGraphErrors(1) for _ in range(3))
    graph_extrapfactor_powlaw_ds, graph_extrapfactor_powlaw_dp, graph_extrapfactor_powlaw_ds_over_dp = (ROOT.TGraphErrors(1) for _ in range(3))
    graph_extrapfactor_tsallissim_ds, graph_extrapfactor_tsallissim_dp, graph_extrapfactor_tsallissim_ds_over_dp = (ROOT.TGraphErrors(1) for _ in range(3))
    graph_extrapfactor_bwsim_ds, graph_extrapfactor_bwsim_dp, graph_extrapfactor_bwsim_ds_over_dp = (ROOT.TGraphErrors(1) for _ in range(3))
    graph_extrapfactor_tsallis_ds.SetNameTitle("graph_extrapfactor_tsallis_ds", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus}) = D_{s}^{#plus}(#it{p}_{T} integrated)/D_{s}^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_tsallis_dp.SetNameTitle("graph_extrapfactor_tsallis_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D^{#plus}) = D^{#plus}(#it{p}_{T} integrated)/D^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_tsallis_ds_over_dp.SetNameTitle("graph_extrapfactor_tsallis_ds_over_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus})/#alpha(D^{#plus})")
    graph_extrapfactor_powlaw_ds.SetNameTitle("graph_extrapfactor_powlaw_ds", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus}) = D_{s}^{#plus}(#it{p}_{T} integrated)/D_{s}^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_powlaw_dp.SetNameTitle("graph_extrapfactor_powlaw_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D^{#plus}) = D^{#plus}(#it{p}_{T} integrated)/D^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_powlaw_ds_over_dp.SetNameTitle("graph_extrapfactor_powlaw_ds_over_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus})/#alpha(D^{#plus})")
    graph_extrapfactor_tsallissim_ds.SetNameTitle("graph_extrapfactor_tsallissim_ds", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus}) = D_{s}^{#plus}(#it{p}_{T} integrated)/D_{s}^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_tsallissim_dp.SetNameTitle("graph_extrapfactor_tsallissim_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D^{#plus}) = D^{#plus}(#it{p}_{T} integrated)/D^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_tsallissim_ds_over_dp.SetNameTitle("graph_extrapfactor_tsallissim_ds_over_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus})/#alpha(D^{#plus})")
    graph_extrapfactor_bwsim_ds.SetNameTitle("graph_extrapfactor_bwsim_ds", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus}) = D_{s}^{#plus}(#it{p}_{T} integrated)/D_{s}^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_bwsim_dp.SetNameTitle("graph_extrapfactor_bwsim_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D^{#plus}) = D^{#plus}(#it{p}_{T} integrated)/D^{#plus}(#it{p}_{T} > 1 GeV/#it{c})")
    graph_extrapfactor_bwsim_ds_over_dp.SetNameTitle("graph_extrapfactor_bwsim_ds_over_dp", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5}; #alpha(D_{s}^{#plus})/#alpha(D^{#plus})")
    set_obj_style(graph_extrapfactor_tsallis_ds, ROOT.kBlack)
    set_obj_style(graph_extrapfactor_tsallis_dp, ROOT.kBlack)
    set_obj_style(graph_extrapfactor_tsallis_ds_over_dp, ROOT.kBlack)
    set_obj_style(graph_extrapfactor_powlaw_ds, ROOT.kGray+1)
    set_obj_style(graph_extrapfactor_powlaw_dp, ROOT.kGray+1)
    set_obj_style(graph_extrapfactor_powlaw_ds_over_dp, ROOT.kGray+1)
    set_obj_style(graph_extrapfactor_tsallissim_ds, ROOT.kGray+2)
    set_obj_style(graph_extrapfactor_tsallissim_dp, ROOT.kGray+2)
    set_obj_style(graph_extrapfactor_tsallissim_ds_over_dp, ROOT.kGray+2)
    set_obj_style(graph_extrapfactor_bwsim_ds, ROOT.kGray+3)
    set_obj_style(graph_extrapfactor_bwsim_dp, ROOT.kGray+3)
    set_obj_style(graph_extrapfactor_bwsim_ds_over_dp, ROOT.kGray+3)

    leg_cent = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
    leg_cent.SetTextSize(0.04)
    leg_cent.SetBorderSize(0)
    leg_cent.SetFillStyle(0)

    # extrapolation factors from models
    leg_model = ROOT.TLegend(0.2, 0.2, 0.5, 0.65)
    leg_model.SetTextSize(0.04)
    leg_model.SetBorderSize(0)
    leg_model.SetFillStyle(0)
    leg_model.AddEntry(graph_extrapfactor_powlaw_ds, "Power law")
    leg_model.AddEntry(graph_extrapfactor_tsallis_ds, "Tsallis")
    leg_model.AddEntry(graph_extrapfactor_tsallissim_ds, "Tsallis (sim fit)")
    leg_model.AddEntry(graph_extrapfactor_bwsim_ds, "BW (sim fit)")
    graph_model_ds, graph_model_dp, graph_model_ratio = ([] for _ in range(3))
    for imodel, input_model in enumerate(config["models"]["inputs"]):
        infile_model = ROOT.TFile.Open(input_model)
        graph_model_ds.append(infile_model.Get("graph_extrapfactor_ds"))
        graph_model_dp.append(infile_model.Get("graph_extrapfactor_dp"))
        graph_model_ratio.append(infile_model.Get("graph_extrapfactor_ds_over_dp"))
        set_obj_style(graph_model_ds[imodel], colors_models[imodel])
        set_obj_style(graph_model_dp[imodel], colors_models[imodel])
        set_obj_style(graph_model_ratio[imodel], colors_models[imodel])
        for imult in range(graph_model_ds[imodel].GetN()):
            graph_model_ds[imodel].SetPointError(imult, 0., graph_model_ds[imodel].GetErrorYlow(imult))
            graph_model_dp[imodel].SetPointError(imult, 0., graph_model_ds[imodel].GetErrorYlow(imult))
            graph_model_ratio[imodel].SetPointError(imult, 0., graph_model_ds[imodel].GetErrorYlow(imult))
            if imult == 0:
                leg_model.AddEntry(graph_model_ds[imodel], config["models"]["labels"][imodel], "lp")

    gstat_ds_over_dp_vis = ROOT.TGraphAsymmErrors(1)
    gstat_ds_over_dp_ptint = ROOT.TGraphAsymmErrors(1)
    gstat_ds_over_dp_vis.SetNameTitle("gstat_ds_over_dp_vis", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
    gstat_ds_over_dp_ptint.SetNameTitle("gstat_ds_over_dp_ptint", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
    gsyst_ds_over_dp_vis = ROOT.TGraphAsymmErrors(1)
    gsyst_ds_over_dp_ptint = ROOT.TGraphAsymmErrors(1)
    gsyst_ds_over_dp_vis.SetNameTitle("gsyst_ds_over_dp_vis", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
    gsyst_ds_over_dp_ptint.SetNameTitle("gsyst_ds_over_dp_ptint", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
    gsyst_br_ds_over_dp_vis = ROOT.TGraphAsymmErrors(1)
    gsyst_br_ds_over_dp_ptint = ROOT.TGraphAsymmErrors(1)
    gsyst_br_ds_over_dp_vis.SetNameTitle("gsyst_br_ds_over_dp_vis", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
    gsyst_br_ds_over_dp_ptint.SetNameTitle("gsyst_br_ds_over_dp_ptint", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
    gsyst_extrap_ds_over_dp_ptint = ROOT.TGraphAsymmErrors(1)
    gsyst_extrap_ds_over_dp_ptint.SetNameTitle("gsyst_extrap_ds_over_dp_ptint", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
    gsyst_tot_nobr_ds_over_dp_ptint = ROOT.TGraphAsymmErrors(1)
    gsyst_tot_nobr_ds_over_dp_ptint.SetNameTitle("gsyst_tot_nobr_ds_over_dp_ptint", ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
    set_obj_style(gstat_ds_over_dp_vis, ROOT.kBlack)
    set_obj_style(gsyst_ds_over_dp_vis, ROOT.kBlack, 0, 0)
    set_obj_style(gsyst_br_ds_over_dp_vis, ROOT.kBlack, 0, 0)
    set_obj_style(gstat_ds_over_dp_ptint, ROOT.kRed+1, ROOT.kOpenCircle)
    set_obj_style(gsyst_ds_over_dp_ptint, ROOT.kRed+1, 0, 0)
    set_obj_style(gsyst_br_ds_over_dp_ptint, ROOT.kRed+1, 0, 0)
    set_obj_style(gsyst_extrap_ds_over_dp_ptint, ROOT.kRed+1, 0, 1000, 0.3)
    set_obj_style(gsyst_tot_nobr_ds_over_dp_ptint, ROOT.kRed+1, 0, 0)

    gsyst_source_ds_over_dp_vis, gsyst_source_ds_over_dp_ptint = ({} for _ in range(2))
    for source in SYST_SOURCES:
        gsyst_source_ds_over_dp_vis[source] = ROOT.TGraphAsymmErrors(1)
        gsyst_source_ds_over_dp_vis[source].SetNameTitle(
            f"gsyst_{source}_ds_over_dp_vis",
            ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
        set_obj_style(gsyst_source_ds_over_dp_vis[source], ROOT.kBlack, 0, 0)
        gsyst_source_ds_over_dp_ptint[source] = ROOT.TGraphAsymmErrors(1)
        gsyst_source_ds_over_dp_ptint[source].SetNameTitle(
            f"gsyst_{source}_ds_over_dp_ptint",
            ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
        set_obj_style(gsyst_source_ds_over_dp_ptint[source], ROOT.kRed+1, 0, 0)

    max_mult = 1.
    for icent, cent in enumerate(config["inputs"]):
        cent_min = cent.split(sep="_")[1]
        cent_max = cent.split(sep="_")[2]
        cent_width = (int(cent_max) - int(cent_min))
        infile = ROOT.TFile.Open(config["inputs"][cent]["file"])
        hist_rawy_ds = infile.Get(f"centrality_{cent_min}_{cent_max}/h_raw_yields_ds_{cent_min}_{cent_max}")
        hist_rawy_dp = infile.Get(f"centrality_{cent_min}_{cent_max}/h_raw_yields_dplus_{cent_min}_{cent_max}")
        hist_eff_ds = infile.Get(f"centrality_{cent_min}_{cent_max}/eff_DsPrompt_cent_{cent_min}_{cent_max}")
        hist_eff_dp = infile.Get(f"centrality_{cent_min}_{cent_max}/eff_DplusPrompt_cent_{cent_min}_{cent_max}")
        hist_frac_ds = infile.Get(f"centrality_{cent_min}_{cent_max}/hRawFracPrompt_cent_{cent_min}_{cent_max};1")
        if not isinstance(hist_frac_ds, ROOT.TH1):
            hist_frac_ds = infile.Get(f"centrality_{cent_min}_{cent_max}/hRawFracDsPrompt_cent_{cent_min}_{cent_max}")
        hist_frac_dp = infile.Get(f"centrality_{cent_min}_{cent_max}/hRawFracPrompt_cent_{cent_min}_{cent_max};2")
        if not isinstance(hist_frac_dp, ROOT.TH1):
            hist_frac_dp = infile.Get(f"centrality_{cent_min}_{cent_max}/hRawFracDplusPrompt_cent_{cent_min}_{cent_max}")
        hist_corry_ds[cent] = hist_rawy_ds.Clone(f"h_corr_yields_ds_{cent_min}_{cent_max}")
        hist_corry_dp[cent] = hist_rawy_dp.Clone(f"h_corr_yields_dplus_{cent_min}_{cent_max}")
        set_obj_style(hist_corry_ds[cent], colors[icent])
        set_obj_style(hist_corry_dp[cent], colors[icent])
        hist_corry_ds[cent].Scale(1./config["br"]["ds_to_phipi_to_kkpi"]/cent_width, "width")
        hist_corry_dp[cent].Scale(1./config["br"]["dplus_to_phipi_to_kkpi"]/cent_width, "width")
        hist_corry_ds[cent].Divide(hist_eff_ds)
        hist_corry_dp[cent].Divide(hist_eff_dp)
        hist_corry_ds[cent].Multiply(hist_frac_ds)
        hist_corry_dp[cent].Multiply(hist_frac_dp)

        leg_cent.AddEntry(hist_corry_ds[cent], f"{cent_min}#minus{cent_max}%", "lp")

        # get multiplicity value
        g_ratio_4mult = infile.Get("g_ratio_dndeta_20_40") # 2-4 GeV/c bin always there for sure
        mult_cent = g_ratio_4mult.GetPointX(config["inputs"][cent]["mult_bin"])
        mult_unc_low = db_sys["systematics"][f"cent_{cent_min}_{cent_max}"]["multiplicity"][0] * mult_cent
        mult_unc_high = db_sys["systematics"][f"cent_{cent_min}_{cent_max}"]["multiplicity"][1] * mult_cent
        if max_mult < mult_cent:
            max_mult = mult_cent

        # pt extrapolation
        pt_min = 0.
        pt_max_ds = config["extrap"]["fit"]["pt_max_ds"]
        pt_max_dplus = config["extrap"]["fit"]["pt_max_dplus"]

        func_tsallis_corry_ds[cent] = ROOT.TF1(f"func_tsallis_corry_ds_{cent_min}_{cent_max}", tsallis, pt_min, pt_max_ds, 4)
        func_tsallis_corry_dp[cent] = ROOT.TF1(f"func_tsallis_corry_dp_{cent_min}_{cent_max}", tsallis, pt_min, pt_max_dplus, 4)
        func_powlaw_corry_ds[cent] = ROOT.TF1(f"func_powlaw_corry_ds_{cent_min}_{cent_max}", power_law, pt_min, pt_max_ds, 4)
        func_powlaw_corry_dp[cent] = ROOT.TF1(f"func_powlaw_corry_dp_{cent_min}_{cent_max}", power_law, pt_min, pt_max_dplus, 4)
        func_tsallissim_corry_ds[cent] = ROOT.TF1(f"func_tsallissim_corry_ds_{cent_min}_{cent_max}", tsallis, pt_min, pt_max_ds, 4)
        func_tsallissim_corry_dp[cent] = ROOT.TF1(f"func_tsallissim_corry_dp_{cent_min}_{cent_max}", tsallis, pt_min, pt_max_dplus, 4)
        func_bwsim_corry_ds[cent] = ROOT.TF1(f"func_bwsim_corry_ds_{cent_min}_{cent_max}", bw_func, pt_min, config["extrap"]["fit"]["pt_max_bw"], 5)
        func_bwsim_corry_dp[cent] = ROOT.TF1(f"func_bwsim_corry_dp_{cent_min}_{cent_max}", bw_func, pt_min, config["extrap"]["fit"]["pt_max_bw"], 5)
        set_obj_style(func_tsallis_corry_ds[cent], colors[icent])
        set_obj_style(func_tsallis_corry_dp[cent], colors[icent])
        set_obj_style(func_powlaw_corry_ds[cent], colors[icent])
        set_obj_style(func_powlaw_corry_dp[cent], colors[icent])
        set_obj_style(func_bwsim_corry_ds[cent], colors[icent])
        set_obj_style(func_bwsim_corry_dp[cent], colors[icent])

        # indpendent fits with power law
        func_powlaw_corry_ds[cent].SetParameters(hist_corry_ds[cent].Integral(), 2.6, 2.8, 2.)
        func_powlaw_corry_dp[cent].SetParameters(hist_corry_dp[cent].Integral(), 2.6, 2.8, 2.)
        fitres_ds = hist_corry_ds[cent].Fit(func_powlaw_corry_ds[cent], "LIMEQ0S")
        fitres_dp = hist_corry_dp[cent].Fit(func_powlaw_corry_dp[cent], "LIMEQ0S")
        # We need to retrieve the covariance matrices to compute the uncertainties on the integrals
        cov_mat_ds = fitres_ds.GetCovarianceMatrix()
        cov_mat_dp = fitres_dp.GetCovarianceMatrix()
        yield_ds = func_powlaw_corry_ds[cent].Integral(0., 1000)
        yield_vis_ds = func_powlaw_corry_ds[cent].Integral(1., 1000)
        yield_dp = func_powlaw_corry_dp[cent].Integral(0., 1000)
        yield_vis_dp = func_powlaw_corry_dp[cent].Integral(1., 1000)
        unc_yield_ds = func_powlaw_corry_ds[cent].IntegralError(0., 1000, fitres_ds.GetParams(), cov_mat_ds.GetMatrixArray())
        unc_yield_vis_ds = func_powlaw_corry_ds[cent].IntegralError(1., 1000, fitres_ds.GetParams(), cov_mat_ds.GetMatrixArray())
        unc_yield_dp = func_powlaw_corry_dp[cent].IntegralError(0., 1000, fitres_dp.GetParams(), cov_mat_dp.GetMatrixArray())
        unc_yield_vis_dp = func_powlaw_corry_dp[cent].IntegralError(1., 1000, fitres_dp.GetParams(), cov_mat_dp.GetMatrixArray())
        extrap_factor_powlaw_ds = compute_extrap_factor(yield_ds, yield_vis_ds, unc_yield_ds, unc_yield_vis_ds)
        extrap_factor_powlaw_dp = compute_extrap_factor(yield_dp, yield_vis_dp, unc_yield_dp, unc_yield_vis_dp)

        # indpendent fits with Tsallis
        func_tsallis_corry_ds[cent].SetParameters(hist_corry_ds[cent].Integral(), 10, 0.5, 1.968)
        func_tsallis_corry_dp[cent].SetParameters(hist_corry_dp[cent].Integral(), 10, 0.5, 1.870)
        func_tsallis_corry_ds[cent].FixParameter(3, 1.968)
        func_tsallis_corry_dp[cent].FixParameter(3, 1.870)
        # func_tsallis_corry_ds[cent].SetParLimits(1, 1.5, 15.)
        # func_tsallis_corry_dp[cent].SetParLimits(1, 1.5, 15.)
        fitres_ds = hist_corry_ds[cent].Fit(func_tsallis_corry_ds[cent], "IMEQ0S")
        fitres_dp = hist_corry_dp[cent].Fit(func_tsallis_corry_dp[cent], "IMEQ0S")
        # We need to retrieve the covariance matrices to compute the uncertainties on the integrals
        cov_mat_ds = fitres_ds.GetCovarianceMatrix()
        cov_mat_dp = fitres_dp.GetCovarianceMatrix()
        yield_ds = func_tsallis_corry_ds[cent].Integral(0., 1000)
        yield_vis_ds = func_tsallis_corry_ds[cent].Integral(1., 1000)
        yield_dp = func_tsallis_corry_dp[cent].Integral(0., 1000)
        yield_vis_dp = func_tsallis_corry_dp[cent].Integral(1., 1000)
        unc_yield_ds = func_tsallis_corry_ds[cent].IntegralError(0., 1000, fitres_ds.GetParams(), cov_mat_ds.GetMatrixArray())
        unc_yield_vis_ds = func_tsallis_corry_ds[cent].IntegralError(1., 1000, fitres_ds.GetParams(), cov_mat_ds.GetMatrixArray())
        unc_yield_dp = func_tsallis_corry_dp[cent].IntegralError(0., 1000, fitres_dp.GetParams(), cov_mat_dp.GetMatrixArray())
        unc_yield_vis_dp = func_tsallis_corry_dp[cent].IntegralError(1., 1000, fitres_dp.GetParams(), cov_mat_dp.GetMatrixArray())
        extrap_factor_tsallis_ds = compute_extrap_factor(yield_ds, yield_vis_ds, unc_yield_ds, unc_yield_vis_ds)
        extrap_factor_tsallis_dp = compute_extrap_factor(yield_dp, yield_vis_dp, unc_yield_dp, unc_yield_vis_dp)

        # let's also do a simultaneous fit with Tsallis and BW
        wf_corrected_yields_ds = ROOT.Math.WrappedMultiTF1(func_tsallissim_corry_ds[cent], 1)
        wf_corrected_yields_dplus = ROOT.Math.WrappedMultiTF1(func_tsallissim_corry_dp[cent], 1)

        opt = ROOT.Fit.DataOptions()
        range_ds = ROOT.Fit.DataRange()
        range_ds.SetRange(pt_min, pt_max_ds)
        data_ds = ROOT.Fit.BinData(opt, range_ds)
        ROOT.Fit.FillData(data_ds, hist_corry_ds[cent])
        range_dplus = ROOT.Fit.DataRange()
        range_dplus.SetRange(pt_min, pt_max_ds)
        data_dplus = ROOT.Fit.BinData(opt, range_dplus)
        ROOT.Fit.FillData(data_dplus, hist_corry_dp[cent])

        chi2_ds = ROOT.Fit.Chi2Function(data_ds, wf_corrected_yields_ds)
        chi2_dplus = ROOT.Fit.Chi2Function(data_dplus, wf_corrected_yields_dplus)
        global_chi2 = GlobalChi2Tsallis(chi2_ds, chi2_dplus)
        fitter_tsallis = ROOT.Fit.Fitter()

        pars_all = np.array(
            [hist_corry_ds[cent].Integral(), 10, 0.5, 1.968, hist_corry_dp[cent].Integral(), 1.870]
        )
        fitter_tsallis.Config().SetParamsSettings(6, pars_all)
        # fix the masses
        fitter_tsallis.Config().ParSettings(3).Fix()
        fitter_tsallis.Config().ParSettings(5).Fix()
        # set limits
        #fitter_tsallis.Config().ParSettings(1).SetLimits(1.5, 1000.)
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
        func_tsallissim_corry_ds[cent].SetFitResult(result_tsallis, pars_tsallis_ds)
        func_tsallissim_corry_dp[cent].SetFitResult(result_tsallis, pars_tsallis_dplus)

        # We need to retrieve the covariance matrices to compute the uncertainties on the integrals
        cov_mat_sim = ROOT.TMatrixDSym(6)
        result_tsallis.GetCovarianceMatrix(cov_mat_sim)
        yield_ds = func_tsallissim_corry_ds[cent].Integral(0., 1000)
        yield_vis_ds = func_tsallissim_corry_ds[cent].Integral(1., 1000)
        yield_dp = func_tsallissim_corry_dp[cent].Integral(0., 1000)
        yield_vis_dp = func_tsallissim_corry_dp[cent].Integral(1., 1000)
        unc_yield_ds = func_tsallissim_corry_ds[cent].IntegralError(0., 1000, result_tsallis.GetParams(), cov_mat_sim.GetMatrixArray())
        unc_yield_vis_ds = func_tsallissim_corry_ds[cent].IntegralError(1., 1000, result_tsallis.GetParams(), cov_mat_sim.GetMatrixArray())
        unc_yield_dp = func_tsallissim_corry_dp[cent].IntegralError(0., 1000, result_tsallis.GetParams(), cov_mat_sim.GetMatrixArray())
        unc_yield_vis_dp = func_tsallissim_corry_dp[cent].IntegralError(1., 1000, result_tsallis.GetParams(), cov_mat_sim.GetMatrixArray())
        extrap_factor_tsallissim_ds = compute_extrap_factor(yield_ds, yield_vis_ds, unc_yield_ds, unc_yield_vis_ds)
        extrap_factor_tsallissim_dp = compute_extrap_factor(yield_dp, yield_vis_dp, unc_yield_dp, unc_yield_vis_dp)

        wf_bw_corrected_yields_ds = ROOT.Math.WrappedMultiTF1(func_bwsim_corry_ds[cent], 1)
        wf_bw_corrected_yields_dplus = ROOT.Math.WrappedMultiTF1(func_bwsim_corry_dp[cent], 1)

        opt = ROOT.Fit.DataOptions()
        range_ds = ROOT.Fit.DataRange()
        range_ds.SetRange(pt_min, config["extrap"]["fit"]["pt_max_bw"])
        data_ds = ROOT.Fit.BinData(opt, range_ds)
        ROOT.Fit.FillData(data_ds, hist_corry_ds[cent])
        range_dplus = ROOT.Fit.DataRange()
        range_dplus.SetRange(pt_min, config["extrap"]["fit"]["pt_max_bw"])
        data_dplus = ROOT.Fit.BinData(opt, range_dplus)
        ROOT.Fit.FillData(data_dplus, hist_corry_dp[cent])

        chi2_ds = ROOT.Fit.Chi2Function(data_ds, wf_bw_corrected_yields_ds)
        chi2_dplus = ROOT.Fit.Chi2Function(data_dplus, wf_bw_corrected_yields_dplus)
        global_chi2 = GlobalChi2BW(chi2_ds, chi2_dplus)
        fitter_bw = ROOT.Fit.Fitter()

        pars_all_bw = np.array(
            [func_powlaw_corry_ds[cent].Integral(0, 4)*1.e5, 1.968, config["extrap"]["fit"]["fix_bw_pars_from_lf"]["beta"][icent],
            config["extrap"]["fit"]["fix_bw_pars_from_lf"]["T"][icent], config["extrap"]["fit"]["fix_bw_pars_from_lf"]["n"][icent],
            func_powlaw_corry_dp[cent].Integral(0, 4)*1.e5, 1.870]
        )
        fitter_bw.Config().SetParamsSettings(7, pars_all_bw)
        # set limits
        fitter_bw.Config().ParSettings(0).SetLimits(0., func_powlaw_corry_ds[cent].Integral(0, 4)*1.e8)
        fitter_bw.Config().ParSettings(2).SetLimits(0.3, 0.99)
        fitter_bw.Config().ParSettings(3).SetLimits(0., 0.5)
        fitter_bw.Config().ParSettings(4).SetLimits(0.1, 5.)
        fitter_bw.Config().ParSettings(5).SetLimits(0., func_powlaw_corry_dp[cent].Integral(0, 4)*1.e8)
        fitter_bw.Config().MinimizerOptions().SetPrintLevel(0)
        fitter_bw.Config().SetMinimizer("Minuit2", "Migrad")
        # fix parameters
        fitter_bw.Config().ParSettings(1).Fix()
        fitter_bw.Config().ParSettings(6).Fix()
        if config["extrap"]["fit"]["fix_bw_pars_from_lf"]["enable"]:
            fitter_bw.Config().ParSettings(2).Fix()
            fitter_bw.Config().ParSettings(3).Fix()
            fitter_bw.Config().ParSettings(4).Fix()

        # we can't pass the Python object global_chi2 directly to FitFCN.
        # It needs to be wrapped in a ROOT::Math::Functor.
        global_chi2_functor = ROOT.Math.Functor(global_chi2, 7)

        # fit FCN function
        # (specify optionally data size and flag to indicate that is a chi2 fit)
        fitter_bw.FitFCN(global_chi2_functor, 0, data_ds.Size() + data_dplus.Size(), 1)
        result_bw = fitter_bw.Result()
        result_bw.Print(ROOT.std.cout)
        func_bwsim_corry_ds[cent].SetFitResult(result_bw, pars_bw_ds)
        func_bwsim_corry_dp[cent].SetFitResult(result_bw, pars_bw_dplus)

        # func_bwsim_corry_ds[cent].SetParameter(0, func_bwsim_corry_ds[cent].GetParameter(0) * hist_corry_ds[cent].Integral(1, 4) / func_bwsim_corry_ds[cent].Integral(1, 4))
        # func_bwsim_corry_dp[cent].SetParameter(0, func_bwsim_corry_dp[cent].GetParameter(0) * hist_corry_dp[cent].Integral(1, 4) / func_bwsim_corry_dp[cent].Integral(1, 4))

        # We need to retrieve the covariance matrices to compute the uncertainties on the integrals
        cov_mat_sim = ROOT.TMatrixDSym(6)
        result_bw.GetCovarianceMatrix(cov_mat_sim)
        yield_ds = func_bwsim_corry_ds[cent].Integral(0., 1000)
        yield_vis_ds = func_bwsim_corry_ds[cent].Integral(1., 1000)
        yield_dp = func_bwsim_corry_dp[cent].Integral(0., 1000)
        yield_vis_dp = func_bwsim_corry_dp[cent].Integral(1., 1000)
        unc_yield_ds = func_bwsim_corry_ds[cent].IntegralError(0., 1000, result_bw.GetParams(), cov_mat_sim.GetMatrixArray())
        unc_yield_vis_ds = func_bwsim_corry_ds[cent].IntegralError(1., 1000, result_bw.GetParams(), cov_mat_sim.GetMatrixArray())
        unc_yield_dp = func_bwsim_corry_dp[cent].IntegralError(0., 1000, result_bw.GetParams(), cov_mat_sim.GetMatrixArray())
        unc_yield_vis_dp = func_bwsim_corry_dp[cent].IntegralError(1., 1000, result_bw.GetParams(), cov_mat_sim.GetMatrixArray())
        extrap_factor_bwsim_ds = compute_extrap_factor(yield_ds, yield_vis_ds, unc_yield_ds, unc_yield_vis_ds)
        extrap_factor_bwsim_dp = compute_extrap_factor(yield_dp, yield_vis_dp, unc_yield_dp, unc_yield_vis_dp)

        # fill graphs of extrapolation factors
        graph_extrapfactor_tsallis_ds.SetPoint(icent, mult_cent, extrap_factor_tsallis_ds[0])
        graph_extrapfactor_tsallis_dp.SetPoint(icent, mult_cent, extrap_factor_tsallis_dp[0])
        graph_extrapfactor_powlaw_ds.SetPoint(icent, mult_cent, extrap_factor_powlaw_ds[0])
        graph_extrapfactor_powlaw_dp.SetPoint(icent, mult_cent, extrap_factor_powlaw_dp[0])
        graph_extrapfactor_tsallissim_dp.SetPoint(icent, mult_cent, extrap_factor_tsallissim_dp[0])
        graph_extrapfactor_tsallissim_ds.SetPoint(icent, mult_cent, extrap_factor_tsallissim_ds[0])
        graph_extrapfactor_bwsim_dp.SetPoint(icent, mult_cent, extrap_factor_bwsim_dp[0])
        graph_extrapfactor_bwsim_ds.SetPoint(icent, mult_cent, extrap_factor_bwsim_ds[0])

        graph_extrapfactor_tsallis_ds_over_dp.SetPoint(icent, mult_cent, extrap_factor_tsallis_ds[0]/extrap_factor_tsallis_dp[0])
        graph_extrapfactor_powlaw_ds_over_dp.SetPoint(icent, mult_cent, extrap_factor_powlaw_ds[0]/extrap_factor_powlaw_dp[0])
        graph_extrapfactor_tsallissim_ds_over_dp.SetPoint(icent, mult_cent, extrap_factor_tsallissim_ds[0]/extrap_factor_tsallissim_dp[0])
        graph_extrapfactor_bwsim_ds_over_dp.SetPoint(icent, mult_cent, extrap_factor_bwsim_ds[0]/extrap_factor_bwsim_dp[0])

        main_extrap_factor = 1.
        if config["extrap"]["main"] == "tsallis":
            main_extrap_factor = extrap_factor_tsallis_ds[0]/extrap_factor_tsallis_dp[0]
        elif config["extrap"]["main"] == "tsallis_sim":
            main_extrap_factor = extrap_factor_tsallissim_ds[0]/extrap_factor_tsallissim_dp[0]
        elif config["extrap"]["main"] == "powlaw":
            main_extrap_factor = extrap_factor_powlaw_ds[0]/extrap_factor_powlaw_dp[0]
        elif config["extrap"]["main"] == "bwsim":
            main_extrap_factor = extrap_factor_bwsim_ds[0]/extrap_factor_bwsim_dp[0]
        elif "model" in config["extrap"]["main"]:
            model_idx = config["extrap"]["main"].split(sep="model")[0]
            main_extrap_factor = graph_model_ds[model_idx].Eval(mult_cent) / graph_model_dp[model_idx].Eval(mult_cent)

        # compute pT-integrated ratio
        yield_meas_ds, yield_meas_dp, unc_yield_meas_ds, unc_yield_meas_dp, syst_yield_low_meas, syst_yield_high_meas = (0. for _ in range(6))
        syst_yield_low_by_source_vis, syst_yield_low_by_source_extrap, syst_yield_high_by_source_vis, syst_yield_high_by_source_extrap = (
            {"rawyield": 0., "np_frac": 0., "bdt_eff": 0., "pt_shape": 0.} for _ in range(4)
        )
        sys_cent = db_sys["systematics"][f"cent_{cent_min}_{cent_max}"]
        for ipt in range(1, hist_corry_ds[cent].GetNbinsX()+1):
            yield_meas_ds += hist_corry_ds[cent].GetBinContent(ipt) * hist_corry_ds[cent].GetBinWidth(ipt)
            unc_yield_meas_ds += hist_corry_ds[cent].GetBinError(ipt)**2 * hist_corry_ds[cent].GetBinWidth(ipt)**2
            yield_meas_dp += hist_corry_dp[cent].GetBinContent(ipt) * hist_corry_dp[cent].GetBinWidth(ipt)
            unc_yield_meas_dp += hist_corry_dp[cent].GetBinError(ipt)**2 * hist_corry_dp[cent].GetBinWidth(ipt)**2
            rel_sys_tot_low = np.sqrt(sys_cent["rawyield"][ipt-1]**2 + sys_cent["np_frac"][ipt-1][0]**2 + sys_cent["bdt_eff"][ipt-1]**2 + sys_cent["pt_shape"][ipt-1]**2)
            rel_sys_tot_high = np.sqrt(sys_cent["rawyield"][ipt-1]**2 + sys_cent["np_frac"][ipt-1][1]**2 + sys_cent["bdt_eff"][ipt-1]**2 + sys_cent["pt_shape"][ipt-1]**2)
            syst_yield_low_meas += hist_corry_ds[cent].GetBinContent(ipt) * rel_sys_tot_low * hist_corry_ds[cent].GetBinWidth(ipt)
            syst_yield_high_meas += hist_corry_ds[cent].GetBinContent(ipt) * rel_sys_tot_high * hist_corry_ds[cent].GetBinWidth(ipt)
            for source in syst_yield_low_by_source_vis:
                if source == "np_frac":
                    syst_yield_low_by_source_vis[source] += hist_corry_ds[cent].GetBinContent(ipt) * sys_cent[source][ipt-1][0] * hist_corry_ds[cent].GetBinWidth(ipt)
                    syst_yield_high_by_source_vis[source] += hist_corry_ds[cent].GetBinContent(ipt) * sys_cent[source][ipt-1][1] * hist_corry_ds[cent].GetBinWidth(ipt)
                else:
                    syst_yield_low_by_source_vis[source] += hist_corry_ds[cent].GetBinContent(ipt) * sys_cent[source][ipt-1] * hist_corry_ds[cent].GetBinWidth(ipt)
                    syst_yield_high_by_source_vis[source] += hist_corry_ds[cent].GetBinContent(ipt) * sys_cent[source][ipt-1] * hist_corry_ds[cent].GetBinWidth(ipt)
        unc_yield_meas_ds = np.sqrt(unc_yield_meas_ds)
        unc_yield_meas_dp = np.sqrt(unc_yield_meas_dp)

        ratio_ptint_vis = yield_meas_ds/yield_meas_dp
        unc_ratio_ptint_vis = np.sqrt(unc_yield_meas_ds**2/yield_meas_ds**2 + unc_yield_meas_dp**2/yield_meas_dp**2) * ratio_ptint_vis
        sys_ratio_low_ptint_vis = syst_yield_low_meas/yield_meas_ds * ratio_ptint_vis
        sys_ratio_high_ptint_vis = syst_yield_high_meas/yield_meas_ds * ratio_ptint_vis
        ratio_ptint_extrap = main_extrap_factor * ratio_ptint_vis
        unc_ratio_ptint_extrap = main_extrap_factor * unc_ratio_ptint_vis
        sys_ratio_low_ptint_extrap = main_extrap_factor * sys_ratio_low_ptint_vis
        sys_ratio_high_ptint_extrap = main_extrap_factor * sys_ratio_high_ptint_vis
        for source in syst_yield_low_by_source_vis:
            syst_yield_low_by_source_extrap[source] += syst_yield_low_by_source_vis[source] * main_extrap_factor
            syst_yield_high_by_source_extrap[source] += syst_yield_high_by_source_vis[source] * main_extrap_factor

        gstat_ds_over_dp_vis.SetPoint(icent, mult_cent, ratio_ptint_vis)
        gstat_ds_over_dp_vis.SetPointError(icent, 0., 0., unc_ratio_ptint_vis, unc_ratio_ptint_vis)
        gstat_ds_over_dp_ptint.SetPoint(icent, mult_cent, ratio_ptint_extrap)
        gstat_ds_over_dp_ptint.SetPointError(icent, 0., 0., unc_ratio_ptint_extrap, unc_ratio_ptint_extrap)

        gsyst_ds_over_dp_vis.SetPoint(icent, mult_cent, ratio_ptint_vis)
        gsyst_ds_over_dp_vis.SetPointError(icent, mult_unc_low, mult_unc_high, sys_ratio_low_ptint_vis, sys_ratio_high_ptint_vis)
        gsyst_ds_over_dp_ptint.SetPoint(icent, mult_cent, ratio_ptint_extrap)
        gsyst_ds_over_dp_ptint.SetPointError(icent, mult_unc_low, mult_unc_high, sys_ratio_low_ptint_extrap, sys_ratio_high_ptint_extrap)

        # same as above, but keeping the sources separate
        for source in SYST_SOURCES:
            sys_ratio_low_source_vis = syst_yield_low_by_source_vis[source]/yield_meas_ds * ratio_ptint_vis
            sys_ratio_high_source_vis = syst_yield_high_by_source_vis[source]/yield_meas_ds * ratio_ptint_vis
            sys_ratio_low_source_extrap = syst_yield_low_by_source_extrap[source]/yield_meas_ds * ratio_ptint_vis
            sys_ratio_high_source_extrap = syst_yield_high_by_source_extrap[source]/yield_meas_ds * ratio_ptint_vis
            gsyst_source_ds_over_dp_vis[source].SetPoint(icent, mult_cent, ratio_ptint_vis)
            gsyst_source_ds_over_dp_vis[source].SetPointError(icent, mult_unc_low, mult_unc_high,
                                                              sys_ratio_low_source_vis, sys_ratio_high_source_vis)
            gsyst_source_ds_over_dp_ptint[source].SetPoint(icent, mult_cent, ratio_ptint_extrap)
            gsyst_source_ds_over_dp_ptint[source].SetPointError(icent, mult_unc_low, mult_unc_high,
                                                                sys_ratio_low_source_extrap, sys_ratio_high_source_extrap)

        gsyst_extrap_ds_over_dp_ptint.SetPoint(icent, mult_cent, ratio_ptint_extrap)
        gsyst_extrap_ds_over_dp_ptint.SetPointError(icent, mult_unc_low*2, mult_unc_high*2, ratio_ptint_extrap * config["extrap"]["unc"][icent],
                                                    ratio_ptint_extrap * config["extrap"]["unc"][icent])

        sys_tot_low = np.sqrt(sys_ratio_low_ptint_extrap**2 + ratio_ptint_extrap**2 * config["extrap"]["unc"][icent]**2)
        sys_tot_high = np.sqrt(sys_ratio_high_ptint_extrap**2 + ratio_ptint_extrap**2 * config["extrap"]["unc"][icent]**2)
        gsyst_tot_nobr_ds_over_dp_ptint.SetPoint(icent, mult_cent, ratio_ptint_extrap)
        gsyst_tot_nobr_ds_over_dp_ptint.SetPointError(icent, mult_unc_low, mult_unc_high, sys_tot_low, sys_tot_high)

        # BR uncertainty overall normalisation, we put it separately
        br_unc_low = np.sqrt(db_sys["systematics"]["br"]["ds"]["high"]**2 + db_sys["systematics"]["br"]["dplus"]["low"]**2)
        br_unc_high = np.sqrt(db_sys["systematics"]["br"]["ds"]["low"]**2 + db_sys["systematics"]["br"]["dplus"]["high"]**2)
        gsyst_br_ds_over_dp_vis.SetPoint(icent, mult_cent, ratio_ptint_vis)
        gsyst_br_ds_over_dp_vis.SetPointError(icent, mult_unc_low, mult_unc_high, ratio_ptint_vis * br_unc_low, ratio_ptint_vis * br_unc_high)
        gsyst_br_ds_over_dp_ptint.SetPoint(icent, mult_cent, ratio_ptint_extrap)
        gsyst_br_ds_over_dp_ptint.SetPointError(icent, mult_unc_low, mult_unc_high, ratio_ptint_extrap * br_unc_low, ratio_ptint_extrap * br_unc_high)

    canv_corry = ROOT.TCanvas("canv_corry", "", 1500, 750)
    canv_corry.Divide(2, 1)
    canv_corry.cd(1).DrawFrame(0., 1.e2, 24., 1.e10, ";#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}(D_{s}^{#plus}) (#it{c} GeV^{#minus1})")
    canv_corry.cd(1).SetLogy()
    for cent in hist_corry_ds:
        hist_corry_ds[cent].Draw("same")
        func_tsallis_corry_ds[cent].Draw("same")
    leg_cent.Draw()
    canv_corry.cd(2).DrawFrame(0., 1.e2, 24., 1.e10, ";#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}(D^{#plus}) (#it{c} GeV^{#minus1})")
    canv_corry.cd(2).SetLogy()
    for cent in hist_corry_dp:
        hist_corry_dp[cent].Draw("same")
        func_tsallis_corry_dp[cent].Draw("same")
    leg_cent.Draw()
    canv_corry.SaveAs(f"corr_yields{config['output']['suffix']}.pdf")

    canv_corry_powlaw = ROOT.TCanvas("canv_corry_powlaw", "", 1500, 750)
    canv_corry_powlaw.Divide(2, 1)
    canv_corry_powlaw.cd(1).DrawFrame(0., 1.e2, 24., 1.e10, ";#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}(D_{s}^{#plus}) (#it{c} GeV^{#minus1})")
    canv_corry_powlaw.cd(1).SetLogy()
    for cent in hist_corry_ds:
        hist_corry_ds[cent].Draw("same")
        func_powlaw_corry_ds[cent].Draw("same")
    leg_cent.Draw()
    canv_corry_powlaw.cd(2).DrawFrame(0., 1.e2, 24., 1.e10, ";#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}(D^{#plus}) (#it{c} GeV^{#minus1})")
    canv_corry_powlaw.cd(2).SetLogy()
    for cent in hist_corry_dp:
        hist_corry_dp[cent].Draw("same")
        func_powlaw_corry_dp[cent].Draw("same")
    leg_cent.Draw()
    canv_corry_powlaw.SaveAs(f"corr_yields_powlaw{config['output']['suffix']}.pdf")

    canv_corry_bwsim = ROOT.TCanvas("canv_corry_bwsim", "", 1500, 750)
    canv_corry_bwsim.Divide(2, 1)
    canv_corry_bwsim.cd(1).DrawFrame(0., 1.e2, 24., 1.e10, ";#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}(D_{s}^{#plus}) (#it{c} GeV^{#minus1})")
    canv_corry_bwsim.cd(1).SetLogy()
    for cent in hist_corry_ds:
        hist_corry_ds[cent].Draw("same")
        func_bwsim_corry_ds[cent].Draw("same")
    leg_cent.Draw()
    canv_corry_bwsim.cd(2).DrawFrame(0., 1.e2, 24., 1.e10, ";#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}(D^{#plus}) (#it{c} GeV^{#minus1})")
    canv_corry_bwsim.cd(2).SetLogy()
    for cent in hist_corry_dp:
        hist_corry_dp[cent].Draw("same")
        func_bwsim_corry_dp[cent].Draw("same")
    leg_cent.Draw()
    canv_corry_bwsim.SaveAs(f"corr_yields_bwsim{config['output']['suffix']}.pdf")

    canv_extrap_fact = ROOT.TCanvas("canv_extrap_fact", "", 1500, 500)
    canv_extrap_fact.Divide(3, 1)
    canv_extrap_fact.cd(1)
    graph_extrapfactor_tsallis_ds.GetXaxis().SetNdivisions(505)
    graph_extrapfactor_tsallis_ds.GetYaxis().SetDecimals()
    graph_extrapfactor_tsallis_ds.GetYaxis().SetRangeUser(0., 2.)
    graph_extrapfactor_tsallis_ds.Draw("acpz")
    graph_extrapfactor_tsallissim_ds.Draw("cpz")
    graph_extrapfactor_bwsim_ds.Draw("cpz")
    graph_extrapfactor_powlaw_ds.Draw("cpz")
    for graph in graph_model_ds:
        graph.Draw("cpz")
    leg_model.Draw()
    canv_extrap_fact.cd(2)
    graph_extrapfactor_tsallis_dp.GetYaxis().SetRangeUser(0., 2.)
    graph_extrapfactor_tsallis_dp.GetYaxis().SetDecimals()
    graph_extrapfactor_tsallis_dp.GetXaxis().SetNdivisions(505)
    graph_extrapfactor_tsallis_dp.Draw("acpz")
    graph_extrapfactor_tsallissim_dp.Draw("cpz")
    graph_extrapfactor_bwsim_dp.Draw("cpz")
    graph_extrapfactor_powlaw_dp.Draw("cpz")
    for graph in graph_model_dp:
        graph.Draw("cpz")
    canv_extrap_fact.cd(3)
    graph_extrapfactor_tsallis_ds_over_dp.GetYaxis().SetRangeUser(0.8, 1.2)
    graph_extrapfactor_tsallis_ds_over_dp.GetYaxis().SetDecimals()
    graph_extrapfactor_tsallis_ds_over_dp.GetXaxis().SetNdivisions(505)
    graph_extrapfactor_tsallis_ds_over_dp.Draw("acpz")
    graph_extrapfactor_tsallissim_ds_over_dp.Draw("cpz")
    graph_extrapfactor_bwsim_ds_over_dp.Draw("cpz")
    graph_extrapfactor_powlaw_ds_over_dp.Draw("cpz")
    for graph in graph_model_ratio:
        graph.Draw("cpz")
    canv_extrap_fact.SaveAs(f"extrap_factor{config['output']['suffix']}.pdf")

    leg = ROOT.TLegend(0.2, 0.75, 0.4, 0.9)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(gstat_ds_over_dp_vis, "Measured #it{p}_{T}", "lp")
    leg.AddEntry(gstat_ds_over_dp_ptint, "Integrated #it{p}_{T}", "lp")
    leg_syst = ROOT.TLegend(0.2, 0.2, 0.4, 0.35)
    leg_syst.SetTextSize(0.04)
    leg_syst.SetBorderSize(0)
    leg_syst.SetFillStyle(0)
    leg_syst.AddEntry(gsyst_ds_over_dp_vis, "Data syst.", "f")
    leg_syst.AddEntry(gsyst_extrap_ds_over_dp_ptint, "Extrapolation syst.", "f")
    canv_ratio = ROOT.TCanvas("canv_ratio", "", 500, 500)
    frame = canv_ratio.DrawFrame(1., 0., max_mult*2, 1., ";d#it{N}_{ch}/d#it{#eta}_{|#it{#eta}|<0.5};#sigma(D_{s}^{#plus})/#sigma(D^{#plus})")
    canv_ratio.SetLogx()
    frame.GetYaxis().SetRangeUser(0., 1.)
    frame.GetYaxis().SetDecimals()
    gsyst_ds_over_dp_vis.Draw("2")
    gstat_ds_over_dp_vis.Draw("pz")
    gsyst_ds_over_dp_ptint.Draw("2")
    gsyst_extrap_ds_over_dp_ptint.Draw("2")
    gstat_ds_over_dp_ptint.Draw("pz")
    leg_syst.Draw()
    leg.Draw()
    canv_ratio.SaveAs(f"ds_over_dp_ratio_ptint{config['output']['suffix']}.pdf")

    outfile = ROOT.TFile(f"ds_over_dp_ratio_ptint{config['output']['suffix']}.root", "recreate")
    canv_ratio.Write()
    gsyst_ds_over_dp_vis.Write()
    gstat_ds_over_dp_vis.Write()
    gsyst_ds_over_dp_ptint.Write()
    gsyst_extrap_ds_over_dp_ptint.Write()
    gstat_ds_over_dp_ptint.Write()
    gsyst_tot_nobr_ds_over_dp_ptint.Write()
    gsyst_br_ds_over_dp_vis.Write()
    gsyst_br_ds_over_dp_ptint.Write()
    for source in SYST_SOURCES:
        gsyst_source_ds_over_dp_vis[source].Write()
        gsyst_source_ds_over_dp_ptint[source].Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute extrapolation factor for Ds and D+ vs imppar')
    parser.add_argument('--config_file', '-c', metavar='text', default="config_extrap_pbpb.yml", help='Input config file')
    args = parser.parse_args()

    extapolate(args.config_file)
