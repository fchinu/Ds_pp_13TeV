"""
Utils functions for model RAA interpolation
"""

import pandas as pd
import ROOT

def read_tamu(infile, graph_name, is_beauty=False):
    """
    Function to read TAMU files and convert predictions to a graph
    """

    df = pd.read_csv(infile, sep=" ", comment="#", header=0)
    graph = ROOT.TGraph(1)
    graph.SetNameTitle(graph_name, ";#it{p}_{T} (GeV/#it{c}); #it{R}_{AA}")
    if is_beauty:
        for ipt, (pt, raa) in enumerate(zip(df["PtCent"].to_numpy(), df["R_AA"].to_numpy())):
            graph.SetPoint(ipt, pt, raa)
            if ipt == len(df)-1 and pt < 50.:
                graph.SetPoint(ipt, 50., raa)
    else:
        for ipt, (pt, raa_min, raa_max) in enumerate(zip(df["PtCent"].to_numpy(),
                                            df["R_AA_min"].to_numpy(),
                                            df["R_AA_max"].to_numpy())):
            raa_cent = (raa_min + raa_max) / 2
            graph.SetPoint(ipt, pt, raa_cent)
            if ipt == len(df)-1 and pt < 50.:
                graph.SetPoint(ipt, 50., raa_cent)

    return graph

def read_lido(infile, graph_name):
    """
    Function to read LIDO files and convert predictions to a graph
    """

    df = pd.read_csv(infile, sep=" ", comment="#", header=0)
    graph = ROOT.TGraph(1)
    graph.SetNameTitle(graph_name, ";#it{p}_{T} (GeV/#it{c}); #it{R}_{AA}")
    for ipt, (pt, raa) in enumerate(zip(df["pT"].to_numpy(), df["Raa"].to_numpy())):
        graph.SetPoint(ipt, pt, raa)
        if ipt == len(df)-1 and pt < 50.:
            graph.SetPoint(ipt, 50., raa)

    return graph

def read_fonll(infile, hist_name, which="central"):
    """
    Function to read FONLL files and convert predictions to a histogram
    """

    columns = ["pt", "central", "min", "max", "min_sc", "max_sc", "min_mass", "max_mass", "min_pdf",
               "max_pdf", "fr=.5.5", "fr=22", "fr=21", "fr=12", "fr=1.5", "fr=.51"]
    df = pd.read_csv(infile, sep=" ", comment="#", names=columns)
    pt_cents = df["pt"].to_numpy()
    delta_pt = pt_cents[1]-pt_cents[0]
    pt_min = pt_cents[0]-delta_pt/2
    pt_max = pt_cents[-1]+delta_pt/2
    hist = ROOT.TH1D(hist_name, ";#it{p}_{T} (GeV/#it{c}); d#sigma/d#it{p}_{T} (a.u.)",
                    len(df), pt_min, pt_max)
    for ipt, xsec in enumerate(df[which].to_numpy()):
        hist.SetBinContent(ipt, xsec)
    return hist

def read_langevin(infile, graph_name):
    """
    Function to read Langevin files and convert predictions to a graph
    """

    df = pd.read_csv(infile, sep="\s+", comment="#", names=["pt", "raa", "uknown_1", "uknown_2"])
    graph = ROOT.TGraph(1)
    graph.SetNameTitle(graph_name, ";#it{p}_{T} (GeV/#it{c}); #it{R}_{AA}")
    for ipt, (pt, raa) in enumerate(zip(df["pt"].to_numpy(), df["raa"].to_numpy())):
        graph.SetPoint(ipt, pt, raa)
        if ipt == len(df)-1 and pt < 50.:
            graph.SetPoint(ipt, 50., raa)
    return graph

def get_centrality_interpolation(graphs, base_name):
    """
    Function to compute centrality interpolation
    """

    cents = [f"{icent*10}_{icent*10+10}" for icent in range(9)] #0-10 to 80-90
    # extend cents with missing centralities if needed
    for key in graphs.keys():
        if key not in cents:
            cents.append(key)    
    model_exists = {}
    for cent_key in cents:
        model_exists[cent_key] = True
        if cent_key not in graphs:
            graphs[cent_key] = ROOT.TGraph(1)
            graphs[cent_key].SetNameTitle(f"{base_name}_{cent_key}",
                                          ";#it{p}_{T} (GeV/#it{c}); #it{R}_{AA}")
            model_exists[cent_key] = False

    for ipt in range(graphs["0_10"].GetN()): # assuming 0-10% always there
        graph_interp = ROOT.TGraph(1)
        func_interp = ROOT.TF1("func_interp", "pol1", 0., 90., 2)
        icent = 0
        for cent_key, cent_exists in model_exists.items():
            if cent_exists:
                cent_min = float(cent_key.split(sep="_")[0])
                cent_max = float(cent_key.split(sep="_")[1])
                cent_mid = (cent_min+cent_max)/2
                graph_interp.SetPoint(
                    icent, cent_mid, graphs[cent_key].GetPointY(ipt))
                icent += 1
        graph_interp.Fit(func_interp, "Q0")
        for icent, cent_key in enumerate(model_exists.keys()):
            if not model_exists[cent_key]:
                cent_min = float(cent_key.split(sep="_")[0])
                cent_max = float(cent_key.split(sep="_")[1])
                cent_mid = (cent_min+cent_max)/2
                graphs[cent_key].SetPoint(
                    ipt, graphs["0_10"].GetPointX(ipt), func_interp.Eval(cent_mid))

def get_fonll_times_raa(hist_fonll, graphs_raa, base_name):
    """
    Function to compute FONLL x RAA
    """
    cents = [0 + 10 * icent for icent in range(10)] #0-10 to 80-90
    hists = {}
    for cent_min, cent_max in zip(cents[:-1], cents[1:]):
        cent_key = f"{cent_min}_{cent_max}"
        hists[cent_key] = hist_fonll.Clone(f"{base_name}_{cent_key}")
        hists[cent_key].SetDirectory(0)
        for ipt in range(1, hists[cent_key].GetNbinsX()+1):
            pt_cent = hists[cent_key].GetBinCenter(ipt)
            hists[cent_key].SetBinContent(
                ipt, hist_fonll.GetBinContent(ipt) * graphs_raa[cent_key].Eval(pt_cent))
        hists[cent_key].Scale(1./hists[cent_key].Integral())

    return hists
