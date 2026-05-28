import argparse
import enum
import ROOT
import numpy as np
import matplotlib as mpl
import ctypes

def __convert_ntracks_to_dn_deta_backward(ntracks):
    return {
        0.507467: 34.6,
        1.014934: 54.6,
        1.304915: 70.1,
        1.594896: 86.0,
        1.884877: 102.3,
        2.319849: 126.2,
        3.117297: 164.7
    }.get(ntracks, 0) / 2.8 # eta gap: 2 < eta < 4.8

def __convert_ntracks_to_dn_deta_error_backward(ntracks):
    return {
        0.507467: 10.5,
        1.014934: 8.8,
        1.304915: 10.2,
        1.594896: 11.5,
        1.884877: 12.7,
        2.319849: 16.7,
        3.117297: 22.1
    }.get(ntracks, 0) / 2.8 # eta gap: 2 < eta < 4.8

def __convert_ntracks_to_dn_deta_forward(ntracks):
    return {
        0.5809129: 32.8,
        1.161826: 49.9,
        1.493776: 61.9,
        1.825726: 74.5,
        2.157676: 87.5,
        2.821577: 108.4
    }.get(ntracks, 0) / 2.8 # eta gap: 2 < eta < 4.8

def __convert_ntracks_to_dn_deta_error_forward(ntracks):
    return {
        0.5809129: 9.8,
        1.161826: 8.6,
        1.493776: 10.0,
        1.825726: 11.4,
        2.157676: 12.5,
        2.821577: 17.0
    }.get(ntracks, 0) / 2.8 # eta gap: 2 < eta < 4.8

def get_discrete_matplotlib_palette(paletteName):
    cmap = mpl.colormaps[paletteName]
    colors = cmap.colors
    ROOTColorIndices = []
    ROOTColors = []
    for color in colors:
        idx = ROOT.TColor.GetFreeColorIndex()
        ROOTColors.append(ROOT.TColor(idx, color[0], color[1], color[2],"color%i" % idx))
        ROOTColorIndices.append(idx)
        
    return ROOTColorIndices, ROOTColors

def draw_graphs(graph_stat, graph_syst, c=ROOT.kBlack, s=1., m=ROOT.kFullCircle):
    graph_stat.SetMarkerStyle(m)
    graph_stat.SetMarkerSize(s)
    graph_stat.SetMarkerColor(c)
    graph_stat.SetLineColor(c)
    graph_stat.SetLineWidth(1)
    graph_syst.SetMarkerStyle(m)
    graph_syst.SetMarkerSize(s)
    graph_syst.SetMarkerColor(c)
    graph_syst.SetLineColor(c)
    graph_syst.SetLineWidth(1)
    graph_syst.SetFillStyle(0)
    graph_stat.Draw("pz, same")
    graph_syst.Draw("pz2, same")


def draw_alice_pbpb(pt_min, pt_max, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Ratios/ratio_w_syst.root") as infile:
        graph_stat_0_20 = infile.Get(f"g_stat_{pt_min*10:.0f}_{pt_max*10:.0f}")
        graph_syst_0_20 = infile.Get(f"g_syst_no_br_{pt_min*10:.0f}_{pt_max*10:.0f}")
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/Data/ds_over_dplus_ratios_2050.root") as infile:
        graph_stat_20_50 = infile.Get(f"g_stat_{pt_min*10:.0f}_{pt_max*10:.0f}")
        graph_syst_20_50 = infile.Get(f"g_syst_no_br_{pt_min*10:.0f}_{pt_max*10:.0f}")
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/Data/ds_over_dplus_ratios_5090.root") as infile:
        graph_stat_50_90 = infile.Get(f"g_stat_{pt_min*10:.0f}_{pt_max*10:.0f}")
        graph_syst_50_90 = infile.Get(f"g_syst_no_br_{pt_min*10:.0f}_{pt_max*10:.0f}")

    if (pt_min, pt_max) == (12, 24):
        #remove first two points for 50-90%
        graph_stat_50_90.SetPoint(2, 0, 1.e9)
        graph_stat_50_90.SetPoint(3, 0, 1.e9)
        graph_syst_50_90.SetPoint(2, 0, 1.e9)
        graph_syst_50_90.SetPoint(3, 0, 1.e9)

    graph_stat_0_20.SetMarkerStyle(m)
    graph_stat_0_20.SetMarkerSize(s)
    graph_stat_0_20.SetMarkerColor(c)
    graph_stat_0_20.SetLineColor(c)
    graph_stat_0_20.SetLineWidth(1)
    graph_stat_0_20.Draw("pz, same")
    graph_syst_0_20.SetMarkerStyle(m)
    graph_syst_0_20.SetMarkerSize(s)
    graph_syst_0_20.SetMarkerColor(c)
    graph_syst_0_20.SetLineColor(c)
    graph_syst_0_20.SetLineWidth(1)
    graph_syst_0_20.SetFillStyle(0)
    graph_syst_0_20.Draw("pz2, same")
    graph_stat_20_50.SetMarkerStyle(m)
    graph_stat_20_50.SetMarkerSize(s)
    graph_stat_20_50.SetMarkerColor(c)
    graph_stat_20_50.SetLineColor(c)
    graph_stat_20_50.SetLineWidth(1)
    graph_stat_20_50.Draw("pz, same")
    graph_syst_20_50.SetMarkerStyle(m)
    graph_syst_20_50.SetMarkerSize(s)
    graph_syst_20_50.SetMarkerColor(c)
    graph_syst_20_50.SetLineColor(c)
    graph_syst_20_50.SetLineWidth(1)
    graph_syst_20_50.SetFillStyle(0)
    graph_syst_20_50.Draw("pz2, same")
    graph_stat_50_90.SetMarkerStyle(m)
    graph_stat_50_90.SetMarkerSize(s)
    graph_stat_50_90.SetMarkerColor(c)
    graph_stat_50_90.SetLineColor(c)
    graph_stat_50_90.SetLineWidth(1)
    graph_stat_50_90.Draw("pz, same")
    graph_syst_50_90.SetMarkerStyle(m)
    graph_syst_50_90.SetMarkerSize(s)
    graph_syst_50_90.SetMarkerColor(c)
    graph_syst_50_90.SetLineColor(c)
    graph_syst_50_90.SetLineWidth(1)
    graph_syst_50_90.SetFillStyle(0)
    graph_syst_50_90.Draw("pz2, same")
    return graph_stat_0_20#, graph_stat_50_90



def draw_alice_pp(pt_min, pt_max, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Ratios/VsMult/w_syst/ratio_w_syst.root") as infile:
        graph_stat = infile.Get(f"g_stat_{pt_min*10:.0f}_{pt_max*10:.0f}")
        graph_syst = infile.Get(f"g_syst_no_br_{pt_min*10:.0f}_{pt_max*10:.0f}")
    
    #set x unc. for the systematic box
    for i in range(graph_syst.GetN()):
        graph_syst.SetPointEXlow(i, graph_syst.GetErrorXlow(i)*5)
        graph_syst.SetPointEXhigh(i, graph_syst.GetErrorXhigh(i)*5)

    graph_stat.SetMarkerStyle(m)
    graph_stat.SetMarkerSize(s)
    graph_stat.SetMarkerColor(c)
    graph_stat.SetLineColor(c)
    graph_stat.SetLineWidth(1)
    graph_syst.SetMarkerStyle(m)
    graph_syst.SetMarkerSize(s)
    graph_syst.SetMarkerColor(c)
    graph_syst.SetLineColor(c)
    graph_syst.SetLineWidth(1)
    graph_syst.SetFillStyle(0)
    graph_stat.DrawClone("pz, same")
    graph_syst.DrawClone("pz2, same")
    return graph_stat #, graph_syst


def draw_pythia_pp_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_Monash")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(1)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3L, same")
    return graph_mid

def draw_pythia_pp_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_Monash")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(1)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3L, same")
    return graph_ft0

def draw_pythia_pp_mod_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_MonashModDstar")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(1)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3L, same")
    return graph_mid

def draw_pythia_pp_mod_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_MonashModDstar")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(1)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3L, same")
    return graph_ft0

def draw_pythia_pp_mode2_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_Mode2")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(1)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3L, same")
    return graph_mid

def draw_pythia_pp_mode2_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_Mode2")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(1)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3L, same")
    return graph_ft0

def draw_pythia_pp_mod_mode2_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_Mode2ModDstar")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(1)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3L, same")
    return graph_mid

def draw_pythia_pp_mod_mode2_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_Mode2ModDstar")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(1)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3L, same")
    return graph_ft0

def draw_pythia_pp_sccr_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_SRRC")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(1)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3L, same")
    return graph_mid

def draw_pythia_pp_sccr_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_SRRC")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(1)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3L, same")
    return graph_ft0

def draw_pythia_pbpb_dstar_tune_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/PYTHIA_Simulations/pbpb/CharmHadRatiosVsMult_angantyr_dstar_tune.root") as infile:
        hist_ft0 = infile.Get(f"gDsDpRatioMidMultInFwdBinsPt{pt_min:.0f}{pt_max:.0f}")
    for i in range(hist_ft0.GetN()):
        hist_ft0.SetPointEXlow(i, 0)
        hist_ft0.SetPointEXhigh(i, 0)
    hist_ft0.SetMarkerStyle(m)
    hist_ft0.SetMarkerSize(s)
    hist_ft0.SetMarkerColor(c)
    hist_ft0.SetLineColor(c)
    hist_ft0.SetLineWidth(1)
    hist_ft0.SetFillStyle(1001)
    hist_ft0.SetFillColorAlpha(c, a)
    hist_ft0.Draw("3L, same")
    return hist_ft0

def draw_pythia_pbpb_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/PYTHIA_Simulations/pbpb/CharmHadRatiosVsMult_angantyr.root") as infile:
        hist_ft0 = infile.Get(f"gDsDpRatioMidMultInFwdBinsPt{pt_min:.0f}{pt_max:.0f}")
    for i in range(hist_ft0.GetN()):
        hist_ft0.SetPointEXlow(i, 0)
        hist_ft0.SetPointEXhigh(i, 0)
    hist_ft0.SetMarkerStyle(m)
    hist_ft0.SetMarkerSize(s)
    hist_ft0.SetMarkerColor(c)
    hist_ft0.SetLineColor(c)
    hist_ft0.SetLineWidth(1)
    hist_ft0.SetFillStyle(1001)
    hist_ft0.SetFillColorAlpha(c, a)
    hist_ft0.Draw("3L, same")
    return hist_ft0

def draw_pythia_pbpb_sccr_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/pbpb/CharmHadRatiosVsMult_merge_cut.root") as infile:
        hist_ft0 = infile.Get(f"gDsDpRatioMidMultInFwdBinsPt{pt_min:.0f}{pt_max:.0f}")
    deta = (3.3-2.1) + (4.9-3.5) 
    for i in range(hist_ft0.GetN()):
        hist_ft0.SetPointX(i, hist_ft0.GetX()[i] / deta)
        hist_ft0.SetPointEXlow(i, 0)
        hist_ft0.SetPointEXhigh(i, 0)
    hist_ft0.SetMarkerStyle(m)
    hist_ft0.SetMarkerSize(s)
    hist_ft0.SetMarkerColor(c)
    hist_ft0.SetLineColor(c)
    hist_ft0.SetLineWidth(1)
    hist_ft0.SetFillStyle(1001)
    hist_ft0.SetFillColorAlpha(c, a)
    hist_ft0.Draw("3L, same")
    return hist_ft0

def draw_epos4hq_pp_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/EPOS4HQ_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_Frag+coal")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(1)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3L, same")
    return graph_mid

def draw_epos4hq_pp_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/EPOS4HQ_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_Frag+coal")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(1)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3L, same")
    return graph_ft0

def draw_epos4hq_pbpb_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/EPOS4HQ_Ds_over_Dplus_PbPb5dot02TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_Frag+coal")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(1)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3L, same")
    return graph_mid

def draw_epos4hq_pbpb_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/EPOS4HQ_Ds_over_Dplus_PbPb5dot02TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_Frag+coal")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(1)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3L, same")
    return graph_ft0

def to_pad_coordinates(x=None, y=None):
    xl, yl, xu, yu = ctypes.c_double(), ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    ROOT.gPad.GetPadPar(xl, yl, xu, yu)
    xl, yl, xu, yu = xl.value, yl.value, xu.value, yu.value
    pw, ph = xu - xl, yu - yl
    lm, rm, tm, bm = ROOT.gPad.GetLeftMargin(), ROOT.gPad.GetRightMargin(), ROOT.gPad.GetTopMargin(), ROOT.gPad.GetBottomMargin()
    fw, fh = pw - pw * lm - pw * rm, ph - ph * bm - ph * tm

    x_pad = (x * fw + pw * lm) / pw if x is not None else None
    y_pad = (y * fh + bm * ph) / ph if y is not None else None

    return x_pad, y_pad

def main(do_pythia, do_epos):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadRightMargin(0.03)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadTopMargin(0.02)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetOptLogx(1)
    ROOT.gStyle.SetLabelFont(43, "XYZ")
    ROOT.gStyle.SetLabelSize(45, "XYZ")
    ROOT.gStyle.SetTitleFont(43, "XYZ")
    ROOT.gStyle.SetTitleSize(55, "XYZ")
    ROOT.gStyle.SetTitleOffset(1., "X")
    ROOT.gStyle.SetLabelOffset(-0.01, "X")
    colors, cols = get_discrete_matplotlib_palette('tab20')

    ORDER = {
        #'PYTHIA_PP_FT0M': 0,
        #'PYTHIA_PP_MOD_FT0M': 1,
        #'PYTHIA_PP_MODE2_FT0M': 2,
        'PYTHIA_PP_MOD_MODE2_FT0M': 3,
        'PYTHIA_PP_SCCR_FT0M': 4,
        'PYTHIA_PBPB_SCCR_FT0M': 5,
        'PYTHIA_PBPB_MOD_FT0M': 6,
        #'PYTHIA_PBPB_FT0M': 7,
        'EPOS4HQ_PP_FT0M': 8,
        'EPOS4HQ_PBPB_FT0M': 9,
        'ALICE_PP': 10,
        'ALICE_PBPB': 11
    }


    functions = [
        draw_pythia_pp_ft0m,
        draw_pythia_pp_mod_ft0m,
        draw_pythia_pp_mode2_ft0m,
        draw_pythia_pp_mod_mode2_ft0m,
        draw_pythia_pp_sccr_ft0m,
        draw_pythia_pbpb_sccr_ft0m,
        draw_pythia_pbpb_dstar_tune_ft0m,
        draw_pythia_pbpb_ft0m,
        draw_epos4hq_pp_ft0m,
        draw_epos4hq_pbpb_ft0m,
        draw_alice_pp,
        draw_alice_pbpb
    ]

    colors_pythia = [
        colors[2],
        colors[0],
        colors[4],
        colors[2],
        colors[0],
        colors[1],
        colors[3],
        colors[18],
        ROOT.kGray,
        ROOT.kGray,
        ROOT.kBlack,
        ROOT.kBlack
    ]

    colors_epos = [
        colors[2],
        colors[0],
        colors[4],
        ROOT.kGray, #colors[2],
        ROOT.kGray, #colors[0],
        ROOT.kGray, #colors[1],
        ROOT.kGray, #colors[3],
        colors[18],
        colors[6],
        colors[7],
        ROOT.kBlack,
        ROOT.kBlack
    ]

    colors_combined = [
        colors[2],
        colors[0],
        colors[4],
        colors[2],
        colors[0],
        colors[1],
        colors[3],
        colors[18],
        colors[6],
        colors[7],
        ROOT.kBlack,
        ROOT.kBlack
    ]

    if do_pythia and do_epos:
        colors_used = colors_combined
    else:
        colors_used = colors_pythia if do_pythia else colors_epos

    markers_used = [
        ROOT.kFullDoubleDiamond,
        ROOT.kFullDiamond,
        ROOT.kFullSquare,
        ROOT.kFullCross,
        ROOT.kFullFourTrianglesPlus,
        ROOT.kFullFourTrianglesPlus,
        ROOT.kFullFourTrianglesPlus,
        ROOT.kFullFourTrianglesPlus,
        ROOT.kFullFourTrianglesPlus,
        ROOT.kFullFourTrianglesPlus,
        ROOT.kOpenDiamond,
        ROOT.kFullCircle
    ]

    sizes_used = [
        3,
        3,
        2,
        2.5,
        2,
        0,
        0,
        0,
        0,
        0,
        3,
        2.5
    ]

    y_text = ROOT.TLatex(0.185, 0.83, '|#it{y}| < 0.5')
    y_text.SetNDC()
    y_text.SetTextFont(42)
    y_text.SetTextSize(0.04)

    pt_text = ROOT.TLatex(0.6, 0.88, '')
    pt_text.SetNDC()
    pt_text.SetTextFont(43)
    pt_text.SetTextSize(40)

    pp_unc_text = ROOT.TLatex(0.6, 0.15, '')
    pp_unc_text.SetNDC()
    pp_unc_text.SetTextFont(43)
    pp_unc_text.SetTextSize(40)

    br_unc_text = ROOT.TLatex(0.6, 0.15, '')
    br_unc_text.SetNDC()
    br_unc_text.SetTextFont(43)
    br_unc_text.SetTextSize(40)

    c = ROOT.TCanvas("canvas", "canvas", 2400, 1600)

    lm = 0.2
    rm = 0.03
    tm = 0.03
    bm = 0.21
    l = 1 / (1 / (1-lm) + 1 + 1/(1-rm))
    h = 1 / (1 / (1-tm) + 1/(1-bm))

    c = ROOT.TCanvas("canvas", "canvas", 2400, 1600)
    pads = []
    # c.Divide(3, 2, 0.000, 0.000)
    # Set margins
    pads.append(ROOT.TPad(f"pad0", f"pad0", 0, h / (1-bm), l/(1-lm), 1.))
    pads.append(ROOT.TPad(f"pad1", f"pad1", l/(1-lm), h / (1-bm), l/(1-lm) + l, 1.))
    pads.append(ROOT.TPad(f"pad2", f"pad2", l/(1-lm) + l, h / (1-bm), l/(1-lm) + l + l/(1-rm), 1.))
    pads.append(ROOT.TPad(f"pad3", f"pad3", 0, 0, l/(1-lm), h / (1-bm)))
    pads.append(ROOT.TPad(f"pad4", f"pad4", l/(1-lm), 0, l/(1-lm) + l, h / (1-bm)))
    pads.append(ROOT.TPad(f"pad5", f"pad5", l/(1-lm) + l, 0, l/(1-lm) + l + l/(1-rm), h / (1-bm)))
    for i_pad in range(3):
        pads[i_pad].SetTopMargin(tm)
        pads[i_pad].SetBottomMargin(0.)
        pads[i_pad].Draw()
    for i_pad in range(3):
        pads[i_pad + 3].SetBottomMargin(bm)
        pads[i_pad + 3].SetTopMargin(0.)
        pads[i_pad + 3].Draw()

    for i in (0, 3):
        pads[i].SetRightMargin(0.)
        pads[i].SetLeftMargin(lm)
    for i in (1, 4):
        pads[i].SetLeftMargin(0.)
        pads[i].SetRightMargin(0.)
    for i in (2, 5):
        pads[i].SetLeftMargin(0.)
        pads[i].SetRightMargin(rm)

    pads[0].cd()
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}| < 0.5};#it{#sigma}(D_{s}^{+})/#it{#sigma}(D^{+})")
    h_frame.GetYaxis().ChangeLabel(1, 1, 0)

    results = [None] * (ORDER["ALICE_PBPB"] + 1)
    for key, i in ORDER.items():
        results[i] = functions[i](1, 2, c=colors_used[i], s=sizes_used[i], m=markers_used[i])

    x, y = to_pad_coordinates(0.05, 0.9)
    alice_text = ROOT.TLatex(x, y, 'ALICE Preliminary')
    alice_text.SetNDC()
    alice_text.SetTextFont(43)
    alice_text.SetTextSize(50)
    alice_text.Draw()

    x, y = to_pad_coordinates(0.05, 0.85)
    pt_text.DrawLatexNDC(x, y, '1#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{2} GeV/#it{c}')


    pythia_header = ROOT.TLatex()
    pythia_header.SetNDC()
    pythia_header.SetTextFont(43)
    pythia_header.SetTextSize(40)
    x, y = to_pad_coordinates(0.04, 0.235)
    pythia_header.DrawLatex(x, y, 'PYTHIA 8')
    x, y = to_pad_coordinates(0.04, 0.18)
    pythia_header.DrawLatex(x, y, 'pp,#kern[0.07]{#sqrt{#it{s}} = 13.6 TeV}')
    x, y = to_pad_coordinates(0.04, 0.12)
    pythia_header.DrawLatex(x, y, 'FT0M multiplicity estimator')

    x_min, y_min = to_pad_coordinates(0.03, 0.02)
    x_max, y_max = to_pad_coordinates(0.99, 0.12)
    
    legend_pythia_no_mod = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_pythia_no_mod.SetBorderSize(0)
    legend_pythia_no_mod.SetFillStyle(0)
    legend_pythia_no_mod.SetTextFont(43)
    legend_pythia_no_mod.SetNColumns(3)
    legend_pythia_no_mod.SetTextSize(35)
    if 'PYTHIA_PP_FT0M' in ORDER:
        legend_pythia_no_mod.AddEntry(results[ORDER['PYTHIA_PP_FT0M']], 'Monash', 'fl')
    if 'PYTHIA_PP_SCCR_FT0M' in ORDER:
        legend_pythia_no_mod.AddEntry(results[ORDER['PYTHIA_PP_SCCR_FT0M']], 'SC#minusCR', 'fl')
    if 'PYTHIA_PP_MOD_MODE2_FT0M' in ORDER:
        legend_pythia_no_mod.AddEntry(results[ORDER['PYTHIA_PP_MOD_MODE2_FT0M']], '#splitline{StringFlav:mesonCvector=1.75}{#lower[-0.2]{CR-BLC Mode 2}}', 'fl')
    if 'PYTHIA_PP_MODE2_FT0M' in ORDER:
        legend_pythia_no_mod.AddEntry(results[ORDER['PYTHIA_PP_MODE2_FT0M']], 'CR-BLC Mode 2', 'fl')
    legend_pythia_no_mod.Draw()

    ROOT.gPad.RedrawAxis()

    pads[1].cd()
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}| < 0.5};#it{#sigma}(D_{s}^{+})/#it{#sigma}(D^{+})")
    for key, i in ORDER.items():
        functions[i](2, 4, c=colors_used[i], s=sizes_used[i], m=markers_used[i])
    #alice_text.Draw()
    #y_text.Draw()
    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '2#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{4} GeV/#it{c}')


    pythia_pbpb_header = ROOT.TLatex()
    pythia_pbpb_header.SetNDC()
    pythia_pbpb_header.SetTextFont(43)
    pythia_pbpb_header.SetTextSize(40)
    x, y = to_pad_coordinates(0.04, 0.235)
    pythia_pbpb_header.DrawLatex(x, y, 'PYTHIA 8 Angantyr')
    x, y = to_pad_coordinates(0.04, 0.18)
    pythia_pbpb_header.DrawLatex(x, y, 'Pb#minusPb,#kern[0.07]{#sqrt{#it{s}_{NN}} = 5.36 TeV}')
    x, y = to_pad_coordinates(0.04, 0.12)
    pythia_pbpb_header.DrawLatex(x, y, 'FT0M multiplicity estimator')

    x_min, y_min = to_pad_coordinates(0.03, 0.02)
    x_max, y_max = to_pad_coordinates(0.97, 0.12)
    legend_pythia_pbpb = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_pythia_pbpb.SetBorderSize(0)
    legend_pythia_pbpb.SetFillStyle(0)
    legend_pythia_pbpb.SetTextFont(43)
    legend_pythia_pbpb.SetNColumns(3)
    legend_pythia_pbpb.SetTextSize(35)
    # legend_pythia_pbpb.SetMargin(0.40)
    if 'PYTHIA_PBPB_SCCR_FT0M' in ORDER:
        legend_pythia_pbpb.AddEntry(results[ORDER['PYTHIA_PBPB_SCCR_FT0M']], 'SC#minusCR', 'fl')
    if 'PYTHIA_PBPB_FT0M' in ORDER:
        legend_pythia_pbpb.AddEntry(results[ORDER['PYTHIA_PBPB_FT0M']], 'Angantyr', 'fl')
    if 'PYTHIA_PBPB_MOD_FT0M' in ORDER:
        legend_pythia_pbpb.AddEntry(results[ORDER['PYTHIA_PBPB_MOD_FT0M']], 'StringFlav:mesonCvector=1.75', 'fl')
    legend_pythia_pbpb.Draw()

    ROOT.gPad.RedrawAxis()

    pads[2].cd()
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}| < 0.5};#it{#sigma}(D_{s}^{+})/#it{#sigma}(D^{+})")
    for key, i in ORDER.items():
        functions[i](4, 6, c=colors_used[i], s=sizes_used[i], m=markers_used[i])

    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '4#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{6} GeV/#it{c}')

    epos4hq_header = ROOT.TLatex()
    epos4hq_header.SetNDC()
    epos4hq_header.SetTextFont(43)
    epos4hq_header.SetTextSize(40)
    x, y = to_pad_coordinates(0.04, 0.23)
    epos4hq_header.DrawLatex(x, y, 'EPOS4HQ')
    x, y = to_pad_coordinates(0.04, 0.17)
    epos4hq_header.DrawLatex(x, y, 'FT0M multiplicity estimator')

    x_min, y_min = to_pad_coordinates(0.03, 0.02)
    x_max, y_max = to_pad_coordinates(0.97, 0.15)
    legend_epos4hq = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_epos4hq.SetBorderSize(0)
    legend_epos4hq.SetFillStyle(0)
    legend_epos4hq.SetTextFont(43)
    legend_epos4hq.SetNColumns(1)
    legend_epos4hq.SetTextSize(40)
    legend_epos4hq.SetMargin(0.1)
    if 'EPOS4HQ_PP_FT0M' in ORDER:
        legend_epos4hq.AddEntry(results[ORDER["EPOS4HQ_PP_FT0M"]], 'pp,#kern[0.07]{#sqrt{#it{s}} = 13.6 TeV}', 'fl')
    if 'EPOS4HQ_PBPB_FT0M' in ORDER:
        legend_epos4hq.AddEntry(results[ORDER["EPOS4HQ_PBPB_FT0M"]], 'Pb#minusPb#kern[0.07]{#sqrt{#it{s}_{NN}} = 5.02 TeV}', 'fl')
    legend_epos4hq.Draw()

    ROOT.gPad.RedrawAxis()

    pads[3].cd()
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}| < 0.5};#it{#sigma}(D_{s}^{+})/#it{#sigma}(D^{+})")
    h_frame.GetYaxis().ChangeLabel(10, 1, 0)
    h_frame.GetXaxis().CenterTitle(True)
    for key, i in ORDER.items():
        functions[i](6, 8, c=colors_used[i], s=sizes_used[i], m=markers_used[i])
    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '6#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{8} GeV/#it{c}')

    x_min, y_min = to_pad_coordinates(0.03, 0.07)
    x_max, y_max = to_pad_coordinates(0.5, 0.17)
    legend_alice_pp = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_alice_pp.SetBorderSize(0)
    legend_alice_pp.SetFillStyle(0)
    legend_alice_pp.SetTextFont(43)
    legend_alice_pp.SetTextSize(40)
    legend_alice_pp.AddEntry(results[ORDER['ALICE_PP']], '#splitline{pp,#kern[0.2]{#sqrt{#it{s}}} = 13.6 TeV, |#it{y}| < 0.5}{FT0M multiplicity estimator}', 'p')
    legend_alice_pp.Draw()

    ROOT.gPad.RedrawAxis()

    pads[4].cd()
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}| < 0.5};#it{#sigma}(D_{s}^{+})/#it{#sigma}(D^{+})")
    h_frame.GetXaxis().ChangeLabel(1, 1, 0)
    h_frame.GetXaxis().CenterTitle(True)
    for key, i in ORDER.items():
        functions[i](8, 12, c=colors_used[i], s=sizes_used[i], m=markers_used[i])
    #alice_text.Draw()
    #y_text.Draw()
    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '8#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[.5]{12} GeV/#it{c}')

    x_min, y_min = to_pad_coordinates(0.03, 0.07)
    x_max, y_max = to_pad_coordinates(0.5, 0.17)
    legend_alice_pbpb = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_alice_pbpb.SetBorderSize(0)
    legend_alice_pbpb.SetFillStyle(0)
    legend_alice_pbpb.SetTextFont(43)
    legend_alice_pbpb.SetTextSize(40)
    legend_alice_pbpb.AddEntry(results[ORDER['ALICE_PBPB']], '#splitline{Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.36 TeV, |#it{y}| < 0.5}{FT0C multiplicity estimator}', 'p')
    legend_alice_pbpb.Draw()

    ROOT.gPad.RedrawAxis()

    pads[5].cd()
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}| < 0.5};#it{#sigma}(D_{s}^{+})/#it{#sigma}(D^{+})")
    h_frame.GetXaxis().ChangeLabel(1, 1, 0)
    h_frame.GetXaxis().CenterTitle(True)
    for key, i in ORDER.items():
        functions[i](12, 24, c=colors_used[i], s=sizes_used[i], m=markers_used[i])
    #alice_text.Draw()
    #y_text.Draw()
    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '12#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[.5]{24} GeV/#it{c}')

    x, y = to_pad_coordinates(0.05, 0.145)
    pp_unc_text.DrawLatexNDC(x, y, 'pp uncertainty on #kern[-0.05]{#it{x}-axis scaled by 5}')

    x, y = to_pad_coordinates(0.05, 0.065)
    br_unc_text.DrawLatexNDC(x, y, '#lower[-0.03]{^{+4.0}}')
    br_unc_text.DrawLatexNDC(x, y, '_{#minus3.8}')
    br_unc_text.DrawLatexNDC(x, y, '#kern[.12]{% BR uncertainty not shown}')

    ROOT.gPad.RedrawAxis()

    if do_pythia and do_epos:
        c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_vs_pred_combined.pdf")
        c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_vs_pred_combined.root")
    elif do_pythia:
        c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_vs_pred_pythia.pdf")
        c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_vs_pred_pythia.root")
    elif do_epos:
        c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_vs_pred_epos.pdf")
        c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_vs_pred_epos.root")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw Ds/D+ ratio vs multiplicity with model predictions')
    parser.add_argument('--pythia', '-p', action='store_true', default=False, help='Input ROOT file')
    parser.add_argument('--epos', '-e', action='store_true', default=False, help='Output PDF file')
    args = parser.parse_args()

    if not args.pythia and not args.epos:
        # Run all combinations
        main(True, True)
        main(True, False)
        main(False, True)

    main(args.pythia, args.epos)