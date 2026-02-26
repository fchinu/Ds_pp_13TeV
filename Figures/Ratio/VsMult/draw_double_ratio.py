import argparse
import enum
import ROOT
import numpy as np
import matplotlib as mpl
import ctypes

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

def get_ratio_pbpb_pp(h_stat_pbpb, h_stat_pp, g_syst_pbpb, g_syst_pp):
    h_stat_ratio = h_stat_pbpb.Clone("h_stat_ratio")
    g_syst_ratio = g_syst_pbpb.Clone("g_syst_ratio")

    h_stat_ratio.Divide(h_stat_pp)

    for i in range(g_syst_pbpb.GetN()):
        g_syst_ratio.SetPointY(i, g_syst_pbpb.GetY()[i] / g_syst_pp.GetY()[i])
        err_pbpb_high = g_syst_pbpb.GetErrorYhigh(i)
        err_pbpb_low = g_syst_pbpb.GetErrorYlow(i)
        err_pp_high = g_syst_pp.GetErrorYhigh(i)
        err_pp_low = g_syst_pp.GetErrorYlow(i)

        err_ratio_high = (g_syst_ratio.GetY()[i]) * np.sqrt( (err_pbpb_high / g_syst_pbpb.GetY()[i])**2 + (err_pp_low / g_syst_pp.GetY()[i])**2 )
        err_ratio_low = (g_syst_ratio.GetY()[i]) * np.sqrt( (err_pbpb_low / g_syst_pbpb.GetY()[i])**2 + (err_pp_high / g_syst_pp.GetY()[i])**2 )
    
        g_syst_ratio.SetPointEYhigh(i, err_ratio_high)
        g_syst_ratio.SetPointEYlow(i, err_ratio_low)

    return h_stat_ratio, g_syst_ratio

def draw_alice_pbpb_pp_ratio(cent_min_pbpb, cent_max_pbpb, cent_min_pp, cent_max_pp, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle):
    g_syst_pbpb = None
    if cent_min_pbpb < 20:
        with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Ratios/ratio_w_syst.root") as infile:
            h_stat_pbpb = infile.Get(f"h_stat_vs_pt_{cent_min_pbpb:.0f}_{cent_max_pbpb:.0f}")
            h_stat_pbpb.SetDirectory(0)
            g_syst_pbpb = infile.Get(f"g_syst_vs_pt_no_br_{cent_min_pbpb:.0f}_{cent_max_pbpb:.0f}")
    elif cent_min_pbpb < 50:
        with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/Data/ds_over_dplus_ratios_2050.root") as infile:
            h_stat_pbpb = infile.Get(f"h_stat_vs_pt_{cent_min_pbpb:.0f}_{cent_max_pbpb:.0f}")
            h_stat_pbpb.SetDirectory(0)
            g_syst_pbpb = infile.Get(f"g_syst_vs_pt_no_br_{cent_min_pbpb:.0f}_{cent_max_pbpb:.0f}")
    elif cent_min_pbpb < 90:
        with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/Data/ds_over_dplus_ratios_5090.root") as infile:
            h_stat_pbpb = infile.Get(f"h_stat_vs_pt_{cent_min_pbpb:.0f}_{cent_max_pbpb:.0f}")
            h_stat_pbpb.SetDirectory(0)
            g_syst_pbpb = infile.Get(f"g_syst_vs_pt_no_br_{cent_min_pbpb:.0f}_{cent_max_pbpb:.0f}")
    
    if cent_min_pbpb >= 70:
        #remove last two points for 70-80% and 80-90%
        h_stat_pbpb.SetBinContent(6, 1.e9)
        if g_syst_pbpb:
            g_syst_pbpb.SetPoint(5, 0, 1.e9)

    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Ratios/VsMult/w_syst/ratio_w_syst.root") as infile:
        graph_stat = infile.Get(f"h_stat_vs_pt_{cent_min_pp:.0f}_{cent_max_pp:.0f}")
        graph_stat.SetDirectory(0)
        graph_syst = infile.Get(f"g_syst_vs_pt_no_br_{cent_min_pp:.0f}_{cent_max_pp:.0f}")

    h_stat_ratio, g_syst_ratio = get_ratio_pbpb_pp(h_stat_pbpb, graph_stat, g_syst_pbpb, graph_syst)

    h_stat_ratio.SetMarkerStyle(m)
    h_stat_ratio.SetMarkerSize(s)
    h_stat_ratio.SetMarkerColor(c)
    h_stat_ratio.SetLineColor(c)
    h_stat_ratio.SetLineWidth(2)
    h_stat_ratio.Draw("pz, same")
    if g_syst_ratio:
        #set x unc. for the systematic box
        for i in range(g_syst_ratio.GetN()):
            g_syst_ratio.SetPointEXlow(i, g_syst_ratio.GetErrorXlow(i)/2)
            g_syst_ratio.SetPointEXhigh(i, g_syst_ratio.GetErrorXhigh(i)/2)
        g_syst_ratio.SetMarkerStyle(m)
        g_syst_ratio.SetMarkerSize(s)
        g_syst_ratio.SetMarkerColor(c)
        g_syst_ratio.SetLineColor(c)
        g_syst_ratio.SetLineWidth(2)
        g_syst_ratio.SetFillStyle(0)
        g_syst_ratio.Draw("E2, same")
    return h_stat_ratio#, graph_stat_50_90


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

if __name__ == '__main__':
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadRightMargin(0.03)
    ROOT.gStyle.SetPadLeftMargin(0.1)
    ROOT.gStyle.SetPadTopMargin(0.02)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    # ROOT.gStyle.SetOptLogx(1)
    # ROOT.gStyle.SetLabelFont(43, "XYZ")
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetTitleOffset(0.9, "Y")
    # ROOT.gStyle.SetTitleFont(43, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    # ROOT.gStyle.SetLabelOffset(0.005, "X")
    colors, cols = get_discrete_matplotlib_palette('tab10')

    y_text = ROOT.TLatex(0.185, 0.83, '|#it{y}| < 0.5')
    y_text.SetNDC()
    y_text.SetTextFont(42)
    y_text.SetTextSize(0.04)


    legend_text = ROOT.TLatex(0.6, 0.25, '')
    legend_text.SetNDC()
    legend_text.SetTextFont(43)
    legend_text.SetTextSize(25)

    br_unc_text = ROOT.TLatex(0.6, 0.25, '')
    br_unc_text.SetNDC()
    br_unc_text.SetTextFont(43)
    br_unc_text.SetTextSize(20)

    c = ROOT.TCanvas("canvas", "canvas", 700, 600)
    h_frame = c.DrawFrame(0., 0., 24., 2.5, ";#it{p}_{T} (GeV/#it{c});(#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}})_{Pb#minusPb}/(#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}})_{pp}")
    h_frame.GetYaxis().ChangeLabel(1, 1, 0)

    alice_pbpb_10_20 = draw_alice_pbpb_pp_ratio(0, 10, 70, 100, c=colors[1], s=2.5, m=ROOT.kFullDiamond)
    alice_pbpb_40_50 = draw_alice_pbpb_pp_ratio(40, 50, 70, 100, c=colors[2], s=2., m=ROOT.kFullSquare)
    alice_pbpb_80_90 = draw_alice_pbpb_pp_ratio(80, 90, 70, 100, c=colors[0], s=2., m=ROOT.kFullCrossX)

    x, y = to_pad_coordinates(0.05, 0.9)
    alice_text = ROOT.TLatex(x, y, 'ALICE Preliminary')
    alice_text.SetNDC()
    alice_text.SetTextFont(43)
    alice_text.SetTextSize(33)
    alice_text.Draw()

    x, y = to_pad_coordinates(0.05, 0.25)
    legend_text.DrawLatexNDC(x, y, 'Pb#minusPb, #kern[-0.1]{#sqrt{#it{s}_{NN}} = 5.36 TeV}')

    x_min, y_min = to_pad_coordinates(0.03, 0.03)
    x_max, y_max = to_pad_coordinates(0.5, 0.23)
    legend = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(43)
    legend.SetTextSize(25)
    legend.AddEntry(alice_pbpb_10_20, '0#minus10% FT0C', 'pl')
    legend.AddEntry(alice_pbpb_40_50, '40#minus50% FT0C', 'pl')
    legend.AddEntry(alice_pbpb_80_90, '80#minus90% FT0C', 'pl')
    legend.Draw()

    x, y = to_pad_coordinates(0.52, 0.82)
    legend_text.DrawLatexNDC(x, y, 'pp, #kern[-0.1]{#sqrt{#it{s}} = 13.6 TeV}')
    x, y = to_pad_coordinates(0.52, 0.77)
    legend_text.DrawLatexNDC(x, y, '70#minus100% FT0M')

    # x, y = to_pad_coordinates(0.5, 0.05)
    # br_unc_text.DrawLatexNDC(x, y, '#lower[-0.03]{^{+4.0}}')
    # br_unc_text.DrawLatexNDC(x, y, '_{#minus3.8}')
    # br_unc_text.DrawLatexNDC(x, y, '#kern[.12]{% BR uncertainty not shown}')

    ROOT.gPad.RedrawAxis()

    c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/double_ratio_alice_vs_pt.pdf")
    c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/double_ratio_alice_vs_pt.root")