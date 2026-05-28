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

def draw_alice_pbpb(cent_min, cent_max, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle):
    g_syst = None
    if cent_min < 20:
        with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Ratios/ratio_w_syst.root") as infile:
            h_stat = infile.Get(f"h_stat_vs_pt_{cent_min:.0f}_{cent_max:.0f}")
            h_stat.SetDirectory(0)
            g_syst = infile.Get(f"g_syst_vs_pt_no_br_{cent_min:.0f}_{cent_max:.0f}")
    elif cent_min < 50:
        with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/Data/ds_over_dplus_ratios_2050.root") as infile:
            h_stat = infile.Get(f"h_stat_vs_pt_{cent_min:.0f}_{cent_max:.0f}")
            h_stat.SetDirectory(0)
            g_syst = infile.Get(f"g_syst_vs_pt_no_br_{cent_min:.0f}_{cent_max:.0f}")
    elif cent_min < 90:
        with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/Data/ds_over_dplus_ratios_5090.root") as infile:
            h_stat = infile.Get(f"h_stat_vs_pt_{cent_min:.0f}_{cent_max:.0f}")
            h_stat.SetDirectory(0)
            g_syst = infile.Get(f"g_syst_vs_pt_no_br_{cent_min:.0f}_{cent_max:.0f}")
    
    if cent_min >= 70:
        #remove last two points for 70-80% and 80-90%
        h_stat.SetBinContent(6, 1.e9)
        if g_syst:
            g_syst.SetPoint(5, 0, 1.e9)

    h_stat.SetMarkerStyle(m)
    h_stat.SetMarkerSize(s)
    h_stat.SetMarkerColor(c)
    h_stat.SetLineColor(c)
    h_stat.SetLineWidth(2)
    h_stat.Draw("pz, same")
    if g_syst:
        #set x unc. for the systematic box
        for i in range(g_syst.GetN()):
            g_syst.SetPointEXlow(i, g_syst.GetErrorXlow(i)/2)
            g_syst.SetPointEXhigh(i, g_syst.GetErrorXhigh(i)/2)
        g_syst.SetMarkerStyle(m)
        g_syst.SetMarkerSize(s)
        g_syst.SetMarkerColor(c)
        g_syst.SetLineColor(c)
        g_syst.SetLineWidth(2)
        g_syst.SetFillStyle(0)
        g_syst.Draw("E2, same")
    return h_stat#, graph_stat_50_90



def draw_alice_pp(cent_min, cent_max, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Ratios/VsMult/w_syst/ratio_w_syst.root") as infile:
        graph_stat = infile.Get(f"h_stat_vs_pt_{cent_min:.0f}_{cent_max:.0f}")
        graph_stat.SetDirectory(0)
        graph_syst = infile.Get(f"g_syst_vs_pt_no_br_{cent_min:.0f}_{cent_max:.0f}")
    
    #set x unc. for the systematic box
    for i in range(graph_syst.GetN()):
        graph_syst.SetPointEXlow(i, graph_syst.GetErrorXlow(i)/2)
        graph_syst.SetPointEXhigh(i, graph_syst.GetErrorXhigh(i)/2)

    graph_stat.SetMarkerStyle(m)
    graph_stat.SetMarkerSize(s)
    graph_stat.SetMarkerColor(c)
    graph_stat.SetLineColor(c)
    graph_stat.SetLineWidth(2)
    graph_syst.SetMarkerStyle(m)
    graph_syst.SetMarkerSize(s)
    graph_syst.SetMarkerColor(c)
    graph_syst.SetLineColor(c)
    graph_syst.SetLineWidth(2)
    graph_syst.SetFillStyle(0)
    graph_stat.DrawClone("pz, same")
    graph_syst.DrawClone("pz2, same")
    return graph_stat #, graph_syst

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
    h_frame = c.DrawFrame(0., 0., 24., 1., ";#it{p}_{T} (GeV/#it{c});#it{#sigma}(D_{s}^{+})/#it{#sigma}(D^{+})")
    h_frame.GetYaxis().ChangeLabel(1, 1, 0)

    alice_pp = draw_alice_pp(70, 100, c=ROOT.kBlack, s=1.5, m=ROOT.kFullCircle)
    alice_pbpb_10_20 = draw_alice_pbpb(0, 10, c=ROOT.kRed+1, s=2., m=ROOT.kFullDiamond)
    alice_pbpb_40_50 = draw_alice_pbpb(40, 50, c=colors[1], s=1.3, m=ROOT.kFullSquare)
    alice_pbpb_80_90 = draw_alice_pbpb(80, 90, c=colors[0], s=1.5, m=ROOT.kFullCrossX)

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
    legend.SetMargin(0.15)
    legend.AddEntry(alice_pbpb_10_20, '0#minus10% FT0C', 'p')
    legend.AddEntry(alice_pbpb_40_50, '40#minus50% FT0C', 'p')
    legend.AddEntry(alice_pbpb_80_90, '80#minus90% FT0C', 'p')
    legend.Draw()

    x, y = to_pad_coordinates(0.52, 0.82)
    legend_text.DrawLatexNDC(x, y, 'pp, #kern[-0.1]{#sqrt{#it{s}} = 13.6 TeV}')

    x_min, y_min = to_pad_coordinates(0.5, 0.75)
    x_max, y_max = to_pad_coordinates(1., 0.8)
    legend_pp = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_pp.SetBorderSize(0)
    legend_pp.SetFillStyle(0)
    legend_pp.SetTextFont(43)
    legend_pp.SetTextSize(25)
    legend_pp.SetMargin(0.15)
    legend_pp.AddEntry(alice_pp, '70#minus100% FT0M', 'p')
    legend_pp.Draw()


    x, y = to_pad_coordinates(0.5, 0.05)
    br_unc_text.DrawLatexNDC(x, y, '#lower[-0.03]{^{+4.0}}')
    br_unc_text.DrawLatexNDC(x, y, '_{#minus3.8}')
    br_unc_text.DrawLatexNDC(x, y, '#kern[.12]{% BR uncertainty not shown}')

    ROOT.gPad.RedrawAxis()

    c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_alice_vs_pt.pdf")
    c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_alice_vs_pt.root")