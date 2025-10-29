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

def draw_alice_pbpb(pt_min, pt_max, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Ratios/ds_over_dplus_ratios_020_norm.root") as infile:
        graph_stat_0_20 = infile.Get(f"g_ratio_dndeta_{pt_min*10:.0f}_{pt_max*10:.0f}")
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Ratios/ds_over_dplus_ratios_5090_norm.root") as infile:
        graph_stat_50_90 = infile.Get(f"g_ratio_dndeta_{pt_min*10:.0f}_{pt_max*10:.0f}")

    if (pt_min, pt_max) == (12, 24):
        #remove first two points for 50-90%
        graph_stat_50_90.SetPoint(2, 0, 1.e9)
        graph_stat_50_90.SetPoint(3, 0, 1.e9)

    graph_stat_0_20.SetMarkerStyle(m)
    graph_stat_0_20.SetMarkerSize(s)
    graph_stat_0_20.SetMarkerColor(c)
    graph_stat_0_20.SetLineColor(c)
    graph_stat_0_20.SetLineWidth(1)
    graph_stat_0_20.Draw("pz, same")
    graph_stat_50_90.SetMarkerStyle(m)
    graph_stat_50_90.SetMarkerSize(s)
    graph_stat_50_90.SetMarkerColor(c)
    graph_stat_50_90.SetLineColor(c)
    graph_stat_50_90.SetLineWidth(1)
    graph_stat_50_90.Draw("pz, same")
    return graph_stat_0_20#, graph_stat_50_90


def draw_alice_pp(pt_min, pt_max, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Ratios/VsMult/w_syst/ratio_with_syst.root") as infile:
        graph_stat = infile.Get(f"g_stat_{pt_min*10:.0f}_{pt_max*10:.0f}")
        graph_syst = infile.Get(f"g_syst_{pt_min*10:.0f}_{pt_max*10:.0f}")
    
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

    y_text = ROOT.TLatex(0.185, 0.83, '|#it{y}| < 0.5')
    y_text.SetNDC()
    y_text.SetTextFont(42)
    y_text.SetTextSize(0.04)

    pt_text = ROOT.TLatex(0.6, 0.88, '')
    pt_text.SetNDC()
    pt_text.SetTextFont(43)
    pt_text.SetTextSize(40)

    c = ROOT.TCanvas("canvas", "canvas", 2400, 1600)
    c.Divide(3, 2, 0.000, 0.000)
    # Set margins
    for i_pad in range(1, 4):
        pad = c.cd(i_pad)
        pad.SetTopMargin(0.02)
    for i_pad in range(4, 7):
        pad = c.cd(i_pad)
        pad.SetBottomMargin(0.16)

    c.cd(1)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{y}|<0.5} ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    h_frame.GetYaxis().ChangeLabel(1, 1, 0)

    alice_pp = draw_alice_pp(1, 2, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    alice_pbpb = draw_alice_pbpb(1, 2, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)

    x, y = to_pad_coordinates(0.05, 0.9)
    alice_text = ROOT.TLatex(x, y, 'ALICE Work in Progress')
    alice_text.SetNDC()
    alice_text.SetTextFont(43)
    alice_text.SetTextSize(50)
    alice_text.Draw()

    x, y = to_pad_coordinates(0.05, 0.85)
    pt_text.DrawLatexNDC(x, y, '1#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{2} GeV/#it{c}')

    ROOT.gPad.RedrawAxis()

    c.cd(2)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{y}|<0.5} ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    draw_alice_pp(2, 4, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    draw_alice_pbpb(2, 4, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)
    #alice_text.Draw()
    #y_text.Draw()
    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '2#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{4} GeV/#it{c}')

    ROOT.gPad.RedrawAxis()

    c.cd(3)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{y}|<0.5} ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    draw_alice_pp(4, 6, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    draw_alice_pbpb(4, 6, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)

    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '4#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{6} GeV/#it{c}')

    ROOT.gPad.RedrawAxis()

    c.cd(4)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{y}|<0.5} ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    h_frame.GetYaxis().ChangeLabel(10, 1, 0)
    draw_alice_pp(6, 8, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    draw_alice_pbpb(6, 8, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)
    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '6#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{8} GeV/#it{c}')

    x_min, y_min = to_pad_coordinates(0.03, 0.07)
    x_max, y_max = to_pad_coordinates(0.5, 0.17)
    legend_alice_pp = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_alice_pp.SetBorderSize(0)
    legend_alice_pp.SetFillStyle(0)
    legend_alice_pp.SetTextFont(43)
    legend_alice_pp.SetTextSize(40)
    legend_alice_pp.AddEntry(alice_pp, '#splitline{pp,#kern[0.2]{#sqrt{#it{s}}} = 13.6 TeV, |#it{y}| < 0.5}{FT0M multiplicity estimator}', 'pl')
    legend_alice_pp.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(5)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{y}|<0.5} ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    h_frame.GetXaxis().ChangeLabel(1, 1, 0)
    draw_alice_pp(8, 12, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    draw_alice_pbpb(8, 12, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)
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
    legend_alice_pbpb.AddEntry(alice_pbpb, '#splitline{Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.36 TeV, |#it{y}| < 0.5}{FT0C multiplicity estimator}', 'pl')
    legend_alice_pbpb.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(6)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{y}|<0.5} ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    h_frame.GetXaxis().ChangeLabel(1, 1, 0)
    draw_alice_pp(12, 24, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    draw_alice_pbpb(12, 24, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)
    #alice_text.Draw()
    #y_text.Draw()
    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '12#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[.5]{24} GeV/#it{c}')

    ROOT.gPad.RedrawAxis()

    c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_alice.pdf")
    c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_alice.root")