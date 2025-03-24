import ROOT
import numpy as np
import matplotlib as mpl

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
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/ds_over_dplus_vs_mult_run2.root") as infile:
        graph_stat = infile.Get(f"pt_{pt_min}_{pt_max}/graph_stat_pbpb2018_pt{pt_min}_{pt_max}")
        graph_syst = infile.Get(f"pt_{pt_min}_{pt_max}/graph_syst_pbpb2018_pt{pt_min}_{pt_max}")
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
    return graph_stat, graph_syst


def draw_alice_pp(pt_min, pt_max, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Ratios/VsMult/w_syst/ratio_with_syst.root") as infile:
        graph_stat = infile.Get(f"g_stat_{pt_min*10:.0f}_{pt_max*10:.0f}")
        graph_syst = infile.Get(f"g_syst_{pt_min*10:.0f}_{pt_max*10:.0f}")
    
    #set x unc. for the systematic box
    for i in range(graph_syst.GetN()):
        graph_syst.SetPointEXlow(i, graph_syst.GetPointX(i)*0.075)
        graph_syst.SetPointEXhigh(i, graph_syst.GetPointX(i)*0.075)

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

def draw_pythia_pbpb_sccr_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/pbpb/CharmHadRatiosVsMult_merge_cut.root") as infile:
        hist_mid = infile.Get(f"hDsDpRatioMidMultPt{pt_min:.0f}{pt_max:.0f}")
    hist_mid.SetMarkerStyle(m)
    hist_mid.SetMarkerSize(s)
    hist_mid.SetMarkerColor(c)
    hist_mid.SetLineColor(c)
    hist_mid.SetLineWidth(1)
    hist_mid.SetFillStyle(1001)
    hist_mid.SetFillColorAlpha(c, a)
    hist_mid.Draw("E3, same")
    hist_mid.Draw("L, same")
    return hist_mid

def draw_pythia_pbpb_sccr_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/pbpb/CharmHadRatiosVsMult_merge_cut.root") as infile:
        hist_ft0 = infile.Get(f"hDsDpRatioFwdMultPt{pt_min:.0f}{pt_max:.0f}")
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

if __name__ == '__main__':
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadRightMargin(0.03)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetOptLogx(1)
    ROOT.gStyle.SetLabelSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleOffset(1.1, "X")
    colors, _ = get_discrete_matplotlib_palette('tab20')
    alice_text = ROOT.TLatex(0.18, 0.88, 'ALICE Preliminary')
    alice_text.SetNDC()
    alice_text.SetTextFont(42)
    alice_text.SetTextSize(0.05)

    y_text = ROOT.TLatex(0.185, 0.83, '|#it{y}| < 0.5')
    y_text.SetNDC()
    y_text.SetTextFont(42)
    y_text.SetTextSize(0.04)

    pt_text = ROOT.TLatex(0.6, 0.88, '')
    pt_text.SetNDC()
    pt_text.SetTextFont(42)
    pt_text.SetTextSize(0.04)

    c = ROOT.TCanvas("canvas", "canvas", 2400, 1600)
    c.Divide(3, 2, 0.0001, 0.0001)
    c.cd(1)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    pythia_ft0m = draw_pythia_pp_ft0m(1, 2, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    pythia_ft0m_mod = draw_pythia_pp_mod_ft0m(1, 2, c=colors[0], s=3, m=ROOT.kFullDiamond)
    pythia_ft0m_mode2 = draw_pythia_pp_mode2_ft0m(1, 2, c=colors[4], s=2, m=ROOT.kFullSquare)
    pythia_ft0m_mod_mode2 = draw_pythia_pp_mod_mode2_ft0m(1, 2, c=colors[8], s=2.5, m=ROOT.kFullCross)
    pythia_ft0m_sccr = draw_pythia_pp_sccr_ft0m(1, 2, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    epos_pp = draw_epos4hq_pp_ft0m(1, 2, c=colors[12], s=0, m=ROOT.kFullFourTrianglesPlus)
    pythia_pbpb_ft0m_sccr = draw_pythia_pbpb_sccr_ft0m(1, 2, c=colors[14], s=0, m=ROOT.kFullFourTrianglesPlus)
    g_alice_pp_1_2 = draw_alice_pp(1, 2, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '1#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{2} GeV/#it{c}')

    legend_alice_pp = ROOT.TLegend(0.175, 0.76, 0.485, 0.86)
    legend_alice_pp.SetBorderSize(0)
    legend_alice_pp.SetFillStyle(0)
    legend_alice_pp.SetTextFont(42)
    legend_alice_pp.SetTextSize(0.04)
    legend_alice_pp.AddEntry(g_alice_pp_1_2, '#splitline{pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5}{Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}}', 'pl')
    # legend_alice_pp.AddEntry('', '#splitline{pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5}{Mult. estim.: |#it{#eta}| < 0.8}', 'pl')
    legend_alice_pp.Draw()

    pythia_header = ROOT.TLatex()
    pythia_header.SetNDC()
    pythia_header.SetTextFont(42)
    pythia_header.SetTextSize(0.04)
    pythia_header.DrawLatex(0.185, 0.27, 'PYTHIA 8, pp')
    pythia_header.DrawLatex(0.185, 0.23, 'Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}')

    legend_pythia_no_mod = ROOT.TLegend(0.175, 0.14, 0.95, 0.22)
    legend_pythia_no_mod.SetBorderSize(0)
    legend_pythia_no_mod.SetFillStyle(0)
    legend_pythia_no_mod.SetTextFont(42)
    legend_pythia_no_mod.SetNColumns(3)
    legend_pythia_no_mod.SetTextSize(0.04)
    legend_pythia_no_mod.AddEntry(pythia_ft0m, 'Monash', 'fl')
    legend_pythia_no_mod.AddEntry(pythia_ft0m_mode2, 'CL-BLC Mode 2', 'fl')
    legend_pythia_no_mod.AddEntry(pythia_ft0m_sccr, 'SC#minusCR', 'fl')
    legend_pythia_no_mod.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(2)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    graph_stat_alice_2_4, graph_syst_alice_2_4 = draw_alice_pbpb(2, 4, c=colors[6], s=2.5, m=ROOT.kOpenCircle)
    draw_pythia_pp_ft0m(2, 4, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mod_ft0m(2, 4, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mode2_ft0m(2, 4, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mod_mode2_ft0m(2, 4, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_sccr_ft0m(2, 4, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(2, 4, c=colors[14], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_epos4hq_pp_ft0m(2, 4, c=colors[12], s=0, m=ROOT.kFullFourTrianglesPlus)
    g_alice_pp_2_4 = draw_alice_pp(2, 4, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '2#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{4} GeV/#it{c}')

    pythia_mod_header = ROOT.TLatex()
    pythia_mod_header.SetNDC()
    pythia_mod_header.SetTextFont(42)
    pythia_mod_header.SetTextSize(0.04)
    pythia_mod_header.DrawLatex(0.185, 0.31, 'PYTHIA 8, pp')
    pythia_mod_header.DrawLatex(0.185, 0.27, 'StringFlav:mesonCvector = 1.75')
    pythia_mod_header.DrawLatex(0.185, 0.23, 'Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}')

    legend_pythia_mod = ROOT.TLegend(0.175, 0.14, 0.95, 0.22)
    legend_pythia_mod.SetBorderSize(0)
    legend_pythia_mod.SetFillStyle(0)
    legend_pythia_mod.SetTextFont(42)
    legend_pythia_mod.SetNColumns(3)
    legend_pythia_mod.SetTextSize(0.04)
    legend_pythia_mod.AddEntry(pythia_ft0m_mod, 'Monash', 'fl')
    legend_pythia_mod.AddEntry(pythia_ft0m_mod_mode2, 'CL-BLC Mode 2', 'fl')
    legend_pythia_mod.Draw("same")

    ROOT.gPad.RedrawAxis()

    c.cd(3)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    draw_pythia_pp_ft0m(4, 6, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mod_ft0m(4, 6, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mode2_ft0m(4, 6, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mod_mode2_ft0m(4, 6, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_sccr_ft0m(4, 6, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(4, 6, c=colors[14], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_epos4hq_pp_ft0m(4, 6, c=colors[12], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_alice_pbpb(4, 6, c=colors[6], s=2.5, m=ROOT.kOpenCircle)
    draw_alice_pp(4, 6, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '4#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{6} GeV/#it{c}')

    legend_alice_pbpb = ROOT.TLegend(0.175, 0.76, 0.485, 0.86)
    legend_alice_pbpb.SetBorderSize(0)
    legend_alice_pbpb.SetFillStyle(0)
    legend_alice_pbpb.SetTextFont(42)
    legend_alice_pbpb.SetTextSize(0.04)
    legend_alice_pbpb.AddEntry(graph_stat_alice_2_4, '#splitline{ALICE, Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{y}| < 0.5}{Mult. estim.: #kern[-0.3]{#font[122]{-}3.7} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}1.7} #kern[-0.9]{#vee} 2.8 < #kern[-0.2]{#it{#eta} < 5.1}}', 'pl')
    legend_alice_pbpb.Draw()

    pythia_pbpb_header = ROOT.TLatex()
    pythia_pbpb_header.SetNDC()
    pythia_pbpb_header.SetTextFont(42)
    pythia_pbpb_header.SetTextSize(0.04)
    pythia_pbpb_header.DrawLatex(0.185, 0.27, 'PYTHIA 8 Angantyr, Pb#minusPb')
    pythia_pbpb_header.DrawLatex(0.185, 0.23, 'Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}')

    legend_pythia_pbpb = ROOT.TLegend(0.175, 0.14, 0.95, 0.22)
    legend_pythia_pbpb.SetBorderSize(0)
    legend_pythia_pbpb.SetFillStyle(0)
    legend_pythia_pbpb.SetTextFont(42)
    legend_pythia_pbpb.SetNColumns(3)
    legend_pythia_pbpb.SetTextSize(0.04)
    legend_pythia_pbpb.AddEntry(pythia_pbpb_ft0m_sccr, 'SC#minusCR', 'fl')
    legend_pythia_pbpb.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(4)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    draw_pythia_pp_ft0m(6, 8, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mod_ft0m(6, 8, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mode2_ft0m(6, 8, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mod_mode2_ft0m(6, 8, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_sccr_ft0m(6, 8, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(6, 8, c=colors[14], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_epos4hq_pp_ft0m(6, 8, c=colors[12], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_alice_pbpb(6, 8, c=colors[6], s=2.5, m=ROOT.kOpenCircle)
    draw_alice_pp(6, 8, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '6#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{8} GeV/#it{c}')

    epos4hq_header = ROOT.TLatex()
    epos4hq_header.SetNDC()
    epos4hq_header.SetTextFont(42)
    epos4hq_header.SetTextSize(0.04)
    epos4hq_header.DrawLatex(0.185, 0.23, 'Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}')

    legend_epos4hq = ROOT.TLegend(0.175, 0.14, 0.95, 0.22)
    legend_epos4hq.SetBorderSize(0)
    legend_epos4hq.SetFillStyle(0)
    legend_epos4hq.SetTextFont(42)
    legend_epos4hq.SetNColumns(3)
    legend_epos4hq.SetTextSize(0.04)
    legend_epos4hq.AddEntry(epos_pp, 'EPOS4HQ, pp', 'fl')
    legend_epos4hq.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(5)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    draw_pythia_pp_ft0m(8, 12, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mod_ft0m(8, 12, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mode2_ft0m(8, 12, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mod_mode2_ft0m(8, 12, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_sccr_ft0m(8, 12, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(8, 12, c=colors[14], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_epos4hq_pp_ft0m(8, 12, c=colors[12], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_alice_pbpb(8, 12, c=colors[6], s=2.5, m=ROOT.kOpenCircle)
    draw_alice_pp(8, 12, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '8#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[.5]{12} GeV/#it{c}')

    ROOT.gPad.RedrawAxis()

    c.cd(6)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    draw_pythia_pp_ft0m(12, 24, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mod_ft0m(12, 24, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mode2_ft0m(12, 24, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mod_mode2_ft0m(12, 24, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_sccr_ft0m(12, 24, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(12, 24, c=colors[14], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_epos4hq_pp_ft0m(12, 24, c=colors[12], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_alice_pp(12, 24, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '12#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[.5]{24} GeV/#it{c}')

    ROOT.gPad.RedrawAxis()

    c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/ratio_vs_pred_2023.pdf")
