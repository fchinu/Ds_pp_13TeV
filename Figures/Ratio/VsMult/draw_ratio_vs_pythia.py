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
    graph_mid.SetLineWidth(0)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3, same")
    return graph_mid

def draw_pythia_pp_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_Monash")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(0)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3, same")
    return graph_ft0

def draw_pythia_pp_mod_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_MonashModDstar")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(0)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3, same")
    return graph_mid

def draw_pythia_pp_mod_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_MonashModDstar")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(0)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3, same")
    return graph_ft0

def draw_pythia_pp_mode2_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_Mode2")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(0)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3, same")
    return graph_mid

def draw_pythia_pp_mode2_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_Mode2")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(0)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3, same")
    return graph_ft0

def draw_pythia_pp_mod_mode2_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_Mode2ModDstar")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(0)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3, same")
    return graph_mid

def draw_pythia_pp_mod_mode2_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_Mode2ModDstar")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(0)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3, same")
    return graph_ft0

def draw_pythia_pp_sccr_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_mid = infile.Get(f"graph_ds_over_dp_p_mid_pt{pt_min:.1f}_{pt_max:.1f}_SRRC")
    graph_mid.SetMarkerStyle(m)
    graph_mid.SetMarkerSize(s)
    graph_mid.SetMarkerColor(c)
    graph_mid.SetLineColor(c)
    graph_mid.SetLineWidth(0)
    graph_mid.SetFillStyle(1001)
    graph_mid.SetFillColorAlpha(c, a)
    graph_mid.Draw("3, same")
    return graph_mid

def draw_pythia_pp_sccr_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/PYTHIA8_Ds_over_Dplus_pp13dot6TeV_vsMult.root") as infile:
        graph_ft0 = infile.Get(f"graph_ds_over_dp_p_ft0m_pt{pt_min:.1f}_{pt_max:.1f}_SRRC")
    graph_ft0.SetMarkerStyle(m)
    graph_ft0.SetMarkerSize(s)
    graph_ft0.SetMarkerColor(c)
    graph_ft0.SetLineColor(c)
    graph_ft0.SetLineWidth(0)
    graph_ft0.SetFillStyle(1001)
    graph_ft0.SetFillColorAlpha(c, a)
    graph_ft0.Draw("3, same")
    return graph_ft0

def draw_pythia_pbpb_sccr_mid(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kOpenDiamond):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/pbpb/CharmHadRatiosVsMult_merge_cut.root") as infile:
        hist_mid = infile.Get(f"hDsDpRatioMidMultPt{pt_min:.0f}{pt_max:.0f}")
        hist_mid.SetDirectory(0)
    # Divide x values by eta gap
    deta = (3.3-2.1) + (4.9-3.5) 
    for i in range(hist_mid.GetN()):
        x, y = ROOT.Double(0), ROOT.Double(0)
        hist_mid.GetPoint(i, x, y)
        hist_mid.SetPoint(i, x/deta, y)
        hist_mid.SetPointEXlow(i, 0)
        hist_mid.SetPointEXhigh(i, 0)
    hist_mid.SetMarkerStyle(m)
    hist_mid.SetMarkerSize(s)
    hist_mid.SetMarkerColor(c)
    hist_mid.SetLineColor(c)
    hist_mid.SetLineWidth(0)
    hist_mid.SetFillStyle(1001)
    hist_mid.SetFillColorAlpha(c, a)
    hist_mid.Draw("E3, same")
    return hist_mid

def draw_pythia_pbpb_sccr_ft0m(pt_min, pt_max, c=ROOT.kBlack, s=3, a=0.5, m=ROOT.kFullCircle):
    with ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/pbpb/CharmHadRatiosVsMult_merge_cut.root") as infile:
        hist_ft0 = infile.Get(f"hDsDpRatioFwdMultPt{pt_min:.0f}{pt_max:.0f}")
        hist_ft0.SetDirectory(0)
    deta = (3.3-2.1) + (4.9-3.5) 
    for i in range(hist_ft0.GetN()):
        x, y = ROOT.Double(0), ROOT.Double(0)
        hist_ft0.GetPoint(i, x, y)
        hist_ft0.SetPoint(i, x/deta, y)
        hist_ft0.SetPointEXlow(i, 0)
        hist_ft0.SetPointEXhigh(i, 0)
    hist_ft0.SetMarkerStyle(m)
    hist_ft0.SetMarkerSize(s)
    hist_ft0.SetMarkerColor(c)
    hist_ft0.SetLineColor(c)
    hist_ft0.SetLineWidth(0)
    hist_ft0.SetFillStyle(1001)
    hist_ft0.SetFillColorAlpha(c, a)
    hist_ft0.Draw("E3, same")
    return hist_ft0

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
    alice_text = ROOT.TLatex(0.18, 0.88, 'This work')
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
    pythia_mid = draw_pythia_pp_mid(1, 2, c=colors[3], s=3, m=ROOT.kOpenDoubleDiamond)
    pythia_ft0m_mod = draw_pythia_pp_mod_ft0m(1, 2, c=colors[0], s=3, m=ROOT.kFullDiamond)
    pythia_mid_mod = draw_pythia_pp_mod_mid(1, 2, c=colors[1], s=3, m=ROOT.kOpenDiamond)
    pythia_ft0m_mode2 = draw_pythia_pp_mode2_ft0m(1, 2, c=colors[4], s=2, m=ROOT.kFullSquare)
    pythia_mid_mode2 = draw_pythia_pp_mode2_mid(1, 2, c=colors[5], s=2, m=ROOT.kOpenSquare)
    pythia_ft0m_mod_mode2 = draw_pythia_pp_mod_mode2_ft0m(1, 2, c=colors[8], s=2.5, m=ROOT.kFullCross)
    pythia_mid_mod_mode2 = draw_pythia_pp_mod_mode2_mid(1, 2, c=colors[9], s=2, m=ROOT.kOpenCross)
    pythia_ft0m_sccr = draw_pythia_pp_sccr_ft0m(1, 2, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    pythia_mid_sccr = draw_pythia_pp_sccr_mid(1, 2, c=colors[11], s=2, m=ROOT.kOpenFourTrianglesPlus)
    pythia_pbpb_ft0m_sccr = draw_pythia_pbpb_sccr_ft0m(1, 2, c=colors[11], s=0, m=ROOT.kFullFourTrianglesPlus)
    pythia_pbpb_mid_sccr = draw_pythia_pbpb_sccr_mid(1, 2, c=colors[10], s=0, m=ROOT.kFullFourTrianglesPlus)
    g_alice_pp_1_2 = draw_alice_pp(1, 2, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '1 < #it{p}_{T} < 2 GeV/#it{c}')

    legend_alice_pp = ROOT.TLegend(0.175, 0.76, 0.485, 0.86)
    legend_alice_pp.SetBorderSize(0)
    legend_alice_pp.SetFillStyle(0)
    legend_alice_pp.SetTextFont(42)
    legend_alice_pp.SetTextSize(0.04)
    legend_alice_pp.AddEntry(g_alice_pp_1_2, '#splitline{pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5}{Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}}', 'pl')
    # legend_alice_pp.AddEntry('', '#splitline{pp, #sqrt{#it{s}} = 13.6 TeV, |#it{y}| < 0.5}{Mult. estim.: |#it{#eta}| < 0.8}', 'pl')
    legend_alice_pp.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(2)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    graph_stat_alice_2_4, graph_syst_alice_2_4 = draw_alice_pbpb(2, 4, c=colors[6], s=2.5, m=ROOT.kOpenCircle)
    draw_pythia_pp_ft0m(2, 4, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mid(2, 4, c=colors[3], s=3, m=ROOT.kOpenDoubleDiamond)
    draw_pythia_pp_mod_ft0m(2, 4, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mod_mid(2, 4, c=colors[1], s=3, m=ROOT.kOpenDiamond)
    draw_pythia_pp_mode2_ft0m(2, 4, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mode2_mid(2, 4, c=colors[5], s=2, m=ROOT.kOpenSquare)
    draw_pythia_pp_mod_mode2_ft0m(2, 4, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_mod_mode2_mid(2, 4, c=colors[9], s=2, m=ROOT.kOpenCross)
    draw_pythia_pp_sccr_ft0m(2, 4, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pp_sccr_mid(2, 4, c=colors[11], s=2, m=ROOT.kOpenFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(2, 4, c=colors[11], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_mid(2, 4, c=colors[10], s=0, m=ROOT.kFullFourTrianglesPlus)
    g_alice_pp_2_4 = draw_alice_pp(2, 4, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '2 < #it{p}_{T} < 4 GeV/#it{c}')

    legend_pythia_sccr = ROOT.TLegend(0.175, 0.14, 0.485, 0.34)
    legend_pythia_sccr.SetBorderSize(0)
    legend_pythia_sccr.SetFillStyle(0)
    legend_pythia_sccr.SetTextFont(42)
    legend_pythia_sccr.SetTextSize(0.04)
    legend_pythia_sccr.SetHeader("PYTHIA 8 with SC-CR")
    legend_pythia_sccr.AddEntry(pythia_ft0m_sccr, 'Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}', 'f')
    legend_pythia_sccr.AddEntry(pythia_mid_sccr, 'Mult. estim.: |#it{#eta}| < 0.8', 'f')
    legend_pythia_sccr.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(3)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    draw_pythia_pp_ft0m(4, 6, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mid(4, 6, c=colors[3], s=3, m=ROOT.kOpenDoubleDiamond)
    draw_pythia_pp_mod_ft0m(4, 6, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mod_mid(4, 6, c=colors[1], s=3, m=ROOT.kOpenDiamond)
    draw_pythia_pp_mode2_ft0m(4, 6, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mode2_mid(4, 6, c=colors[5], s=2, m=ROOT.kOpenSquare)
    draw_pythia_pp_mod_mode2_ft0m(4, 6, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_mod_mode2_mid(4, 6, c=colors[9], s=2, m=ROOT.kOpenCross)
    draw_pythia_pp_sccr_ft0m(4, 6, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pp_sccr_mid(4, 6, c=colors[11], s=2, m=ROOT.kOpenFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(4, 6, c=colors[11], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_mid(4, 6, c=colors[10], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_alice_pbpb(4, 6, c=colors[6], s=2.5, m=ROOT.kOpenCircle)
    draw_alice_pp(4, 6, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '4 < #it{p}_{T} < 6 GeV/#it{c}')

    legend_alice_pbpb = ROOT.TLegend(0.175, 0.76, 0.485, 0.86)
    legend_alice_pbpb.SetBorderSize(0)
    legend_alice_pbpb.SetFillStyle(0)
    legend_alice_pbpb.SetTextFont(42)
    legend_alice_pbpb.SetTextSize(0.04)
    legend_alice_pbpb.AddEntry(graph_stat_alice_2_4, '#splitline{ALICE, Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{y}| < 0.5}{Mult. estim.: #kern[-0.3]{#font[122]{-}3.7} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}1.7} #kern[-0.9]{#vee} 2.8 < #kern[-0.2]{#it{#eta} < 5.1}}', 'pl')
    legend_alice_pbpb.Draw()

    legend_pythia_mod_mode2 = ROOT.TLegend(0.175, 0.14, 0.485, 0.34)
    legend_pythia_mod_mode2.SetBorderSize(0)
    legend_pythia_mod_mode2.SetFillStyle(0)
    legend_pythia_mod_mode2.SetTextFont(42)
    legend_pythia_mod_mode2.SetTextSize(0.04)
    legend_pythia_mod_mode2.SetHeader("#lower[-0.1]{#splitline{PYTHIA 8 with CL-BLC Mode 2}{StringFlav:mesonCvector = 1.75}}")
    legend_pythia_mod_mode2.AddEntry(pythia_ft0m_mod_mode2, 'Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}', 'f')
    legend_pythia_mod_mode2.AddEntry(pythia_mid_mod_mode2, 'Mult. estim.: |#it{#eta}| < 0.8', 'f')
    legend_pythia_mod_mode2.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(4)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    draw_pythia_pp_ft0m(6, 8, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mid(6, 8, c=colors[3], s=3, m=ROOT.kOpenDoubleDiamond)
    draw_pythia_pp_mod_ft0m(6, 8, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mod_mid(6, 8, c=colors[1], s=3, m=ROOT.kOpenDiamond)
    draw_pythia_pp_mode2_ft0m(6, 8, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mode2_mid(6, 8, c=colors[5], s=2, m=ROOT.kOpenSquare)
    draw_pythia_pp_mod_mode2_ft0m(6, 8, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_mod_mode2_mid(6, 8, c=colors[9], s=2, m=ROOT.kOpenCross)
    draw_pythia_pp_sccr_ft0m(6, 8, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pp_sccr_mid(6, 8, c=colors[11], s=2, m=ROOT.kOpenFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(6, 8, c=colors[11], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_mid(6, 8, c=colors[10], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_alice_pbpb(6, 8, c=colors[6], s=2.5, m=ROOT.kOpenCircle)
    draw_alice_pp(6, 8, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '6 < #it{p}_{T} < 8 GeV/#it{c}')

    legend_pythia = ROOT.TLegend(0.175, 0.14, 0.485, 0.34)
    legend_pythia.SetBorderSize(0)
    legend_pythia.SetFillStyle(0)
    legend_pythia.SetTextFont(42)
    legend_pythia.SetTextSize(0.04)
    legend_pythia.SetHeader("PYTHIA 8 Monash")
    legend_pythia.AddEntry(pythia_ft0m, 'Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}', 'f')
    legend_pythia.AddEntry(pythia_mid, 'Mult. estim.: |#it{#eta}| < 0.8', 'f')
    legend_pythia.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(5)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    draw_pythia_pp_ft0m(8, 12, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mid(8, 12, c=colors[3], s=3, m=ROOT.kOpenDoubleDiamond)
    draw_pythia_pp_mod_ft0m(8, 12, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mod_mid(8, 12, c=colors[1], s=3, m=ROOT.kOpenDiamond)
    draw_pythia_pp_mode2_ft0m(8, 12, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mode2_mid(8, 12, c=colors[5], s=2, m=ROOT.kOpenSquare)
    draw_pythia_pp_mod_mode2_ft0m(8, 12, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_mod_mode2_mid(8, 12, c=colors[9], s=2, m=ROOT.kOpenCross)
    draw_pythia_pp_sccr_ft0m(8, 12, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pp_sccr_mid(8, 12, c=colors[11], s=2, m=ROOT.kOpenFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(8, 12, c=colors[11], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_mid(8, 12, c=colors[10], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_alice_pbpb(8, 12, c=colors[6], s=2.5, m=ROOT.kOpenCircle)
    draw_alice_pp(8, 12, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '8 < #it{p}_{T} < 12 GeV/#it{c}')

    legend_pythia_mod = ROOT.TLegend(0.175, 0.14, 0.485, 0.34)
    legend_pythia_mod.SetBorderSize(0)
    legend_pythia_mod.SetFillStyle(0)
    legend_pythia_mod.SetTextFont(42)
    legend_pythia_mod.SetTextSize(0.04)
    legend_pythia_mod.SetHeader("#lower[-0.1]{#splitline{PYTHIA 8 Monash}{StringFlav:mesonCvector = 1.75}}")
    legend_pythia_mod.AddEntry(pythia_ft0m_mod, 'Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}', 'f')
    legend_pythia_mod.AddEntry(pythia_mid_mod, 'Mult. estim.: |#it{#eta}| < 0.8', 'f')
    legend_pythia_mod.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(6)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 1., ";d#it{N}_{ch}/d#it{#eta};#sigma_{D_{s}^{+}}/#sigma_{D^{+}}")
    draw_pythia_pp_ft0m(12, 24, c=colors[2], s=3, m=ROOT.kFullDoubleDiamond)
    draw_pythia_pp_mid(12, 24, c=colors[3], s=3, m=ROOT.kOpenDoubleDiamond)
    draw_pythia_pp_mod_ft0m(12, 24, c=colors[0], s=3, m=ROOT.kFullDiamond)
    draw_pythia_pp_mod_mid(12, 24, c=colors[1], s=3, m=ROOT.kOpenDiamond)
    draw_pythia_pp_mode2_ft0m(12, 24, c=colors[4], s=2, m=ROOT.kFullSquare)
    draw_pythia_pp_mode2_mid(12, 24, c=colors[5], s=2, m=ROOT.kOpenSquare)
    draw_pythia_pp_mod_mode2_ft0m(12, 24, c=colors[8], s=2.5, m=ROOT.kFullCross)
    draw_pythia_pp_mod_mode2_mid(12, 24, c=colors[9], s=2, m=ROOT.kOpenCross)
    draw_pythia_pp_sccr_ft0m(12, 24, c=colors[10], s=2, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pp_sccr_mid(12, 24, c=colors[11], s=2, m=ROOT.kOpenFourTrianglesPlus)
    draw_pythia_pbpb_sccr_ft0m(12, 24, c=colors[11], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_pythia_pbpb_sccr_mid(12, 24, c=colors[10], s=0, m=ROOT.kFullFourTrianglesPlus)
    draw_alice_pp(12, 24, c=colors[6], s=2.5, m=ROOT.kFullCircle)
    #alice_text.Draw()
    #y_text.Draw()
    pt_text.DrawLatexNDC(0.6, 0.88, '12 < #it{p}_{T} < 24 GeV/#it{c}')

    legend_pythia_mode2 = ROOT.TLegend(0.175, 0.14, 0.485, 0.34)
    legend_pythia_mode2.SetBorderSize(0)
    legend_pythia_mode2.SetFillStyle(0)
    legend_pythia_mode2.SetTextFont(42)
    legend_pythia_mode2.SetTextSize(0.04)
    legend_pythia_mode2.SetHeader("PYTHIA 8 with CL-BLC Mode 2")
    legend_pythia_mode2.AddEntry(pythia_ft0m_mode2, 'Mult. estim.: #kern[-0.3]{#font[122]{-}3.3} < #kern[-0.9]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.1} #kern[-0.9]{#vee} 3.5 < #kern[-0.2]{#it{#eta} < 4.9}', 'f')
    legend_pythia_mode2.AddEntry(pythia_mid_mode2, 'Mult. estim.: |#it{#eta}| < 0.8', 'f')
    legend_pythia_mode2.Draw()

    ROOT.gPad.RedrawAxis()

    c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Ratio/VsMult/ratio_vs_pythia_mod_2023.pdf")
