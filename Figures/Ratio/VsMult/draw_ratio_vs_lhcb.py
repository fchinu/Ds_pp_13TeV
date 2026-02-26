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

def get_graphs_1():
    x = [0.507467, 1.014934, 1.304915, 1.594896, 1.884877, 2.319849, 3.117297]
    y = [0.5037849, 0.5021154, 0.5206722, 0.540023, 0.5788887, 0.636868, 0.7932779]
    ex = [__convert_ntracks_to_dn_deta_error_backward(ntracks) for ntracks in x]
    x = [__convert_ntracks_to_dn_deta_backward(ntracks) for ntracks in x]
    ey_stat = [0.007068203, 0.009820427, 0.007475892, 0.01087115, 0.01055111, 0.01054322, 0.03467543]
    ey_syst = [0.0190297, 0.02275265, 0.01829271, 0.01770099, 0.02003852, 0.02052603, 0.02703643]
    graph_stat = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_stat))
    ex = (np.array(ex)/2.).tolist()
    graph_syst = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_syst))
    return graph_stat, graph_syst

def get_graphs_2():
    x = [0.5809129, 1.161826, 1.493776, 1.825726, 2.157676, 2.821577]
    y = [0.4600249, 0.4990606, 0.4858854, 0.5110171, 0.520232, 0.5888317]
    ex = [__convert_ntracks_to_dn_deta_error_forward(ntracks) for ntracks in x]
    x = [__convert_ntracks_to_dn_deta_forward(ntracks) for ntracks in x]
    ey_stat = [0.01378267, 0.01297882, 0.01083654, 0.009990128, 0.0132344, 0.01192524]
    ey_syst = [0.01479059, 0.01543968, 0.01504221, 0.01539608, 0.01436491, 0.01854724]
    graph_stat = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_stat))
    ex = (np.array(ex)/2.).tolist()
    graph_syst = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_syst))
    return graph_stat, graph_syst

def get_graphs_3():
    x = [0.507467, 1.014934, 1.304915, 1.594896, 1.884877, 2.319849, 3.117297]
    y = [0.4676848, 0.5024395, 0.5402505, 0.5726697, 0.5958571, 0.6476449, 0.772239]
    ex = [__convert_ntracks_to_dn_deta_error_backward(ntracks) for ntracks in x]
    x = [__convert_ntracks_to_dn_deta_backward(ntracks) for ntracks in x]
    ey_stat = [0.009040637, 0.009738218, 0.009419727, 0.01282185, 0.01527444, 0.01191757, 0.0231139]
    ey_syst = [0.02035057, 0.02136169, 0.02013609, 0.01970749, 0.01754084, 0.01803121, 0.02203735]
    graph_stat = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_stat))
    ex = (np.array(ex)/2.).tolist()
    graph_syst = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_syst))
    return graph_stat, graph_syst

def get_graphs_4():
    x = [0.5809129, 1.161826, 1.493776, 1.825726, 2.157676, 2.821577]
    y = [0.4629038, 0.47614, 0.4972322, 0.5333656, 0.5037638, 0.6000951]
    ex = [__convert_ntracks_to_dn_deta_error_forward(ntracks) for ntracks in x]
    x = [__convert_ntracks_to_dn_deta_forward(ntracks) for ntracks in x]
    ey_stat = [0.009689121, 0.01365575, 0.01046424, 0.01317811, 0.01137415, 0.01448644]
    ey_syst = [0.01421114, 0.01479545, 0.01277154, 0.01485013, 0.01381297, 0.01588298]
    graph_stat = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_stat))
    ex = (np.array(ex)/2.).tolist()
    graph_syst = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_syst))
    return graph_stat, graph_syst

def get_graphs_5():
    x = [0.507467, 1.014934, 1.304915, 1.594896, 1.884877, 2.319849, 3.117297]
    y = [0.5026822, 0.4896216, 0.4662211, 0.5468596, 0.5657158, 0.6337634, 0.7081699]
    ex = [__convert_ntracks_to_dn_deta_error_backward(ntracks) for ntracks in x]
    x = [__convert_ntracks_to_dn_deta_backward(ntracks) for ntracks in x]
    ey_stat = [0.01711132, 0.01710068, 0.0153787, 0.0173467, 0.01988859, 0.01986312, 0.03870117]
    ey_syst = [0.03349001, 0.03404179, 0.02638654, 0.02720649, 0.02607376, 0.02191929, 0.0477596]
    graph_stat = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_stat))
    ex = (np.array(ex)/2.).tolist()
    graph_syst = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_syst))
    return graph_stat, graph_syst

def get_graphs_6():
    x = [0.5809129, 1.161826, 1.493776, 1.825726, 2.157676, 2.821577]
    y = [0.4676977, 0.4771303, 0.494122, 0.5045162, 0.5301832, 0.5727374]
    ex = [__convert_ntracks_to_dn_deta_error_forward(ntracks) for ntracks in x]
    x = [__convert_ntracks_to_dn_deta_forward(ntracks) for ntracks in x]
    ey_stat = [0.01667869, 0.0147616, 0.01438878, 0.0178378, 0.02040139, 0.0321286]
    ey_syst = [0.01827252, 0.02060807, 0.0175403, 0.0176477, 0.0180787, 0.02639767]
    graph_stat = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_stat))
    ex = (np.array(ex)/2.).tolist()
    graph_syst = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_syst))
    return graph_stat, graph_syst

def get_graphs_7():
    x = [0.507467, 1.014934, 1.304915, 1.594896, 1.884877, 2.319849, 3.117297]
    y = [0.5339488, 0.4969878, 0.5277335, 0.6015078, 0.4943096, 0.5615348, 0.684751]
    ex = [__convert_ntracks_to_dn_deta_error_backward(ntracks) for ntracks in x]
    x = [__convert_ntracks_to_dn_deta_backward(ntracks) for ntracks in x]
    ey_stat = [0.02995252, 0.0246888, 0.02913669, 0.04254611, 0.03037722, 0.03127429, 0.06118703]
    ey_syst = [0.05469621, 0.04566076, 0.04435726, 0.0460332, 0.03392438, 0.03103332, 0.03932987]
    graph_stat = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_stat))
    ex = (np.array(ex)/2.).tolist()
    graph_syst = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_syst))
    return graph_stat, graph_syst

def get_graphs_8():
    x = [0.5809129, 1.161826, 1.493776, 1.825726, 2.157676, 2.821577]
    y = [0.4082578, 0.5174834, 0.4760793, 0.5002983, 0.551419, 0.5872638]
    ex = [__convert_ntracks_to_dn_deta_error_forward(ntracks) for ntracks in x]
    x = [__convert_ntracks_to_dn_deta_forward(ntracks) for ntracks in x]
    ey_stat = [0.02116116, 0.0232252, 0.02675465, 0.02952689, 0.02915208, 0.03201709]
    ey_syst = [0.02249501, 0.02602947, 0.02099161, 0.02439124, 0.02518392, 0.02266324]
    graph_stat = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_stat))
    ex = (np.array(ex)/2.).tolist()
    graph_syst = ROOT.TGraphErrors(len(x), np.array(x), np.array(y), np.array(ex), np.array(ey_syst))
    return graph_stat, graph_syst

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

def draw_halo(graph, s=1., m=ROOT.kFullCircle, c=ROOT.kBlack):
    graph.SetMarkerStyle(m)
    graph.SetMarkerSize(s)
    graph.SetMarkerColor(c)
    graph.SetLineWidth(1)
    graph.Draw("pz, same")

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

    pp_unc_text = ROOT.TLatex(0.6, 0.15, '')
    pp_unc_text.SetNDC()
    pp_unc_text.SetTextFont(43)
    pp_unc_text.SetTextSize(40)

    br_unc_text = ROOT.TLatex(0.6, 0.15, '')
    br_unc_text.SetNDC()
    br_unc_text.SetTextFont(43)
    br_unc_text.SetTextSize(40)

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
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    h_frame.GetYaxis().ChangeLabel(1, 1, 0)

    alice_pp = draw_alice_pp(1, 2, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    alice_pbpb = draw_alice_pbpb(1, 2, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)

    x, y = to_pad_coordinates(0.05, 0.9)
    alice_text = ROOT.TLatex(x, y, 'ALICE Preliminary')
    alice_text.SetNDC()
    alice_text.SetTextFont(43)
    alice_text.SetTextSize(50)
    alice_text.Draw()

    x, y = to_pad_coordinates(0.05, 0.85)
    pt_text.DrawLatexNDC(x, y, '1#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{2} GeV/#it{c}')

    graph_stat_lhcb_2_4_b, graph_syst_lhcb_2_4_b = get_graphs_1()

    header_lhcb_b = ROOT.TLatex()
    header_lhcb_b.SetNDC()
    header_lhcb_b.SetTextFont(43)
    header_lhcb_b.SetTextSize(35)
    x, y = to_pad_coordinates(0.06, 0.20)
    header_lhcb_b.DrawLatex(x, y, 'LHCb, #sqrt{#it{s}_{NN}} = 8.16 TeV')

    x_min, y_min = to_pad_coordinates(0.05, 0.07)
    x_max, y_max = to_pad_coordinates(0.5, 0.17)

    legend_lhcb_b = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_lhcb_b.SetBorderSize(0)
    legend_lhcb_b.SetFillStyle(0)
    legend_lhcb_b.SetTextFont(43)
    legend_lhcb_b.SetTextSize(35)
    legend_lhcb_b.AddEntry(graph_stat_lhcb_2_4_b, '#splitline{Pb#font[122]{-}p, #kern[-0.3]{#font[122]{-}4.3} < #kern[-0.3]{#it{y}} < #kern[-0.3]{#font[122]{-}2.8}}{Mult. estim.: #kern[-0.3]{#font[122]{-}4.8} < #kern[-0.3]{#it{#eta}} < #kern[-0.3]{#font[122]{-}2.0}}', 'pl')
    legend_lhcb_b.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(2)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    graph_stat_lhcb_2_4_f, graph_syst_lhcb_2_4_f = get_graphs_2()
    draw_graphs(graph_stat_lhcb_2_4_f, graph_syst_lhcb_2_4_f, c=colors[7], s=4., m=ROOT.kFullDiamond)
    graph_stat_lhcb_2_4_f_halo = graph_stat_lhcb_2_4_f.Clone("graph_stat_lhcb_2_4_f_halo")
    draw_halo(graph_stat_lhcb_2_4_f_halo, s=4.3, m=ROOT.kOpenDiamond, c=ROOT.kRed+2)
    draw_graphs(graph_stat_lhcb_2_4_b, graph_syst_lhcb_2_4_b, c=colors[6], s=4., m=ROOT.kOpenDiamond)
    draw_alice_pp(2, 4, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    draw_alice_pbpb(2, 4, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)
    #alice_text.Draw()
    #y_text.Draw()
    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '2#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{4} GeV/#it{c}')

    header_lhcb_f = ROOT.TLatex()
    header_lhcb_f.SetNDC()
    header_lhcb_f.SetTextFont(43)
    header_lhcb_f.SetTextSize(35)
    x, y = to_pad_coordinates(0.06, 0.20)
    header_lhcb_f.DrawLatex(x, y, 'LHCb, #sqrt{#it{s}_{NN}} = 8.16 TeV')

    x_min, y_min = to_pad_coordinates(0.05, 0.07)
    x_max, y_max = to_pad_coordinates(0.5, 0.17)
    legend_lhcb_f = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_lhcb_f.SetBorderSize(0)
    legend_lhcb_f.SetFillStyle(0)
    legend_lhcb_f.SetTextFont(43)
    legend_lhcb_f.SetTextSize(35)
    legend_lhcb_f.AddEntry(graph_stat_lhcb_2_4_f, '#splitline{p#font[122]{-}Pb, 1.8 < #kern[-0.3]{#it{y}} < 3.3}{Mult. estim.: 2.0 < #kern[-0.4]{#it{#eta}} < 4.8}', 'pl')
    legend_lhcb_f.Draw()

    legend_lhcb_f_halo = ROOT.TLegend(x_min, y_min, x_max, y_max)
    legend_lhcb_f_halo.SetBorderSize(0)
    legend_lhcb_f_halo.SetFillStyle(0)
    legend_lhcb_f_halo.SetTextFont(43)
    legend_lhcb_f_halo.SetTextSize(35)
    legend_lhcb_f_halo.AddEntry(graph_stat_lhcb_2_4_f_halo, '#splitline{p#font[122]{-}Pb, 1.8 < #kern[-0.3]{#it{y}} < 3.3}{Mult. estim.: 2.0 < #kern[-0.4]{#it{#eta}} < 4.8}', 'pl')
    legend_lhcb_f_halo.Draw()

    ROOT.gPad.RedrawAxis()

    c.cd(3)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    graph_stat_lhcb_4_6_f, graph_syst_lhcb_4_6_f = get_graphs_4()
    draw_graphs(graph_stat_lhcb_4_6_f, graph_syst_lhcb_4_6_f, c=colors[7], s=4., m=ROOT.kFullDiamond)
    graph_stat_lhcb_4_6_f_halo = graph_stat_lhcb_4_6_f.Clone("graph_stat_lhcb_4_6_f_halo")
    draw_halo(graph_stat_lhcb_4_6_f_halo, s=4.3, m=ROOT.kOpenDiamond, c=ROOT.kRed+2)
    graph_stat_lhcb_4_6_b, graph_syst_lhcb_4_6_b = get_graphs_3()
    draw_graphs(graph_stat_lhcb_4_6_b, graph_syst_lhcb_4_6_b, c=colors[6], s=4., m=ROOT.kOpenDiamond)
    draw_alice_pp(4, 6, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    draw_alice_pbpb(4, 6, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)

    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '4#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[1.]{6} GeV/#it{c}')

    x, y = to_pad_coordinates(0.05, 0.1)
    br_unc_text.DrawLatexNDC(x, y, '#lower[-0.03]{^{+4.0}}')
    br_unc_text.DrawLatexNDC(x, y, '_{#minus3.8}')
    br_unc_text.DrawLatexNDC(x, y, '#kern[.12]{% BR uncertainty not shown}')

    ROOT.gPad.RedrawAxis()

    c.cd(4)
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    h_frame.GetYaxis().ChangeLabel(10, 1, 0)
    graph_stat_lhcb_6_8_f, graph_syst_lhcb_6_8_f = get_graphs_6()
    draw_graphs(graph_stat_lhcb_6_8_f, graph_syst_lhcb_6_8_f, c=colors[7], s=4., m=ROOT.kFullDiamond)
    graph_stat_lhcb_6_8_f_halo = graph_stat_lhcb_6_8_f.Clone("graph_stat_lhcb_6_8_f_halo")
    draw_halo(graph_stat_lhcb_6_8_f_halo, s=4.3, m=ROOT.kOpenDiamond, c=ROOT.kRed+2)
    graph_stat_lhcb_6_8_b, graph_syst_lhcb_6_8_b = get_graphs_5()
    draw_graphs(graph_stat_lhcb_6_8_b, graph_syst_lhcb_6_8_b, c=colors[6], s=4., m=ROOT.kOpenDiamond)
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
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    h_frame.GetXaxis().ChangeLabel(1, 1, 0)
    graph_stat_lhcb_8_12_f, graph_syst_lhcb_8_12_f = get_graphs_8()
    draw_graphs(graph_stat_lhcb_8_12_f, graph_syst_lhcb_8_12_f, c=colors[7], s=4., m=ROOT.kFullDiamond)
    graph_stat_lhcb_8_12_f_halo = graph_stat_lhcb_8_12_f.Clone("graph_stat_lhcb_8_12_f_halo")
    draw_halo(graph_stat_lhcb_8_12_f_halo, s=4.3, m=ROOT.kOpenDiamond, c=ROOT.kRed+2)
    graph_stat_lhcb_8_12_b, graph_syst_lhcb_8_12_b = get_graphs_7()
    draw_graphs(graph_stat_lhcb_8_12_b, graph_syst_lhcb_8_12_b, c=colors[6], s=4., m=ROOT.kOpenDiamond)
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
    h_frame = ROOT.gPad.DrawFrame(1., 0., 5000., 0.9, ";#LTd#it{N}_{ch}/d#it{#eta}#GT ;#it{#sigma}_{D_{s}^{+}}/#it{#sigma}_{D^{+}}")
    h_frame.GetXaxis().ChangeLabel(1, 1, 0)
    draw_alice_pp(12, 24, c=ROOT.kBlack, s=2.5, m=ROOT.kFullCircle)
    draw_alice_pbpb(12, 24, c=ROOT.kBlack, s=2.5, m=ROOT.kOpenCircle)
    #alice_text.Draw()
    #y_text.Draw()
    x, y = to_pad_coordinates(0.05, 0.9)
    pt_text.DrawLatexNDC(x, y, '12#kern[1.]{<}#kern[.5]{#it{p}_{T}}#kern[.5]{<}#kern[.5]{24} GeV/#it{c}')

    x, y = to_pad_coordinates(0.05, 0.1)
    pp_unc_text.DrawLatexNDC(x, y, 'pp uncertainty on #kern[-0.05]{#it{x}-axis scaled by 5}')

    ROOT.gPad.RedrawAxis()

    c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_lhcb.pdf")
    c.SaveAs("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Figures/Ratio/VsMult/ratio_lhcb.root")
