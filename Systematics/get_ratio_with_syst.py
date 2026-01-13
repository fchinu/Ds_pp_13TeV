import argparse
import numpy as np
import yaml
import ROOT

def get_unc_dn_deta_pbpb(centmin, centmax):
    """"
    mult_unc_low = [0.05, 0.07, np.sqrt(0.11**2 + 0.09**2)/2, np.sqrt((0.21/4)**2 + (0.15/4)**2 + (0.13/2)**2), np.sqrt(0.21**2 + 0.18**2)/2, 0.36]
    mult_unc_high = [0.07, 0.10, np.sqrt(0.16**2 + 0.13**2)/2, np.sqrt((0.18/4)**2 + (0.19/4)**2 + (0.19/2)**2), np.sqrt(0.27**2 + 0.24**2)/2, 0.42]
    for a in zip(mult_unc_low[::-1], mult_unc_high[::-1]):
         print(a)
    """
    mapping = {
        (0, 10): (34.2, 34.2),
        (10, 20): (33, 33),
        (20, 30): (25, 25),
        (30, 40): (19, 19),
        (40, 50): (14, 14),
        (50, 60): (11, 11),
        (60, 70): (8, 8),
        (70, 80): (5, 5),
        (80, 90): (2.8, 2.8)}
    return mapping[(centmin, centmax)]

def get_unc_dn_deta_pp(centmin, centmax):
    """"
    mult_unc_low = [0.05, 0.07, np.sqrt(0.11**2 + 0.09**2)/2, np.sqrt((0.21/4)**2 + (0.15/4)**2 + (0.13/2)**2), np.sqrt(0.21**2 + 0.18**2)/2, 0.36]
    mult_unc_high = [0.07, 0.10, np.sqrt(0.16**2 + 0.13**2)/2, np.sqrt((0.18/4)**2 + (0.19/4)**2 + (0.19/2)**2), np.sqrt(0.27**2 + 0.24**2)/2, 0.42]
    for a in zip(mult_unc_low[::-1], mult_unc_high[::-1]):
         print(a)
    """
    mapping = {
        (0, 1): (0.36, 0.42),
        (1, 10): (0.1382931668593933, 0.18062391868188443),
        (10, 30): (0.091583295420071, 0.11535271995059328),
        (30, 50): (0.07106335201775947, 0.10307764064044152),
        (50, 70): (0.07, 0.1),
        (70, 100): (0.05, 0.07),
        (0, 100): (0.09, 0.13)}
    return mapping[(centmin, centmax)]

mult_unc_low = [0.05, 0.07, np.sqrt(0.11**2 + 0.09**2)/2, np.sqrt((0.21/4)**2 + (0.15/4)**2 + (0.13/2)**2), np.sqrt(0.21**2 + 0.18**2)/2, 0.36]
mult_unc_high = [0.07, 0.10, np.sqrt(0.16**2 + 0.13**2)/2, np.sqrt((0.18/4)**2 + (0.19/4)**2 + (0.19/2)**2), np.sqrt(0.27**2 + 0.24**2)/2, 0.42]

def get_cross_sec_with_syst(config_file_name):
    with open(config_file_name, 'r') as f:
        config = yaml.safe_load(f)
    with open(config['inputs']['cutset'], 'r') as f:
        cutset = yaml.safe_load(f)
    
    pt_mins = cutset['pt']['min']
    pt_maxs = cutset['pt']['max']

    cent_mins = cutset['cent']['min']
    cent_maxs = cutset['cent']['max']

    f = ROOT.TFile.Open(config['output_name'], 'recreate')
    f.Close()

    # remove 0-100 bin
    cent_tuples = [(cent_mins[i], cent_maxs[i]) for i in range(len(cent_mins))]
    if (0, 100) in cent_tuples:
        idx = cent_tuples.index((0, 100))
        cent_mins.pop(idx)
        cent_maxs.pop(idx)

    uncs_low = {}
    uncs_high = {}
    uncs_no_br_low = {}
    uncs_no_br_high = {}
    for pt_min, pt_max in zip(pt_mins, pt_maxs):
        uncs_low[(pt_min, pt_max)] = []
        uncs_high[(pt_min, pt_max)] = []
        uncs_no_br_low[(pt_min, pt_max)] = []
        uncs_no_br_high[(pt_min, pt_max)] = []

        with ROOT.TFile.Open(config['inputs']['ratio_file']) as f:
            g_ratio = f.Get(f'g_ratio_dndeta_{pt_min*10:.0f}_{pt_max*10:.0f}')

        # set dN/deta uncertainties
        for i_cent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
            unc_dn_deta_low, unc_dn_deta_high = get_unc_dn_deta_pp(cent_min, cent_max) if config['is_pp'] else get_unc_dn_deta_pbpb(cent_min, cent_max)
            g_ratio.SetPointEXlow(i_cent, unc_dn_deta_low)
            g_ratio.SetPointEXhigh(i_cent, unc_dn_deta_high)

        # Get systematics graphs (pt dependent)
        g_systs_rel_lower = {}
        g_systs_rel_upper = {}
        g_systs_lower = {}
        g_systs_upper = {}
        for syst_name, file_name in config['inputs']['syst_files'].items():

            with ROOT.TFile.Open(file_name) as f:
                try:
                    try:
                        h_syst = f.Get('assigned_syst')
                        g_systs_rel_lower[syst_name] = g_ratio.Clone()
                        g_systs_rel_upper[syst_name] = g_ratio.Clone()
                        for i_cent in range(g_systs_rel_lower[syst_name].GetN()):
                            g_systs_rel_lower[syst_name].SetPointY(
                                i_cent,
                                h_syst.GetBinContent(h_syst.FindBin((pt_min+pt_max)/2))
                            )
                            g_systs_rel_upper[syst_name].SetPointY(
                                i_cent,
                                h_syst.GetBinContent(h_syst.FindBin((pt_min+pt_max)/2))
                            )

                    except:
                        g_systs_rel_lower[syst_name] = g_ratio.Clone()
                        g_systs_rel_upper[syst_name] = g_ratio.Clone()
                        for i_cent in range(g_systs_rel_lower[syst_name].GetN()):
                            syst_vs_pt = f.Get(f'assigned_syst_{cent_mins[i_cent]}_{cent_maxs[i_cent]}')
                            g_systs_rel_lower[syst_name].SetPointY(
                                i_cent,
                                syst_vs_pt.GetBinContent(syst_vs_pt.FindBin((pt_min+pt_max)/2))
                            )
                            g_systs_rel_upper[syst_name].SetPointY(
                                i_cent,
                                syst_vs_pt.GetBinContent(syst_vs_pt.FindBin((pt_min+pt_max)/2))
                            )
                    g_systs_lower[syst_name] = g_systs_rel_lower[syst_name].Clone()
                    g_systs_upper[syst_name] = g_systs_rel_upper[syst_name].Clone()
                    for i_cent in range(g_systs_lower[syst_name].GetN()):
                        g_systs_lower[syst_name].SetPointY(
                            i_cent,
                            g_systs_rel_lower[syst_name].GetPointY(i_cent)*g_ratio.GetPointY(i_cent)
                        )
                        g_systs_upper[syst_name].SetPointY(
                            i_cent,
                            g_systs_rel_upper[syst_name].GetPointY(i_cent)*g_ratio.GetPointY(i_cent)
                        )
                except:
                    g_systs_rel_lower[syst_name] = g_ratio.Clone()
                    g_systs_rel_upper[syst_name] = g_ratio.Clone()
                    for i_cent in range(g_systs_rel_lower[syst_name].GetN()):
                        syst_vs_pt = f.Get(f'assigned_syst_lower_{cent_mins[i_cent]}_{cent_maxs[i_cent]}')
                        g_systs_rel_lower[syst_name].SetPointY(
                            i_cent,
                            syst_vs_pt.GetBinContent(syst_vs_pt.FindBin((pt_min+pt_max)/2))
                        )
                        syst_vs_pt = f.Get(f'assigned_syst_upper_{cent_mins[i_cent]}_{cent_maxs[i_cent]}')
                        g_systs_rel_upper[syst_name].SetPointY(
                            i_cent,
                            syst_vs_pt.GetBinContent(syst_vs_pt.FindBin((pt_min+pt_max)/2))
                        )
                    g_systs_lower[syst_name] = g_systs_rel_lower[syst_name].Clone()
                    g_systs_upper[syst_name] = g_systs_rel_upper[syst_name].Clone()
                    for i_cent in range(g_systs_lower[syst_name].GetN()):
                        g_systs_lower[syst_name].SetPointY(
                            i_cent,
                            g_systs_rel_lower[syst_name].GetPointY(i_cent)*g_ratio.GetPointY(i_cent)
                        )
                        g_systs_upper[syst_name].SetPointY(
                            i_cent,
                            g_systs_rel_upper[syst_name].GetPointY(i_cent)*g_ratio.GetPointY(i_cent)
                        )

        # Get systematics graphs (cent independent)
        unc_br_ds_low_rel = config["br"][0]["unc_low"] / config["br"][0]["value"]
        unc_br_ds_high_rel = config["br"][0]["unc_high"] / config["br"][0]["value"]
        unc_br_dplus_low_rel = config["br"][1]["unc_low"] / config["br"][1]["value"]
        unc_br_dplus_high_rel = config["br"][1]["unc_high"] / config["br"][1]["value"]

        unc_br_ratio_low_rel = np.sqrt(
            unc_br_ds_low_rel**2 + unc_br_dplus_high_rel**2
        )
        unc_br_ratio_high_rel = np.sqrt(
            unc_br_ds_high_rel**2 + unc_br_dplus_low_rel**2
        )

        g_systs_rel_lower['br'] = g_ratio.Clone()
        for i_cent in range(g_systs_rel_lower['br'].GetN()):
            g_systs_rel_lower['br'].SetPointY(i_cent, unc_br_ratio_low_rel)

        g_systs_rel_upper['br'] = g_ratio.Clone()
        for i_cent in range(g_systs_rel_upper['br'].GetN()):
            g_systs_rel_upper['br'].SetPointY(i_cent, unc_br_ratio_high_rel)

        with ROOT.TFile.Open(config['output_name'], 'update') as f:
            for syst_name, g_syst in g_systs_rel_lower.items():
                g_syst.Write(f'{syst_name}_lower_{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            for syst_name, g_syst in g_systs_rel_upper.items():
                g_syst.Write(f'{syst_name}_upper_{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            for syst_name, g_syst in g_systs_lower.items():
                g_syst.Write(f'{syst_name}_lower_{pt_min*10:.0f}_{pt_max*10:.0f}')
            for syst_name, g_syst in g_systs_upper.items():
                g_syst.Write(f'{syst_name}_upper_{pt_min*10:.0f}_{pt_max*10:.0f}')

        g_total_syst_lower = g_ratio.Clone("g_total_syst_lower")
        g_total_syst_upper = g_ratio.Clone("g_total_syst_upper")
        g_total_syst_rel_lower = g_ratio.Clone("g_total_syst_rel_lower")
        g_total_syst_rel_upper = g_ratio.Clone("g_total_syst_rel_upper")

        g_total_syst_no_br_lower = g_ratio.Clone("g_total_syst_no_br_lower")
        g_total_syst_no_br_upper = g_ratio.Clone("g_total_syst_no_br_upper")
        g_total_syst_no_br_rel_lower = g_ratio.Clone("g_total_syst_no_br_rel_lower")
        g_total_syst_no_br_rel_upper = g_ratio.Clone("g_total_syst_no_br_rel_upper")

        # Evaluate total systematics
        for i in range(g_ratio.GetN()):
            unc_lower = 0
            unc_upper = 0
            for syst_name, g_syst in g_systs_rel_lower.items():
                if "br" not in syst_name:
                    unc_lower += g_syst.GetPointY(i)**2
            for syst_name, g_syst in g_systs_rel_upper.items():
                if "br" not in syst_name:
                    unc_upper += g_syst.GetPointY(i)**2

            g_total_syst_no_br_rel_lower.SetPointY(i, np.sqrt(unc_lower))
            g_total_syst_no_br_rel_upper.SetPointY(i, np.sqrt(unc_upper))
            g_total_syst_no_br_lower.SetPointY(i, g_ratio.GetPointY(i)*np.sqrt(unc_lower))
            g_total_syst_no_br_upper.SetPointY(i, g_ratio.GetPointY(i)*np.sqrt(unc_upper))

            uncs_no_br_low[(pt_min, pt_max)].append(g_ratio.GetPointY(i)*np.sqrt(unc_lower))
            uncs_no_br_high[(pt_min, pt_max)].append(g_ratio.GetPointY(i)*np.sqrt(unc_upper))

            unc_lower = np.sqrt(unc_lower + g_systs_rel_lower['br'].GetPointY(i)**2)
            unc_upper = np.sqrt(unc_upper + g_systs_rel_upper['br'].GetPointY(i)**2)

            g_total_syst_rel_lower.SetPointY(i, unc_lower)
            g_total_syst_rel_upper.SetPointY(i, unc_upper)
            g_total_syst_lower.SetPointY(i, g_ratio.GetPointY(i)*unc_lower)
            g_total_syst_upper.SetPointY(i, g_ratio.GetPointY(i)*unc_upper)

            uncs_low[(pt_min, pt_max)].append(g_ratio.GetPointY(i)*unc_lower)
            uncs_high[(pt_min, pt_max)].append(g_ratio.GetPointY(i)*unc_upper)

        g_syst = g_ratio.Clone("g_syst")
        for i in range(g_ratio.GetN()):
            g_syst.SetPointEYlow(i, g_total_syst_lower.GetPointY(i))
            g_syst.SetPointEYhigh(i, g_total_syst_upper.GetPointY(i))

        g_syst_no_br = g_ratio.Clone("g_syst_no_br")
        for i in range(g_ratio.GetN()):
            g_syst_no_br.SetPointEYlow(i, g_total_syst_no_br_lower.GetPointY(i))
            g_syst_no_br.SetPointEYhigh(i, g_total_syst_no_br_upper.GetPointY(i))
        
        with ROOT.TFile.Open(config['output_name'], 'update') as f:
            g_total_syst_upper.Write(f'g_total_syst_upper_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_total_syst_lower.Write(f'g_total_syst_lower_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_total_syst_no_br_upper.Write(f'g_total_syst_no_br_upper_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_total_syst_no_br_lower.Write(f'g_total_syst_no_br_lower_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_total_syst_rel_upper.Write(f'g_total_syst_upper_{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            g_total_syst_rel_lower.Write(f'g_total_syst_lower_{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            g_total_syst_no_br_rel_lower.Write(f'g_total_syst_no_br_lower{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            g_total_syst_no_br_rel_upper.Write(f'g_total_syst_no_br_upper_{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            g_syst.Write(f'g_syst_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_syst_no_br.Write(f'g_syst_no_br_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_ratio.Write(f'g_stat_{pt_min*10:.0f}_{pt_max*10:.0f}')

    for i_cent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
        with ROOT.TFile.Open(config['inputs']['ratio_file']) as f:
            h_stat_vs_pt = f.Get(f'centrality_{cent_min}_{cent_max}/h_ratio')
            g_syst_vs_pt = ROOT.TGraphAsymmErrors(h_stat_vs_pt)
            g_syst_vs_pt_no_br = ROOT.TGraphAsymmErrors(h_stat_vs_pt)
            for i in range(g_syst_vs_pt.GetN()):
                g_syst_vs_pt.SetPointEYlow(i, uncs_low[(pt_min, pt_max)][i_cent])
                g_syst_vs_pt.SetPointEYhigh(i, uncs_high[(pt_min, pt_max)][i_cent])
                g_syst_vs_pt_no_br.SetPointEYlow(i, uncs_no_br_low[(pt_min, pt_max)][i_cent])
                g_syst_vs_pt_no_br.SetPointEYhigh(i, uncs_no_br_high[(pt_min, pt_max)][i_cent])
        with ROOT.TFile.Open(config['output_name'], 'update') as f:
            h_stat_vs_pt.Write(f'h_stat_vs_pt_{cent_min}_{cent_max}')
            g_syst_vs_pt.Write(f'g_syst_vs_pt_{cent_min}_{cent_max}')
            g_syst_vs_pt_no_br.Write(f'g_syst_vs_pt_no_br_{cent_min}_{cent_max}')

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Get cross section with systematics')
    parser.add_argument('config_file', help='Configuration file')
    args = parser.parse_args()

    ROOT.TH1.AddDirectory(False)

    get_cross_sec_with_syst(args.config_file)
