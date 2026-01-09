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
    for pt_min, pt_max in zip(pt_mins, pt_maxs):
        uncs_low[(pt_min, pt_max)] = []
        uncs_high[(pt_min, pt_max)] = []

        with ROOT.TFile.Open(config['inputs']['ratio_file']) as f:
            g_ratio = f.Get(f'g_ratio_dndeta_{pt_min*10:.0f}_{pt_max*10:.0f}')

        # set dN/deta uncertainties
        for i_cent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
            unc_dn_deta_low, unc_dn_deta_high = get_unc_dn_deta_pp(cent_min, cent_max) if config['is_pp'] else get_unc_dn_deta_pbpb(cent_min, cent_max)
            g_ratio.SetPointEXlow(i_cent, unc_dn_deta_low)
            g_ratio.SetPointEXhigh(i_cent, unc_dn_deta_high)

        # Get systematics graphs (pt dependent)
        g_systs_rel = {}
        g_systs = {}
        for syst_name, file_name in config['inputs']['syst_files'].items():
            print(syst_name)
            with ROOT.TFile.Open(file_name) as f:
                try:
                    h_syst = f.Get('assigned_syst')
                    g_systs_rel[syst_name] = g_ratio.Clone()
                    for i_cent in range(g_systs_rel[syst_name].GetN()):
                        g_systs_rel[syst_name].SetPointY(
                            i_cent,
                            h_syst.GetBinContent(h_syst.FindBin((pt_min+pt_max)/2))
                        )

                except:
                    g_systs_rel[syst_name] = g_ratio.Clone()
                    current_g_syst = g_systs_rel[syst_name]
                    for i_cent in range(current_g_syst.GetN()):
                        print(f"assigned_syst_{cent_mins[i_cent]}_{cent_maxs[i_cent]}")
                        syst_vs_pt = f.Get(f'assigned_syst_{cent_mins[i_cent]}_{cent_maxs[i_cent]}')
                        current_g_syst.SetPointY(
                            i_cent,
                            syst_vs_pt.GetBinContent(syst_vs_pt.FindBin((pt_min+pt_max)/2))
                        )
                g_systs[syst_name] = g_systs_rel[syst_name].Clone()
                for i_cent in range(g_systs[syst_name].GetN()):
                    g_systs[syst_name].SetPointY(
                        i_cent,
                        g_systs_rel[syst_name].GetPointY(i_cent)*g_ratio.GetPointY(i_cent)
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

        g_systs_rel['br_low'] = g_ratio.Clone()
        for i_cent in range(g_systs_rel['br_low'].GetN()):
            g_systs_rel['br_low'].SetPointY(i_cent, unc_br_ratio_low_rel)

        g_systs_rel['br_high'] = g_ratio.Clone()
        for i_cent in range(g_systs_rel['br_high'].GetN()):
            g_systs_rel['br_high'].SetPointY(i_cent, unc_br_ratio_high_rel)

        with ROOT.TFile.Open(config['output_name'], 'update') as f:
            for syst_name, g_syst in g_systs_rel.items():
                g_syst.Write(f'{syst_name}_{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            for syst_name, g_syst in g_systs.items():
                g_syst.Write(f'{syst_name}_{pt_min*10:.0f}_{pt_max*10:.0f}')
        

        g_total_syst_low = g_ratio.Clone("g_total_syst_low")
        g_total_syst_high = g_ratio.Clone("g_total_syst_high")
        g_total_syst_low_rel = g_ratio.Clone("g_total_syst_low_rel")
        g_total_syst_high_rel = g_ratio.Clone("g_total_syst_high_rel")

        g_total_syst_no_br = g_ratio.Clone("g_total_syst_no_br")
        g_total_syst_no_br_rel = g_ratio.Clone("g_total_syst_no_br_rel")
        # Evaluate total systematics
        for i in range(g_ratio.GetN()):
            unc = 0
            for syst_name, g_syst in g_systs_rel.items():
                if "br" not in syst_name:
                    unc += g_syst.GetPointY(i)**2

            unc_low = np.sqrt(unc + g_systs_rel['br_low'].GetPointY(i)**2)
            unc_high = np.sqrt(unc + g_systs_rel['br_high'].GetPointY(i)**2)

            g_total_syst_low_rel.SetPointY(i, unc_low)
            g_total_syst_high_rel.SetPointY(i, unc_high)
            g_total_syst_no_br_rel.SetPointY(i, np.sqrt(unc))

            g_total_syst_low.SetPointY(i, g_ratio.GetPointY(i)*unc_low)
            g_total_syst_high.SetPointY(i, g_ratio.GetPointY(i)*unc_high)
            g_total_syst_no_br.SetPointY(i, g_ratio.GetPointY(i)*np.sqrt(unc))
            uncs_low[(pt_min, pt_max)].append(g_ratio.GetPointY(i)*unc_low)
            uncs_high[(pt_min, pt_max)].append(g_ratio.GetPointY(i)*unc_high)

        g_syst = g_ratio.Clone("g_syst")
        for i in range(g_ratio.GetN()):
            g_syst.SetPointEYlow(i, g_total_syst_low.GetPointY(i))
            g_syst.SetPointEYhigh(i, g_total_syst_high.GetPointY(i))
        
        with ROOT.TFile.Open(config['output_name'], 'update') as f:
            g_total_syst_high.Write(f'g_total_syst_high_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_total_syst_low.Write(f'g_total_syst_low_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_total_syst_no_br.Write(f'g_total_syst_no_br_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_total_syst_high_rel.Write(f'g_total_syst_high_{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            g_total_syst_low_rel.Write(f'g_total_syst_low_{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            g_total_syst_no_br_rel.Write(f'g_total_syst_no_br_{pt_min*10:.0f}_{pt_max*10:.0f}_rel')
            g_syst.Write(f'g_syst_{pt_min*10:.0f}_{pt_max*10:.0f}')
            g_ratio.Write(f'g_stat_{pt_min*10:.0f}_{pt_max*10:.0f}')

    for i_cent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
        with ROOT.TFile.Open(config['inputs']['ratio_file']) as f:
            h_stat_vs_pt = f.Get(f'centrality_{cent_min}_{cent_max}/h_ratio')
            g_syst_vs_pt = ROOT.TGraphAsymmErrors(h_stat_vs_pt)
            for i in range(g_syst_vs_pt.GetN()):
                g_syst_vs_pt.SetPointEYlow(i, uncs_low[(pt_min, pt_max)][i_cent])
                g_syst_vs_pt.SetPointEYhigh(i, uncs_high[(pt_min, pt_max)][i_cent])
        with ROOT.TFile.Open(config['output_name'], 'update') as f:
            g_syst_vs_pt.Write(f'g_syst_vs_pt_{cent_min}_{cent_max}')
            h_stat_vs_pt.Write(f'h_stat_vs_pt_{cent_min}_{cent_max}')

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Get cross section with systematics')
    parser.add_argument('config_file', help='Configuration file')
    args = parser.parse_args()

    ROOT.TH1.AddDirectory(False)

    get_cross_sec_with_syst(args.config_file)
