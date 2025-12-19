import ROOT
import numpy as np

ROOT.TH1.AddDirectory(False)

CONFIG_PBPB = {
    "D0": {
        "Prompt": {
            (0, 10): ["Table 1", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2plus", "Hist1D_y1_e2minus"],
            (30, 50): ["Table 2", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2plus", "Hist1D_y1_e2minus"],
        },
        "NonPrompt": {
            (0, 10): ["Table 1a", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2"],
            (30, 50): ["Table 1a", "Hist1D_y2", "Hist1D_y2_e1", "Hist1D_y2_e2"],
        }
    },
    "Ds": {
        "Prompt": {
            (0, 10): ["Table 1", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2plus", "Hist1D_y1_e2minus"],
            (30, 50): ["Table 2", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2plus", "Hist1D_y1_e2minus"],
        },
        "NonPrompt": {
            (0, 10): ["Table 0", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2"],
            (30, 50): ["Table 1", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2"],
        }
    }
}

CONFIG_PP = {
    "Dplus": {
        "Prompt":  ["Table 5", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2plus", "Hist1D_y1_e2minus"],
        "NonPrompt":  ["Table 2", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2"],
    },
    "Ds": {
        "Prompt":  ["Table 6", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2plus", "Hist1D_y1_e2minus"],
        "NonPrompt":  ["Table 3", "Hist1D_y1", "Hist1D_y1_e1", "Hist1D_y1_e2"],
    }
}

def get_stat_syst_combined_from_hepdata_asymm(h_central, h_stat_unc, h_syst_unc_plus, h_syst_unc_minus):
    """
    Combine statistical and systematic uncertainties from HepData histograms

    Parameters:
    - h_central (TH1): central values histogram
    - h_stat_unc (TH1): statistical uncertainties histogram
    - h_syst_unc_plus (TH1): systematic uncertainties (plus) histogram
    - h_syst_unc_minus (TH1): systematic uncertainties (minus) histogram

    Returns:
    - TH1: statistical uncertainties histogram
    - TGraphAsymmErrors: systematic uncertainties graph
    - TH1: combined uncertainties histogram
    """
    h_stat = h_central.Clone("h_stat")
    h_stat.SetDirectory(0)
    for i in range(1, h_stat.GetNbinsX() + 1):
        h_stat.SetBinError(i, h_stat_unc.GetBinContent(i))
        
    g_syst = ROOT.TGraphAsymmErrors(h_central.GetNbinsX())
    g_syst.SetName("g_syst")
    for i in range(1, h_central.GetNbinsX() + 1):
        g_syst.SetPoint(i - 1, h_central.GetBinCenter(i), h_central.GetBinContent(i))
        if h_syst_unc_minus.GetBinContent(i) < 0:
            g_syst.SetPointError(i - 1, h_central.GetBinWidth(i) / 2, h_central.GetBinWidth(i) / 2, -h_syst_unc_minus.GetBinContent(i), h_syst_unc_plus.GetBinContent(i))
        else:
            g_syst.SetPointError(i - 1, h_central.GetBinWidth(i) / 2, h_central.GetBinWidth(i) / 2, h_syst_unc_minus.GetBinContent(i), h_syst_unc_plus.GetBinContent(i))

    h_combined = ROOT.TH1F(h_stat)
    h_combined.SetName("h_combined")
    h_combined.SetDirectory(0)
    for i in range(1, h_stat.GetNbinsX() + 1):
        h_combined.SetBinError(i, np.sqrt(g_syst.GetErrorY(i-1)**2 + h_stat_unc.GetBinContent(i)**2))

    return h_stat, g_syst, h_combined

def get_stat_syst_combined_from_hepdata(h_central, h_stat_unc, h_syst_unc):
    """
    Combine statistical and systematic uncertainties from HepData histograms

    Parameters:
    - h_central (TH1): central values histogram
    - h_stat_unc (TH1): statistical uncertainties histogram
    - h_syst_unc (TH1): systematic uncertainties histogram

    Returns:
    - TH1: statistical uncertainties histogram
    - TGraphAsymmErrors: systematic uncertainties graph
    - TH1: combined uncertainties histogram
    """
    h_stat = h_central.Clone("h_stat")
    h_stat.SetDirectory(0)
    for i in range(1, h_stat.GetNbinsX() + 1):
        h_stat.SetBinError(i, h_stat_unc.GetBinContent(i))

    g_syst = ROOT.TGraphAsymmErrors(h_central.GetNbinsX())
    g_syst.SetName("g_syst")
    for i in range(1, h_central.GetNbinsX() + 1):
        g_syst.SetPoint(i - 1, h_central.GetBinCenter(i), h_central.GetBinContent(i))
        g_syst.SetPointError(i - 1, h_central.GetBinWidth(i) / 2, h_central.GetBinWidth(i) / 2, h_syst_unc.GetBinContent(i), h_syst_unc.GetBinContent(i))

    h_combined = ROOT.TH1F(h_stat)
    h_combined.SetName("h_combined")
    h_combined.SetDirectory(0)
    for i in range(1, h_stat.GetNbinsX() + 1):
        h_combined.SetBinError(i, np.sqrt(g_syst.GetErrorY(i-1)**2 + h_stat_unc.GetBinContent(i)**2))

    return h_stat, g_syst, h_combined


def get_yield(particle, is_pbpb, isprompt, cent_min=None, cent_max=None):
    """
    Get Ds non-prompt yield from ROOT file
    """
    if is_pbpb:
        filename = f"data/{particle}{'Prompt' if isprompt else 'NonPrompt'}Yield_PbPb5TeV_Run2_{cent_min}{cent_max}.root"
        config = CONFIG_PBPB[particle]["Prompt" if isprompt else "NonPrompt"][(cent_min, cent_max)]
    else:
        filename = f"data/{particle}{'Prompt' if isprompt else 'NonPrompt'}Yield_pp5TeV_Run2.root"
        config = CONFIG_PP[particle]["Prompt" if isprompt else "NonPrompt"]

    with ROOT.TFile.Open(filename) as f:
        table_name, hist_central_name, hist_stat_name, *hist_syst_names = config

        dir_table = f.Get(table_name)
        h_central = dir_table.Get(hist_central_name)
        h_stat_unc = dir_table.Get(hist_stat_name)
        if len(hist_syst_names) == 2:
            h_syst_unc_plus = dir_table.Get(hist_syst_names[0])
            h_syst_unc_minus = dir_table.Get(hist_syst_names[1])
            return get_stat_syst_combined_from_hepdata_asymm(h_central, h_stat_unc, h_syst_unc_plus, h_syst_unc_minus)
        else:
            h_syst_unc = dir_table.Get(hist_syst_names[0])
            return get_stat_syst_combined_from_hepdata(h_central, h_stat_unc, h_syst_unc)

if __name__ == "__main__":
    cent_classes = [(0, 10), (30, 50)]
    with ROOT.TFile.Open(f"data/Data_PbPb5TeV_Run2.root", "RECREATE") as outfile:
        ds_dir = outfile.mkdir("Ds")
        d0_dir = outfile.mkdir("D0")
        ds_dir_prompt = ds_dir.mkdir("Prompt")
        d0_dir_prompt = d0_dir.mkdir("Prompt")
        ds_dir_nonprompt = ds_dir.mkdir("NonPrompt")
        d0_dir_nonprompt = d0_dir.mkdir("NonPrompt")
        for cent_min, cent_max in cent_classes:
            h_d0_prompt, g_d0_prompt_syst, h_d0_prompt_combined = get_yield("D0", True, True, cent_min, cent_max)
            h_d0_nonprompt, g_d0_nonprompt_syst, h_d0_nonprompt_combined = get_yield("D0", True, False, cent_min, cent_max)
            h_ds_prompt, g_ds_prompt_syst, h_ds_prompt_combined = get_yield("Ds", True, True, cent_min, cent_max)
            h_ds_nonprompt, g_ds_nonprompt_syst, h_ds_nonprompt_combined = get_yield("Ds", True, False, cent_min, cent_max)

            ds_prompt_dir_cent = ds_dir_prompt.mkdir(f"{cent_min}{cent_max}")
            ds_nonprompt_dir_cent = ds_dir_nonprompt.mkdir(f"{cent_min}{cent_max}")
            d0_prompt_dir_cent = d0_dir_prompt.mkdir(f"{cent_min}{cent_max}")
            d0_nonprompt_dir_cent = d0_dir_nonprompt.mkdir(f"{cent_min}{cent_max}")
            ds_prompt_dir_cent.cd()
            for dsobj in [h_ds_prompt, g_ds_prompt_syst, h_ds_prompt_combined]:
                dsobj.Write()
            d0_prompt_dir_cent.cd()
            for d0obj in [h_d0_prompt, g_d0_prompt_syst, h_d0_prompt_combined]:
                d0obj.Write()
            ds_nonprompt_dir_cent.cd()
            for dsobj in [h_ds_nonprompt, g_ds_nonprompt_syst, h_ds_nonprompt_combined]:
                dsobj.Write()
            d0_nonprompt_dir_cent.cd()
            for d0obj in [h_d0_nonprompt, g_d0_nonprompt_syst, h_d0_nonprompt_combined]:
                d0obj.Write()

    with ROOT.TFile.Open(f"data/Data_pp5TeV_Run2.root", "RECREATE") as outfile:
        ds_dir = outfile.mkdir("Ds")
        dplus_dir = outfile.mkdir("Dplus")
        ds_dir_prompt = ds_dir.mkdir("Prompt")
        dplus_dir_prompt = dplus_dir.mkdir("Prompt")
        ds_dir_nonprompt = ds_dir.mkdir("NonPrompt")
        dplus_dir_nonprompt = dplus_dir.mkdir("NonPrompt")
        h_dplus_prompt, g_dplus_prompt_syst, h_dplus_prompt_combined = get_yield("Dplus", False, True)
        h_dplus_nonprompt, g_dplus_nonprompt_syst, h_dplus_nonprompt_combined = get_yield("Dplus", False, False)
        h_ds_prompt, g_ds_prompt_syst, h_ds_prompt_combined = get_yield("Ds", False, True)
        h_ds_nonprompt, g_ds_nonprompt_syst, h_ds_nonprompt_combined = get_yield("Ds", False, False)

        ds_dir_prompt.cd()
        for dsobj in [h_ds_prompt, g_ds_prompt_syst, h_ds_prompt_combined]:
            dsobj.Write()
        dplus_dir_prompt.cd()
        for dplusobj in [h_dplus_prompt, g_dplus_prompt_syst, h_dplus_prompt_combined]:
            dplusobj.Write()
        ds_dir_nonprompt.cd()
        for dsobj in [h_ds_nonprompt, g_ds_nonprompt_syst, h_ds_nonprompt_combined]:
            dsobj.Write()
        dplus_dir_nonprompt.cd()
        for dplusobj in [h_dplus_nonprompt, g_dplus_nonprompt_syst, h_dplus_nonprompt_combined]:
            dplusobj.Write()