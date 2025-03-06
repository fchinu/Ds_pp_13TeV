import yaml
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib.patches import Rectangle
'''
    Macro for the evaluation of cut-variation systematic uncertainty
    Main Function: PlotCutVariationsOnePtBin                        
'''

def get_raw_prompt_fraction(corry_p, corry_np, effacc_p, effacc_np):
    """
    Helper function to get the raw prompt fraction given the efficiencies

    Parameters
    -----------------------------------------------------
    - corry_p: float
        corrected yield for prompt signal
    - corry_np: float
        corrected yield for non-prompt signal
    - effacc_p: str
        eff x acc for prompt signal
    - effacc_np: str
        eff x acc for non-prompt signal

    Returns
    -----------------------------------------------------
    - f_p: (float, float)
        raw prompt fraction with its uncertainty
    """

    rawy_p = effacc_p * corry_p
    rawy_np = effacc_np * corry_np
    f_p = rawy_p / (rawy_p + rawy_np)

    return f_p

def get_ratio(rawy_ds, rawy_dplus, eff_ds, eff_dplus, frac_ds, frac_dplus, config):
    """
    Helper function to get the ratio of Ds and Dplus raw yields

    Parameters
    -----------------------------------------------------
    - rawy_ds: float
        raw yield for Ds
    - rawy_dplus: float
        raw yield for Dplus
    - eff_ds: float
        efficiency for Ds
    - eff_dplus: float
        efficiency for Dplus
    - frac_ds: float
        prompt fraction for Ds
    - frac_dplus: float
        prompt fraction for Dplus

    Returns
    -----------------------------------------------------
    - ratio: (float, float)
        ratio of Ds and Dplus raw yields with its uncertainty
    """

    ratio = (rawy_ds * frac_ds / eff_ds / config["br"]["ds_to_phipi_to_kkpi"]) / (rawy_dplus * frac_dplus / eff_dplus / config["br"]["dplus_to_phipi_to_kkpi"])

    return ratio


def PlotCutVariationsOnePtBin(cfgFileName):
  
    # load inputs
    with open(cfgFileName, 'r') as file:
        config = yaml.safe_load(file)

    with open(config["cutset_file"], 'r') as file:
        cutset = yaml.safe_load(file)

    cent_mins = cutset["cent"]["min"][1:] # do not consider 0-100%
    cent_maxs = cutset["cent"]["max"][1:] # do not consider 0-100%

    inDirName = config["inputs"]["directory"]
    inCommonFileNameRawY = config["inputs"]["commonfilenames"]["rawyield"]
    inCommonFileNameEff = config["inputs"]["commonfilenames"]["efficiency"]
    inFileNamePromptFracDs = config["inputs"]["prompt_frac_ds"]
    inFileNamePromptFracDplus = config["inputs"]["prompt_frac_dplus"]
    cutSetSuffix = config["inputs"]["cutsets"]
    nFiles = len(cutSetSuffix)

    outFileName = config["outfilename"]
    outFile = ROOT.TFile(outFileName, "RECREATE")
    outFile.Close()

    maxChi2 = config["quality"]["maxchisquare"]
    minSignif = config["quality"]["minsignif"]
    minRelSignif = config["quality"]["minrelsignif"]
    minRelEff = config["quality"]["minreleff"]
    maxRelEff = config["quality"]["maxreleff"]
    fillThrRelEff = config["quality"]["fillthrreleff"]
    fRelativeVariation = bool(config["plots"]["plotrelativevar"])
    relAssignedSyst = config["plots"]["relassignedsyst"] if "relassignedsyst" in config["plots"] else []


    for i_cent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):

        hRatios = []
        hRawYieldsDs = []
        hRawYieldsDplus = []
        hSignificancesDs = []
        hSignificancesDplus = []
        hSoverBsDs = []
        hSoverBsDplus = []
        hchi2s = []
        hEffPromptsDs = []
        hEffPromptsDplus = []
        hEffFDsDs = []
        hEffFDsDplus = []
        for iFile in range(nFiles):

            infile_rawyield = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNameRawY}{cutSetSuffix[iFile]}.root")
            hRawYieldsDs.append(infile_rawyield.Get(f"h_raw_yields_ds_{cent_min}_{cent_max}"))
            hSignificancesDs.append(infile_rawyield.Get(f"h_significance_ds_{cent_min}_{cent_max}"))
            hSoverBsDs.append(infile_rawyield.Get(f"h_s_over_b_ds_{cent_min}_{cent_max}"))
            hchi2s.append(infile_rawyield.Get(f"h_chi2_{cent_min}_{cent_max}"))
            hRawYieldsDs[iFile].SetDirectory(0)
            hSignificancesDs[iFile].SetDirectory(0)
            hSoverBsDs[iFile].SetDirectory(0)
            hchi2s[iFile].SetDirectory(0)

            hRawYieldsDplus.append(infile_rawyield.Get(f"h_raw_yields_dplus_{cent_min}_{cent_max}"))
            hSignificancesDplus.append(infile_rawyield.Get(f"h_significance_dplus_{cent_min}_{cent_max}"))
            hSoverBsDplus.append(infile_rawyield.Get(f"h_s_over_b_dplus_{cent_min}_{cent_max}"))
            hRawYieldsDplus[iFile].SetDirectory(0)
            hSignificancesDplus[iFile].SetDirectory(0)
            hSoverBsDplus[iFile].SetDirectory(0)
            infile_rawyield.Close()

            infile_eff = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNameEff}{cutSetSuffix[iFile]}.root")
            hEffPromptsDs.append(infile_eff.Get(f"eff_DsPrompt_cent_{cent_min}_{cent_max}"))
            hEffFDsDs.append(infile_eff.Get(f"eff_DsNonPrompt_cent_{cent_min}_{cent_max}"))
            hEffPromptsDs[iFile].SetDirectory(0)
            hEffFDsDs[iFile].SetDirectory(0)

            hEffPromptsDplus.append(infile_eff.Get(f"eff_DplusPrompt_cent_{cent_min}_{cent_max}"))
            hEffFDsDplus.append(infile_eff.Get(f"eff_DplusNonPrompt_cent_{cent_min}_{cent_max}"))
            hEffPromptsDplus[iFile].SetDirectory(0)
            hEffFDsDplus[iFile].SetDirectory(0)
            infile_eff.Close()

        inFileNamePromptFracDs_cent = inFileNamePromptFracDs.replace("MB.root", f"{cent_min}_{cent_max}.root")
        infile_DsPromptFrac = ROOT.TFile.Open(inFileNamePromptFracDs_cent)
        hNCorrPromptDs = infile_DsPromptFrac.Get("hCorrYieldsPrompt")
        hNCorrNonPromptDs = infile_DsPromptFrac.Get("hCorrYieldsNonPrompt")
        hNCorrPromptDs.SetDirectory(0)
        hNCorrNonPromptDs.SetDirectory(0)
        infile_DsPromptFrac.Close()

        inFileNamePromptFracDplus_cent = inFileNamePromptFracDplus.replace("MB.root", f"{cent_min}_{cent_max}.root")
        infile_DplusPromptFrac = ROOT.TFile.Open(inFileNamePromptFracDplus_cent)
        hNCorrPromptDplus = infile_DplusPromptFrac.Get("hCorrYieldsPrompt")
        hNCorrNonPromptDplus = infile_DplusPromptFrac.Get("hCorrYieldsNonPrompt")
        hNCorrPromptDplus.SetDirectory(0)
        hNCorrNonPromptDplus.SetDirectory(0)
        infile_DplusPromptFrac.Close()

        hSyst = hRawYieldsDs[0].Clone("hSyst")
        hAssignedSyst = hRawYieldsDs[0].Clone("hAssignedSyst")

        for iPt in range(hRawYieldsDs[0].GetNbinsX()):
            
            Ratios = []
            RawYieldsDs = []
            RawYieldsDplus = []
            SignificancesDs = []
            SignificancesDplus = []
            SoverBsDs = []
            SoverBsDplus = []
            chi2s = []
            EffPromptsDs = []
            EffPromptsDplus = []
            EffFDsDs = []
            EffFDsDplus = []
            PromptFracsDs = []
            PromptFracsDplus = []
            FDFracsDs = []
            FDFracsDplus = []

            RatiosErr = []
            RawYieldsDsErr = []
            RawYieldsDplusErr = []
            SignificancesDsErr = []
            SignificancesDplusErr = []
            SoverBsDsErr = []
            SoverBsDplusErr = []
            EffPromptsDsErr = []
            EffPromptsDplusErr = []
            EffFDsDsErr = []
            EffFDsDplusErr = []
            PromptFracsDsErr = []
            PromptFracsDplusErr = []
            FDFracsDsErr = []
            FDFracsDplusErr = []

            pass_counter = 0
            for iFile in range(nFiles):
                
                # chi2 quality check
                if hchi2s[iFile].GetBinContent(iPt+1) > maxChi2:
                    continue

                # Significance quality check
                if hSignificancesDs[iFile].GetBinContent(iPt+1) < minSignif or hSignificancesDplus[iFile].GetBinContent(iPt+1) < minSignif:
                    continue
                if hSignificancesDs[iFile].GetBinContent(iPt+1) < hSignificancesDs[0].GetBinContent(iPt+1) * minRelSignif or hSignificancesDplus[iFile].GetBinContent(iPt+1) < hSignificancesDplus[0].GetBinContent(iPt+1) * minRelSignif:
                    continue
                
                # Efficiency quality check
                if hEffPromptsDs[iFile].GetBinContent(iPt+1) < hEffPromptsDs[0].GetBinContent(iPt+1) * minRelEff or hEffPromptsDplus[iFile].GetBinContent(iPt+1) < hEffPromptsDplus[iFile].GetBinContent(iPt+1) * minRelEff:
                    continue


                RawYieldsDs.append(hRawYieldsDs[iFile].GetBinContent(iPt+1))
                RawYieldsDsErr.append(hRawYieldsDs[iFile].GetBinError(iPt+1))
                SignificancesDs.append(hSignificancesDs[iFile].GetBinContent(iPt+1))
                SignificancesDsErr.append(hSignificancesDs[iFile].GetBinError(iPt+1))
                SoverBsDs.append(hSoverBsDs[iFile].GetBinContent(iPt+1))
                SoverBsDsErr.append(hSoverBsDs[iFile].GetBinError(iPt+1))
                EffPromptsDs.append(hEffPromptsDs[iFile].GetBinContent(iPt+1))
                EffPromptsDsErr.append(hEffPromptsDs[iFile].GetBinError(iPt+1))
                EffFDsDs.append(hEffFDsDs[iFile].GetBinContent(iPt+1))
                EffFDsDsErr.append(hEffFDsDs[iFile].GetBinError(iPt+1))
                PromptFracsDs.append(get_raw_prompt_fraction(
                    hNCorrPromptDs.GetBinContent(iPt+1),
                    hNCorrNonPromptDs.GetBinContent(iPt+1),
                    hEffPromptsDs[iFile].GetBinContent(iPt+1),
                    hEffFDsDs[iFile].GetBinContent(iPt+1)
                ))

                RawYieldsDplus.append(hRawYieldsDplus[iFile].GetBinContent(iPt+1))
                RawYieldsDplusErr.append(hRawYieldsDplus[iFile].GetBinError(iPt+1))
                SignificancesDplus.append(hSignificancesDplus[iFile].GetBinContent(iPt+1))
                SignificancesDplusErr.append(hSignificancesDplus[iFile].GetBinError(iPt+1))
                SoverBsDplus.append(hSoverBsDplus[iFile].GetBinContent(iPt+1))
                SoverBsDplusErr.append(hSoverBsDplus[iFile].GetBinError(iPt+1))
                EffPromptsDplus.append(hEffPromptsDplus[iFile].GetBinContent(iPt+1))
                EffPromptsDplusErr.append(hEffPromptsDplus[iFile].GetBinError(iPt+1))
                EffFDsDplus.append(hEffFDsDplus[iFile].GetBinContent(iPt+1))
                EffFDsDplusErr.append(hEffFDsDplus[iFile].GetBinError(iPt+1))
                PromptFracsDplus.append(get_raw_prompt_fraction(
                    hNCorrPromptDplus.GetBinContent(iPt+1),
                    hNCorrNonPromptDplus.GetBinContent(iPt+1),
                    hEffPromptsDplus[iFile].GetBinContent(iPt+1),
                    hEffFDsDplus[iFile].GetBinContent(iPt+1)
                ))

                Ratios.append(get_ratio(
                    hRawYieldsDs[iFile].GetBinContent(iPt+1),
                    hRawYieldsDplus[iFile].GetBinContent(iPt+1),
                    hEffPromptsDs[iFile].GetBinContent(iPt+1),
                    hEffPromptsDplus[iFile].GetBinContent(iPt+1),
                    PromptFracsDs[pass_counter],
                    PromptFracsDplus[pass_counter],
                    config
                ))
                ratio_err = Ratios[-1] * np.sqrt(
                    (RawYieldsDsErr[-1]/RawYieldsDs[-1])**2 +
                    (RawYieldsDplusErr[-1]/RawYieldsDplus[-1])**2 +
                    (EffPromptsDsErr[-1]/EffPromptsDs[-1])**2 +
                    (EffPromptsDplusErr[-1]/EffPromptsDplus[-1])**2
                )
                RatiosErr.append(ratio_err)
                chi2s.append(hchi2s[iFile].GetBinContent(iPt+1))
                pass_counter += 1
            
            # Create figures
            fig, axs = plt.subplots(2, 4, figsize=(32, 20))

            # Raw yields
            yerrDs = np.array(RawYieldsDs)/RawYieldsDs[0]*np.sqrt(np.square(np.array(RawYieldsDsErr)/np.array(RawYieldsDs)) + (RawYieldsDsErr[0]/RawYieldsDs[0])**2)
            yerrDplus = np.array(RawYieldsDplus)/RawYieldsDplus[0]*np.sqrt(np.square(np.array(RawYieldsDplusErr)/np.array(RawYieldsDplus)) + (RawYieldsDplusErr[0]/RawYieldsDplus[0])**2)
            axs[0, 0].scatter(range(len(RawYieldsDs)), np.array(RawYieldsDs)/RawYieldsDs[0], c='r',  label='$D_s^+$', s=100)  # Increase marker size to 100
            axs[0, 0].scatter(range(len(RawYieldsDs)), np.array(RawYieldsDplus)/RawYieldsDplus[0], c='dodgerblue', label='$D^+$', s=100)  # Increase marker size to 100
            axs[0, 0].set_title('Raw Yields', fontsize=20)
            axs[0, 0].set_ylabel('Raw Yield/Raw Yield (central)', fontsize=16)
            axs[0, 0].set_xlabel('Cut set', fontsize=16)
            axs[0, 0].legend(fontsize=14)

            # Efficiencies
            yerrDsPrompt = np.array(EffPromptsDs)/EffPromptsDs[0]*np.sqrt(np.square(np.array(EffPromptsDsErr)/np.array(EffPromptsDs)) + (EffPromptsDsErr[0]/EffPromptsDs[0])**2)
            yerrDsFD = np.array(EffFDsDs)/EffFDsDs[0]*np.sqrt(np.square(np.array(EffFDsDsErr)/np.array(EffFDsDs)) + (EffFDsDsErr[0]/EffFDsDs[0])**2)
            axs[0, 1].scatter(range(len(RawYieldsDs)), np.array(EffPromptsDs)/EffPromptsDs[0], c='r', label='Prompt $D_s^+$', s=100)  # Increase marker size to 100
            axs[0, 1].scatter(range(len(RawYieldsDs)), np.array(EffFDsDs)/EffFDsDs[0], c='dodgerblue', label='FD $D_s^+$', s=100)  # Increase marker size to 100
            axs[0, 1].set_title('Efficiencies $D_s^+$', fontsize=20)
            axs[0, 1].set_ylabel('Efficiency/Efficiency (central)', fontsize=16)
            axs[0, 1].set_xlabel('Cut set', fontsize=16)
            axs[0, 1].legend(fontsize=14)

            # Efficiencies
            yerrDplusPrompt = np.array(EffPromptsDplus)/EffPromptsDplus[0]*np.sqrt(np.square(np.array(EffPromptsDplusErr)/np.array(EffPromptsDplus)) + (EffPromptsDplusErr[0]/EffPromptsDplus[0])**2)
            yerrDplusFD = np.array(EffFDsDplus)/EffFDsDplus[0]*np.sqrt(np.square(np.array(EffFDsDplusErr)/np.array(EffFDsDplus)) + (EffFDsDplusErr[0]/EffPromptsDplus[0])**2)
            axs[0, 2].scatter(range(len(RawYieldsDs)), np.array(EffPromptsDplus)/EffPromptsDplus[0], c='r', label='Prompt $D^+$', s=100)  # Increase marker size to 100
            axs[0, 2].scatter(range(len(RawYieldsDs)), np.array(EffFDsDplus)/EffFDsDplus[0], c='dodgerblue', label='FD $D^+$', s=100)  # Increase marker size to 100
            axs[0, 2].set_title('Efficiencies $D^+$', fontsize=20)
            axs[0, 2].set_ylabel('Efficiency/Efficiency (central)', fontsize=16)
            axs[0, 2].set_xlabel('Cut set', fontsize=16)
            axs[0, 2].legend(fontsize=14)

            # fprompt
            axs[0, 3].scatter(range(len(RawYieldsDs)), np.array(PromptFracsDs)/PromptFracsDs[0], c='r', label='$D_s^+$', s=100)  # Increase marker size to 100
            axs[0, 3].scatter(range(len(RawYieldsDs)), np.array(PromptFracsDplus)/PromptFracsDplus[0], c='dodgerblue', label='$D^+$', s=100)  # Increase marker size to 100
            axs[0, 3].set_title('Prompt fraction', fontsize=20)
            axs[0, 3].set_ylabel('Prompt fraction/Prompt fraction (central)', fontsize=16)
            axs[0, 3].set_xlabel('Cut set', fontsize=16)
            axs[0, 3].legend(fontsize=14)

            # Significance
            axs[1, 0].errorbar(range(len(RawYieldsDs)), np.array(SignificancesDs), yerr=np.array(SignificancesDsErr), fmt='o', c='r', label='$D_s^+$', markersize=10)  # Increase marker size to 10
            axs[1, 0].errorbar(range(len(RawYieldsDs)), np.array(SignificancesDplus), yerr=np.array(SignificancesDplusErr), fmt='o', c='dodgerblue', label='$D^+$', markersize=10)  # Increase marker size to 10
            axs[1, 0].set_title('Significance', fontsize=20)
            axs[1, 0].set_ylabel('Significance', fontsize=16)
            axs[1, 0].set_xlabel('Cut set', fontsize=16)
            axs[1, 0].legend(fontsize=14)

            # S/B
            axs[1, 1].errorbar(range(len(RawYieldsDs)), np.array(SoverBsDs), yerr=np.array(SoverBsDsErr), fmt='o', c='r', label='$D_s^+$', markersize=10)  # Increase marker size to 10
            axs[1, 1].errorbar(range(len(RawYieldsDs)), np.array(SoverBsDplus), yerr=np.array(SoverBsDplusErr), fmt='o', c='dodgerblue', label='$D^+$', markersize=10)  # Increase marker size to 10
            axs[1, 1].set_title('S/B', fontsize=20)
            axs[1, 1].set_ylabel('S/B', fontsize=16)
            axs[1, 1].set_xlabel('Cut set', fontsize=16)
            axs[1, 1].legend(fontsize=14)

            # Ratio
            axs[1, 2].errorbar(range(len(RawYieldsDs)), np.array(Ratios), yerr=np.array(RatiosErr), fmt='o', c='k', markersize=10)  # Increase marker size to 10
            axs[1, 2].set_title('Ratio', fontsize=20)
            axs[1, 2].set_ylabel('Ratio', fontsize=16)
            axs[1, 2].set_xlabel('Cut set', fontsize=16)

            # Ratio histogram
            axs[1, 3].hist(np.array(Ratios)/Ratios[0], bins=20, histtype='step', color='k')
            axs[1, 3].set_title('Ratio', fontsize=20)
            axs[1, 3].set_ylabel('Counts', fontsize=16)
            axs[1, 3].set_xlabel('Ratio/Ratio (central)', fontsize=16)

            rms = np.std(Ratios)
            mean = np.mean(Ratios)
            rmsAndShift = np.sqrt((rms/Ratios[0])**2 + (mean/Ratios[0] - 1)**2)

            print(f"{hRawYieldsDs[0].GetBinLowEdge(iPt+1)} < pT < {hRawYieldsDs[0].GetBinLowEdge(iPt+1)+hRawYieldsDs[0].GetBinWidth(iPt+1)}: RMS+shift = {rmsAndShift}")
            hSyst.SetBinContent(iPt+1, rmsAndShift)
            hSyst.SetBinError(iPt+1, 0)

            axs[1, 3].add_patch(Rectangle((1 - rmsAndShift, 0), 2*rmsAndShift, axs[1, 3].get_ylim()[1], facecolor="dodgerblue", alpha=0.4, ec='darkblue', lw=5, label=f'$\sqrt{{RMS^2 + shift^2}}$ = {rmsAndShift*100:.2f}%'))
            syst = config['plots']['relassignedsyst'][iPt][i_cent]
            axs[1, 3].add_patch(Rectangle((1 - syst, 0), 2*syst, axs[1, 3].get_ylim()[1], facecolor="r", alpha=0.5, ec='firebrick', lw=5, label='Assigned syst'))
            hAssignedSyst.SetBinContent(iPt+1, syst)
            hAssignedSyst.SetBinError(iPt+1, 0)

            xMax = max(1+rmsAndShift, 1+syst, max(Ratios)/Ratios[0])
            xMin = min(1-rmsAndShift, 1-syst, min(Ratios)/Ratios[0])
            axs[1, 3].set_xlim(xMin*0.95, xMax*1.05)

            axs[1, 3].legend(fontsize=14)

            # Save figure and ROOT file
            plt.savefig(f"{outFileName}_pt_{hRawYieldsDs[0].GetBinLowEdge(iPt+1)}_{hRawYieldsDs[0].GetBinLowEdge(iPt+1)+hRawYieldsDs[0].GetBinWidth(iPt+1)}_{cent_min}_{cent_max}.png", bbox_inches='tight')

        outFile = ROOT.TFile(outFileName, "UPDATE")
        hSyst.Write(f"rms_shifts_sum_quadrature_{cent_min}_{cent_max}")
        hAssignedSyst.Write(f"assigned_syst_{cent_min}_{cent_max}")
        outFile.Close()

if __name__ == "__main__":
    args = argparse.ArgumentParser(description='Plot cut variations for one pt bin')
    args.add_argument('config', type=str, help='Path to the configuration file')
    args = args.parse_args()

    PlotCutVariationsOnePtBin(args.config)