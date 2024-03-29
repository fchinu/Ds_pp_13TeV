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

def PlotCutVariationsOnePtBin(cfgFileName):
  
    # load inputs
    with open(cfgFileName, 'r') as file:
        config = yaml.safe_load(file)
    inDirName = config["inputs"]["directory"]
    inCommonFileNameRawY = config["inputs"]["commonfilenames"]["rawyield"]
    inCommonFileNameEff = config["inputs"]["commonfilenames"]["efficiency"]
    inCommonFileNameRatios = config["inputs"]["commonfilenames"]["ratios"]
    cutSetSuffix = config["inputs"]["cutsets"]
    nFiles = len(cutSetSuffix)

    outFileName = config["outfilename"]

    maxChi2 = config["quality"]["maxchisquare"]
    minSignif = config["quality"]["minsignif"]
    minRelSignif = config["quality"]["minrelsignif"]
    minRelEff = config["quality"]["minreleff"]
    maxRelEff = config["quality"]["maxreleff"]
    fillThrRelEff = config["quality"]["fillthrreleff"]
    fRelativeVariation = bool(config["plots"]["plotrelativevar"])
    relAssignedSyst = config["plots"]["relassignedsyst"] if "relassignedsyst" in config["plots"] else []

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
    hPromptFracsDs = []
    hPromptFracsDplus = []
    hFDFracsDs = []
    hFDFracsDplus = []

    crossSectionTitle = "dN/dpt"

    for iFile in range(nFiles):

        infile_ratios = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNameRatios}{cutSetSuffix[iFile]}.root")
        hRatios.append(infile_ratios.Get("hRatio"))
        hRatios[iFile].SetDirectory(0)
        infile_ratios.Close()

        infile_rawyield = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNameRawY}{cutSetSuffix[iFile]}.root")
        hRawYieldsDs.append(infile_rawyield.Get("hRawYields"))
        hSignificancesDs.append(infile_rawyield.Get("hRawYieldsSignificance"))
        hSoverBsDs.append(infile_rawyield.Get("hRawYieldsSoverB"))
        hchi2s.append(infile_rawyield.Get("hRawYieldsChiSquare"))
        hRawYieldsDs[iFile].SetDirectory(0)
        hSignificancesDs[iFile].SetDirectory(0)
        hSoverBsDs[iFile].SetDirectory(0)
        hchi2s[iFile].SetDirectory(0)

        hRawYieldsDplus.append(infile_rawyield.Get("hRawYieldsSecondPeak"))
        hSignificancesDplus.append(infile_rawyield.Get("hRawYieldsSignificanceSecondPeak"))
        hSoverBsDplus.append(infile_rawyield.Get("hRawYieldsSoverBSecondPeak"))
        hRawYieldsDplus[iFile].SetDirectory(0)
        hSignificancesDplus[iFile].SetDirectory(0)
        hSoverBsDplus[iFile].SetDirectory(0)
        infile_rawyield.Close()

        infile_eff = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNameEff}{cutSetSuffix[iFile]}.root")
        hEffPromptsDs.append(infile_eff.Get("Eff_DsPrompt"))
        hEffFDsDs.append(infile_eff.Get("Eff_DsFD"))
        hEffPromptsDs[iFile].SetDirectory(0)
        hEffFDsDs[iFile].SetDirectory(0)

        hEffPromptsDplus.append(infile_eff.Get("Eff_DplusPrompt"))
        hEffFDsDplus.append(infile_eff.Get("Eff_DplusFD"))
        hEffPromptsDplus[iFile].SetDirectory(0)
        hEffFDsDplus[iFile].SetDirectory(0)
        infile_eff.Close()


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

        for iFile in range(nFiles):

            Ratios.append(hRatios[iFile].GetBinContent(iPt+1))
            RatiosErr.append(hRatios[iFile].GetBinError(iPt+1))
            chi2s.append(hchi2s[iFile].GetBinContent(iPt+1))

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
            #PromptFracsDs.append(hPromptFracsDs[iFile].GetBinContent(iPt+1))
            #PromptFracsDsErr.append(hPromptFracsDs[iFile].GetBinError(iPt+1))
            #FDFracsDs.append(hFDFracsDs[iFile].GetBinContent(iPt+1))
            #FDFracsDsErr.append(hFDFracsDs[iFile].GetBinError(iPt+1))

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
            #PromptFracsDplus.append(hPromptFracsDplus[iFile].GetBinContent(iPt+1))
            #PromptFracsDplusErr.append(hPromptFracsDplus[iFile].GetBinError(iPt+1))
            #FDFracsDplus.append(hFDFracsDplus[iFile].GetBinContent(iPt+1))
            #FDFracsDplusErr.append(hFDFracsDplus[iFile].GetBinError(iPt+1))
        
        # Create figures
        fig, axs = plt.subplots(2, 4, figsize=(32, 20))

        # Raw yields
        yerrDs = np.array(RawYieldsDs)/RawYieldsDs[0]*np.sqrt(np.square(np.array(RawYieldsDsErr)/np.array(RawYieldsDs)) + (RawYieldsDsErr[0]/RawYieldsDs[0])**2)
        yerrDplus = np.array(RawYieldsDplus)/RawYieldsDplus[0]*np.sqrt(np.square(np.array(RawYieldsDplusErr)/np.array(RawYieldsDplus)) + (RawYieldsDplusErr[0]/RawYieldsDplus[0])**2)
        axs[0, 0].scatter(range(nFiles), np.array(RawYieldsDs)/RawYieldsDs[0], c='r',  label='$D_s^+$', s=100)  # Increase marker size to 100
        axs[0, 0].scatter(range(nFiles), np.array(RawYieldsDplus)/RawYieldsDplus[0], c='dodgerblue', label='$D^+$', s=100)  # Increase marker size to 100
        axs[0, 0].set_title('Raw Yields', fontsize=20)
        axs[0, 0].set_ylabel('Raw Yield/Raw Yield (central)', fontsize=16)
        axs[0, 0].set_xlabel('Cut set', fontsize=16)
        axs[0, 0].legend(fontsize=14)

        # Efficiencies
        yerrDsPrompt = np.array(EffPromptsDs)/EffPromptsDs[0]*np.sqrt(np.square(np.array(EffPromptsDsErr)/np.array(EffPromptsDs)) + (EffPromptsDsErr[0]/EffPromptsDs[0])**2)
        yerrDsFD = np.array(EffFDsDs)/EffFDsDs[0]*np.sqrt(np.square(np.array(EffFDsDsErr)/np.array(EffFDsDs)) + (EffFDsDsErr[0]/EffFDsDs[0])**2)
        axs[0, 1].scatter(range(nFiles), np.array(EffPromptsDs)/EffPromptsDs[0], c='r', label='Prompt $D_s^+$', s=100)  # Increase marker size to 100
        axs[0, 1].scatter(range(nFiles), np.array(EffFDsDs)/EffFDsDs[0], c='dodgerblue', label='FD $D_s^+$', s=100)  # Increase marker size to 100
        axs[0, 1].set_title('Efficiencies $D_s^+$', fontsize=20)
        axs[0, 1].set_ylabel('Efficiency/Efficiency (central)', fontsize=16)
        axs[0, 1].set_xlabel('Cut set', fontsize=16)
        axs[0, 1].legend(fontsize=14)

        # Efficiencies
        yerrDplusPrompt = np.array(EffPromptsDplus)/EffPromptsDplus[0]*np.sqrt(np.square(np.array(EffPromptsDplusErr)/np.array(EffPromptsDplus)) + (EffPromptsDplusErr[0]/EffPromptsDplus[0])**2)
        yerrDplusFD = np.array(EffFDsDplus)/EffFDsDplus[0]*np.sqrt(np.square(np.array(EffFDsDplusErr)/np.array(EffFDsDplus)) + (EffFDsDplusErr[0]/EffPromptsDplus[0])**2)
        axs[0, 2].scatter(range(nFiles), np.array(EffPromptsDplus)/EffPromptsDplus[0], c='r', label='Prompt $D^+$', s=100)  # Increase marker size to 100
        axs[0, 2].scatter(range(nFiles), np.array(EffFDsDplus)/EffFDsDplus[0], c='dodgerblue', label='FD $D^+$', s=100)  # Increase marker size to 100
        axs[0, 2].set_title('Efficiencies $D^+$', fontsize=20)
        axs[0, 2].set_ylabel('Efficiency/Efficiency (central)', fontsize=16)
        axs[0, 2].set_xlabel('Cut set', fontsize=16)
        axs[0, 2].legend(fontsize=14)

        # fprompt
        axs[0, 3].set_title('Prompt fraction', fontsize=20)
        axs[0, 3].set_ylabel('Prompt fraction/Prompt fraction (central)', fontsize=16)
        axs[0, 3].set_xlabel('Cut set', fontsize=16)
        # Skip for now

        # Significance
        axs[1, 0].errorbar(range(nFiles), np.array(SignificancesDs), yerr=np.array(SignificancesDsErr), fmt='o', c='r', label='$D_s^+$', markersize=10)  # Increase marker size to 10
        axs[1, 0].errorbar(range(nFiles), np.array(SignificancesDplus), yerr=np.array(SignificancesDplusErr), fmt='o', c='dodgerblue', label='$D^+$', markersize=10)  # Increase marker size to 10
        axs[1, 0].set_title('Significance', fontsize=20)
        axs[1, 0].set_ylabel('Significance', fontsize=16)
        axs[1, 0].set_xlabel('Cut set', fontsize=16)
        axs[1, 0].legend(fontsize=14)

        # S/B
        axs[1, 1].errorbar(range(nFiles), np.array(SoverBsDs), yerr=np.array(SoverBsDsErr), fmt='o', c='r', label='$D_s^+$', markersize=10)  # Increase marker size to 10
        axs[1, 1].errorbar(range(nFiles), np.array(SoverBsDplus), yerr=np.array(SoverBsDplusErr), fmt='o', c='dodgerblue', label='$D^+$', markersize=10)  # Increase marker size to 10
        axs[1, 1].set_title('S/B', fontsize=20)
        axs[1, 1].set_ylabel('S/B', fontsize=16)
        axs[1, 1].set_xlabel('Cut set', fontsize=16)
        axs[1, 1].legend(fontsize=14)

        # Ratio
        axs[1, 2].errorbar(range(nFiles), np.array(Ratios), yerr=np.array(RatiosErr), fmt='o', c='k', markersize=10)  # Increase marker size to 10
        axs[1, 2].set_title('Ratio', fontsize=20)
        axs[1, 2].set_ylabel('Ratio', fontsize=16)
        axs[1, 2].set_xlabel('Cut set', fontsize=16)

        # Ratio histogram
        axs[1, 3].hist(np.array(Ratios)/Ratios[0], bins=20, histtype='step', color='k')
        axs[1, 3].set_title('Ratio', fontsize=20)
        axs[1, 3].set_ylabel('Counts', fontsize=16)
        axs[1, 3].set_xlabel('Ratio/Ratio (central)', fontsize=16)

        rms = np.std(Ratios)
        axs[1, 3].add_patch(Rectangle((1 - rms, 0), 2*rms, axs[1, 3].get_ylim()[1], facecolor="dodgerblue", alpha=0.4, ec='darkblue', lw=5, label='RMS'))
        syst = config['plots']['relassignedsyst'][iPt]
        axs[1, 3].add_patch(Rectangle((1 - syst, 0), 2*syst, axs[1, 3].get_ylim()[1], facecolor="r", alpha=0.5, ec='firebrick', lw=5, label='Assigned syst'))

        axs[1, 3].legend(fontsize=14)
        # Save figure
        plt.savefig(f"{outFileName}_pt_{hRawYieldsDs[0].GetBinLowEdge(iPt+1)}_{hRawYieldsDs[0].GetBinLowEdge(iPt+1)+hRawYieldsDs[0].GetBinWidth(iPt+1)}.png", bbox_inches='tight')

if __name__ == "__main__":
    args = argparse.ArgumentParser(description='Plot cut variations for one pt bin')
    args.add_argument('config', type=str, help='Path to the configuration file')
    args = args.parse_args()

    PlotCutVariationsOnePtBin(args.config)