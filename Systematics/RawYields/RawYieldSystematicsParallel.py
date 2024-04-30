import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
import argparse
import yaml
import itertools
import zfit
import ROOT
from particle import Particle
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
import pickle
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import pandas as pd 
import seaborn as sns


def fit(input_file, input_file_bkgtempl, output_dir, pt_min, pt_max, sigmaMultFactorSecPeak, **kwargs):
    """
    Method for fitting
    """
    mass_min, mass_max, rebin, bkgfunc, sigmaDs, sigmaDplus = kwargs.get("trial", (1.74, 2.10, 2, "chebpol2", "kFree", "kFree"))
    idx = kwargs.get("idx", 0)

    data_hdl = DataHandler(data=input_file, var_name="fM",
                           histoname=f'hMass_{pt_min*10:.0f}_{pt_max*10:.0f}',
                           limits=[mass_min,mass_max], rebin=rebin)
    data_corr_bkg = DataHandler(data=input_file_bkgtempl, var_name="fM",
                                histoname=f'hDplusTemplate_{pt_min*10:.0f}_{pt_max*10:.0f}',
                                limits=[mass_min,mass_max], rebin=rebin)

    fitter = F2MassFitter(data_hdl, name_signal_pdf=["gaussian", "gaussian"],
                          name_background_pdf=["hist", bkgfunc],
                          name=f"ds_pt{pt_min*10:.0f}_{pt_max*10:.0f}_{idx}", chi2_loss=False,
                          verbosity=1, tol=1.e-1)

    # bkg initialisation
    if bkgfunc == "expo":
        fitter.set_background_initpar(1, "lam", -2)
    elif bkgfunc == "chebpol2":
        fitter.set_background_initpar(1, "c0", 0.6)
        fitter.set_background_initpar(1, "c1", -0.2)
        fitter.set_background_initpar(1, "c2", 0.01)
    elif bkgfunc == "chebpol3":
        fitter.set_background_initpar(1, "c0", 0.4)
        fitter.set_background_initpar(1, "c1", -0.2)
        fitter.set_background_initpar(1, "c2", -0.01)
        fitter.set_background_initpar(1, "c3", 0.01)

    fitter.set_background_template(0, data_corr_bkg)
    fitter.set_background_initpar(0, "frac", 0.01, limits=[0., 1.])

    # signals initialisation
    fitter.set_particle_mass(0, pdg_id=431, limits=[Particle.from_pdgid(431).mass*0.99e-3,
                                                    Particle.from_pdgid(431).mass*1.01e-3])
    fitter.set_signal_initpar(0, "sigma", 0.008, limits=[0.002, 0.030])
    fitter.set_signal_initpar(0, "frac", 0.1, limits=[0., 1.])
    fitter.set_particle_mass(1, pdg_id=411,
                             limits=[Particle.from_pdgid(411).mass*0.99e-3,
                                     Particle.from_pdgid(411).mass*1.01e-3])
    fitter.set_signal_initpar(1, "sigma", 0.006, limits=[0.002, 0.030])
    fitter.set_signal_initpar(1, "frac", 0.1, limits=[0., 1.])

    if sigmaDplus == "kFixed":
        if sigmaMultFactorSecPeak:
            fit_result = fitter.mass_zfit() # First fit to get the sigma
            fitter.set_signal_initpar(1, "sigma", fitter.get_sigma(0)[0]*sigmaMultFactorSecPeak, fix=True)
        else:
            print("ERROR: No sigma factor provided, please check your config file")
            sys.exit()
    elif sigmaDplus == "kFixedPlus10Perc":
        if sigmaMultFactorSecPeak:
            fit_result = fitter.mass_zfit() # First fit to get the sigma
            fitter.set_signal_initpar(1, "sigma", fitter.get_sigma(0)[0]*sigmaMultFactorSecPeak*1.1, fix=True)
        else:
            print("ERROR: No sigma factor provided, please check your config file")
            sys.exit()
    elif sigmaDplus == "kFixedMinus10Perc":
        if sigmaMultFactorSecPeak:
            fit_result = fitter.mass_zfit() # First fit to get the sigma
            fitter.set_signal_initpar(1, "sigma", fitter.get_sigma(0)[0]*sigmaMultFactorSecPeak*0.9, fix=True)
        else:
            print("ERROR: No sigma factor provided, please check your config file")
            sys.exit()


    try:
        fit_result = fitter.mass_zfit()
        if fit_result.converged and kwargs.get("save", False):
            loc = ["lower left", "upper left"]
            ax_title = r"$M(\mathrm{KK\pi})$ GeV$/c^2$"
            fig = fitter.plot_mass_fit(style="ATLAS", show_extra_info=True,
                                    figsize=(8, 8), extra_info_loc=loc,
                                    axis_title=ax_title)
            figres = fitter.plot_raw_residuals(figsize=(8, 8), style="ATLAS",
                                            extra_info_loc=loc, axis_title=ax_title)
            fig.savefig(f"{output_dir}/ds_mass_pt{pt_min:.1f}_{pt_max:.1f}.pdf")
            figres.savefig(f"{output_dir}/ds_massres_pt{pt_min:.1f}_{pt_max:.1f}.pdf")

        outDict = {"rawyields": [fitter.get_raw_yield(i) for i in range(2)],
                "rawyields_bincounting": [[fitter.get_raw_yield_bincounting(i, nsigma=j) for j in kwargs.get("nsigma", [3, 5])] for i in range(2)],
                "sigma": [fitter.get_sigma(i) for i in range(2)],
                "mean": [fitter.get_mass(i) for i in range(2)],
                "chi2": fitter.get_chi2_ndf(),
                "significance": [fitter.get_significance(i) for i in range(2)],
                "signal": [fitter.get_signal(i) for i in range(2)],
                "background": [fitter.get_background(i) for i in range(2)],
                "converged": fit_result.converged}

    except Exception as e:
        print("An error occurred during mass fitting:", e)
        outDict = {"rawyields": -999.,
                "rawyields_bincounting": -999.,
                "sigma": -999.,
                "mean": -999.,
                "chi2": -999.,
                "significance": -999.,
                "signal": -999.,
                "background": -999.,
                "converged": False } 

    return outDict


def ProduceFigure(multiTrialDict, multiTrialCfg, ptMin, ptMax):
    fig, axs = plt.subplots(2, 2, figsize=(20, 15))

    colors = iter([plt.cm.tab10(i) for i in range(10)])
    colorsRatios = iter([plt.cm.tab10(i) for i in range(10)])

    # Plot the results
    axs[0, 0].errorbar(x=range(1,len(multiTrialDict["rawyieldsDs"])+1), y=multiTrialDict["rawyieldsDs"], yerr=multiTrialDict["rawyieldsDs_err"], fmt='o', c=next(colors), label=f'$D_s^+$')
    axs[0, 0].errorbar(x=range(1,len(multiTrialDict["rawyieldsDs"])+1), y=multiTrialDict["rawyieldsDplus"], yerr=multiTrialDict["rawyieldsDplus_err"], fmt='o', c=next(colors), label=f'$D^+$')
    for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
        axs[0, 0].errorbar(x=list(range(len(multiTrialDict["binCountDs"][i])*(i+1) + 1,len(multiTrialDict["binCountDs"][i])*(i+2) + 1)), y=multiTrialDict["binCountDs"][i], yerr=multiTrialDict["binCountDs_err"][i], fmt='o', c=next(colors), label=f'$D_s^+$, Bin counting {nsigma} $\sigma$')
        axs[0, 0].errorbar(x=list(range(len(multiTrialDict["binCountDs"][i])*(i+1) + 1,len(multiTrialDict["binCountDs"][i])*(i+2) + 1)), y=multiTrialDict["binCountDplus"][i], yerr=multiTrialDict["binCountDs_err"][i], fmt='o', c=next(colors), label=f'$D^+$, Bin counting {nsigma} $\sigma$')
    axs[0, 0].set_xlim(0, len(multiTrialDict["binCountDs"][i])*(len(multiTrialCfg['bincounting']['nsigma'])+1)+1)
    axs[0, 0].set_xlabel('Trial', fontsize=14)
    axs[0, 0].set_ylabel('Raw yield', fontsize=14)
    axs[0, 0].legend(fontsize=12)

    # Draw the central values
    axs[0, 0].axhline(y=multiTrialDict["hRawYieldsDsCentral"].GetBinContent(multiTrialDict["hRawYieldsDsCentral"].FindBin(ptMin+0.05)), color='r', linestyle='--')
    axs[0, 0].axhline(y=multiTrialDict["hRawYieldsDplusCentral"].GetBinContent(multiTrialDict["hRawYieldsDplusCentral"].FindBin(ptMin+0.05)), color='b', linestyle='--')

    axs[0, 1].hist(multiTrialDict["ratios"], bins=30, color=next(colorsRatios), alpha=0.7, label='Fit', histtype='stepfilled', ec="k", linewidth=2)

    for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
        axs[0, 1].hist(multiTrialDict["binCountRatios"][i], bins=30,  color=next(colorsRatios), alpha=0.3, label=f'Bin Counting {nsigma} $\sigma$', histtype='stepfilled', ec="k")
    # draw information
    axs[0, 1].set_xlabel('$D_s^+/D^+$ Ratio', fontsize=14)
    axs[0, 1].set_ylabel('Counts', fontsize=14)
    axs[0, 1].legend(fontsize=12, loc='upper right')
    info = 'Fit:\n'
    info += fr'$\mu =$ {np.mean(multiTrialDict["ratios"]):.3f}''\n'
    info += fr'$\sigma =$ {np.std(multiTrialDict["ratios"]):.3f}''\n'
    for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
        info += f'{nsigma}$\sigma$ Bin counting:\n'
        info += fr'$\mu =$ {np.mean(multiTrialDict["binCountRatios"][i]):.3f}''\n'
        info += fr'$\sigma =$ {np.std(multiTrialDict["binCountRatios"][i]):.3f}''\n'
    anchored_text_fit = AnchoredText(info,
                                     loc = 'upper left',
                                     frameon=False)
    axs[0, 1].add_artist(anchored_text_fit)
    # Draw the central values
    axs[0, 1].axvline(x=multiTrialDict["hRawYieldsDsCentral"].GetBinContent(multiTrialDict["hRawYieldsDsCentral"].FindBin(ptMin+0.05))/multiTrialDict["hRawYieldsDplusCentral"].GetBinContent(multiTrialDict["hRawYieldsDplusCentral"].FindBin(ptMin+0.05)), color='r', linestyle='--')

    axs[1, 0].errorbar(x=range(1,len(multiTrialDict["rawyieldsDs"])+1), y=multiTrialDict["sigmasDs"], yerr=multiTrialDict["sigmasDs_err"], fmt='o', c='r', label=f'$D_s^+$')
    axs[1, 0].errorbar(x=range(1,len(multiTrialDict["rawyieldsDs"])+1), y=multiTrialDict["sigmasDplus"], yerr=multiTrialDict["sigmasDplus_err"], fmt='o', c='b', label=f'$D^+$')
    axs[1, 0].set_xlim(0, len(multiTrialDict["binCountDs"][i])*(len(multiTrialCfg['bincounting']['nsigma'])+1)+1)
    axs[1, 0].set_xlabel('Trial', fontsize=14)
    axs[1, 0].set_ylabel('Width ($GeV/c^2$)', fontsize=14)
    axs[1, 0].legend(fontsize=12)
    
    # Draw the central values
    axs[1, 0].axhline(y=multiTrialDict["hSigmaDsCentral"].GetBinContent(multiTrialDict["hSigmaDsCentral"].FindBin(ptMin+0.05)), color='g', linestyle='--')
    axs[1, 0].axhline(y=multiTrialDict["hSigmaDplusCentral"].GetBinContent(multiTrialDict["hSigmaDplusCentral"].FindBin(ptMin+0.05)), color='orange', linestyle='--')

    axs[1, 1].scatter(x=range(1,len(multiTrialDict["rawyieldsDs"])+1), y=multiTrialDict["chi2s"])
    axs[1, 1].set_xlim(0, len(multiTrialDict["binCountDs"][i])*(len(multiTrialCfg['bincounting']['nsigma'])+1)+1)
    axs[1, 1].set_xlabel('Trial', fontsize=14)
    axs[1, 1].set_ylabel('$\chi^2/ndf$', fontsize=14)

    plt.show()
    fig.savefig(f'/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/MultiTrial/pt{ptMin*10}_{ptMax*10}.png',bbox_inches='tight')

    Trials = np.array(multiTrialDict["convergedTrials"])
    Trials = Trials.transpose()

    min_mass = Trials[0]
    max_mass = Trials[1]
    rebin = Trials[2]
    bkg_func = Trials[3]

    df = pd.DataFrame(multiTrialDict, columns=["rawyieldsDs","rawyieldsDplus","ratios","sigmasDs","sigmasDplus","chi2s","convergedTrials"])
    df["min_mass"] = min_mass
    df["max_mass"] = max_mass
    df["rebin"] = rebin
    df["bkg_func"] = bkg_func

    Variations = ["min_mass", "max_mass", "rebin", "bkg_func"]
    combinations = set(itertools.combinations(Variations, 2))

    plt.figure(figsize=(7, 7))
    for comb in combinations:

        ### Ratio plottings
        print(f"Plotting {comb[0]} vs {comb[1]}")
        sns.stripplot(
        data=df, x=comb[0], y="ratios", hue=comb[1],
        dodge=0.5, alpha=.5, legend=False,
        )
        sns.pointplot(
            data=df, x=comb[0], y="ratios", hue=comb[1],
            dodge=0.5, linestyle="none", errorbar=None,
            marker="_", markersize=20, markeredgewidth=3,
        )

        # if folder does not exist, create it
        if not os.path.exists(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/Ratio/{comb[0]}_{comb[1]}"):
            os.makedirs(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/Ratio/{comb[0]}_{comb[1]}")
        
        plt.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/Ratio/{comb[0]}_{comb[1]}/pt{ptMin*10:.1f}_{ptMax*10:.1f}.png")
        #Clear figure
        plt.clf()

        
        sns.stripplot(
        data=df, x=comb[1], y="ratios", hue=comb[0],
        dodge=0.5, alpha=.5, legend=False,
        )
        sns.pointplot(
            data=df, x=comb[1], y="ratios", hue=comb[0],
            dodge=0.5, linestyle="none", errorbar=None,
            marker="_", markersize=20, markeredgewidth=3,
        )

        # if folder does not exist, create it
        if not os.path.exists(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/Ratio/{comb[1]}_{comb[0]}"):
            os.makedirs(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/Ratio/{comb[1]}_{comb[0]}")
        
        plt.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/Ratio/{comb[1]}_{comb[0]}/pt{ptMin*10:.1f}_{ptMax*10:.1f}.png")
        #Clear figure
        plt.clf()

        ### Ds raw yields plottings
        sns.stripplot(
        data=df, x=comb[0], y="rawyieldsDs", hue=comb[1],
        dodge=0.5, alpha=.5, legend=False,
        )
        sns.pointplot(
            data=df, x=comb[0], y="rawyieldsDs", hue=comb[1],
            dodge=0.5, linestyle="none", errorbar=None,
            marker="_", markersize=20, markeredgewidth=3,
        )

        # if folder does not exist, create it
        if not os.path.exists(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DsRawYields/{comb[0]}_{comb[1]}"):
            os.makedirs(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DsRawYields/{comb[0]}_{comb[1]}")
        
        plt.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DsRawYields/{comb[0]}_{comb[1]}/pt{ptMin*10:.1f}_{ptMax*10:.1f}.png")
        #Clear figure
        plt.clf()

        
        sns.stripplot(
        data=df, x=comb[1], y="rawyieldsDs", hue=comb[0],
        dodge=0.5, alpha=.5, legend=False,
        )
        sns.pointplot(
            data=df, x=comb[1], y="rawyieldsDs", hue=comb[0],
            dodge=0.5, linestyle="none", errorbar=None,
            marker="_", markersize=20, markeredgewidth=3,
        )

        # if folder does not exist, create it
        if not os.path.exists(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DsRawYields/{comb[1]}_{comb[0]}"):
            os.makedirs(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DsRawYields/{comb[1]}_{comb[0]}")
        
        plt.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DsRawYields/{comb[1]}_{comb[0]}/pt{ptMin*10:.1f}_{ptMax*10:.1f}.png")
        #Clear figure

        ### Dplus raw yields plottings
        sns.stripplot(
        data=df, x=comb[0], y="rawyieldsDplus", hue=comb[1],
        dodge=0.5, alpha=.5, legend=False,
        )
        sns.pointplot(
            data=df, x=comb[0], y="rawyieldsDplus", hue=comb[1],
            dodge=0.5, linestyle="none", errorbar=None,
            marker="_", markersize=20, markeredgewidth=3,
        )

        # if folder does not exist, create it
        if not os.path.exists(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DplusRawYields/{comb[0]}_{comb[1]}"):
            os.makedirs(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DplusRawYields/{comb[0]}_{comb[1]}")
        
        plt.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DplusRawYields/{comb[0]}_{comb[1]}/pt{ptMin*10:.1f}_{ptMax*10:.1f}.png")
        #Clear figure
        plt.clf()

        
        sns.stripplot(
        data=df, x=comb[1], y="rawyieldsDplus", hue=comb[0],
        dodge=0.5, alpha=.5, legend=False,
        )
        sns.pointplot(
            data=df, x=comb[1], y="rawyieldsDplus", hue=comb[0],
            dodge=0.5, linestyle="none", errorbar=None,
            marker="_", markersize=20, markeredgewidth=3,
        )

        # if folder does not exist, create it
        if not os.path.exists(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DplusRawYields/{comb[1]}_{comb[0]}"):
            os.makedirs(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DplusRawYields/{comb[1]}_{comb[0]}")
        
        plt.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/DplusRawYields/{comb[1]}_{comb[0]}/pt{ptMin*10:.1f}_{ptMax*10:.1f}.png")
        #Clear figure
        plt.clf()



def doMultiTrial(config: dict, fitConfig, ptMin, ptMax):
    """
    Method for doing multiTrials evaluation of raw yields systematics
    """

    refFileName = config['reffilenames']['data']
    refFile = ROOT.TFile.Open(refFileName)
    hRawYieldsDsCentral = refFile.Get('hRawYields')
    hRawYieldsDsCentral.SetDirectory(0)
    hRawYieldsDplusCentral = refFile.Get('hRawYieldsSecondPeak')
    hRawYieldsDplusCentral.SetDirectory(0)
    hSigmaDsCentral = refFile.Get('hRawYieldsSigma')
    hSigmaDsCentral.SetDirectory(0)
    hSigmaDplusCentral = refFile.Get('hRawYieldsSigmaSecondPeak')
    hSigmaDplusCentral.SetDirectory(0)
    refFile.Close()

    # Do the combination of the trials
    multiTrialCfg = config['multitrial']
    mins = multiTrialCfg['mins']
    maxs = multiTrialCfg['maxs']
    rebins = multiTrialCfg['rebins']
    bkgfuncs = multiTrialCfg['bkgfuncs']
    ptbins = multiTrialCfg['ptbins']
    sigmasDs = multiTrialCfg['sigmaDs']
    sigmasDplus = multiTrialCfg['sigmaDplus']
  
    # Do the combination of the trials
    trials = list(itertools.product(mins, maxs, rebins, bkgfuncs, sigmasDs, sigmasDplus))

    rawyieldsDs = []
    rawyieldsDs_err = []
    binCountDs = [[] for _ in multiTrialCfg['bincounting']['nsigma']]
    binCountDs_err = [[] for _ in multiTrialCfg['bincounting']['nsigma']]

    rawyieldsDplus = []
    rawyieldsDplus_err = []
    binCountDplus = [[] for _ in multiTrialCfg['bincounting']['nsigma']]
    binCountDplus_err = [[] for _ in multiTrialCfg['bincounting']['nsigma']]

    sigmasDs = []
    sigmasDs_err = []

    sigmasDplus = []
    sigmasDplus_err = []

    chi2s = []
    ratios = []
    binCountRatios = [[] for _ in multiTrialCfg['bincounting']['nsigma']]

    convergedTrials = []
    
    results = []
    with ProcessPoolExecutor(max_workers=10) as executor:
        for idx, trial in enumerate(trials):
            results.append((executor.submit(fit, config['reffilenames']['data'], fitConfig["pp13TeVPrompt"]['TemplateFile'],
                "", ptMin, ptMax, nsigma=multiTrialCfg['bincounting']['nsigma'], save=False, idx=idx, sigmaMultFactorSecPeak=fitConfig["pp13TeVPrompt"]["SigmaMultFactorSecPeak"],
                trial = trial), idx))
    
    for result, idx in results:
        if result.result()['converged'] and result.result()['chi2'] < 2:
            convergedTrials.append(trials[idx])
            rawyieldsDs.append(result.result()['rawyields'][0][0])
            rawyieldsDs_err.append(result.result()['rawyields'][0][1])

            rawyieldsDplus.append(result.result()['rawyields'][1][0])
            rawyieldsDplus_err.append(result.result()['rawyields'][1][1])

            sigmasDs.append(result.result()['sigma'][0][0])
            sigmasDs_err.append(result.result()['sigma'][0][1])
            sigmasDplus.append(result.result()['sigma'][1][0])
            sigmasDplus_err.append(result.result()['sigma'][1][1])

            chi2s.append(result.result()['chi2'])
            ratios.append(rawyieldsDs[-1]/rawyieldsDplus[-1])

            for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
                binCountDs[i].append(result.result()['rawyields_bincounting'][0][i][0])
                binCountDs_err[i].append(result.result()['rawyields_bincounting'][0][i][1])
                binCountDplus[i].append(result.result()['rawyields_bincounting'][1][i][0])
                binCountDplus_err[i].append(result.result()['rawyields_bincounting'][1][i][1])
                binCountRatios[i].append(binCountDs[i][-1]/binCountDplus[i][-1])

    del results

    multiTrialDict = {  'rawyieldsDs': rawyieldsDs, 'rawyieldsDs_err': rawyieldsDs_err, \
                        'rawyieldsDplus': rawyieldsDplus, 'rawyieldsDplus_err': rawyieldsDplus_err, \
                        'binCountDs': binCountDs, 'binCountDs_err': binCountDs_err, \
                        'binCountDplus': binCountDplus, 'binCountDplus_err': binCountDplus_err, \
                        'sigmasDs': sigmasDs, 'sigmasDs_err': sigmasDs_err, \
                        'sigmasDplus': sigmasDplus, 'sigmasDplus_err': sigmasDplus_err, \
                        'chi2s': chi2s, 'ratios': ratios, 'binCountRatios': binCountRatios,\
                        'hRawYieldsDsCentral': hRawYieldsDsCentral, 'hRawYieldsDplusCentral': hRawYieldsDplusCentral,\
                        'hSigmaDsCentral': hSigmaDsCentral, 'hSigmaDplusCentral': hSigmaDplusCentral,\
                        'trials': trials, 'convergedTrials': convergedTrials}
    
    # Save the results
    with open(f'/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/results/pt{ptMin*10}_{ptMax*10}.pkl', 'wb') as f:
        pickle.dump(multiTrialDict, f)

    ProduceFigure(multiTrialDict, multiTrialCfg, ptMin, ptMax)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fit Ds')
    parser.add_argument('configFile', type=str, help='Path to the configuration file')
    parser.add_argument('FitConfig', type=str, help='Path to the fit configuration file')
    parser.add_argument('--ptmin', '-pmi', type=float, help='minimum pT', default=2.5)
    parser.add_argument('--ptmax', '-pma', type=float, help='maximum pT', default=3.0)
    args = parser.parse_args()

    # Load the configuration files
    with open(args.configFile, 'r') as f:
        config = yaml.safe_load(f)

    with open(args.FitConfig, 'r') as f:
        FitConfig = yaml.safe_load(f)

    zfit.run.set_cpus_explicit(intra=3, inter=3)

    doMultiTrial(config, FitConfig, args.ptmin, args.ptmax)   