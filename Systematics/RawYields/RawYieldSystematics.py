import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
import argparse
import yaml
from itertools import product
import zfit
import ROOT
from particle import Particle
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
import concurrent.futures
import pickle
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""

def fit(dataHandl, data_corr_bkg, bkgfunc, idx):

    fitter = F2MassFitter(dataHandl, name_signal_pdf=["gaussian", "gaussian"],
                          name_background_pdf=[bkgfunc, "hist"],
                          name=f"ds_pt{idx}", chi2_loss=False,
                          verbosity=7, tol=1.e-1)

    # bkg initialisation
    fitter.set_background_initpar(0, "c0", 0.6)
    fitter.set_background_initpar(0, "c1", -0.2)
    fitter.set_background_initpar(0, "c2", 0.008)
    fitter.set_background_initpar(0, "frac", 0.7)
    fitter.set_background_template(1, data_corr_bkg)

    # signals initialisation
    fitter.set_particle_mass(0, pdg_id=431,
                             limits=[Particle.from_pdgid(431).mass*0.97e-3,
                                     Particle.from_pdgid(431).mass*1.03e-3])
    fitter.set_signal_initpar(0, "sigma", 0.01, limits=[0.005, 0.03])
    fitter.set_signal_initpar(0, "frac", 0.1, limits=[0., 1.])
    fitter.set_particle_mass(1, pdg_id=411,
                             limits=[Particle.from_pdgid(411).mass*0.97e-3,
                                     Particle.from_pdgid(411).mass*1.03e-3])
    fitter.set_signal_initpar(1, "sigma", 0.01, limits=[0.005, 0.03])
    fitter.set_signal_initpar(1, "frac", 0.1, limits=[0., 1.])

    fit_result = fitter.mass_zfit()

    return fitter, fit_result

def ProcessThisTrial(trial, config, ptMin, ptMax, idx):
    mass_min, mass_max, rebin, bkgfunc = trial
    print(f"mass_min: {mass_min}, mass_max: {mass_max}, rebin: {rebin}, bkgfunc: {bkgfunc}")
    
    dataHandl = DataHandler(data=config['reffilenames']['data'], var_name="fM",
                            histoname=f'hMass_{ptMin*10:.0f}_{ptMax*10:.0f}',
                            limits=[mass_min,mass_max], rebin=rebin)
    
    data_corr_bkg = DataHandler(data=config['reffilenames']['template'], var_name="fM",
                                    histoname=f'hDplusTemplate_{ptMin*10:.0f}_{ptMax*10:.0f}',
                                    limits=[mass_min,mass_max], rebin=rebin)
    
    # Do the fit
    return fit(dataHandl, data_corr_bkg, bkgfunc, idx)


def doMultiTrial(config: dict, ptMin, ptMax):
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

  
    # Do the combination of the trials
    trials = list(product(mins, maxs, rebins, bkgfuncs))

    fig, axs = plt.subplots(2, 2, figsize=(20, 15))

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
    
    for idx, trial in enumerate(trials):
        fitter, fit_result = ProcessThisTrial(trials[idx], config, ptMin, ptMax, idx)
        
        if fit_result.converged:
            rawyieldsDs.append(fitter.get_raw_yield(0)[0])
            rawyieldsDs_err.append(fitter.get_raw_yield(0)[1])

            rawyieldsDplus.append(fitter.get_raw_yield(1)[0])
            rawyieldsDplus_err.append(fitter.get_raw_yield(1)[1])

            sigmasDs.append(fitter.get_sigma(0)[0])
            sigmasDs_err.append(fitter.get_sigma(0)[1])
            sigmasDplus.append(fitter.get_sigma(1)[0])
            sigmasDplus_err.append(fitter.get_sigma(1)[1])

            chi2s.append(fitter.get_chi2_ndf())
            ratios.append(rawyieldsDs[-1]/rawyieldsDplus[-1])

            for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
                binCountDs[i].append(fitter.get_raw_yield_bincounting(0, nsigma=nsigma)[0])
                binCountDs_err[i].append(fitter.get_raw_yield_bincounting(0, nsigma=nsigma)[1])
                binCountDplus[i].append(fitter.get_raw_yield_bincounting(1, nsigma=nsigma)[0])
                binCountDplus_err[i].append(fitter.get_raw_yield_bincounting(1, nsigma=nsigma)[1])
                binCountRatios[i].append(binCountDs[i][-1]/binCountDplus[i][-1])

        del fitter

    multiTrialDict = {  'rawyieldsDs': rawyieldsDs, 'rawyieldsDs_err': rawyieldsDs_err, \
                        'rawyieldsDplus': rawyieldsDplus, 'rawyieldsDplus_err': rawyieldsDplus_err, \
                        'binCountDs': binCountDs, 'binCountDs_err': binCountDs_err, \
                        'binCountDplus': binCountDplus, 'binCountDplus_err': binCountDplus_err, \
                        'sigmasDs': sigmasDs, 'sigmasDs_err': sigmasDs_err, \
                        'sigmasDplus': sigmasDplus, 'sigmasDplus_err': sigmasDplus_err, \
                        'chi2s': chi2s, 'ratios': ratios, 'binCountRatios': binCountRatios}
    
    # Save the results
    with open(f'/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/results/pt{ptMin*10}_{ptMax*10}.pkl', 'wb') as f:
        pickle.dump(multiTrialDict, f)
        
    colors = iter([plt.cm.tab10(i) for i in range(10)])
    colorsRatios = iter([plt.cm.tab10(i) for i in range(10)])

    # Plot the results
    axs[0, 0].errorbar(x=range(1,len(rawyieldsDs)+1), y=rawyieldsDs, yerr=rawyieldsDs_err, fmt='o', c=next(colors), label=f'$D_s^+$')
    axs[0, 0].errorbar(x=range(1,len(rawyieldsDs)+1), y=rawyieldsDplus, yerr=rawyieldsDplus_err, fmt='o', c=next(colors), label=f'$D^+$')
    for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
        axs[0, 0].errorbar(x=list(range(len(trials)*(i+1) + 1,len(trials)*(i+2) + 1)), y=binCountDs[i], yerr=binCountDs_err[i], fmt='o', c=next(colors), label=f'$D_s^+$, Bin counting {nsigma} $\sigma$')
        axs[0, 0].errorbar(x=list(range(len(trials)*(i+1) + 1,len(trials)*(i+2) + 1)), y=binCountDplus[i], yerr=binCountDs_err[i], fmt='o', c=next(colors), label=f'$D^+$, Bin counting {nsigma} $\sigma$')
    axs[0, 0].set_xlim(0, len(trials)*(len(multiTrialCfg['bincounting']['nsigma'])+1)+1)
    axs[0, 0].set_xlabel('Trial', fontsize=14)
    axs[0, 0].set_ylabel('Raw yield', fontsize=14)
    axs[0, 0].legend(fontsize=12)

    # Draw the central values
    axs[0, 0].axhline(y=hRawYieldsDsCentral.GetBinContent(hRawYieldsDsCentral.FindBin(ptMin+0.05)), color='r', linestyle='--')
    axs[0, 0].axhline(y=hRawYieldsDplusCentral.GetBinContent(hRawYieldsDplusCentral.FindBin(ptMin+0.05)), color='b', linestyle='--')

    axs[0, 1].hist(ratios, bins=30, color=next(colorsRatios), alpha=0.3, label='Fit', histtype='stepfilled', ec="k")

    for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
        axs[0, 1].hist(binCountRatios[i], bins=30,  color=next(colorsRatios), alpha=0.3, label=f'Bin Counting {nsigma} $\sigma$', histtype='stepfilled', ec="k")
    # draw information
    axs[0, 1].set_xlabel('$D_s^+/D^+$ Ratio', fontsize=14)
    axs[0, 1].set_ylabel('Counts', fontsize=14)
    axs[0, 1].legend(fontsize=12, loc='upper right')
    info = 'Fit:\n'
    info += fr'$\mu =$ {np.mean(ratios):.3f}''\n'
    info += fr'$\sigma =$ {np.std(ratios):.3f}''\n'
    for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
        info += f'{nsigma}$\sigma$ Bin counting:\n'
        info += fr'$\mu =$ {np.mean(binCountRatios[i]):.3f}''\n'
        info += fr'$\sigma =$ {np.std(binCountRatios[i]):.3f}''\n'
    anchored_text_fit = AnchoredText(info,
                                     loc = 'upper left',
                                     frameon=False)
    axs[0, 1].add_artist(anchored_text_fit)
    # Draw the central values
    axs[0, 1].axvline(x=hRawYieldsDsCentral.GetBinContent(hRawYieldsDsCentral.FindBin(ptMin+0.05))/hRawYieldsDplusCentral.GetBinContent(hRawYieldsDplusCentral.FindBin(ptMin+0.05)), color='r', linestyle='--')

    axs[1, 0].errorbar(x=range(1,len(rawyieldsDs)+1), y=sigmasDs, yerr=sigmasDs_err, fmt='o', c='r', label=f'$D_s^+$')
    axs[1, 0].errorbar(x=range(1,len(rawyieldsDs)+1), y=sigmasDplus, yerr=sigmasDplus_err, fmt='o', c='b', label=f'$D^+$')
    axs[1, 0].set_xlim(0, len(trials)*(len(multiTrialCfg['bincounting']['nsigma'])+1)+1)
    axs[1, 0].set_xlabel('Trial', fontsize=14)
    axs[1, 0].set_ylabel('Width ($Mev/c^2$)', fontsize=14)
    axs[1, 0].legend(fontsize=12)
    
    # Draw the central values
    axs[1, 0].axhline(y=hSigmaDsCentral.GetBinContent(hSigmaDsCentral.FindBin(ptMin+0.05)), color='g', linestyle='--')
    axs[1, 0].axhline(y=hSigmaDplusCentral.GetBinContent(hSigmaDplusCentral.FindBin(ptMin+0.05)), color='orange', linestyle='--')

    axs[1, 1].scatter(x=range(1,len(rawyieldsDs)+1), y=chi2s)
    axs[1, 1].set_xlim(0, len(trials)*(len(multiTrialCfg['bincounting']['nsigma'])+1)+1)
    axs[1, 1].set_xlabel('Trial', fontsize=14)
    axs[1, 1].set_ylabel('$\chi^2/ndf$', fontsize=14)

    plt.show()
    fig.savefig(f'/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/pt{ptMin*10}_{ptMax*10}.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fit Ds')
    parser.add_argument('configFile', type=str, help='Path to the configuration file')
    parser.add_argument('--ptmin', '-pmi', type=float, help='minimum pT', default=2.5)
    parser.add_argument('--ptmax', '-pma', type=float, help='maximum pT', default=3.0)
    args = parser.parse_args()

    # Load the configuration file
    with open(args.configFile, 'r') as f:
        config = yaml.safe_load(f)

    zfit.run.set_cpus_explicit(intra=30, inter=30)

    doMultiTrial(config, args.ptmin, args.ptmax)