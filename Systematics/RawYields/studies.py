import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns
import argparse
import itertools
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

outdir = "/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots"

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('ptmin', metavar=float, default=2,
                    help='minimum pT')
parser.add_argument('ptmax', metavar=float, default=2.5,
                    help='maximum pT')
args = parser.parse_args()

with open(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/results/pt{float(args.ptmin)*10:.1f}_{float(args.ptmax)*10:.1f}.pkl", "rb") as f:
    data = pickle.load(f)

Trials = np.array(data["convergedTrials"])
Trials = Trials.transpose()

min_mass = Trials[0]
max_mass = Trials[1]
rebin = Trials[2]
bkg_func = Trials[3]

df = pd.DataFrame(data, columns=["rawyieldsDs","rawyieldsDplus","ratios","sigmasDs","sigmasDplus","chi2s","convergedTrials"])
df["min_mass"] = min_mass
df["max_mass"] = max_mass
df["rebin"] = rebin
df["bkg_func"] = bkg_func

Variations = ["min_mass", "max_mass", "rebin", "bkg_func"]
combinations = set(itertools.combinations(Variations, 2))

for comb in combinations:
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
    if not os.path.exists(f"{outdir}/Distributions_Figures/{comb[0]}_{comb[1]}"):
        os.makedirs(f"{outdir}/Distributions_Figures/{comb[0]}_{comb[1]}")
    
    plt.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/{comb[0]}_{comb[1]}/pt{float(args.ptmin)*10:.1f}_{float(args.ptmax)*10:.1f}.png")
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
    if not os.path.exists(f"{outdir}/Distributions_Figures/{comb[1]}_{comb[0]}"):
        os.makedirs(f"{outdir}/Distributions_Figures/{comb[1]}_{comb[0]}")
    
    plt.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/plots/Distributions_Figures/{comb[1]}_{comb[0]}/pt{float(args.ptmin)*10:.1f}_{float(args.ptmax)*10:.1f}.png")
    #Clear figure
    plt.clf()


