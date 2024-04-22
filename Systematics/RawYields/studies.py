import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns

with open("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/results/pt120.0_240.0.pkl", "rb") as f:
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


sns.stripplot(
    data=df, x="rebin", y="ratios", hue="bkg_func",
    dodge=True, alpha=.5, legend=False,
)
sns.pointplot(
    data=df, x="rebin", y="ratios", hue="bkg_func",
    dodge=.4, linestyle="none", errorbar=None,
    marker="_", markersize=20, markeredgewidth=3,
)

plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/figure.png")