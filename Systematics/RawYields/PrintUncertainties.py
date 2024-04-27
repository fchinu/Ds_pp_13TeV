import numpy as np
import pickle

ptMins = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,8,12] 
ptMaxs = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,8,12,24]

for ptMin,ptMax in zip(ptMins,ptMaxs):
    with open(f"/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/results/pt{ptMin*10:.1f}_{ptMax*10:.1f}.pkl", "rb") as f:
        results = pickle.load(f)
    central_value = results["hRawYieldsDsCentral"].GetBinContent(results["hRawYieldsDsCentral"].FindBin(ptMin+0.05))
    central_value = central_value/results["hRawYieldsDplusCentral"].GetBinContent(results["hRawYieldsDplusCentral"].FindBin(ptMin+0.05))
    mean = np.mean(results["ratios"])
    rms = np.std(results["ratios"])
    shift = mean - central_value
    syst = np.sqrt(rms**2 + shift**2)
    syst = syst/central_value
    #print(syst, central_value, mean, rms)
    print(f"{ptMin:.1f}-{ptMax:.1f} & {syst:.2f}")
