import numpy as np
import matplotlib.pyplot as plt
import uproot
import pandas as pd
import yaml
import flarefly
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
configCutsName = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/cutset_pp13TeV_binary.yml"
projectionsName = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/Data/Projections_Data.root"
configFitName = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/config_Ds_Fit_pp13TeV.yml"
# Load the configs
with open(configCutsName, 'r') as stream:
    configCuts = yaml.safe_load(stream)
with open(configFitName, 'r') as stream:
    configFit = yaml.safe_load(stream)
# define PDFs lists
signal_pdfs = ["gaussian", "gaussian"]
background_pdfs = ["Pol2"]

# define the ids
DplusId = 0 # because signal_pdfs[gaussian_id] = "gaussian"
DsId = 1 # because signal_pdfs[gaussian_id] = "gaussian"
expo_id = 0     # because background_pdfs[expo_id] = "expo"
fitters = []
start = 6

for iPt, (ptMin, ptMax) in enumerate(zip(configFit["pp13TeVPrompt"]["PtMin"][5:6], configFit["pp13TeVPrompt"]["PtMax"][5:6])):

    dataHandler = DataHandler(data=projectionsName, var_name="fM", histoname=f'hMass_{ptMin*10:.0f}_{ptMax*10:.0f}', rebin = 4, limits=[1.71,2.1], use_zfit=True)
    # Fit the data
    fitters.append(F2MassFitter(data_handler=dataHandler,
                                name_signal_pdf=signal_pdfs,
                                name_background_pdf=background_pdfs,
                                name=f"{iPt+start}{background_pdfs[expo_id]}_{signal_pdfs[DplusId]}_{signal_pdfs[DsId]}"))
    fitters[iPt].set_particle_mass(DsId, mass = 1.968)
    fitters[iPt].set_particle_mass(DplusId, mass = 1.868)
    fitters[iPt].set_signal_initpar(DsId, "sigma", 0.01, limits=[0.001, 0.03])
    fitters[iPt].set_signal_initpar(DplusId, "sigma", 0.01, limits=[0.001, 0.03])

    fitters[iPt].set_background_initpar(expo_id, "a", -2, fix=False)
    fitters[iPt].set_background_initpar(expo_id, "b", 2, fix=False)
    fitters[iPt].set_background_initpar(expo_id, "c", -20, fix=False)

    fitters[iPt].mass_zfit()
    # Print the results
    print(f"Pt: {ptMin} - {ptMax}")
    print(f"Signal yields: {fitters[iPt].get_raw_yield(DplusId)}")
    print(f"Signal yields: {fitters[iPt].get_raw_yield(DsId)}")
    #print(f"Background yields: {fitters[iPt].backgroundYields[expo_id]}")
    #print(f"Signal significance: {fitters[iPt].signalSignificance[DplusId]}")
    #print(f"Signal significance: {fitters[iPt].signalSignificance[DsId]}")
    #print(f"Chi2: {fitters[iPt].chi2}")
    #print(f"Chi2/ndf: {fitters[iPt].chi2/fitters[iPt].ndf}")
fig = fitters[0].plot_mass_fit(style="ATLAS",
                    show_extra_info=False,
                    extra_info_loc=['upper left', 'center right'])
fig.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/test.png")

for massFitter in fitters:
    massFitter.dump_to_root("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/test.root", option='recreate')


