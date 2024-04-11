"""
Test for binned fit with flarefly.F2MassFitter
"""

import os
# import sys
import argparse
# import numpy as np
import time
from particle import Particle
import zfit
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter


def fit(input_file, input_file_bkgtempl, output_dir, pt_min, pt_max):
    """
    Method for fitting
    """
    zfit.run.set_cpus_explicit(intra=30, inter=30)

    start_data_init = time.time()

    data_hdl = DataHandler(data=input_file, var_name="fM",
                           histoname=f'hMass_{pt_min*10:.0f}_{pt_max*10:.0f}',
                           limits=[1.7,2.10], rebin=1)
    data_corr_bkg = DataHandler(data=input_file_bkgtempl, var_name="fM",
                                histoname=f'hDplusTemplate_{pt_min*10:.0f}_{pt_max*10:.0f}',
                                limits=[1.7,2.10], rebin=1)
    stop_data_init = time.time()

    start_fit_init = time.time()
    fitter = F2MassFitter(data_hdl, name_signal_pdf=["gaussian", "gaussian"],
                          name_background_pdf=["expo", "hist"],
                          name=f"ds_pt{pt_min*10:.0f}_{pt_max*10:.0f}", chi2_loss=False,
                          verbosity=7, tol=1.e-1)

    # bkg initialisation
    fitter.set_background_initpar(0, "lam", -2)
    fitter.set_background_initpar(0, "frac", 0.7)
    fitter.set_background_template(1, data_corr_bkg)


    # signals initialisation
    fitter.set_particle_mass(0, pdg_id=431, limits=[Particle.from_pdgid(431).mass*0.99e-3,
                                                    Particle.from_pdgid(431).mass*1.01e-3])
    fitter.set_signal_initpar(0, "sigma", 0.008, limits=[0.005, 0.030])
    fitter.set_signal_initpar(0, "frac", 0.1, limits=[0., 1.])
    fitter.set_particle_mass(1, pdg_id=411,
                             limits=[Particle.from_pdgid(411).mass*0.99e-3,
                                     Particle.from_pdgid(411).mass*1.01e-3])
    fitter.set_signal_initpar(1, "sigma", 0.008, limits=[0.005, 0.030])
    fitter.set_signal_initpar(1, "frac", 0.1, limits=[0., 1.])
    stop_fit_init = time.time()

    start_fit = time.time()
    fit_result =fitter.mass_zfit()
    stop_fit = time.time()

    stop_save, start_save = 0., 0.
    stop_getrawy, start_getrawy = 0., 0.
    rawy_fit, rawy_binc = 0., 0.
    if fit_result.converged:
        start_save = time.time()
        loc = ["lower left", "upper left"]
        ax_title = r"$M(\mathrm{KK\pi})$ GeV$/c^2$"
        fig = fitter.plot_mass_fit(style="ATLAS", show_extra_info=True,
                                figsize=(8, 8), extra_info_loc=loc,
                                axis_title=ax_title)
        figres = fitter.plot_raw_residuals(figsize=(8, 8), style="ATLAS",
                                        extra_info_loc=loc, axis_title=ax_title)
        fig.savefig(f"{output_dir}/ds_mass_pt{pt_min:.0f}_{pt_max:.0f}.pdf")
        figres.savefig(f"{output_dir}/ds_massres_pt{pt_min:.0f}_{pt_max:.0f}.pdf")
        fitter.dump_to_root(os.path.join(output_dir, "ds_fit.root"), option="recreate",
                            suffix=f"_ds_pt{pt_min:.0f}_{pt_max:.0f}")
        stop_save = time.time()

        start_getrawy = time.time()
        rawy_fit = fitter.get_raw_yield(0)
        rawy_binc = fitter.get_raw_yield_bincounting(0, nsigma=3)
        stop_getrawy = time.time()

    print("\n\nTime data handler initialisation: ", stop_data_init - start_data_init)
    print("Time fitter initialisation: ", stop_fit_init - start_fit_init)
    print("Time fit: ", stop_fit - start_fit)
    print("Time save outputs: ", stop_save - start_save)
    print("Time getting raw yields: ", stop_getrawy - start_getrawy)
    print(f"\n\nDs raw-yield fit = {rawy_fit}, Ds raw-yield bin counting = {rawy_binc}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_file", "-i", metavar="text", default="/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/Data/Projections_Data.root",
                        help="input root file", required=True)
    parser.add_argument("--input_file_bkgtempl", "-b", metavar="text",
                        default="/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/DplusForTemplateHistos_Train165702.root",
                        help="input root file with D+ template", required=False)
    parser.add_argument("--output_dir", "-o", metavar="text", default=".",
                        help="output directory", required=False)
    parser.add_argument("--pt_min", "-pmi", type=float, default=2.,
                        help="min pt", required=False)
    parser.add_argument("--pt_max", "-pma", type=float, default=2.5,
                        help="max pt", required=False)
    args = parser.parse_args()

    fit(args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min, args.pt_max)
