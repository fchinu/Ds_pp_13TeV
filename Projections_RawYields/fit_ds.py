"""
Test for binned fit with flarefly.F2MassFitter
"""

import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
# import sys
import argparse
# import numpy as np
import time
from particle import Particle
import zfit
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter

from concurrent.futures import ProcessPoolExecutor

def fit(input_file, input_file_bkgtempl, output_dir, pt_min, pt_max):
    """
    Method for fitting
    """
    

    start_data_init = time.time()

    data_hdl = DataHandler(data=input_file, var_name="fM",
                           histoname=f'hMass_{pt_min*10:.0f}_{pt_max*10:.0f}',
                           limits=[1.7,2.10], rebin=8)
    data_corr_bkg = DataHandler(data=input_file_bkgtempl, var_name="fM",
                                histoname=f'hDplusTemplate_{pt_min*10:.0f}_{pt_max*10:.0f}',
                                limits=[1.7,2.10], rebin=8)
    stop_data_init = time.time()

    start_fit_init = time.time()
    fitter = F2MassFitter(data_hdl, name_signal_pdf=["gaussian", "gaussian"],
                          name_background_pdf=["chebpol2", "hist"],
                          name=f"ds_pt{pt_min*10:.0f}_{pt_max*10:.0f}", chi2_loss=False,
                          verbosity=7, tol=1.e-1)

    # bkg initialisation
    fitter.set_background_initpar(0, "c0", 0.6)
    fitter.set_background_initpar(0, "c1", -0.2)
    fitter.set_background_initpar(0, "c2", 0.01)
    fitter.set_background_initpar(0, "frac", 0.7)
    fitter.set_background_template(1, data_corr_bkg)
    fitter.set_background_initpar(1, "frac", 0.1, limits=[0., 1.])

    # signals initialisation
    fitter.set_particle_mass(0, pdg_id=431, limits=[Particle.from_pdgid(431).mass*0.99e-3,
                                                    Particle.from_pdgid(431).mass*1.01e-3])
    fitter.set_signal_initpar(0, "sigma", 0.008, limits=[0.005, 0.035])
    fitter.set_signal_initpar(0, "frac", 0.1, limits=[0., 1.])
    fitter.set_particle_mass(1, pdg_id=411,
                             limits=[Particle.from_pdgid(411).mass*0.99e-3,
                                     Particle.from_pdgid(411).mass*1.01e-3])
    fitter.set_signal_initpar(1, "sigma", 0.008, limits=[0.005, 0.035])
    fitter.set_signal_initpar(1, "frac", 0.1, limits=[0., 1.])
    stop_fit_init = time.time()

    start_fit = time.time()
    fit_result =fitter.mass_zfit()
    stop_fit = time.time()

    stop_save, start_save = 0., 0.
    if fit_result.converged:
        start_save = time.time()
        loc = ["lower left", "upper left"]
        ax_title = r"$M(\mathrm{KK\pi})$ GeV$/c^2$"
        fig = fitter.plot_mass_fit(style="ATLAS", show_extra_info=True,
                                figsize=(8, 8), extra_info_loc=loc,
                                axis_title=ax_title)
        figres = fitter.plot_raw_residuals(figsize=(8, 8), style="ATLAS",
                                        extra_info_loc=loc, axis_title=ax_title)
        fig.savefig(f"{output_dir}/ds_mass_pt{pt_min:.1f}_{pt_max:.1f}.pdf")
        figres.savefig(f"{output_dir}/ds_massres_pt{pt_min:.1f}_{pt_max:.1f}.pdf")
        fitter.dump_to_root(os.path.join(output_dir, "ds_fit.root"), option="recreate",
                            suffix=f"_ds_pt{pt_min:.1f}_{pt_max:.1f}")
        stop_save = time.time()

    print("\n\nTime data handler initialisation: ", stop_data_init - start_data_init)
    print("Time fitter initialisation: ", stop_fit_init - start_fit_init)
    print("Time fit: ", stop_fit - start_fit)
    print("Time save outputs: ", stop_save - start_save)
    return {"rawyields": [fitter.get_raw_yield(i) for i in range(2)],
            "sigma": [fitter.get_sigma(i) for i in range(2)],
            "mean": [fitter.get_mass(i) for i in range(2)],
            "chi2": fitter.get_chi2_ndf(),
            "significance": [fitter.get_significance(i) for i in range(2)],
            "signal": [fitter.get_signal(i) for i in range(2)],
            "background": [fitter.get_background(i) for i in range(2)]}

if __name__ == "__main__":  
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_file", "-i", metavar="text", default="AnalysisResults.root",
                        help="input root file", required=True)
    parser.add_argument("--input_file_bkgtempl", "-b", metavar="text",
                        default="DplusForTemplateHistos_Train165702.root",
                        help="input root file with D+ template", required=False)
    parser.add_argument("--output_dir", "-o", metavar="text", default=".",
                        help="output directory", required=False)
    parser.add_argument("--pt_min", "-pmi", type=float, default=3.,
                        help="min pt", required=False)
    parser.add_argument("--pt_max", "-pma", type=float, default=3.5,
                        help="max pt", required=False)
    args = parser.parse_args()

    zfit.run.set_cpus_explicit(intra=30, inter=30)

    futures = []
    start_all = time.time()
    if False:
        with ProcessPoolExecutor(max_workers=6) as executor:
            futures.append(executor.submit(fit, args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min, args.pt_max))
            futures.append(executor.submit(fit, args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min+0.5, args.pt_max+0.5))
            futures.append(executor.submit(fit, args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min+1, args.pt_max+1))
            futures.append(executor.submit(fit, args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min+1.5, args.pt_max+1.5))
            #futures.append(executor.submit(fit, args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min+2, args.pt_max+2))
            #futures.append(executor.submit(fit, args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min+2.5, args.pt_max+2.5))
        print(time.time() - start_all)
        for fut in futures:
            print(fut.result()["sigma"])
    else:
        d = fit(args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min, args.pt_max)
        print(d["sigma"])
        print(d["chi2"])
        #fit(args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min+0.5, args.pt_max+0.5)
        #fit(args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min, args.pt_max)
        #fit(args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min+0.5, args.pt_max+0.5)
        #fit(args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min, args.pt_max)
        #fit(args.input_file, args.input_file_bkgtempl, args.output_dir, args.pt_min+0.5, args.pt_max+0.5)
