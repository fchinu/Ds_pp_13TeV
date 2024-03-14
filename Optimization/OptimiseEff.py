'''
python script for the optimisation of the working point using TTrees or pandas dataframes as input
run: python ScanSelectionsTree.py cfgFileName.yml outFileName.root
'''

import sys
import argparse
import time
import itertools
import pandas as pd
import numpy as np
import yaml
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter


def ProduceQuery(vars, cuts, loweruppers):
    '''
    Produce the query string for the pandas dataframe
    '''
    query = ""
    for var, cut, lowerupper in zip(vars, cuts, loweruppers):
        if lowerupper == "upper":
            query += f"{var} < {cut} and "
        elif lowerupper == "lower":
            query += f"{var} > {cut} and "
        else:
            exit(f"Error: {lowerupper} is not a valid option for lowerupper, use upper or lower")
    query = query[:-5]  #Remove the last "and"
    return query


def ComputeEfficiency(df, vars, cuts, loweruppers):
    '''
    Compute the efficiency of the cuts
    '''
    eff = 1
    query_str = ProduceQuery(vars, cuts, loweruppers)
    eff = len(df.query(query_str))/len(df)
    return eff

def ComputeSignificance(df, vars, cuts, loweruppers):
    # define PDFs lists
    signal_pdfs = ["gaussian", "gaussian"]
    background_pdfs = ["powlaw"]
    DsId = 0
    DplusId = 1
    bkgId = 0
    query_str = ProduceQuery(vars, cuts, loweruppers)
    #for idx, (var, cut, lowerupper) in enumerate(zip(vars, cuts, loweruppers)):
    dfsel = df.query(query_str)
    data = DataHandler(data=dfsel,
                var_name="fM")
    fitter = F2MassFitter(data_handler=data,
                    name_signal_pdf=signal_pdfs,
                    name_background_pdf=background_pdfs,
                    name=f"Fitter_{0}")

    fitter.set_particle_mass(DsId, mass = 1.968)
    fitter.set_particle_mass(DplusId, mass = 1.868)
    fitter.set_signal_initpar(DsId, "sigma", 0.1, limits=[0.05, 0.3])
    fitter.set_signal_initpar(DplusId, "sigma", 0.1, limits=[0.05, 0.3])
    fitter.set_background_initpar(bkgId, "lam", -1.5)
    fitter.mass_zfit()
    plot_mass_fit = fitter.plot_mass_fit(style="ATLAS",
                        show_extra_info=True,
                        extra_info_loc=['upper left', 'center right'])

    plot_mass_fit.savefig(f"/home/fchinu/Run3/Ds_pp_13TeV/Optimization/mass_fit_{0}.png")

if __name__ == '__main__':
    #parser = argparse.ArgumentParser(description='Script for the optimisation of the working point using TTrees or pandas dataframes as input')
    #parser.add_argument('cfg', metavar='cfg', type=str, help='YAML config file')
    #parser.add_argument('out', metavar='out', type=str, help='output file name')
    #args = parser.parse_args()
#
    #with open(args.cfg, 'r') as yaml_file:
    #    config = yaml.safe_load(yaml_file)
#
    ## Load the data
    #if config["input"]["type"] == "parquet":
    #    df = pd.read_parquet(config["input"]["file"])
    #else:
    #    exit(f"Error: {config['input']['type']} is not a valid option for input type, use parquet")
    #df.query(f"fPt > {config['ptmin']} and fPt < {config['ptmax']}", inplace=True)
#
    ## Compute the efficiency
    #eff = ComputeEfficiency(df, config["vars"], config["cuts"], config["loweruppers"])
    #print(f"Efficiency: {eff}")
#
    ## Compute the significance
    #ComputeSignificance(df, config["vars"], config["cuts"], config["loweruppers"])
#
    ## Save the dataframe
    #df.to_parquet(args.out)
    #print(f"Dataframe saved in {args.out}")
    #del df
    #print("Done")

    df = pd.read_parquet("/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt4_4.5/Data_pT_4_4.5_ModelApplied.parquet.gzip")
    ComputeSignificance(df, ["ML_output_Bkg", "ML_output_Prompt"], [0.3, 0.4], ["upper", "lower"])