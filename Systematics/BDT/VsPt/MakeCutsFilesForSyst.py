'''
python script to create yaml files with set of cuts for cut-variation studies
'''

import os
import argparse
import copy
from itertools import product
import yaml

def get_variation_mult(edge, kind): # pylint: disable=too-many-return-statements
    if kind == 'loose_1':
        if edge == 'min':
            return -1.
        if edge == 'max':
            return 1.
    if kind == 'loose_2':
        if edge == 'min':
            return -2.
        if edge == 'max':
            return 2.
    if kind == 'tight_1':
        if edge == 'min':
            return 1.
        if edge == 'max':
            return -1.
    if kind == 'tight_2':
        if edge == 'min':
            return 2.
        if edge == 'max':
            return -2.
    return 0.

def check_value(new_value, cut_lim, histo_lim, edge):
    if edge == 'min':
        return cut_lim <= new_value < histo_lim
    if edge == 'max':
        return histo_lim < new_value <= cut_lim
    return False

def make_cuts():
    var_name = ['cospkphi3', 'cp', 'cpxy', 'd0', 'deltamKK', 'dl', 'dlxy', 'ndlxy', 'sigvtx', 'topo']
    var_key = ['CosPiKPhi3', 'CosP', 'CosPXY', 'd0', 'DeltaMassKK', 'DecL', 'DecLXY', 'NormDecLXY',
               'SigmaVtx', 'd0d0Exp']
    edge_to_vary = ['min', 'min', 'min', 'max', 'max', 'min', 'min', 'min', 'max', 'max']
    step_variation = [0.5, 0.5, 0.5, 10., 1., 10., 10., 1., 5., 0.5]
    histo_lims = [3.5, 100., 100., 0., 0., 105., 105., 10.5, 0., 0.]
    variation_kind = ['loose_1', 'loose_2', 'tight_1', 'tight_2']
    # [-1, -2, +1, +2]*step_variation added to central value

    in_dir = 'configfiles/'
    cut_file_central = 'cutset_3050_central_2018.yml'
    cut_file_loose = 'cutset_3050_loose_2018.yml'
    out_dir = 'configfiles/syst_cuts_Ds_3050/'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open(in_dir + cut_file_central, 'r') as cut_file_yml:
        cutset = yaml.load(cut_file_yml, yaml.FullLoader)
    with open(in_dir + cut_file_loose, 'r') as cut_file_loose_yml:
        cutset_loose = yaml.load(cut_file_loose_yml, yaml.FullLoader)

    cutset_central = copy.deepcopy(cutset)
    with open(out_dir + cut_file_central, 'w') as outfile:
        yaml.dump(cutset_central, outfile, default_flow_style=False)

    for name, key, edge, step, histo_lim in zip(var_name, var_key, edge_to_vary, step_variation, histo_lims):
        loose_values = cutset_loose['cutvars'][key][edge]
        for kind in variation_kind:
            cutset_mod = copy.deepcopy(cutset)
            mult_value = get_variation_mult(edge, kind)
            modified_list = []
            for value, cut_lim in zip(cutset_mod['cutvars'][key][edge], loose_values):
                new_value = value + step * mult_value
                if check_value(new_value, cut_lim, histo_lim, edge):
                    modified_list.append(new_value)
                else:
                    modified_list.append(value)
            cutset_mod['cutvars'][key][edge] = modified_list
            cut_file_mod = cut_file_central.replace('central_2018', name + '_' + kind)
            with open(out_dir + cut_file_mod, 'w') as outfile_mod:
                yaml.dump(cutset_mod, outfile_mod, default_flow_style=False)

def make_cuts_ml():
    var_key = ['ML_output_Bkg', 'ML_output_Prompt']
    var_tag = ['outBkg', 'outPrompt'] # used in file names to reduce length
    step_variation_pos = [
                    {"0.5": 0.005, "1.0": 0.02, "1.5": 0.2, "2.0": 0.2, "2.5": 0.2, "3.0": 0.2, "3.5": 0.2, \
                         "4.0": 0.2, "4.5": 0.2, "5.0": 0.2, "5.5": 0.1, "6.0": 0.1, "8.0": 0.1, "12.0": 0.1},
                    {"0.5": 0.15, "1.0": 0.15, "1.5": 0.15, "2.0": 0.15, "2.5": 0.15, "3.0": 0.15, "3.5": 0.15, \
                         "4.0": 0.15, "4.5": 0.15, "5.0": 0.15, "5.5": 0.15, "6.0": 0.15, "8.0": 0.15, "12.0": 0.15}     
                    ]
    step_variation_neg = [
                    {"0.5": 0.002, "1.0": 0.01, "1.5": 0.04, "2.0": 0.05, "2.5": 0.05, "3.0": 0.05, "3.5": 0.05, \
                         "4.0": 0.05, "4.5": 0.05, "5.0": 0.05, "5.5": 0.1, "6.0": 0.1, "8.0": 0.1, "12.0": 0.1},
                    {"0.5": 0.05, "1.0": 0.05, "1.5": 0.05, "2.0": 0.05, "2.5": 0.05, "3.0": 0.05, "3.5": 0.05, \
                         "4.0": 0.05, "4.5": 0.05, "5.0": 0.05, "5.5": 0.05, "6.0": 0.05, "8.0": 0.05, "12.0": 0.05}     
                    ]

    num_step_pos = [3,3]
    num_step_neg = [3,3]
    edge_to_vary = ['max', 'min']

    in_dir = '/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/configs/'
    cut_file_central = 'cutset_pp13TeV_binary.yml'
    out_dir = '/home/fchinu/Run3/Ds_pp_13TeV/Systematics/BDT/configs/'
    out_file_tag = 'cutset_ML_'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open(in_dir + cut_file_central, 'r') as cut_file_yml:
        cutset = yaml.load(cut_file_yml, yaml.FullLoader)

    #same number of steps for all variables
    neg_steps = [[-i for i in range(1, step_neg + 1)] for step_neg in num_step_neg]
    for idx, _ in enumerate(neg_steps):
        neg_steps[idx] = neg_steps[idx][::-1]
    pos_steps = [list(range(0, step_pos + 1)) for step_pos in num_step_pos]
    steps = [neg_step + pos_step for neg_step, pos_step in zip(neg_steps, pos_steps)]

    n_combinations = len(var_key)

    for prod_steps in product(*steps):
        print(prod_steps)
        cutset_mod = copy.deepcopy(cutset)
        file_tag = str()
        for i, step in enumerate(prod_steps):
            modified_list = []
            cuts = cutset_mod['cutvars'][var_key[i]]
            for min_val, max_val, pt_min in zip(cuts['min'], cuts['max'], cutset_mod['cutvars']['Pt']['min']):
                if edge_to_vary[i] == 'min':
                    if step < 0.:
                        new_value = min_val + step * step_variation_neg[i][f'{pt_min:.1f}']
                    else:
                        new_value = min_val + step * step_variation_pos[i][f'{pt_min:.1f}']
                    if(new_value < 0. or new_value >= max_val):
                        print("Warning: cut is negative or min value is greater then max value")
                        new_value = min_val
                else:
                    if step < 0.:
                        new_value = max_val + step * step_variation_neg[i][f'{pt_min:.1f}']
                    else:
                        new_value = max_val + step * step_variation_pos[i][f'{pt_min:.1f}']
                    if(new_value > 1. or new_value <= min_val):
                        print("Warning: cut is > 1 or min value is greater then max value")
                        new_value = max_val
                modified_list.append(new_value)
            cuts[edge_to_vary[i]] = modified_list
            step_name = 'pos'
            if step < 0.:
                step_name = 'neg'
                step += num_step_neg[i] + 1
            # more than 100 files unlikely
            name = f'_{var_tag[i]}_{step_name}{str(int(step)).zfill(2)}'
            file_tag += name
        cut_file_mod = f'{out_file_tag}{file_tag}.yml'
        with open(out_dir + cut_file_mod, 'w') as outfile_mod:
            yaml.dump(cutset_mod, outfile_mod, default_flow_style=False)

def main():
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument("--ml", help="make cuts for ml", action="store_true")
    args = parser.parse_args()
    if args.ml:
        make_cuts_ml()
    else:
        make_cuts()

main()
