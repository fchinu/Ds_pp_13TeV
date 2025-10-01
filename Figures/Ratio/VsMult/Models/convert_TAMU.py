import pandas as pd
import ROOT

FileName="HeRapp_DsD0_PDG_LHM.dat"

dfpp = pd.read_csv(FileName, sep=' ', header=None, names=['pT', 'DsOverD0LowMult', 'DsOverD0HighMult'])

DpOverD0 = 0.4391 # from https://arxiv.org/abs/1902.08889 (RQM, 170 MeV), assuming independent of everything (pt, mult, etc)

gTAMU_vspt_lowmult = ROOT.TGraphErrors(len(dfpp))
gTAMU_vspt_highmult = ROOT.TGraphErrors(len(dfpp))
for ipt, (pt, ratio_lm, ratio_hm) in enumerate(
    zip(dfpp["pT"].to_numpy(), dfpp["DsOverD0LowMult"].to_numpy(), dfpp["DsOverD0HighMult"].to_numpy())):
    gTAMU_vspt_lowmult.SetPoint(ipt, pt, ratio_lm/DpOverD0)
    gTAMU_vspt_highmult.SetPoint(ipt, pt, ratio_hm/DpOverD0)
    gTAMU_vspt_lowmult.SetPointError(ipt, 0., 0.)
    gTAMU_vspt_highmult.SetPointError(ipt, 0., 0.)

pt_mins = [1., 2., 4., 6., 8.]
pt_maxs = [2., 4., 6., 8., 12.]
mult_values = [3.1, 37.8]
gTAMU_pp = []
for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
    gTAMU_pp.append(ROOT.TGraphErrors(2))
    gTAMU_pp[ipt].SetName(f"graph_ds_over_dp_pdg_T170_pt{pt_min:.1f}_{pt_max:.1f}")
    pt_cent = (pt_max + pt_min) / 2
    gTAMU_pp[ipt].SetPoint(0, mult_values[0], gTAMU_vspt_lowmult.Eval(pt_cent))
    gTAMU_pp[ipt].SetPoint(1, mult_values[1], gTAMU_vspt_highmult.Eval(pt_cent))

outFile = ROOT.TFile("TAMU_pp_13TeV_vsmult_DsOverDp.root", "recreate")
for graph in gTAMU_pp:
    graph.Write()
outFile.Close()
