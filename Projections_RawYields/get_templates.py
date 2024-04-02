import pandas as pd
import ROOT


df = pd.read_parquet("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train165702/DplusForTemplate_Train165702.parquet")
infile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/Data/Projections_Data.root")
hist_templ = []
pt_mins, pt_maxs = [], []
for key in infile.GetListOfKeys():
    if "hMass" in key.GetName():
        hist_name = key.GetName()
        hist = infile.Get(hist_name)
        hist.SetDirectory(0)
        hist_templ.append(hist.Clone(hist_name.replace("hMass", "hDplusTemplate")))
        hist_templ[-1].SetDirectory(0)
        hist_templ[-1].Reset()
        pt_mins.append(float(hist_name.split(sep="_")[1])/10)
        pt_maxs.append(float(hist_name.split(sep="_")[2])/10)

for pt, mass in zip(df["fPt"].to_numpy(), df["fM"].to_numpy()):
    iptbin = -1
    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        if pt_min <= pt < pt_max:
            iptbin = ipt
            break
    if iptbin >= 0:
        hist_templ[iptbin].Fill(mass)

output = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/DplusForTemplateHistos_Train165702.root", "recreate")
for hist in hist_templ:
    hist.Smooth(1000)
    hist.Write()
output.Close()