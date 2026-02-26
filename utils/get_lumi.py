from pathlib import Path
import uproot


indir = Path("/home/fchinu/Run3/Ds_Dp_ratio_PbPb/Datasets/Data/Train530885")
luminometer = "TCE"

files = list(indir.glob("*.root"))

total_lumi = 0
for file in files:
    with uproot.open(file) as f:
        lumi = f["eventselection-run3"]["luminosity"][f"hLumi{luminometer}afterBCcuts"]
        lumi_sum = lumi.values().sum()
        collisions = f["hf-candidate-creator-3prong"]["hCollisions"]
        z_vtx_eff = collisions.values()[-1] / collisions.values()[-2]  #efficiency of the z vertex cut
        lumi_sum *= z_vtx_eff
        total_lumi += lumi_sum

print(f"Total luminosity for {luminometer}: {total_lumi:.2f} ub^-1")
