import json

cents = [0, 1, 10, 30, 50, 70, 100]

with open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/mult_differential/2023/ds/doublecb/Config_ds.json") as f:
    data = json.load(f)

data_new = data.copy()
histo_name = data["rawyields"]["histoname"]
output_file = data["output"]["file"]
histo_eff_names = data["efficiencies"]["histonames"]

for cent_min, cent_max in zip(cents[:-1], cents[1:]):
    data_new["rawyields"]["histoname"] = histo_name.replace("0_100", f"{cent_min}_{cent_max}")
    data_new["output"]["file"] = output_file.replace("MB", f"{cent_min}_{cent_max}")
    data_new["efficiencies"]["histonames"] = {name: histo_eff_names[name].replace("0_100", f"{cent_min}_{cent_max}") for name in histo_eff_names}
    print(data["rawyields"]["histoname"])
    print(data_new["rawyields"]["histoname"])
    with open(f"/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/mult_differential/2023/ds/doublecb/Config_ds_{cent_min}_{cent_max}.json", "w") as f:
        json.dump(data_new, f)

cents = [0, 1, 10, 30, 50, 70, 100]

with open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/mult_differential/2023/dplus/doublecb/Config_ds.json") as f:
    data = json.load(f)

data_new = data.copy()
histo_name = data["rawyields"]["histoname"]
output_file = data["output"]["file"]
histo_eff_names = data["efficiencies"]["histonames"]

for cent_min, cent_max in zip(cents[:-1], cents[1:]):
    data_new["rawyields"]["histoname"] = histo_name.replace("0_100", f"{cent_min}_{cent_max}")
    data_new["output"]["file"] = output_file.replace("MB", f"{cent_min}_{cent_max}")
    data_new["efficiencies"]["histonames"] = {name: histo_eff_names[name].replace("0_100", f"{cent_min}_{cent_max}") for name in histo_eff_names}
    print(data["rawyields"]["histoname"])
    print(data_new["rawyields"]["histoname"])
    with open(f"/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/VsMult/mult_differential/2023/dplus/doublecb/Config_ds_{cent_min}_{cent_max}.json", "w") as f:
        json.dump(data_new, f)