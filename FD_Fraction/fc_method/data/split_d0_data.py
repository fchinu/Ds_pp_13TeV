"""
Script to split D0 meson HEPData based on centrality classes.
"""

import ROOT

objects = {
    "010": [
        "Hist1D_y1",
        "Hist1D_y1_e1",
        "Hist1D_y1_e2",
        "Hist1D_y1_e3",
        "Graph1D_y1"
        ],
    "3050": [
        "Hist1D_y2",
        "Hist1D_y2_e1",
        "Hist1D_y2_e2",
        "Hist1D_y2_e3",
        "Graph1D_y2"
    ]
}


folder_name = "Table 1a"
with ROOT.TFile.Open("D0NonPromptYield_PbPb5TeV_Run2_010_3050.root") as infile:
    for cent_class, obj_names in objects.items():
        with ROOT.TFile(f"D0NonPromptYield_PbPb5TeV_Run2_{cent_class}.root", "RECREATE") as outfile:
            outfile.mkdir(folder_name)
            outfile.cd(folder_name)
            for name in obj_names:
                print(name)
                obj = infile.Get(folder_name + "/" + name)
                obj.Write()
