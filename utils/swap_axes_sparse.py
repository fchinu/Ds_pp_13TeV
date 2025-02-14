import ROOT
import numpy as np

def reorder_axes(hist, axis_map):
    """ Reorder the axes of a THnSparse histogram based on axis_map """
    dim = hist.GetNdimensions()
    if sorted(axis_map.keys()) != list(range(dim)) or sorted(axis_map.values()) != list(range(dim)):
        raise ValueError(f"Invalid axis map for {hist.GetName()}")

    # Extract bin settings
    nbins, xmin, xmax = [0]*dim, [0]*dim, [0]*dim
    titles = [""] * dim
    for i in range(dim):
        nbins[axis_map[i]] = hist.GetAxis(i).GetNbins()
        xmin[axis_map[i]] = hist.GetAxis(i).GetXmin()
        xmax[axis_map[i]] = hist.GetAxis(i).GetXmax()
        titles[axis_map[i]] = hist.GetAxis(i).GetTitle()

    # Create new THnSparse histogram with reordered axes
    new_hist = ROOT.THnSparseF(hist.GetName(), hist.GetName(), dim, np.asarray(nbins, dtype=np.int32), np.asarray(xmin, dtype=np.float64), np.asarray(xmax, dtype=np.float64))
    ROOT.SetOwnership(new_hist, False)  # Avoid memory issues

    # Buffer for bin coordinates
    coord = np.zeros(dim, dtype=np.int32)

    # Iterate over all bins
    for i in range(hist.GetNbins()):
        bin_content = hist.GetBinContent(i, coord)
        
        # Get coordinates indices correctly
        reordered_coord = np.zeros(dim, dtype=np.float64)
        for i_dim in range(dim):
            reordered_coord[axis_map[i_dim]] = hist.GetAxis(i_dim).GetBinCenter(int(coord[i_dim]))

        new_hist.Fill(reordered_coord, bin_content)

    return new_hist


def process_directory(directory, output_directory, path=""):
    """ Recursively process all objects in a directory """
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        name = obj.GetName()
        full_path = f"{path}/{name}" if path else name
        
        if isinstance(obj, ROOT.TDirectory):
            # Create a corresponding directory in the output file
            new_dir = output_directory.mkdir(name)
            new_dir.cd()
            process_directory(obj, new_dir, full_path)  # Recursive call
            output_directory.cd()
        elif isinstance(obj, ROOT.THnSparse) and full_path in axis_reorder_config:
            reordered_hist = reorder_axes(obj, axis_reorder_config[full_path])
            reordered_hist.Write()
        else:
            obj.Write()

# Define axis reordering for specific histograms
axis_reorder_config = {
    "hf-task-ds/MC/Ds/Prompt/hSparseMass": {
        0: 0, #mass
        1: 1, #pt
        2: 2, #centrality
        3: 6, #npvcontributions
        4: 3, #bkg score
        5: 4, #prompt score
        6: 5  #fd score
    },
    "hf-task-ds/MC/Dplus/Prompt/hSparseMass": {
        0: 0, #mass
        1: 1, #pt
        2: 2, #centrality
        3: 6, #npvcontributions
        4: 3, #bkg score
        5: 4, #prompt score
        6: 5  #fd score
    },
    "hf-task-ds/MC/Dplus/NonPrompt/hSparseMass": {
        0: 0, #mass
        1: 1, #pt
        2: 2, #centrality
        3: 6, #npvcontributions
        4: 3, #bkg score
        5: 4, #prompt score
        6: 5  #fd score
    },
    "hf-task-ds/MC/Dplus/Bkg/hSparseMass": {
        0: 0, #mass
        1: 1, #pt
        2: 2, #centrality
        3: 6, #npvcontributions
        4: 3, #bkg score
        5: 4, #prompt score
        6: 5  #fd score
    },
    "hf-task-ds/MC/Ds/NonPrompt/hPtYNPvContribGen": {
        0: 0, #pt
        1: 1, #y
        2: 2, #npvcontributions
        3: 5, #ptB
        4: 4, #flagB
        5: 3  #centrality
    },
    "hf-task-ds/MC/Dplus/NonPrompt/hPtYNPvContribGen": {
        0: 0, #pt
        1: 1, #y
        2: 2, #npvcontributions
        3: 5, #ptB
        4: 4, #flagB
        5: 3  #centrality
    },
    # Add more histograms and their respective reordering here
}

# Load input file
input_file = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train347749/AnalysisResults_LHC24h1_0.root", "READ")
output_file = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train347749/AnalysisResults_LHC24h1.root", "RECREATE")

process_directory(input_file, output_file)

output_file.Close()
input_file.Close()
