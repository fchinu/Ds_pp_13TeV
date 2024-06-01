import ROOT 
import numpy as np
import pandas as pd 
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from DfUtils import read_parquet_in_batches

from sklearn.metrics import confusion_matrix


def set_matplotlib_palette(paletteName):
    n_colors = 255
    cmap = mpl.colormaps[paletteName]
    colors = cmap(np.linspace(0, 1, n_colors))
    stops = array('d', np.linspace(0, 1, n_colors))
    red = array('d', colors[:, 0])
    green = array('d', colors[:, 1])
    blue = array('d', colors[:, 2])
    ROOT.TColor.CreateGradientColorTable(n_colors, stops, red, green, blue, 255)
    ROOT.gStyle.SetNumberContours(255)

def get_label(x):
    if x <0.5:
        return 0
    elif x < 1.5:
        return 1
    else:
        return 2

dfDsPrompt = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train213175/LHC24d3a_PromptDs.parquet", "2<fPt<3")
dfDsFD = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train213175/LHC24d3a_NonPromptDs.parquet", "2<fPt<3")
dfBkg = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train191566/LHC22o.parquet", "2<fPt<3 and (1.7 < fM < 1.75 or  2.1 < fM < 2.15)")


nPrompt = len(dfDsPrompt)
nFD = len(dfDsFD)
nBkg = len(dfBkg)
nCandToKeep = min([nPrompt, nFD, nBkg])

TotDf = pd.concat([dfBkg.iloc[:nCandToKeep], dfDsPrompt.iloc[:nCandToKeep], dfDsFD.iloc[:nCandToKeep]], sort=True)
LabelsArray = np.array([0]*nCandToKeep + [1]*nCandToKeep + [2]*nCandToKeep)

TrainSet, TestSet, yTrain, yTest = train_test_split(TotDf, LabelsArray, test_size=0.2, random_state=42)

# Load the model 
ModelHandl = ModelHandler()
ModelHandl.load_model_handler("/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt2_3/ModelHandler_pT_2_3.pickle")

# Get the confusion matrix
yPredRaw = ModelHandl.predict(TestSet)
yPredLabels = np.argmax(yPredRaw, axis=1)
cm = confusion_matrix(yTest, yPredLabels)
cm = np.flip(cm, axis=1)

# Plot the confusion matrix
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleH(0.08)

c = ROOT.TCanvas("c", "c", 800, 600)
c.SetLeftMargin(0.2)
c.SetRightMargin(0.05)
c.SetBottomMargin(0.05)
c.SetTopMargin(0.15)

hCm = ROOT.TH2F("hCm", ";True labels;Predicted labels", 3, 0, 3, 3, 0, 3)
for i in range(3):
    for j in range(3):
        hCm.SetBinContent(i+1, j+1, cm[i, j])

hCm.GetXaxis().SetBinLabel(1, "Bkg.")
hCm.GetXaxis().SetBinLabel(2, "Prompt D_{#lower[-0.4]{s}}^{+}")
hCm.GetXaxis().SetBinLabel(3, "FD D_{s}^{+}")
hCm.GetYaxis().SetBinLabel(3, "Bkg.")
hCm.GetYaxis().SetBinLabel(2, "Prompt D_{#lower[-0.4]{s}}^{+}")
hCm.GetYaxis().SetBinLabel(1, "FD D_{s}^{+}")

hCm.GetXaxis().SetTitleOffset(1.)
hCm.GetYaxis().SetTitleOffset(1.7)
hCm.GetXaxis().SetLabelOffset(0.005)
hCm.GetXaxis().CenterTitle()
hCm.GetYaxis().CenterTitle()
hCm.GetXaxis().SetTitleSize(0.06)
hCm.GetYaxis().SetTitleSize(0.06)
hCm.GetXaxis().SetLabelSize(0.06)
hCm.GetYaxis().SetLabelSize(0.06)

hCm.Draw("X+")

boxes = []
for i in range(3):
    for j in range(3):
        boxes.append(ROOT.TPaveText(hCm.GetXaxis().GetBinLowEdge(i+1), hCm.GetYaxis().GetBinLowEdge(j+1), hCm.GetXaxis().GetBinUpEdge(i+1), hCm.GetYaxis().GetBinUpEdge(j+1)))
        #boxes.append(ROOT.TPaveText(0,0,1,1))
        boxes[-1].AddText(str(cm[i, j]))
        if (2-i)==j:
            boxes[-1].SetFillColor(ROOT.kTeal+4)
        else:
            boxes[-1].SetFillColor(ROOT.kRed-4)
        boxes[-1].SetTextSize(0.05)
        boxes[-1].SetTextAlign(22)
        boxes[-1].SetBorderSize(1)
        boxes[-1].SetLineWidth(2)
        boxes[-1].Draw("same")

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/confusion_matrix/ConfusionMatrix.pdf")

