import ROOT 
import numpy as np
import pandas as pd 
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from DfUtils import read_parquet_in_batches
from itertools import combinations
from array import array
import matplotlib as mpl

from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve, auc


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

def get_discrete_matplotlib_palette(paletteName):
    cmap = mpl.colormaps[paletteName]
    colors = cmap.colors
    ROOTColorIndices = []
    ROOTColors = []
    for color in colors:
        idx = ROOT.TColor.GetFreeColorIndex()
        ROOTColors.append(ROOT.TColor(idx, color[0], color[1], color[2],"color%i" % idx))
        ROOTColorIndices.append(idx)
        
    return ROOTColorIndices, ROOTColors

dfDsPrompt = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train189890/LHC22b1b_PromptDs_Train.parquet", "2<fPt<3")
dfDsFD = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train189892/LHC22b1a_NonPromptDs_Train.parquet", "2<fPt<3")
dfBkg = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train189072/LHC22o_pass6_small.parquet", "2<fPt<3 and (1.7 < fM < 1.75 or  2.1 < fM < 2.15)")
dfBkg = dfBkg.sample(frac=1, random_state=42)


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

# Get the ROCs
yPredTrainRaw = ModelHandl.predict(TrainSet, False)
yPredTestRaw = ModelHandl.predict(TestSet, False)

# Plot the ROCs
ROOT.gStyle.SetOptStat(0)
colors, ROOTcols = get_discrete_matplotlib_palette("tab10")

#ROOT.gStyle.SetPalette(ROOT.kRedBlue)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleSize(23,"padTitle")
ROOT.gStyle.SetTitleH(0.08)


# Create the canvas
c = ROOT.TCanvas("c", "c", 800, 800)
c.Divide(2,2)

# Prompt Ds
c.cd(1)

ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.13)
ROOT.gPad.SetBottomMargin(0.13)
ROOT.gPad.SetTopMargin(0.1)


hFrame = ROOT.gPad.DrawFrame(0, 1.e-4, 1, 5,"Background score distribution;Background score;Counts (a.u.)")
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetYaxis().SetLabelSize(0.04)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1.1)
hFrame.GetXaxis().SetLabelSize(0.04)

# Create the correlation matrix
# Create the histogram
hDistrBkgDsTrain = []
hDistrBkgDsBorders = []
hDistrBkgDsTest = []

for i, label in enumerate(['Bkg', 'Prompt D_{s}^{+}', 'Prompt D_{s}^{+}']):
    hDistrBkgDsTrain.append(ROOT.TH1D(f"hDistrBkgDs_{label}", "", 50, 0, 1))
    hDistrBkgDsTrain[i].SetFillColorAlpha(colors[i], 0.55)
    hDistrBkgDsTrain[i].SetLineWidth(0)
    for score in np.array(yPredTrainRaw)[yTrain == i]:
        hDistrBkgDsTrain[i].Fill(score[0])


    hDistrBkgDsTest.append(ROOT.TH1D(f"hDistrBkgDs_{label}", "", 50, 0, 1))
    hDistrBkgDsTest[i].SetLineColor(colors[i])
    hDistrBkgDsTest[i].SetLineWidth(2)
    hDistrBkgDsTest[i].SetMarkerColor(colors[i])
    hDistrBkgDsTest[i].SetMarkerStyle(ROOT.kFullCircle)
    hDistrBkgDsTest[i].SetMarkerSize(0.75)
    for score in np.array(yPredTestRaw)[yTest == i]:
        hDistrBkgDsTest[i].Fill(score[0])

for hist in hDistrBkgDsTrain:
    hist.DrawNormalized("hist same")

for hist in hDistrBkgDsTest:
    hist.DrawNormalized("hist,same")
    hist.DrawNormalized("pe same")

legBkg = ROOT.TLegend(0.35, 0.6, 0.6, 0.84)
legBkg.SetBorderSize(0)
legBkg.SetFillStyle(0)
legBkg.SetTextSize(0.03)
legBkg.AddEntry(hDistrBkgDsTrain[0], "Background, Training set", "f")
legBkg.AddEntry(hDistrBkgDsTrain[1], "Prompt D_{s}^{+}, Training set", "f")
legBkg.AddEntry(hDistrBkgDsTrain[2], "Non-prompt D_{s}^{+}, Training set", "f")
legBkg.AddEntry(hDistrBkgDsTest[0], "Background, Test set", "pl")
legBkg.AddEntry(hDistrBkgDsTest[1], "Prompt D_{s}^{+}, Test set", "pl")
legBkg.AddEntry(hDistrBkgDsTest[2], "Non-prompt D_{s}^{+}, Test set", "pl")
legBkg.Draw()


# Draw a grid
#ROOT.gPad.SetGrid()
ROOT.gPad.SetLogy()
ROOT.gPad.SetTickx()
ROOT.gPad.SetTicky()
ROOT.gPad.RedrawAxis()

# FD Ds
c.cd(2)
ROOT.gStyle.SetTitleSize(0.1,"padTitle")
ROOT.gStyle.SetTitleH(0.064)

ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.13)
ROOT.gPad.SetBottomMargin(0.13)
ROOT.gPad.SetTopMargin(0.1)


hFrame = ROOT.gPad.DrawFrame(0, 1.e-4, 1, 5,"Prompt D_{#lower[-0.3]{s}}^{+} score distribution;Prompt D_{#lower[-0.3]{s}}^{+} score;Counts (a.u.)")
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetYaxis().SetLabelSize(0.04)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1.1)
hFrame.GetXaxis().SetLabelSize(0.04)

# Create the correlation matrix
# Create the histogram
hDistrPromptDsTrain = []
hDistrPromptDsBorders = []
hDistrPromptDsTest = []

for i, label in enumerate(['Bkg', 'Prompt D_{s}^{+}', 'Prompt D_{s}^{+}']):
    hDistrPromptDsTrain.append(ROOT.TH1D(f"hDistrPromptDs_{label}", "", 50, 0, 1))
    hDistrPromptDsTrain[i].SetFillColorAlpha(colors[i], 0.55)
    hDistrPromptDsTrain[i].SetLineWidth(0)
    for score in np.array(yPredTrainRaw)[yTrain == i]:
        hDistrPromptDsTrain[i].Fill(score[1])


    hDistrPromptDsTest.append(ROOT.TH1D(f"hDistrPromptDs_{label}", "", 50, 0, 1))
    hDistrPromptDsTest[i].SetLineColor(colors[i])
    hDistrPromptDsTest[i].SetLineWidth(2)
    hDistrPromptDsTest[i].SetMarkerColor(colors[i])
    hDistrPromptDsTest[i].SetMarkerStyle(ROOT.kFullCircle)
    hDistrPromptDsTest[i].SetMarkerSize(0.75)
    for score in np.array(yPredTestRaw)[yTest == i]:
        hDistrPromptDsTest[i].Fill(score[1])

for hist in hDistrPromptDsTrain:
    hist.DrawNormalized("hist same")

for hist in hDistrPromptDsTest:
    hist.DrawNormalized("hist,same")
    hist.DrawNormalized("pe same")

legPrompt = ROOT.TLegend(0.35, 0.6, 0.6, 0.84)
legPrompt.SetBorderSize(0)
legPrompt.SetFillStyle(0)
legPrompt.SetTextSize(0.03)
legPrompt.AddEntry(hDistrPromptDsTrain[0], "Background, Training set", "f")
legPrompt.AddEntry(hDistrPromptDsTrain[1], "Prompt D_{s}^{+}, Training set", "f")
legPrompt.AddEntry(hDistrPromptDsTrain[2], "Non-prompt D_{s}^{+}, Training set", "f")
legPrompt.AddEntry(hDistrPromptDsTest[0], "Background, Test set", "pl")
legPrompt.AddEntry(hDistrPromptDsTest[1], "Prompt D_{s}^{+}, Test set", "pl")
legPrompt.AddEntry(hDistrPromptDsTest[2], "Non-prompt D_{s}^{+}, Test set", "pl")
legPrompt.Draw()


# Draw a grid
#ROOT.gPad.SetGrid()
ROOT.gPad.SetLogy()
ROOT.gPad.SetTickx()
ROOT.gPad.SetTicky()
ROOT.gPad.RedrawAxis()

# Bkg
c.cd(3)
ROOT.gStyle.SetTitleSize(0.1,"padTitle")
ROOT.gStyle.SetTitleH(0.055)

ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.13)
ROOT.gPad.SetBottomMargin(0.13)
ROOT.gPad.SetTopMargin(0.1)


hFrame = ROOT.gPad.DrawFrame(0, 1.e-4, 1, 5,"Non-prompt D_{#lower[-0.3]{s}}^{+} score distribution;Non-prompt D_{#lower[-0.3]{s}}^{+} score;Counts (a.u.)")
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetYaxis().SetLabelSize(0.04)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1.1)
hFrame.GetXaxis().SetLabelSize(0.04)

# Create the correlation matrix
# Create the histogram
hDistrFDDsTrain = []
hDistrFDDsBorders = []
hDistrFDDsTest = []

for i, label in enumerate(['Bkg', 'Prompt D_{s}^{+}', 'FD D_{s}^{+}']):
    hDistrFDDsTrain.append(ROOT.TH1D(f"hDistrFDDs_{label}", "", 50, 0, 1))
    hDistrFDDsTrain[i].SetFillColorAlpha(colors[i], 0.55)
    hDistrFDDsTrain[i].SetLineWidth(0)
    for score in np.array(yPredTrainRaw)[yTrain == i]:
        hDistrFDDsTrain[i].Fill(score[2])


    hDistrFDDsTest.append(ROOT.TH1D(f"hDistrFDDs_{label}", "", 50, 0, 1))
    hDistrFDDsTest[i].SetLineColor(colors[i])
    hDistrFDDsTest[i].SetLineWidth(2)
    hDistrFDDsTest[i].SetMarkerColor(colors[i])
    hDistrFDDsTest[i].SetMarkerStyle(ROOT.kFullCircle)
    hDistrFDDsTest[i].SetMarkerSize(0.75)
    for score in np.array(yPredTestRaw)[yTest == i]:
        hDistrFDDsTest[i].Fill(score[2])

for hist in hDistrFDDsTrain:
    hist.DrawNormalized("hist same")

for hist in hDistrFDDsTest:
    hist.DrawNormalized("hist,same")
    hist.DrawNormalized("pe same")

legFD = ROOT.TLegend(0.35, 0.6, 0.6, 0.84)
legFD.SetBorderSize(0)
legFD.SetFillStyle(0)
legFD.SetTextSize(0.03)
legFD.AddEntry(hDistrFDDsTrain[0], "Background, Training set", "f")
legFD.AddEntry(hDistrFDDsTrain[1], "Prompt D_{s}^{+}, Training set", "f")
legFD.AddEntry(hDistrFDDsTrain[2], "Non-prompt D_{s}^{+}, Training set", "f")
legFD.AddEntry(hDistrFDDsTest[0], "Background, Test set", "pl")
legFD.AddEntry(hDistrFDDsTest[1], "Prompt D_{s}^{+}, Test set", "pl")
legFD.AddEntry(hDistrFDDsTest[2], "Non-prompt D_{s}^{+}, Test set", "pl")
legFD.Draw()


# Draw a grid
#ROOT.gPad.SetGrid()
ROOT.gPad.SetLogy()
ROOT.gPad.SetTickx()
ROOT.gPad.SetTicky()
ROOT.gPad.RedrawAxis()


# Thesis text and collision information
c.cd(4)
thesisText = ROOT.TLatex(0.25, 0.7, 'This Thesis')
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.06)
thesisText.Draw("same")

collisionText = ROOT.TLatex(0.25, 0.64, 'pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.05]{V}')
collisionText.SetNDC()
collisionText.SetTextFont(42)
collisionText.SetTextSize(0.05)
collisionText.Draw("same")

DecayText = ROOT.TLatex(0.25, 0.55, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.25, 0.49, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ptText = ROOT.TLatex(0.25, 0.43, '2 < #it{p}_{T} < 3 Ge#kern[-0.05]{V}/#it{c}')
ptText.SetNDC()
ptText.SetTextFont(42)
ptText.SetTextSize(0.05)
ptText.Draw("same")

outfile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Figures/ML_output_distribution/Distributions.root", "RECREATE")
#hCorrPromptDs.Write()
#hCorrFDDs.Write()
#hCorrBkg.Write()
c.Write()
outfile.Close()

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/ML_output_distribution/Distributions.pdf")
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/ML_output_distribution/Distributions.png")



