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

def get_roc_ovo(y_truth, y_score, n_classes, labels, average):
    """
    Utility function for plot_roc in the multi-class case. Calculate and plot the
    ROC curves with the one-vs-one approach
    """
    outDict = {}
    for i_comb, (aaa, bbb) in enumerate(combinations(range(n_classes), 2)):
        a_mask = y_truth == aaa
        b_mask = y_truth == bbb
        ab_mask = np.logical_or(a_mask, b_mask)
        a_true = a_mask[ab_mask]
        b_true = b_mask[ab_mask]
        fpr_a, tpr_a, _ = roc_curve(a_true, y_score[ab_mask, aaa])
        roc_auc_a = auc(fpr_a, tpr_a)
        outDict[f'{labels[aaa]}__{labels[bbb]}'] = (fpr_a, tpr_a, roc_auc_a)
        fpr_b, tpr_b, _ = roc_curve(b_true, y_score[ab_mask, bbb])
        roc_auc_b = auc(fpr_b, tpr_b)
        outDict[f'{labels[bbb]}__{labels[aaa]}'] = (fpr_b, tpr_b, roc_auc_b)
    global_roc_auc = roc_auc_score(y_truth, y_score, average=average, multi_class='ovo')
    return outDict, global_roc_auc



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

dictTrain_roc, globalTrain_auc = get_roc_ovo(yTrain, yPredTrainRaw, 3, ['Bkg', 'Prompt D_{#lower[-0.5]{s}}^{+}', 'FD D_{#lower[-0.3]{s}}^{+}'], 'macro')
dictTest_roc, globalTest_auc = get_roc_ovo(yTest, yPredTestRaw, 3, ['Bkg', 'Prompt D_{#lower[-0.5]{s}}^{+}', 'FD D_{#lower[-0.3]{s}}^{+}'], 'macro')


# Plot the ROCs
ROOT.gStyle.SetOptStat(0)
colors, ROOTcols = get_discrete_matplotlib_palette("tab10")
c = ROOT.TCanvas("c", "c", 800, 800)
c.SetLeftMargin(0.1)
c.SetRightMargin(0.05)
c.SetBottomMargin(0.1)
c.SetTopMargin(0.05)

hFrame = c.DrawFrame(0, 0, 1, 1.2, ";False Positive Rate;True Positive Rate")
hFrame.GetYaxis().SetNdivisions(-512)
hFrame.GetXaxis().SetTickLength(0.025)

legend = ROOT.TLegend(0.37, 0.13, 0.67, 0.33)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.028)

graphsTrain = []
graphsTest = []
dictTest_roc = {'Bkg__Prompt D_{#lower[-0.5]{s}}^{+}': dictTest_roc['Bkg__Prompt D_{#lower[-0.5]{s}}^{+}'], 
                 'Bkg__FD D_{#lower[-0.3]{s}}^{+}': dictTest_roc['Bkg__FD D_{#lower[-0.3]{s}}^{+}'], 
                 'Prompt D_{#lower[-0.5]{s}}^{+}__FD D_{#lower[-0.3]{s}}^{+}': dictTest_roc['Prompt D_{#lower[-0.5]{s}}^{+}__FD D_{#lower[-0.3]{s}}^{+}']}

for i, comb in enumerate(dictTest_roc.keys()):
    #fpr, tpr, rocauc = dictTrain_roc[comb]
    #graphsTrain.append(ROOT.TGraph(len(fpr), array('d', fpr), array('d', tpr)))
    #graphsTrain[-1].SetTitle()
    #graphsTrain[-1].SetLineColor(colors[i])
    #graphsTrain[-1].SetLineWidth(2)
    #graphsTrain[-1].SetLineStyle(2)
    #graphsTrain[-1].Draw("L,same")
    #legend.AddEntry(graphsTrain[-1], f"Train: {comb.replace('__', ' vs ')} (ROC AUC = {rocauc:.3f})", "l")
    fpr, tpr, rocauc = dictTest_roc[comb]
    graphsTest.append(ROOT.TGraph(len(fpr), array('d', fpr), array('d', tpr)))
    graphsTest[-1].SetTitle()
    graphsTest[-1].SetLineColor(colors[i])
    graphsTest[-1].SetLineWidth(2)
    graphsTest[-1].Draw("L,same")
    legend.AddEntry(graphsTest[-1], f"{comb.replace('__', ' vs ')} (ROC AUC = {rocauc:.3f})", "l")

RandomClassifier = ROOT.TGraph(2, array('d', [0, 1]), array('d', [0, 1]))
RandomClassifier.SetLineColor(ROOT.kGray)
RandomClassifier.SetLineStyle(9)
RandomClassifier.SetLineWidth(2)
RandomClassifier.Draw("L,same")
legend.AddEntry(RandomClassifier, "Luck (ROC AUC = 0.5)", "l")

LineAtOne = ROOT.TLine(0, 1, 1, 1)
LineAtOne.SetLineColor(ROOT.kBlack)
LineAtOne.SetLineStyle(9)
LineAtOne.SetLineWidth(2)
LineAtOne.Draw("same")

legend.AddEntry("", f"Average OvO ROC AUC = {globalTest_auc:.3f}", "")

legend.Draw()
c.Update()
c.RedrawAxis()
c.SetTickx()
c.SetTicky()
#c.SetGrid()

thesisText = ROOT.TLatex(0.15, 0.9, "This Thesis")
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.045)
thesisText.Draw()

ppText = ROOT.TLatex(0.15, 0.86, "pp collisions")
ppText.SetNDC()
ppText.SetTextFont(42)
ppText.SetTextSize(0.04)
ppText.Draw()

EnergyText = ROOT.TLatex(0.15, 0.82, '#sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}')
EnergyText.SetNDC()
EnergyText.SetTextFont(42)
EnergyText.SetTextSize(0.04)
EnergyText.Draw("same")

DecayText = ROOT.TLatex(0.5, 0.9, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.04)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.5, 0.86, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.04)
ConjText.Draw("same")

ptText = ROOT.TLatex(0.5, 0.82, '2 < #it{p}_{#lower[-0.15]{T}} < 3 Ge#kern[-0.03]{V}/#it{c}')
ptText.SetNDC()
ptText.SetTextFont(42)
ptText.SetTextSize(0.04)
ptText.Draw("same")

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/ROC/ROC_Reduced.png")
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/ROC/ROC_Reduced.pdf")

'''

dictTest_roc = {'Bkg__Prompt D_{#lower[-0.5]{s}}^{+}': dictTest_roc['Bkg__Prompt D_{#lower[-0.5]{s}}^{+}'], 
                 'Bkg__FD D_{#lower[-0.3]{s}}^{+}': dictTest_roc['Bkg__FD D_{#lower[-0.3]{s}}^{+}'], 
                 'Prompt D_{#lower[-0.5]{s}}^{+}__FD D_{#lower[-0.3]{s}}^{+}': dictTest_roc['Prompt D_{#lower[-0.5]{s}}^{+}__FD D_{#lower[-0.3]{s}}^{+}']}

for i, comb in enumerate(dictTest_roc.keys()):
    #fpr, tpr, rocauc = dictTrain_roc[comb]
    #graphsTrain.append(ROOT.TGraph(len(fpr), array('d', fpr), array('d', tpr)))
    #graphsTrain[-1].SetTitle()
    #graphsTrain[-1].SetLineColor(colors[i])
    #graphsTrain[-1].SetLineWidth(2)
    #graphsTrain[-1].SetLineStyle(2)
    #graphsTrain[-1].Draw("L,same")
    #legend.AddEntry(graphsTrain[-1], f"Train: {comb.replace('__', ' vs ')} (ROC AUC = {rocauc:.3f})", "l")
    fpr, tpr, rocauc = dictTest_roc[comb]
    graphsTest.append(ROOT.TGraph(len(fpr), array('d', fpr), array('d', tpr)))
    graphsTest[-1].SetTitle()
    graphsTest[-1].SetLineColor(colors[i])
    graphsTest[-1].SetLineWidth(2)
    graphsTest[-1].Draw("L,same")
    reduced_Legend.AddEntry(graphsTest[-1], f"Test: {comb.replace('__', ' vs ')} (ROC AUC = {rocauc:.3f})", "l")

RandomClassifier = ROOT.TGraph(2, array('d', [0, 1]), array('d', [0, 1]))
RandomClassifier.SetLineColor(ROOT.kGray)
RandomClassifier.SetLineStyle(9)
RandomClassifier.SetLineWidth(2)
RandomClassifier.Draw("L,same")
reduced_Legend.AddEntry(RandomClassifier, "Luck (ROC AUC = 0.5)", "l")

LineAtOne = ROOT.TLine(0, 1, 1, 1)
LineAtOne.SetLineColor(ROOT.kBlack)
LineAtOne.SetLineStyle(9)
LineAtOne.SetLineWidth(2)
LineAtOne.Draw("same")

reduced_Legend.AddEntry("", f"Test: Average OvO ROC AUC = {globalTest_auc:.3f}", "")

legend.Draw()
c.Update()
c.RedrawAxis()
c.SetTickx()
c.SetTicky()
c.SetGrid()

thesisText = ROOT.TLatex(0.15, 0.9, "This Thesis")
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.045)
thesisText.Draw()

ppText = ROOT.TLatex(0.15, 0.86, "pp collisions")
ppText.SetNDC()
ppText.SetTextFont(42)
ppText.SetTextSize(0.04)
ppText.Draw()

EnergyText = ROOT.TLatex(0.15, 0.82, '#sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}')
EnergyText.SetNDC()
EnergyText.SetTextFont(42)
EnergyText.SetTextSize(0.04)
EnergyText.Draw("same")

DecayText = ROOT.TLatex(0.5, 0.9, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.04)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.5, 0.86, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.04)
ConjText.Draw("same")

ptText = ROOT.TLatex(0.5, 0.82, '2 < #it{p}_{#lower[-0.15]{T}} < 3 Ge#kern[-0.03]{V}/#it{c}')
ptText.SetNDC()
ptText.SetTextFont(42)
ptText.SetTextSize(0.04)
ptText.Draw("same")

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/ROC/ROC_Reduced.png")
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/ROC/ROC_Reduced.pdf")'''