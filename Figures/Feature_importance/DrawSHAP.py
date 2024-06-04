import ROOT 
import numpy as np
import pandas as pd 
import shap
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from DfUtils import read_parquet_in_batches
from itertools import combinations
from array import array
import matplotlib as mpl
import matplotlib.pyplot as plt

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

VarsToDraw = ['fCpa', 'fCpaXY', 'fDecayLength', 'fDecayLengthXY', 'fImpactParameterXY', 'fAbsCos3PiK', 'fImpactParameter0', 'fImpactParameter1',
                       'fImpactParameter2', 'fNSigTpcTofPi0', 'fNSigTpcTofKa0', 'fNSigTpcTofPi1', 'fNSigTpcTofKa1', 'fNSigTpcTofPi2', 'fNSigTpcTofKa2']
VarsTitles = [r'cos$\theta_p$', r'cos$\theta_p^{xy}$', r'$L$', r'$L^{xy}$', r'$d^{xy}$', r'|cos$^3\theta$(K)|', r'$d^{xy}_{0}$', r'$d^{xy}_{1}$',
             r'$d^{xy}_{2}$', r'n$\sigma_{\mathrm{comb}}^{\pi}(0)$', r'n$\sigma_\mathrm{comb}^{K}(0)$', r'n$\sigma_\mathrm{comb}^{\pi}(1)$', r'n$\sigma_\mathrm{comb}^{K}(1)$', r'n$\sigma_\mathrm{comb}^{\pi}(2)$', 
             r'n$\sigma_\mathrm{comb}^{K}(2)$']


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

# Get the SHAP values
yPredTrainRaw = ModelHandl.predict(TrainSet, False)
yPredTestRaw = ModelHandl.predict(TestSet, False)

class_labels, class_counts = np.unique(yTest, return_counts=True)
n_classes = len(class_labels)
n_sample=10000
for class_count in class_counts:
    n_sample = min(n_sample, class_count)

subs = []
for class_lab in class_labels:
    subs.append(TestSet[yTest == class_lab].sample(n_sample))

df_subs = pd.concat(subs)
df_subs = df_subs[VarsToDraw]


explainer = shap.TreeExplainer(ModelHandl.get_original_model())
shap_values = explainer.shap_values(df_subs, approximate=True)

print(np.array(shap_values).shape)

h = ROOT.TH2D("h", "h", 100, -10, 10, 15, 0, 15)
for i, feat in enumerate(shap_values[0].T):
    for j, sh in enumerate(feat):
        h.Fill(sh, i)
        
outfile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Feature_importance/shap.root", "RECREATE")
h.Write()
outfile.Close()


fig = plt.figure(figsize=(10, 10))
cmap = mpl.cm.get_cmap('coolwarm', 256)
shap.summary_plot(shap_values[1], df_subs, feature_names= VarsTitles, plot_size=(
                18, 9), show=False, cmap=cmap, color_bar=False)
cbar = plt.colorbar(ticks=[0,1]) #Removes ticks (cmap range == 0.888, 0.99)
colbar = fig.axes[-1]
cbar_min, cbar_max = colbar.get_ylim()
colbar.text(2.2, cbar_min, 'Low', ha='right', va='bottom', color='black', size=12)
colbar.text(2.2, cbar_max, 'High', ha='right', va='top', color='black', size=12)
#cbar.ax.set_yticklabels(['Low', 'High'])
#cbar.set_label(label='Feature value', size=15, pad=1)
colbar.set_ylabel('Feature value', fontsize=15)
plt.xlabel(r"SHAP value (impact on model's prompt D$_s^+$ scores)", fontsize=15)
plt.ylabel("Feature", fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Add This Thesis text
plt.text(0.1, 0.5, "This Thesis", fontsize=25, transform=plt.gcf().transFigure)
plt.text(0.1, 0.45, r"pp collisions, $\sqrt{s}=13.6$ TeV", fontsize=20, transform=plt.gcf().transFigure)
plt.text(0.1, 0.3, r"$\mathrm{D_s^+}, \mathrm{D^+} \rightarrow \phi\pi^+ \rightarrow \mathrm{K^+K^-}\pi^+$", fontsize=20, transform=plt.gcf().transFigure)
plt.text(0.1, 0.25, "and charge conjugate", fontsize=20, transform=plt.gcf().transFigure)
plt.text(0.1, 0.2, r"2 < $p_\mathrm{T} < 3$ GeV/$c$", fontsize=20, transform=plt.gcf().transFigure)




plt.tight_layout()


plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Feature_importance/shap.png", bbox_inches='tight')
plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Feature_importance/shap.pdf", bbox_inches='tight')