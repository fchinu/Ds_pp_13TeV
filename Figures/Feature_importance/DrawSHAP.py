import ROOT 
import numpy as np
import pandas as pd 
import shap
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette
from DfUtils import read_parquet_in_batches
from itertools import combinations
from array import array
import matplotlib as mpl
import matplotlib.pyplot as plt

from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve, auc

dfDsPrompt = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train189890/LHC22b1b_PromptDs_Train.parquet", "2<fPt<3")
dfDsFD = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train189892/LHC22b1a_NonPromptDs_Train.parquet", "2<fPt<3")
dfBkg = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train189072/LHC22o_pass6_small.parquet", "2<fPt<3 and (1.7 < fM < 1.75 or  2.1 < fM < 2.15)")
dfBkg = dfBkg.sample(frac=1, random_state=42)

VarsToDraw = ['fCpa', 'fCpaXY', 'fDecayLength', 'fDecayLengthXY', 'fImpactParameterXY', 'fAbsCos3PiK', 'fImpactParameter0', 'fImpactParameter1',
                       'fImpactParameter2', 'fNSigTpcTofPi0', 'fNSigTpcTofKa0', 'fNSigTpcTofPi1', 'fNSigTpcTofKa1', 'fNSigTpcTofPi2', 'fNSigTpcTofKa2']
VarsTitles = [r'cos$\theta_p$', r'cos$\theta_p^{xy}$', r'$L$', r'$L^{xy}$', r'$d^{xy}$', r"|cos$^3\theta'$(K)|", r'$d^{xy}_{0}$', r'$d^{xy}_{1}$',
             r'$d^{xy}_{2}$', r'$n\sigma_{\mathrm{comb}}^{\pi}(0)$', r'$n\sigma_\mathrm{comb}^{K}(0)$', r'$n\sigma_\mathrm{comb}^{\pi}(1)$', r'$n\sigma_\mathrm{comb}^{K}(1)$', r'$n\sigma_\mathrm{comb}^{\pi}(2)$', 
             r'$n\sigma_\mathrm{comb}^{K}(2)$']
VarsTitles_ROOT = ['#it{M}', '#it{p}_{T}', 'cos#it{#theta}_{p}', 'cos#it{#theta}_{p}^{xy}', 'L', 'L^{xy}', 'd^{xy}', "|cos^{3}#it{#theta'}(K)|", 'd^{xy}_{0}', 'd^{xy}_{1}',
             'd^{xy}_{2}', '#it{n}#sigma_{comb}^{#pi}(0)', '#it{n}#sigma_{comb}^{K}(0)', '#it{n}#sigma_{comb}^{#pi}(1)', '#it{n}#sigma_{comb}^{K}(1)', '#it{n}#sigma_{comb}^{#pi}(2)', 
             '#it{n}#sigma_{comb}^{K}(2)']


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
colbar.text(2.2, cbar_min, 'Low', ha='right', va='bottom', color='black', size=15)
colbar.text(2.2, cbar_max, 'High', ha='right', va='top', color='black', size=15)
#cbar.ax.set_yticklabels(['Low', 'High'])
#cbar.set_label(label='Feature value', size=15, pad=1)
colbar.set_ylabel('Feature value', fontsize=18)
plt.xlabel(r"SHAP value (impact on model's prompt D$_s^+$ scores)", fontsize=18)
plt.ylabel("Feature", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

# Add This Thesis text
plt.text(0.12, 0.5, "This Thesis", fontsize=25, transform=plt.gcf().transFigure)
plt.text(0.12, 0.45, r"pp collisions, $\sqrt{s}=13.6$ TeV", fontsize=20, transform=plt.gcf().transFigure)
plt.text(0.12, 0.3, r"$\mathrm{D}_\mathrm{s}^+, \mathrm{D^+} \rightarrow \phi\pi^+ \rightarrow \mathrm{K^+K^-}\pi^+$", fontsize=20, transform=plt.gcf().transFigure)
plt.text(0.12, 0.25, "and charge conjugate", fontsize=20, transform=plt.gcf().transFigure)
plt.text(0.12, 0.2, r"2 < $p_\mathrm{T} < 3$ GeV/$c$", fontsize=20, transform=plt.gcf().transFigure)

plt.tight_layout()


plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Feature_importance/shap.png", bbox_inches='tight')
plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Feature_importance/shap.pdf", bbox_inches='tight')

mean_bkg = np.mean(np.abs(shap_values[0]), axis=0)
mean_prompt = np.mean(np.abs(shap_values[1]), axis=0)
mean_FD = np.mean(np.abs(shap_values[2]), axis=0)

print(mean_bkg.shape)

ordered_index = np.argsort(mean_bkg + mean_prompt + mean_FD)

# Draw summary bar plot

fig, ax = plt.subplots(figsize=(15, 10))
index = np.arange(len(VarsTitles))
bar_width = 0.7
opacity = 0.8

rects1 = plt.barh(index+0.75, mean_bkg[ordered_index], bar_width, alpha=opacity, label='Background')
rects2 = plt.barh(index+0.75, mean_prompt[ordered_index], bar_width, alpha=opacity, left=mean_bkg[ordered_index], label=r'Prompt $\mathrm{D_s^+}$')
rects3 = plt.barh(index+0.75, mean_FD[ordered_index], bar_width, alpha=opacity, left=mean_bkg[ordered_index] + mean_prompt[ordered_index], label=r'Non-prompt $\mathrm{D_s^+}$')


plt.xlabel('Mean |SHAP value| (average impact on model output magnitude)', fontsize=20)
plt.ylabel('Feature', fontsize=20)
#plt.title('Feature importance', fontsize=18)
plt.yticks(index + bar_width, [VarsTitles[i] for i in ordered_index], fontsize=20)
plt.xticks(fontsize=20)
plt.legend(fontsize=20)

plt.set_cmap('tab10')

# Add This Thesis text
plt.text(0.55, 0.6, "This Thesis", fontsize=40, transform=plt.gcf().transFigure)
plt.text(0.55, 0.55, r"pp collisions, $\sqrt{s}=13.6$ TeV", fontsize=30, transform=plt.gcf().transFigure)
plt.text(0.55, 0.4, r"$\mathrm{D}_\mathrm{s}^+, \mathrm{D^+} \rightarrow \phi\pi^+ \rightarrow \mathrm{K^+K^-}\pi^+$", fontsize=30, transform=plt.gcf().transFigure)
plt.text(0.55, 0.35, "and charge conjugate", fontsize=30, transform=plt.gcf().transFigure)
plt.text(0.55, 0.3, r"2 < $p_\mathrm{T} < 3$ GeV/$c$", fontsize=30, transform=plt.gcf().transFigure)

plt.tight_layout()

plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Feature_importance/SHAP_summary2_3.png", bbox_inches='tight')
plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Feature_importance/SHAP_summary2_3.pdf", bbox_inches='tight')


# 8-12

dfDsPrompt = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train189890/LHC22b1b_PromptDs_Train.parquet", "8<fPt<12")
dfDsFD = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train189892/LHC22b1a_NonPromptDs_Train.parquet", "8<fPt<12")
dfBkg = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train189072/LHC22o_pass6_small.parquet", "8<fPt<12 and (1.7 < fM < 1.75 or  2.1 < fM < 2.15)")
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
ModelHandl.load_model_handler("/home/fchinu/Run3/Ds_pp_13TeV/ML/Training/pt8_12/ModelHandler_pT_8_12.pickle")

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

mean_bkg = np.mean(np.abs(shap_values[0]), axis=0)
mean_prompt = np.mean(np.abs(shap_values[1]), axis=0)
mean_FD = np.mean(np.abs(shap_values[2]), axis=0)

print(mean_bkg.shape)

ordered_index = np.argsort(mean_bkg + mean_prompt + mean_FD)

# Draw summary bar plot

fig, ax = plt.subplots(figsize=(15, 10))
index = np.arange(len(VarsTitles))
bar_width = 0.7
opacity = 0.8

rects1 = plt.barh(index+0.75, mean_bkg[ordered_index], bar_width, alpha=opacity, label='Background')
rects2 = plt.barh(index+0.75, mean_prompt[ordered_index], bar_width, alpha=opacity, left=mean_bkg[ordered_index], label=r'Prompt $\mathrm{D_s^+}$')
rects3 = plt.barh(index+0.75, mean_FD[ordered_index], bar_width, alpha=opacity, left=mean_bkg[ordered_index] + mean_prompt[ordered_index], label=r'Non-prompt $\mathrm{D_s^+}$')


plt.xlabel('Mean |SHAP value| (average impact on model output magnitude)', fontsize=20)
plt.ylabel('Feature', fontsize=20)
#plt.title('Feature importance', fontsize=18)
plt.yticks(index + bar_width, [VarsTitles[i] for i in ordered_index], fontsize=20)
plt.xticks(fontsize=20)
plt.legend(fontsize=20)

plt.set_cmap('tab10')

# Add This Thesis text
plt.text(0.55, 0.6, "This Thesis", fontsize=40, transform=plt.gcf().transFigure)
plt.text(0.55, 0.55, r"pp collisions, $\sqrt{s}=13.6$ TeV", fontsize=30, transform=plt.gcf().transFigure)
plt.text(0.55, 0.4, r"$\mathrm{D}_\mathrm{s}^+, \mathrm{D^+} \rightarrow \phi\pi^+ \rightarrow \mathrm{K^+K^-}\pi^+$", fontsize=30, transform=plt.gcf().transFigure)
plt.text(0.55, 0.35, "and charge conjugate", fontsize=30, transform=plt.gcf().transFigure)
plt.text(0.55, 0.3, r"8 < $p_\mathrm{T} < 12$ GeV/$c$", fontsize=30, transform=plt.gcf().transFigure)

plt.tight_layout()

plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Feature_importance/SHAP_summary8_12.png", bbox_inches='tight')
plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Feature_importance/SHAP_summary8_12.pdf", bbox_inches='tight')
