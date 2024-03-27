import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils/")
from EvaluateEfficiency import calculate_efficiencies_with_unc
# Load the data
DplusPromptDf = pd.read_parquet("/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt4_4.5/DplusPrompt_pT_4_4.5_ModelApplied.parquet.gzip")
DsPromptDf = pd.read_parquet("/home/fchinu/Run3/Ds_pp_13TeV/ML/Application/pt4_4.5/DsPrompt_pT_4_4.5_ModelApplied.parquet.gzip")

# Define the cuts
cuts = np.linspace(0, 0.3, 100)
effsDs = []
effsErrDs = []
effsDplus = []
effsErrDplus = []

for cut in cuts:
    eff, error = calculate_efficiencies_with_unc(DsPromptDf, f"ML_output_Bkg < {cut}")
    effsDs.append(eff)
    effsErrDs.append(error)
    eff, error = calculate_efficiencies_with_unc(DplusPromptDf, f"ML_output_Bkg < {cut}")
    effsDplus.append(eff)
    effsErrDplus.append(error)

# Create a DataFrame for Seaborn
data = pd.DataFrame({'Cut': cuts, 'Efficiency Ds': effsDs, 'Error Ds': effsErrDs, 'Efficiency Dplus': effsDplus, 'Error Dplus': effsErrDplus})

# Set the Seaborn style
sns.set(style="whitegrid")

# Plot
plt.figure(figsize=(10, 6))
plt.errorbar(data['Cut'], data['Efficiency Ds'], yerr=data['Error Ds'], label='Ds')
plt.errorbar(data['Cut'], data['Efficiency Dplus'], yerr=data['Error Dplus'], label='Dplus')

# Set plot labels and title
plt.xlabel('BDT Cut on background score')
plt.ylabel('Efficiency')
plt.title('Efficiencies for Ds and Dplus vs. BDT Cut on background score')

# Show the legend
plt.legend()

# Save the plot
plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/Tests/EfficiencyScan/EfficiencyScan.pdf")

# Show the plot
plt.show()
