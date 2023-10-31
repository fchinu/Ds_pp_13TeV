import pandas as pd
import matplotlib.pyplot as plt

#infile = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/o2_Data_AO2D_Merged_Skimmed.parquet"
infile = "/home/fchinu/Run3/Ds_pp_13TeV/Datasets/o2_Data_AO2D_Merged.parquet"
variables = ["fCpa","fDecayLength"]
df = pd.read_parquet(infile)

# Draw scatter plot
plt.scatter(df["fCpa"], df["fDecayLength"], s=1)
plt.xlabel("fCpa")
plt.ylabel("fDecayLength")

# Save scatterplot as png
plt.savefig("/home/fchinu/Run3/Ds_pp_13TeV/DatasetComparison/VariablesScatterplots/fCpa_vs_fDecayLength_for_Data_Ds_Rebecca.png")
#plt.show()

# Calculate correlation
corr = df[variables].corr()

# Print correlation
print("Correlaion Matrix:") 
print(corr)
