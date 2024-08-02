import ROOT 
import pandas as pd
import numpy as np 
from array import array
import matplotlib as mpl
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from DfUtils import read_parquet_in_batches


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

dfDsPrompt = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train213175/LHC24d3a_PromptDs.parquet", "2<fPt<3")
dfDsFD = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train213175/LHC24d3a_NonPromptDs.parquet", "2<fPt<3")
dfBkg = read_parquet_in_batches("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train191566/LHC22o.parquet", "2<fPt<3 and (1.7 < fM < 1.75 or  2.1 < fM < 2.15)")

# Define the variables to plot
VarsToDraw = ['fM', 'fPt', 'fCpa', 'fCpaXY', 'fDecayLength', 'fDecayLengthXY', 'fImpactParameterXY', 'fAbsCos3PiK', 'fImpactParameter0', 'fImpactParameter1',
                       'fImpactParameter2', 'fNSigTpcTofPi0', 'fNSigTpcTofKa0', 'fNSigTpcTofPi1', 'fNSigTpcTofKa1', 'fNSigTpcTofPi2', 'fNSigTpcTofKa2']
VarsTitles = ['#it{M}', '#it{p}_{T}', 'cos#it{#theta}_{p}', 'cos#it{#theta}_{p}^{xy}', 'L', 'L^{xy}', 'd^{xy}', "|cos^{3}#it{#theta'}(K)|", 'd^{xy}_{0}', 'd^{xy}_{1}',
             'd^{xy}_{2}', '#it{n}#sigma_{comb}^{#pi}(0)', '#it{n}#sigma_{comb}^{K}(0)', '#it{n}#sigma_{comb}^{#pi}(1)', '#it{n}#sigma_{comb}^{K}(1)', '#it{n}#sigma_{comb}^{#pi}(2)', 
             '#it{n}#sigma_{comb}^{K}(2)']

# Reverse the axis
VarsToDraw.reverse()
VarsTitles.reverse()

set_matplotlib_palette("coolwarm")
#ROOT.gStyle.SetPalette(ROOT.kRedBlue)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleSize(0.04,"padTitle")
ROOT.gStyle.SetTitleH(0.05)

# Create the canvas
c = ROOT.TCanvas("c", "c", 800, 800)
c.Divide(2,2)

# Prompt Ds
c.cd(1)

# Create the correlation matrix
CorrMatrix = np.zeros((len(VarsToDraw), len(VarsToDraw)))
for i, var1 in enumerate(VarsToDraw):
    for j, var2 in enumerate(VarsToDraw):
        CorrMatrix[i, j] = dfDsPrompt[var1].corr(dfDsPrompt[var2])

ROOT.gPad.SetRightMargin(0.15)
ROOT.gPad.SetLeftMargin(0.16)
ROOT.gPad.SetBottomMargin(0.16)
ROOT.gPad.SetTopMargin(0.08)
# Create the histogram
hCorrPromptDs = ROOT.TH2D("hCorrPromptDs", "Prompt D_{s}^{+} correlation matrix", len(VarsToDraw), 0, len(VarsToDraw), len(VarsToDraw), 0, len(VarsToDraw))
hCorrPromptDs.GetZaxis().SetRangeUser(-1, 1)
hCorrPromptDs.GetXaxis().SetLabelSize(0.05)
hCorrPromptDs.GetYaxis().SetLabelSize(0.05)
for i in range(len(VarsToDraw)):
    hCorrPromptDs.GetXaxis().SetBinLabel(i+1, VarsTitles[i])
    hCorrPromptDs.GetYaxis().SetBinLabel(i+1, VarsTitles[i])
    hCorrPromptDs.GetXaxis().LabelsOption("v")
    for j in range(len(VarsToDraw)):
        hCorrPromptDs.SetBinContent(i+1, j+1, CorrMatrix[i, j])

# Draw the histogram
hCorrPromptDs.Draw("COLZ ")

# Draw a grid
ROOT.gPad.SetGrid()
ROOT.gPad.SetTickx()
ROOT.gPad.SetTicky()
ROOT.gPad.RedrawAxis()

# FD Ds
c.cd(2)

# Create the correlation matrix
CorrMatrix = np.zeros((len(VarsToDraw), len(VarsToDraw)))
for i, var1 in enumerate(VarsToDraw):
    for j, var2 in enumerate(VarsToDraw):
        CorrMatrix[i, j] = dfDsFD[var1].corr(dfDsFD[var2])

ROOT.gPad.SetRightMargin(0.15)
ROOT.gPad.SetLeftMargin(0.16)
ROOT.gPad.SetBottomMargin(0.16)
ROOT.gPad.SetTopMargin(0.08)
# Create the histogram
hCorrFDDs = ROOT.TH2D("hCorrFDDs", "FD D_{s}^{+} correlation matrix", len(VarsToDraw), 0, len(VarsToDraw), len(VarsToDraw), 0, len(VarsToDraw))
hCorrFDDs.GetZaxis().SetRangeUser(-1, 1)
hCorrFDDs.GetXaxis().SetLabelSize(0.05)
hCorrFDDs.GetYaxis().SetLabelSize(0.05)
for i in range(len(VarsToDraw)):
    hCorrFDDs.GetXaxis().SetBinLabel(i+1, VarsTitles[i])
    hCorrFDDs.GetYaxis().SetBinLabel(i+1, VarsTitles[i])
    hCorrFDDs.GetXaxis().LabelsOption("v")
    for j in range(len(VarsToDraw)):
        hCorrFDDs.SetBinContent(i+1, j+1, CorrMatrix[i, j])

# Draw the histogram
hCorrFDDs.Draw("COLZ ")

# Draw a grid
ROOT.gPad.SetGrid()
ROOT.gPad.SetTickx()
ROOT.gPad.SetTicky()
ROOT.gPad.RedrawAxis()

# Bkg
c.cd(3)

# Create the correlation matrix
CorrMatrix = np.zeros((len(VarsToDraw), len(VarsToDraw)))
for i, var1 in enumerate(VarsToDraw):
    for j, var2 in enumerate(VarsToDraw):
        CorrMatrix[i, j] = dfBkg[var1].corr(dfBkg[var2])

ROOT.gPad.SetRightMargin(0.15)
ROOT.gPad.SetLeftMargin(0.16)
ROOT.gPad.SetBottomMargin(0.16)
ROOT.gPad.SetTopMargin(0.08)
# Create the histogram
hCorrBkg = ROOT.TH2D("hCorrBkg", "Bkg. correlation matrix", len(VarsToDraw), 0, len(VarsToDraw), len(VarsToDraw), 0, len(VarsToDraw))
hCorrBkg.GetZaxis().SetRangeUser(-1, 1)
hCorrBkg.GetXaxis().SetLabelSize(0.05)
hCorrBkg.GetYaxis().SetLabelSize(0.05)
for i in range(len(VarsToDraw)):
    hCorrBkg.GetXaxis().SetBinLabel(i+1, VarsTitles[i])
    hCorrBkg.GetYaxis().SetBinLabel(i+1, VarsTitles[i])
    hCorrBkg.GetXaxis().LabelsOption("v")
    for j in range(len(VarsToDraw)):
        hCorrBkg.SetBinContent(i+1, j+1, CorrMatrix[i, j])

# Draw the histogram
hCorrBkg.Draw("COLZ ")

# Draw a grid
ROOT.gPad.SetGrid()
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

outfile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Figures/correlation_matrix/CorrelationMatrix.root", "RECREATE")
hCorrPromptDs.Write()
hCorrFDDs.Write()
hCorrBkg.Write()
c.Write()
outfile.Close()

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/correlation_matrix/CorrelationMatrix.pdf")



