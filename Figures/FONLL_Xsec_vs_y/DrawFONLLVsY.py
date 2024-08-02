import ROOT 
import pandas as pd 
import numpy as np

# Load the data
df = pd.read_csv("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FONLL_Xsec_vs_y/Predictions.txt", sep=" ")
df["min"] = df["central"] - df["min_sc"]
df["max"] = df["max_sc"] - df["central"]
step = df["y"][1] - df["y"][0]
step = [step/2]*len(df["y"])

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

# Create the graph
graph = ROOT.TGraphAsymmErrors(len(df["y"]), np.asarray(df["y"],'d'), np.asarray(df["central"],'d'), np.asarray(step,'d'), np.asarray(step,'d'), np.asarray(df["min"],'d'), np.asarray(df["max"],'d'))
graph.SetFillColorAlpha(ROOT.kRed, 0.5)
graph.SetLineColor(ROOT.kRed)
graph.SetMarkerColor(ROOT.kRed)


# Create the canvas
canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
canvas.SetLeftMargin(0.15)
canvas.SetRightMargin(0.05)
canvas.SetTopMargin(0.06)
canvas.SetBottomMargin(0.12)
hFrame = canvas.DrawFrame(-10., 0, 10, 2.5e9, ";#it{y};d#it{#sigma}/d#it{y} (pb)")
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.3)
hFrame.GetYaxis().SetLabelSize(0.05)
hFrame.GetXaxis().SetLabelSize(0.05)

graph.Draw("2,same")

gLine = graph.Clone("gLine")
gLine.SetLineWidth(2)
for i in range(gLine.GetN()):
    gLine.SetPointEYhigh(i, 0)
    gLine.SetPointEYlow(i, 0)
    
gLine.Draw("Z,same")

#graph.Draw("pe,same")

thesisText = ROOT.TLatex(0.2, 0.85, "This Thesis")
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.07)
thesisText.Draw()

ppText = ROOT.TLatex(0.2, 0.8, "pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
ppText.SetNDC()
ppText.SetTextFont(42)
ppText.SetTextSize(0.05)
ppText.Draw()

FONLLText = ROOT.TLatex(0.2, 0.73, "FONLL calculations")
FONLLText.SetNDC()
FONLLText.SetTextFont(42)
FONLLText.SetTextSize(0.05)
FONLLText.Draw()

CharmText = ROOT.TLatex(0.2, 0.67, "Charm-quark production cross-section")
CharmText.SetNDC()
CharmText.SetTextFont(42)
CharmText.SetTextSize(0.05)
CharmText.Draw()

ROOT.gPad.RedrawAxis()

canvas.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FONLL_Xsec_vs_y/FONLLVsY.pdf")
canvas.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/FONLL_Xsec_vs_y/FONLLVsY.png")