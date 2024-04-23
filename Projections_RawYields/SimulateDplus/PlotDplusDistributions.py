import uproot
import pandas as pd
import numpy as np
import ROOT

MyRed = ROOT.TColor.GetFreeColorIndex()
cMyRed = ROOT.TColor(MyRed, 247./255, 104./255, 104./255, 'MyRed', 1.0)

ROOT.gStyle.SetPalette(ROOT.kLake)
colors = ROOT.TColor.GetPalette()
# Load the data
TreeFileName = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/SimulateDplus/dplus_to_Kpipi_pythia8_seed42.root"

Tree = uproot.open(TreeFileName)["treeD"]
df = Tree.arrays(library="pd")

# Save the data
ptMins = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20]
ptMaxs = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24]

ptEdges = ptMins + [ptMaxs[-1]]

histos = []
for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    df_pt = df.query(f"ptD > {ptMin} and ptD < {ptMax}")
    histos.append(ROOT.TH1F(f"histo_{ptMin}_{ptMax}", f"histo_{ptMin}_{ptMax}", 320, 1.8, 5))
    histos[-1].SetLineWidth(2)
    histos[-1].SetLineColor(ROOT.kBlack)
    histos[-1].SetFillColor(colors[len(colors)//(len(ptMins)+1)*(iPt+1)])
    for mass in df_pt["massKaKaPi"]:
        histos[-1].Fill(mass)

histo2d = ROOT.TH2F("histo2d", "histo2d", 13, np.asarray(ptEdges, "d"), 320, 1.8, 5)
for mass, pt in zip(df["massKaKaPi"], df["ptD"]):
    histo2d.Fill(pt, mass)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadLeftMargin(0.17)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gStyle.SetPadBottomMargin(0.13)


canvas = ROOT.TCanvas("canvas", "canvas", 1200, 900)
canvas.Divide(4, 3, 0.005, 0.0001)
ThesisLabel = ROOT.TLatex(0.4, 0.8, "This Thesis")
ThesisLabel.SetNDC()
ThesisLabel.SetTextFont(42)
ThesisLabel.SetTextSize(0.08)
PYTHIALabel = ROOT.TLatex(0.4, 0.73, "PYTHIA 8 simulation")
PYTHIALabel.SetNDC()
PYTHIALabel.SetTextFont(42)
PYTHIALabel.SetTextSize(0.05)
EventsLabel = ROOT.TLatex(0.4, 0.66, "10^{7} generated D^{+}")
EventsLabel.SetNDC()
EventsLabel.SetTextFont(42)
EventsLabel.SetTextSize(0.05)
DecayLabel = ROOT.TLatex(0.4, 0.59, "D^{+} #rightarrow #pi^{+} K^{#font[122]{-}} #pi^{+}")
DecayLabel.SetNDC()
DecayLabel.SetTextFont(42)
DecayLabel.SetTextSize(0.05)
MisIDLabel = ROOT.TLatex(0.4, 0.52, "#pi^{+} mis-identified as K^{+}")
MisIDLabel.SetNDC()
MisIDLabel.SetTextFont(42)
MisIDLabel.SetTextSize(0.05)
ptLabels = []

for i, h in enumerate(histos):
    hFrame = canvas.cd(i+1).DrawFrame(1.8, 0, 4.5, 0.1, ";Mass (GeV/#it{c}^{2});Counts (a.u.)")
    NormHist = h.DrawNormalized("hist,same")
    hFrame.GetYaxis().SetRangeUser(0, NormHist.GetMaximum()*1.1)
    hFrame.GetYaxis().SetTitleOffset(1.3)
    hFrame.GetYaxis().SetTitleSize(0.06)
    hFrame.GetYaxis().SetLabelSize(0.05)
    hFrame.GetYaxis().SetMaxDigits(3)
    hFrame.GetXaxis().SetTitleOffset(1.)
    hFrame.GetXaxis().SetTitleSize(0.06)
    hFrame.GetXaxis().SetLabelSize(0.05)
    ThesisLabel.Draw("same")
    DecayLabel.Draw("same")
    MisIDLabel.Draw("same")
    PYTHIALabel.Draw("same")
    EventsLabel.Draw("same")
    ptLabels.append(ROOT.TLatex(0.4, 0.45, f"{ptMins[i]} < #it{{p}}_{{T}} (D^{{+}}) < {ptMaxs[i]} GeV/#it{{c}}"))
    ptLabels[-1].SetNDC()
    ptLabels[-1].SetTextFont(42)
    ptLabels[-1].SetTextSize(0.05)
    ptLabels[-1].Draw("same")


canvas.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/SimulateDplus/outhistos.png")
canvas.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/SimulateDplus/outhistos.pdf")


outfileName = "/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/SimulateDplus/outhistos.root"
outfile = ROOT.TFile(outfileName, "RECREATE")
for h in histos:
    h.Write()
histo2d.Write()
canvas.Write()
outfile.Close()


