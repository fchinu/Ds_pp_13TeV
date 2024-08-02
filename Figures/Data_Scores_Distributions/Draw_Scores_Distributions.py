import ROOT
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette


colors, _ = get_discrete_matplotlib_palette("tab10")



# Open the file
file = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train197555/LHC22o_AnalysisResults.root")
# Get the THnSparse
sparse = file.Get("hf-task-ds/Data/hSparseMass")

sparse.GetAxis(1).SetRangeUser(2.0, 2.5)

hBkg = sparse.Projection(3)
hBkg.SetDirectory(0)
hPrompt = sparse.Projection(4)
hPrompt.SetDirectory(0)
hFD = sparse.Projection(5)
hFD.SetDirectory(0)

file.Close()

# Set style
hBkg.SetFillColorAlpha(colors[0], 0.55)
hBkg.SetLineColor(colors[0])
hBkg.SetMarkerColor(colors[0])
hBkg.SetMarkerStyle(ROOT.kFullCircle)
hBkg.SetLineWidth(2)
hPrompt.SetFillColorAlpha(colors[1], 0.55)
hPrompt.SetLineColor(colors[1])
hPrompt.SetMarkerColor(colors[1])
hPrompt.SetMarkerStyle(ROOT.kFullCircle)
hPrompt.SetLineWidth(2)
hFD.SetFillColorAlpha(colors[2], 0.55)
hFD.SetLineColor(colors[2])
hFD.SetMarkerColor(colors[2])
hFD.SetMarkerStyle(ROOT.kFullCircle)
hFD.SetLineWidth(2)

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

hBkg.Rebin(20)
hPrompt.Rebin(2)
hFD.Rebin(2)

c = ROOT.TCanvas("c", "c", 800, 600)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.13)
ROOT.gPad.SetBottomMargin(0.13)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetLogy()

hFrame = ROOT.gPad.DrawFrame(0, 1.e-4, 1, 5,";BDT score;Counts (a.u.)")
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.2)
hFrame.GetYaxis().SetLabelSize(0.04)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1.1)
hFrame.GetXaxis().SetLabelSize(0.04)

hBkg.DrawNormalized("hist,same")
hPrompt.DrawNormalized("hist,same")
#hFD.DrawNormalized("hist,same")
hBkg.DrawNormalized("p,same")
hPrompt.DrawNormalized("p,same")
#hFD.DrawNormalized("p,same")

legend = ROOT.TLegend(0.65, 0.65, 0.9, 0.75)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.045)

legend.AddEntry(hBkg, "Background", "pl")
legend.AddEntry(hPrompt, "Prompt", "pl")
#legend.AddEntry(hFD, "Non#font[122]{-}prompt", "pl")

legend.Draw()

thesisText = ROOT.TLatex(0.2, 0.85, 'This Thesis')
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.06)
thesisText.Draw("same")

collisionText = ROOT.TLatex(0.2, 0.8, 'pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.05]{V}')
collisionText.SetNDC()
collisionText.SetTextFont(42)
collisionText.SetTextSize(0.05)
collisionText.Draw("same")

DecayText = ROOT.TLatex(0.2, 0.72, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.2, 0.67, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ptText = ROOT.TLatex(0.2, 0.6, '2.0 < #it{p}_{T} < 2.5 Ge#kern[-0.05]{V}/#it{c}')
ptText.SetNDC()
ptText.SetTextFont(42)
ptText.SetTextSize(0.05)
ptText.Draw("same")

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Data_Scores_Distributions/Data_ScoresDistributions.png")
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Data_Scores_Distributions/Data_ScoresDistributions.pdf")


