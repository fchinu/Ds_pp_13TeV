import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette

colors, _ = get_discrete_matplotlib_palette("tab20")
colors = colors[15::-1]


# Get inputs
histos = []
infile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Systematics/FD_Fraction/FDSystematic.root")
canvas = infile.Get("canvasCombinations").GetPrimitive("canvasCombinations_1")
for i in range(16):
    histos.append(canvas.GetPrimitive(f"hComb{i}"))
    histos[-1].SetDirectory(0)
    histos[-1].SetMarkerStyle(ROOT.kFullCircle)
    histos[-1].SetMarkerColor(colors[i])
    histos[-1].SetLineColor(colors[i])
    histos[-1].SetMarkerSize(1)
    histos[-1].SetLineWidth(2)
infile.Close()

hRatios = []
for idx, h in enumerate(histos):
    if idx == 0:
        continue
    hRatios.append(h.Clone())
    hRatios[-1].Divide(histos[0])

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

c = ROOT.TCanvas("c", "c", 1000, 600)
c.Divide(2,1,0.0001,0.0001)

c.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.12)
hFrame = ROOT.gPad.DrawFrame(0, 0.75, 24, 1.3, ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c}); #it{f}_{prompt}^{D_{s}^{+}}/f_{prompt}^{D^{+}}")
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.3)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetXaxis().SetLabelSize(0.04)
hFrame.GetYaxis().SetLabelSize(0.04)

for h in histos:
    h.Draw("same")


thesisText = ROOT.TLatex(0.2, 0.85, "This Thesis")
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.07)
thesisText.Draw()

ppText = ROOT.TLatex(0.2, 0.79, "pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
ppText.SetNDC()
ppText.SetTextFont(42)
ppText.SetTextSize(0.05)
ppText.Draw()

DecayText = ROOT.TLatex(0.2, 0.71, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.05)
DecayText.Draw("same")

ConjText = ROOT.TLatex(0.2, 0.66, 'and charge conjugate')
ConjText.SetNDC()
ConjText.SetTextFont(42)
ConjText.SetTextSize(0.05)
ConjText.Draw("same")

ROOT.gPad.RedrawAxis()

c.cd(2)
ROOT.gPad.SetLeftMargin(0.15)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetBottomMargin(0.12)
hFrame = ROOT.gPad.DrawFrame(0, 0.85, 24, 1.3, ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c}); Ratio to central case")
hFrame.GetYaxis().SetTitleSize(0.05)
hFrame.GetYaxis().SetTitleOffset(1.3)
hFrame.GetXaxis().SetTitleSize(0.05)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetXaxis().SetLabelSize(0.04)
hFrame.GetYaxis().SetLabelSize(0.04)

for h in hRatios:
    h.Draw("same")

legend = ROOT.TLegend(0.2, 0.75, 0.9, 0.92)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.04)
legend.SetNColumns(4)

for idx, h in enumerate(histos):
    if idx == 0:
        legend.AddEntry(h, "Central", "lp")
    else:
        legend.AddEntry(h, f"Trial {idx}", "lp")
legend.Draw("same") 


ROOT.gPad.RedrawAxis()

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Systematics/fprompt/PromptFracSyst.pdf")
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Systematics/fprompt/PromptFracSyst.png")
