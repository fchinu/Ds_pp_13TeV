import pandas as pd
from array import array
import ROOT

def set_custom_palette():
    n_colors = 5
    stops = array('d', [0.00, 0.01, 0.3, 0.8, 1])
    red   = array('d', [0.00, 141/255,210/255, 236/255, 241/255])
    green = array('d', [0.00, 3/255,   99/255, 232/255, 242/255])
    blue  = array('d', [0.00, 22/255,   8/255,  10/255, 217/255])
    
    ROOT.TColor.CreateGradientColorTable(n_colors, stops, red, green, blue, 255)
    ROOT.gStyle.SetNumberContours(255)

# Read parquet
df = pd.read_parquet("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train213175/LHC24d3a_PromptDplus.parquet")

dfDs = pd.read_parquet("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/MC/Train213175/LHC24d3a_PromptDs.parquet")
df = pd.concat([df, dfDs])

ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)
set_custom_palette()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.TGaxis.SetMaxDigits(3)

canvas = ROOT.TCanvas("c", "c", 800, 600)
ROOT.gPad.SetRightMargin(0.18)
ROOT.gPad.SetLeftMargin(0.1)
ROOT.gPad.SetBottomMargin(0.1)
ROOT.gPad.SetTopMargin(0.05)
hFrame = canvas.DrawFrame(0, -1.5, 24, 2, ";#it{p}_{T} (GeV/#it{c});#it{y}; Counts (a.u.)")

hFrame.GetZaxis().SetTitleSize(0.05)
hFrame.GetZaxis().SetTitleOffset(0.1)
hFrame.GetXaxis().SetTitleSize(0.04)
hFrame.GetXaxis().SetTitleOffset(1)
hFrame.GetYaxis().SetTitleSize(0.04)
hFrame.GetYaxis().SetTitleOffset(1)

histo = ROOT.TH2F("histo", ";#it{p}_{T} (Ge#kern[-0.05]{V}/#it{c});#it{y}; Normalised counts", 100, 0, 24, 100, -1.5, 1.5)
histo.GetZaxis().SetTitleSize(0.04)
histo.GetZaxis().SetTitleOffset(1.4)
for i, (pt, y) in df[["fPt", "fY"]].iterrows():
    histo.Fill(pt, y)

#histo.Smooth(1)
histo.DrawNormalized("colz;same")

thesisText = ROOT.TLatex(0.15, 0.85, 'This Thesis')
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.055)
thesisText.Draw("same")

collisionText = ROOT.TLatex(0.15, 0.8, 'pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.05]{V}')
collisionText.SetNDC()
collisionText.SetTextFont(42)
collisionText.SetTextSize(0.04)
collisionText.Draw("same")

DecayText = ROOT.TLatex(0.15, 0.74, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+} and charge conj.')
DecayText.SetNDC()
DecayText.SetTextFont(42)
DecayText.SetTextSize(0.04)
DecayText.Draw("same")

PtText = ROOT.TLatex(0.15, 0.15, '#it{p}_{T}^{dau} > 0.1 Ge#kern[-0.05]{V}/#it{c}')
PtText.SetNDC()
PtText.SetTextFont(42)
PtText.SetTextSize(0.04)
PtText.Draw("same")

EtaText = ROOT.TLatex(0.15, 0.20, '|#it{#eta}^{dau}| < 0.8')
EtaText.SetNDC()
EtaText.SetTextFont(42)
EtaText.SetTextSize(0.04)
EtaText.Draw("same")

outfile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/Figures/PtVsY/PtVsY.root", "RECREATE")
canvas.Write()
histo.Write()
outfile.Close()

canvas.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/PtVsY/PtVsY.pdf")
    