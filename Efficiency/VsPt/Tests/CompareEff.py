import ROOT
import matplotlib.pyplot as plt

infileEff = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/Efficiencies.root")
hEffDsPrompt = infileEff.Get("Eff_DsPrompt")
hEffDsPrompt.SetDirectory(0)
hEffDplusPrompt = infileEff.Get("Eff_DplusPrompt")
hEffDplusPrompt.SetDirectory(0)

infileEffNoNorm = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/EfficienciesnoDL.root")
hEffDsPromptNoDL = infileEffNoNorm.Get("Eff_DsPrompt")
hEffDsPromptNoDL.SetDirectory(0)
hEffDplusPromptNoDL = infileEffNoNorm.Get("Eff_DplusPrompt")
hEffDplusPromptNoDL.SetDirectory(0)

ratioDplusOverDs = hEffDplusPrompt.Clone("ratioDplusOverDs")
ratioDplusOverDs.Divide(hEffDsPrompt)
ratioDplusOverDs.SetMarkerStyle(20)
ratioDplusOverDs.SetMarkerSize(2)
ratioDplusOverDs.SetLineColor(ROOT.kRed)
ratioDplusOverDs.SetMarkerColor(ROOT.kRed)


ratioDplusOverDsNoNorm = hEffDplusPromptNoDL.Clone("ratioDplusOverDsNoNorm")
ratioDplusOverDsNoNorm.Divide(hEffDsPromptNoDL)
ratioDplusOverDsNoNorm.SetMarkerStyle(20)
ratioDplusOverDsNoNorm.SetMarkerSize(2)
ratioDplusOverDsNoNorm.SetLineColor(ROOT.kAzure+3)
ratioDplusOverDsNoNorm.SetMarkerColor(ROOT.kAzure+3)

c = ROOT.TCanvas("c", "c", 800, 600)
c.DrawFrame(0, 1, 24, 3.5, ";p_{T} [GeV/c]; D^{+} / D_{s}^{+} efficiency")
ratioDplusOverDs.Draw("same")
ratioDplusOverDsNoNorm.Draw("same")

leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(ratioDplusOverDs, "Prompt, with Norm DL", "l")
leg.AddEntry(ratioDplusOverDsNoNorm, "Prompt, without Norm DL", "l")
leg.Draw()

ROOT.gStyle.SetOptStat(0)

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/ComparisonEff_noDL.pdf")
