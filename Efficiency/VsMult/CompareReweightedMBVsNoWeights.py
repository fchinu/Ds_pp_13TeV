import ROOT 

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadRightMargin(0.05)

if __name__ == "__main__":
    infileWeights = ROOT.TFile.Open('/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/VsMult/Efficiency_LHC24d3b.root')
    effDsPromptWeights = infileWeights.Get('Eff_DsPrompt_Cent_0_100')
    effDsPromptWeights.SetDirectory(0)
    infileWeights.Close()

    infileNoWeights = ROOT.TFile.Open('/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/VsMult/Efficiency_noWeights.root')
    effDsPromptNoWeights = infileNoWeights.Get('Eff_DsPrompt')
    effDsPromptNoWeights.SetDirectory(0)
    infileNoWeights.Close()

    effDsPromptWeights.SetMarkerColor(ROOT.kRed)
    effDsPromptWeights.SetTitle("Efficiency;#it{p}_{T} (GeV/#it{c});Efficiency")
    effDsPromptWeights.SetMarkerStyle(ROOT.kFullCircle)
    effDsPromptWeights.SetMarkerSize(1.5)
    effDsPromptWeights.SetLineColor(ROOT.kRed)

    effDsPromptNoWeights.SetMarkerColor(ROOT.kAzure+3)
    effDsPromptNoWeights.SetMarkerStyle(ROOT.kFullCircle)
    effDsPromptNoWeights.SetMarkerSize(1.5)
    effDsPromptNoWeights.SetLineColor(ROOT.kAzure+3)

    
    leg = ROOT.TLegend(0.35, 0.15, 0.7, 0.25)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(effDsPromptWeights, "Prompt D_{s}^{+} with weights")
    leg.AddEntry(effDsPromptNoWeights, "Prompt D_{s}^{+} without weights")

    c = ROOT.TCanvas("c", "c", 800, 600)
    c.Divide(2,1)

    c.cd(1).SetLogy()
    effDsPromptWeights.Draw()
    effDsPromptNoWeights.Draw("same")
    leg.Draw("same")

    c.cd(2)
    hRatio = effDsPromptWeights.Clone("hRatio")
    hRatio.SetTitle("Ratio;#it{p}_{T} (GeV/#it{c});Ratio")
    hRatio.Divide(effDsPromptNoWeights)
    hRatio.Draw()

    c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/VsMult/CompareEffMB_Weights.pdf")
