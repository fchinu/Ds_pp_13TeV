import ROOT 
import pickle
import yaml
import numpy as np
import argparse
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette


def PlotCutVariationsOnePtBin(cfgFileName):
  
    # load inputs
    with open(cfgFileName, 'r') as file:
        config = yaml.safe_load(file)
    inDirName = config["inputs"]["directory"]
    inCommonFileNameRawY = config["inputs"]["commonfilenames"]["rawyield"]
    inCommonFileNameEff = config["inputs"]["commonfilenames"]["efficiency"]
    inCommonFileNamePromptFracDs = config["inputs"]["commonfilenames"]["prompt_frac_ds"]
    inCommonFileNamePromptFracDplus = config["inputs"]["commonfilenames"]["prompt_frac_dplus"]
    inCommonFileNameRatios = config["inputs"]["commonfilenames"]["ratios"]
    cutSetSuffix = config["inputs"]["cutsets"]
    nFiles = len(cutSetSuffix)

    outFileName = config["outfilename"]

    maxChi2 = config["quality"]["maxchisquare"]
    minSignif = config["quality"]["minsignif"]
    minRelSignif = config["quality"]["minrelsignif"]
    minRelEff = config["quality"]["minreleff"]
    maxRelEff = config["quality"]["maxreleff"]
    fillThrRelEff = config["quality"]["fillthrreleff"]
    fRelativeVariation = bool(config["plots"]["plotrelativevar"])
    relAssignedSyst = config["plots"]["relassignedsyst"] if "relassignedsyst" in config["plots"] else []
    
    colors, _ = get_discrete_matplotlib_palette("tab10")

    hRatios = []
    hRawYieldsDs = []
    hRawYieldsDplus = []
    hSignificancesDs = []
    hSignificancesDplus = []
    hSoverBsDs = []
    hSoverBsDplus = []
    hchi2s = []
    hEffPromptsDs = []
    hEffPromptsDplus = []
    hEffFDsDs = []
    hEffFDsDplus = []
    hPromptFracsDs = []
    hPromptFracsDplus = []
    hFDFracsDs = []
    hFDFracsDplus = []

    crossSectionTitle = "dN/dpt"

    for iFile in range(nFiles):

        infile_ratios = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNameRatios}{cutSetSuffix[iFile]}.root")
        hRatios.append(infile_ratios.Get("hRatio"))
        hRatios[iFile].SetDirectory(0)
        infile_ratios.Close()

        infile_rawyield = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNameRawY}{cutSetSuffix[iFile]}.root")
        hRawYieldsDs.append(infile_rawyield.Get("hRawYields"))
        hSignificancesDs.append(infile_rawyield.Get("hRawYieldsSignificance"))
        hSoverBsDs.append(infile_rawyield.Get("hRawYieldsSoverB"))
        hchi2s.append(infile_rawyield.Get("hRawYieldsChiSquare"))
        hRawYieldsDs[iFile].SetDirectory(0)
        hSignificancesDs[iFile].SetDirectory(0)
        hSoverBsDs[iFile].SetDirectory(0)
        hchi2s[iFile].SetDirectory(0)

        hRawYieldsDplus.append(infile_rawyield.Get("hRawYieldsSecondPeak"))
        hSignificancesDplus.append(infile_rawyield.Get("hRawYieldsSignificanceSecondPeak"))
        hSoverBsDplus.append(infile_rawyield.Get("hRawYieldsSoverBSecondPeak"))
        hRawYieldsDplus[iFile].SetDirectory(0)
        hSignificancesDplus[iFile].SetDirectory(0)
        hSoverBsDplus[iFile].SetDirectory(0)
        infile_rawyield.Close()

        infile_eff = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNameEff}{cutSetSuffix[iFile]}.root")
        hEffPromptsDs.append(infile_eff.Get("Eff_DsPrompt"))
        hEffFDsDs.append(infile_eff.Get("Eff_DsFD"))
        hEffPromptsDs[iFile].SetDirectory(0)
        hEffFDsDs[iFile].SetDirectory(0)

        hEffPromptsDplus.append(infile_eff.Get("Eff_DplusPrompt"))
        hEffFDsDplus.append(infile_eff.Get("Eff_DplusFD"))
        hEffPromptsDplus[iFile].SetDirectory(0)
        hEffFDsDplus[iFile].SetDirectory(0)
        infile_eff.Close()

        infile_DsPromptFrac = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNamePromptFracDs}{cutSetSuffix[iFile]}.root")
        hPromptFracsDs.append(infile_DsPromptFrac.Get("hRawFracPrompt"))
        hPromptFracsDs[iFile].SetDirectory(0)
        infile_DsPromptFrac.Close()

        infile_DplusPromptFrac = ROOT.TFile.Open(f"{inDirName}/{inCommonFileNamePromptFracDplus}{cutSetSuffix[iFile]}.root")
        hPromptFracsDplus.append(infile_DplusPromptFrac.Get("hRawFracPrompt"))
        hPromptFracsDplus[iFile].SetDirectory(0)
        infile_DplusPromptFrac.Close()

    hSyst = hRawYieldsDs[0].Clone("hSyst")

    gRawYieldsDsVsCut = ROOT.TGraphErrors()
    gRawYieldsDplusVsCut = ROOT.TGraphErrors()
    gEffPromptsDsVsCut = ROOT.TGraphErrors()
    gEffPromptsDplusVsCut = ROOT.TGraphErrors()
    gEffFDsDsVsCut = ROOT.TGraphErrors()
    gEffFDsDplusVsCut = ROOT.TGraphErrors()
    gPromptFracsDsVsCut = ROOT.TGraphErrors()
    gPromptFracsDplusVsCut = ROOT.TGraphErrors()
    gSignificancesDsVsCut = ROOT.TGraphErrors()
    gSignificancesDplusVsCut = ROOT.TGraphErrors()
    gSoverBsDsVsCut = ROOT.TGraphErrors()
    gSoverBsDplusVsCut = ROOT.TGraphErrors()
    gRatioVsCut = ROOT.TGraphErrors()
    hRatioDistr = ROOT.TH1D("hRatioDistr", "hRatioDistr", 80, 0.7, 1.1)

    iPt = 8     # 4.5<pt<5 GeV/c

    for iFile in range(nFiles):
        
        # chi2 quality check
        if hchi2s[iFile].GetBinContent(iPt+1) > 2:
            continue

        # Significance quality check
        if hSignificancesDs[iFile].GetBinContent(iPt+1) < minSignif or hSignificancesDplus[iFile].GetBinContent(iPt+1) < minSignif:
            print(f"Significance quality check failed for file {iFile}")
            continue
        if hSignificancesDs[iFile].GetBinContent(iPt+1) < hSignificancesDs[0].GetBinContent(iPt+1) * minRelSignif or hSignificancesDplus[iFile].GetBinContent(iPt+1) < hSignificancesDplus[0].GetBinContent(iPt+1) * minRelSignif:
            print(f"Relative significance quality check failed for file {iFile}")
            continue
        
        # Efficiency quality check
        if hEffPromptsDs[iFile].GetBinContent(iPt+1) < hEffPromptsDs[0].GetBinContent(iPt+1) * minRelEff or hEffPromptsDplus[iFile].GetBinContent(iPt+1) < hEffPromptsDplus[iFile].GetBinContent(iPt+1) * minRelEff:
            print(f"Relative efficiency quality check failed for file {iFile}")
            #continue

        gRawYieldsDsVsCut.SetPoint(iFile, iFile, hRawYieldsDs[iFile].GetBinContent(iPt+1)/hRawYieldsDs[0].GetBinContent(iPt+1))
        gRawYieldsDplusVsCut.SetPoint(iFile, iFile, hRawYieldsDplus[iFile].GetBinContent(iPt+1)/hRawYieldsDplus[0].GetBinContent(iPt+1))
        gEffPromptsDsVsCut.SetPoint(iFile, iFile, hEffPromptsDs[iFile].GetBinContent(iPt+1)/hEffPromptsDs[0].GetBinContent(iPt+1))
        gEffPromptsDplusVsCut.SetPoint(iFile, iFile, hEffPromptsDplus[iFile].GetBinContent(iPt+1)/hEffPromptsDplus[0].GetBinContent(iPt+1))
        gEffFDsDsVsCut.SetPoint(iFile, iFile, hEffFDsDs[iFile].GetBinContent(iPt+1)/hEffFDsDs[0].GetBinContent(iPt+1))
        gEffFDsDplusVsCut.SetPoint(iFile, iFile, hEffFDsDplus[iFile].GetBinContent(iPt+1)/hEffFDsDplus[0].GetBinContent(iPt+1))
        gPromptFracsDsVsCut.SetPoint(iFile, iFile, hPromptFracsDs[iFile].GetBinContent(iPt+1)/hPromptFracsDs[0].GetBinContent(iPt+1))
        gPromptFracsDplusVsCut.SetPoint(iFile, iFile, hPromptFracsDplus[iFile].GetBinContent(iPt+1)/hPromptFracsDplus[0].GetBinContent(iPt+1))
        gSignificancesDsVsCut.SetPoint(iFile, iFile, hSignificancesDs[iFile].GetBinContent(iPt+1))
        gSignificancesDsVsCut.SetPointError(iFile, 0, hSignificancesDs[iFile].GetBinError(iPt+1))
        gSignificancesDplusVsCut.SetPoint(iFile, iFile, hSignificancesDplus[iFile].GetBinContent(iPt+1))
        gSignificancesDplusVsCut.SetPointError(iFile, 0, hSignificancesDplus[iFile].GetBinError(iPt+1))
        gSoverBsDsVsCut.SetPoint(iFile, iFile, hSoverBsDs[iFile].GetBinContent(iPt+1))
        gSoverBsDsVsCut.SetPointError(iFile, 0, hSoverBsDs[iFile].GetBinError(iPt+1))
        gSoverBsDplusVsCut.SetPoint(iFile, iFile, hSoverBsDplus[iFile].GetBinContent(iPt+1))
        gSoverBsDplusVsCut.SetPointError(iFile, 0, hSoverBsDplus[iFile].GetBinError(iPt+1))
        gRatioVsCut.SetPoint(iFile, iFile, hRatios[iFile].GetBinContent(iPt+1))
        gRatioVsCut.SetPointError(iFile, 0, hRatios[iFile].GetBinError(iPt+1))
        hRatioDistr.Fill(hRatios[iFile].GetBinContent(iPt+1)/hRatios[0].GetBinContent(iPt+1))


    # Set style
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    gRawYieldsDsVsCut.SetMarkerStyle(ROOT.kFullCircle)
    gRawYieldsDsVsCut.SetMarkerColor(colors[0])
    gRawYieldsDsVsCut.SetLineColor(colors[0])
    gRawYieldsDsVsCut.SetMarkerSize(1.5)
    gRawYieldsDsVsCut.SetLineWidth(2)

    gRawYieldsDplusVsCut.SetMarkerStyle(ROOT.kFullDiamond)
    gRawYieldsDplusVsCut.SetMarkerColor(colors[1])
    gRawYieldsDplusVsCut.SetLineColor(colors[1])
    gRawYieldsDplusVsCut.SetMarkerSize(2)
    gRawYieldsDplusVsCut.SetLineWidth(2)

    gEffPromptsDsVsCut.SetMarkerStyle(ROOT.kFullCircle)
    gEffPromptsDsVsCut.SetMarkerColor(colors[0])
    gEffPromptsDsVsCut.SetLineColor(colors[0])
    gEffPromptsDsVsCut.SetMarkerSize(1.5)
    gEffPromptsDsVsCut.SetLineWidth(2)

    gEffFDsDsVsCut.SetMarkerStyle(ROOT.kOpenCircle)
    gEffFDsDsVsCut.SetMarkerColor(colors[0])
    gEffFDsDsVsCut.SetLineColor(colors[0])
    gEffFDsDsVsCut.SetMarkerSize(1.5)
    gEffFDsDsVsCut.SetLineWidth(2)

    gEffPromptsDplusVsCut.SetMarkerStyle(ROOT.kFullDiamond)
    gEffPromptsDplusVsCut.SetMarkerColor(colors[1])
    gEffPromptsDplusVsCut.SetLineColor(colors[1])
    gEffPromptsDplusVsCut.SetMarkerSize(2)
    gEffPromptsDplusVsCut.SetLineWidth(2)

    gEffFDsDplusVsCut.SetMarkerStyle(ROOT.kOpenDiamond)
    gEffFDsDplusVsCut.SetMarkerColor(colors[1])
    gEffFDsDplusVsCut.SetLineColor(colors[1])
    gEffFDsDplusVsCut.SetMarkerSize(2)
    gEffFDsDplusVsCut.SetLineWidth(2)

    gPromptFracsDsVsCut.SetMarkerStyle(ROOT.kFullCircle)
    gPromptFracsDsVsCut.SetMarkerColor(colors[0])
    gPromptFracsDsVsCut.SetLineColor(colors[0])
    gPromptFracsDsVsCut.SetMarkerSize(1.5)
    gPromptFracsDsVsCut.SetLineWidth(2)

    gPromptFracsDplusVsCut.SetMarkerStyle(ROOT.kFullDiamond)
    gPromptFracsDplusVsCut.SetMarkerColor(colors[1])
    gPromptFracsDplusVsCut.SetLineColor(colors[1])
    gPromptFracsDplusVsCut.SetMarkerSize(2)
    gPromptFracsDplusVsCut.SetLineWidth(2)
        
    gSignificancesDsVsCut.SetMarkerStyle(ROOT.kFullCircle)
    gSignificancesDsVsCut.SetMarkerColor(colors[0])
    gSignificancesDsVsCut.SetLineColor(colors[0])
    gSignificancesDsVsCut.SetMarkerSize(1.5)
    gSignificancesDsVsCut.SetLineWidth(2)

    gSignificancesDplusVsCut.SetMarkerStyle(ROOT.kFullDiamond)
    gSignificancesDplusVsCut.SetMarkerColor(colors[1])
    gSignificancesDplusVsCut.SetLineColor(colors[1])
    gSignificancesDplusVsCut.SetMarkerSize(2)
    gSignificancesDplusVsCut.SetLineWidth(2)

    gSoverBsDsVsCut.SetMarkerStyle(ROOT.kFullCircle)
    gSoverBsDsVsCut.SetMarkerColor(colors[0])
    gSoverBsDsVsCut.SetLineColor(colors[0])
    gSoverBsDsVsCut.SetMarkerSize(1.5)
    gSoverBsDsVsCut.SetLineWidth(2)

    gSoverBsDplusVsCut.SetMarkerStyle(ROOT.kFullDiamond)
    gSoverBsDplusVsCut.SetMarkerColor(colors[1])
    gSoverBsDplusVsCut.SetLineColor(colors[1])
    gSoverBsDplusVsCut.SetMarkerSize(2)
    gSoverBsDplusVsCut.SetLineWidth(2)

    gRatioVsCut.SetMarkerStyle(ROOT.kFullCircle)
    gRatioVsCut.SetMarkerColor(ROOT.kBlack)
    gRatioVsCut.SetLineColor(ROOT.kBlack)
    gRatioVsCut.SetMarkerSize(1.5)
    gRatioVsCut.SetLineWidth(2)

    hRatioDistr.SetLineColor(ROOT.kBlack)
    hRatioDistr.SetLineWidth(2)
    hRatioDistr.SetFillStyle(0)
    
    
    # Create figures
    c = ROOT.TCanvas("c", "c", 1600, 1200)
    c.Divide(4, 2, 0.001, 0.001)

    # Raw yields
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(-1, 0.35, nFiles+1, 1.65, ";BDT selection; Raw yield/Raw yield (central)")
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleOffset(1.4)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleOffset(0.9)
    hFrame.GetXaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelOffset(0.001)
    hFrame.GetXaxis().SetLabelOffset(0.001)

    ROOT.gStyle.SetLineStyleString(11,"50 15")
    lineAtOneRY = ROOT.TLine(-1, 1, nFiles+1, 1)
    lineAtOneRY.SetLineStyle(2)
    lineAtOneRY.SetLineColor(ROOT.kGray+1)
    lineAtOneRY.SetLineWidth(2)
    lineAtOneRY.Draw("same")
    
    gRawYieldsDsVsCut.Draw("pez same")
    gRawYieldsDplusVsCut.Draw("pez same")

    legRY = ROOT.TLegend(0.2, 0.6, 0.5, 0.7)
    legRY.SetBorderSize(0)
    legRY.SetFillStyle(0)
    legRY.SetTextFont(42)
    legRY.SetTextSize(0.05)
    legRY.AddEntry(gRawYieldsDsVsCut, "D_{s}^{+}", "lp")
    legRY.AddEntry(gRawYieldsDplusVsCut, "D^{+}", "lp")
    legRY.Draw("same")

    thesisText = ROOT.TLatex(0.2, 0.89, "This Thesis")
    thesisText.SetNDC()
    thesisText.SetTextFont(42)
    thesisText.SetTextSize(0.07)
    thesisText.Draw()

    ppText = ROOT.TLatex(0.2, 0.85, "pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
    ppText.SetNDC()
    ppText.SetTextFont(42)
    ppText.SetTextSize(0.05)
    ppText.Draw()

    DecayText = ROOT.TLatex(0.2, 0.8, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
    DecayText.SetNDC()
    DecayText.SetTextFont(42)
    DecayText.SetTextSize(0.05)
    DecayText.Draw("same")

    ConjText = ROOT.TLatex(0.2, 0.765, 'and charge conjugate')
    ConjText.SetNDC()
    ConjText.SetTextFont(42)
    ConjText.SetTextSize(0.05)
    ConjText.Draw("same")

    pTText = ROOT.TLatex(0.2, 0.725, '4.5 < #it{p}_{T} < 5.0 (Ge#kern[-0.05]{V}/#it{c})')
    pTText.SetNDC()
    pTText.SetTextFont(42)
    pTText.SetTextSize(0.05)
    pTText.Draw("same")

    ROOT.gPad.RedrawAxis()


    # Efficiencies
    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(-1, 0.2, nFiles+1, 1.8, ";BDT selection; D_{s}^{+} Efficiency/D_{#lower[-0.3]{s}}^{+} Efficiency (central)")
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleOffset(1.4)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleOffset(0.9)
    hFrame.GetXaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelOffset(0.001)
    hFrame.GetXaxis().SetLabelOffset(0.001)

    ROOT.gStyle.SetLineStyleString(11,"50 15")
    lineAtOneEffDs = ROOT.TLine(-1, 1, nFiles+1, 1)
    lineAtOneEffDs.SetLineStyle(2)
    lineAtOneEffDs.SetLineColor(ROOT.kGray+1)
    lineAtOneEffDs.SetLineWidth(2)
    lineAtOneEffDs.Draw("same")
    
    gEffPromptsDsVsCut.Draw("pez same")
    gEffFDsDsVsCut.Draw("pez same")

    legEffDs = ROOT.TLegend(0.2, 0.8, 0.5, 0.9)
    legEffDs.SetBorderSize(0)
    legEffDs.SetFillStyle(0)
    legEffDs.SetTextFont(42)
    legEffDs.SetTextSize(0.05)
    legEffDs.AddEntry(gEffPromptsDsVsCut, "Prompt D_{#lower[-0.3]{s}}^{+}", "lp")
    legEffDs.AddEntry(gEffFDsDsVsCut, "Non-prompt D_{#lower[-0.3]{s}}^{+}", "lp")
    legEffDs.Draw("same")

    ROOT.gPad.RedrawAxis()


    # Efficiencies
    c.cd(3)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(-1, 0.2, nFiles+1, 1.8, ";BDT selection; D^{+} Efficiency/D^{+} Efficiency (central)")
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleOffset(1.4)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleOffset(0.9)
    hFrame.GetXaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelOffset(0.001)
    hFrame.GetXaxis().SetLabelOffset(0.001)

    ROOT.gStyle.SetLineStyleString(11,"50 15")
    lineAtOneEffDplus = ROOT.TLine(-1, 1, nFiles+1, 1)
    lineAtOneEffDplus.SetLineStyle(2)
    lineAtOneEffDplus.SetLineColor(ROOT.kGray+1)
    lineAtOneEffDplus.SetLineWidth(2)
    lineAtOneEffDplus.Draw("same")
    
    gEffPromptsDplusVsCut.Draw("pez same")
    gEffFDsDplusVsCut.Draw("pez same")

    legEffDplus = ROOT.TLegend(0.2, 0.8, 0.5, 0.9)
    legEffDplus.SetBorderSize(0)
    legEffDplus.SetFillStyle(0)
    legEffDplus.SetTextFont(42)
    legEffDplus.SetTextSize(0.05)
    legEffDplus.AddEntry(gEffPromptsDplusVsCut, "Prompt D^{+}", "lp")
    legEffDplus.AddEntry(gEffFDsDplusVsCut, "Non-prompt D^{+}", "lp")
    legEffDplus.Draw("same")

    ROOT.gPad.RedrawAxis()

    # fprompt
    c.cd(4)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(-1, 0.96, nFiles+1, 1.04, ";BDT selection; Prompt fraction/Prompt fraction (central)")
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleOffset(1.4)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleOffset(0.9)
    hFrame.GetXaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelOffset(0.001)
    hFrame.GetXaxis().SetLabelOffset(0.001)

    ROOT.gStyle.SetLineStyleString(11,"50 15")
    lineAtOnefpropmt = ROOT.TLine(-1, 1, nFiles+1, 1)
    lineAtOnefpropmt.SetLineStyle(2)
    lineAtOnefpropmt.SetLineColor(ROOT.kGray+1)
    lineAtOnefpropmt.SetLineWidth(2)
    lineAtOnefpropmt.Draw("same")
    
    gPromptFracsDsVsCut.Draw("pez same")
    gPromptFracsDplusVsCut.Draw("pez same")

    legfprompt = ROOT.TLegend(0.2, 0.8, 0.5, 0.9)
    legfprompt.SetBorderSize(0)
    legfprompt.SetFillStyle(0)
    legfprompt.SetTextFont(42)
    legfprompt.SetTextSize(0.05)
    legfprompt.AddEntry(gPromptFracsDsVsCut, "D_{s}^{+}", "lp")
    legfprompt.AddEntry(gPromptFracsDplusVsCut, "D^{+}", "lp")
    legfprompt.Draw("same")

    ROOT.gPad.RedrawAxis()


    # Significance
    c.cd(5)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(-1, 0, nFiles+1, 85, ";BDT selection; Significance")
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleOffset(1.4)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleOffset(0.9)
    hFrame.GetXaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelOffset(0.001)
    hFrame.GetXaxis().SetLabelOffset(0.001)
    
    gSignificancesDsVsCut.Draw("pez same")
    gSignificancesDplusVsCut.Draw("pez same")

    legSignificance = ROOT.TLegend(0.2, 0.8, 0.5, 0.9)
    legSignificance.SetBorderSize(0)
    legSignificance.SetFillStyle(0)
    legSignificance.SetTextFont(42)
    legSignificance.SetTextSize(0.05)
    legSignificance.AddEntry(gSignificancesDsVsCut, "D_{s}^{+}", "lp")
    legSignificance.AddEntry(gSignificancesDplusVsCut, "D^{+}", "lp")
    legSignificance.Draw("same")

    ROOT.gPad.RedrawAxis()
 

    # S/B
    c.cd(6)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(-1, 0, nFiles+1, 3., ";BDT selection; S/B")
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleOffset(1.4)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleOffset(0.9)
    hFrame.GetXaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelOffset(0.001)
    hFrame.GetXaxis().SetLabelOffset(0.001)
    
    gSoverBsDsVsCut.Draw("pez same")
    gSoverBsDplusVsCut.Draw("pez same")

    legSoverB = ROOT.TLegend(0.2, 0.8, 0.5, 0.9)
    legSoverB.SetBorderSize(0)
    legSoverB.SetFillStyle(0)
    legSoverB.SetTextFont(42)
    legSoverB.SetTextSize(0.05)
    legSoverB.AddEntry(gSoverBsDsVsCut, "D_{s}^{+}", "lp")
    legSoverB.AddEntry(gSoverBsDplusVsCut, "D^{+}", "lp")
    legSoverB.Draw("same")

    ROOT.gPad.RedrawAxis()


    # Ratio
    c.cd(7)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(-1, 0.42, nFiles+1, 0.54, ";BDT selection; D_{s}^{+}/D^{+} ratio")
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleOffset(1.4)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleOffset(0.9)
    hFrame.GetXaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelOffset(0.001)
    hFrame.GetXaxis().SetLabelOffset(0.001)

    hRMS = ROOT.TH1D("hRMS", "hRMS", 1, -1, nFiles+1)
    hRMS.SetBinContent(1, gRatioVsCut.GetPointY(0))
    hRMS.SetBinError(1, gRatioVsCut.GetRMS(2))
   
    # Convert hRMS to TGraphError
    gRMS = ROOT.TGraphErrors(hRMS)
    gRMS.SetLineColor(colors[3])
    gRMS.SetLineWidth(2)
    gRMS.SetFillColorAlpha(colors[3], 0.5)
    gRMS.SetFillStyle(1001)
    
    gRMS.Draw("5 same")

    ROOT.gStyle.SetLineStyleString(11,"50 15")
    lineRatio = ROOT.TLine(-1, gRatioVsCut.GetPointY(0), nFiles+1, gRatioVsCut.GetPointY(0))
    lineRatio.SetLineStyle(1)
    lineRatio.SetLineColor(colors[3])
    lineRatio.SetLineWidth(2)
    lineRatio.Draw("same")
    
    gRatioVsCut.Draw("pez same")

    legRatioVsCut = ROOT.TLegend(0.2, 0.8, 0.5, 0.9)
    legRatioVsCut.SetBorderSize(0)
    legRatioVsCut.SetFillStyle(0)
    legRatioVsCut.SetTextFont(42)
    legRatioVsCut.SetTextSize(0.05)
    legRatioVsCut.AddEntry(gRatioVsCut, "Data", "lp")
    legRatioVsCut.AddEntry(gRMS, "RMS", "f")
    legRatioVsCut.Draw("same")

    ROOT.gPad.RedrawAxis()



    # Ratio histogram
    c.cd(8)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(0.85, 0, 1.05, 12, ";(D_{s}^{+}/D^{+} ratio) / (D_{s}^{+}/D^{+} ratio (central));Counts")
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetTitleOffset(1.4)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleOffset(0.9)
    hFrame.GetXaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetYaxis().SetLabelOffset(0.001)
    hFrame.GetXaxis().SetLabelOffset(0.001)
   
    gRMS2 = ROOT.TGraphErrors()
    gRMS2.SetPoint(0, 1, 6)
    gRMS2.SetPointError(0, np.sqrt(gRatioVsCut.GetRMS(2)**2+(gRatioVsCut.GetMean(2)-gRatioVsCut.GetPointY(0))**2)/gRatioVsCut.GetPointY(0), 6)
    gRMS2.SetLineColor(colors[3])
    gRMS2.SetLineWidth(2)
    gRMS2.SetFillColorAlpha(colors[3], 0.5)
    gRMS2.SetFillStyle(1001)
    
    gRMS2.Draw("5 same")

    # Convert hRMS to TGraphError
    gSyst = ROOT.TGraphErrors()
    gSyst.SetPoint(0, 1, 6)
    gSyst.SetPointError(0, 0.02, 6)
    gSyst.SetLineColor(colors[0])
    gSyst.SetLineWidth(2)
    gSyst.SetFillColorAlpha(colors[0], 0.3)
    gSyst.SetFillStyle(1001)
    
    gSyst.Draw("5 same")
    
    hRatioDistr.Draw("hist same")

    legRatio = ROOT.TLegend(0.2, 0.75, 0.5, 0.9)
    legRatio.SetBorderSize(0)
    legRatio.SetFillStyle(0)
    legRatio.SetTextFont(42)
    legRatio.SetTextSize(0.05)
    legRatio.AddEntry(hRatioDistr, "Data", "l")
    legRatio.AddEntry(gRMS, "#sqrt{RMS^{2}+#Delta^{2}}", "f")
    legRatio.AddEntry(gSyst, "Assigned syst.", "f")
    legRatio.Draw("same")

    ROOT.gPad.RedrawAxis()


    c.SaveAs(f"/home/fchinu/Run3/Ds_pp_13TeV/Figures/Systematics/BDT/BDTsyst.pdf")
    c.SaveAs(f"/home/fchinu/Run3/Ds_pp_13TeV/Figures/Systematics/BDT/BDTsyst.png")

if __name__ == "__main__":
    args = argparse.ArgumentParser(description='Plot cut variations for one pt bin')
    args.add_argument('config', type=str, help='Path to the configuration file')
    args = args.parse_args()

    PlotCutVariationsOnePtBin(args.config)