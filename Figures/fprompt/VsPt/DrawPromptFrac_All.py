import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette, set_matplotlib_palette

colors, ROOTcols = get_discrete_matplotlib_palette('tab10')
set_matplotlib_palette("coolwarm")

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

ptMins = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 12] 
ptMaxs = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 8, 12, 24]

chi2Ds = [20.71,25.28,10.33,23.49,34.38,30.13,95.12,30.94,45.62,32.90,26.56,63.20,64.99,172.59]

# Get Ds histos
hCovs = []
hEffPropmts = []
hEffNonPrompts = []
hPromptFracs = []
hNonPromptFracs = []
hRawDs = []
hRawPrompts = []
hRawNonPrompts = []
hRawSums = []

infileDs = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root")

for ptMin, ptMax in zip(ptMins, ptMaxs):
    # Correlation matrix
    hCovs.append(infileDs.Get("hCorrMatrixCutSets_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hCovs[-1].SetDirectory(0)

    # Efficiency
    hEffPropmts.append(infileDs.Get("hEffPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hEffPropmts[-1].SetDirectory(0)
    hEffNonPrompts.append(infileDs.Get("hEffNonPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hEffNonPrompts[-1].SetDirectory(0)


    # Prompt fraction
    hPromptFracs.append(infileDs.Get("hFracPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hPromptFracs[-1].SetDirectory(0)
    hNonPromptFracs.append(infileDs.Get("hFracNonPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hNonPromptFracs[-1].SetDirectory(0)

    # Raw yields
    hRawDs.append(infileDs.Get("hRawYieldVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hRawDs[-1].SetDirectory(0)
    hRawPrompts.append(infileDs.Get("hRawYieldPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hRawPrompts[-1].SetDirectory(0)
    hRawNonPrompts.append(infileDs.Get("hRawYieldNonPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hRawNonPrompts[-1].SetDirectory(0)
    hRawSums.append(infileDs.Get("hRawYieldSumVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hRawSums[-1].SetDirectory(0)

    # Set style

    # Efficiency
    hEffPropmts[-1].SetLineColor(colors[3])
    hEffPropmts[-1].SetLineWidth(2)
    hEffPropmts[-1].SetMarkerColor(colors[3])
    hEffPropmts[-1].SetMarkerStyle(ROOT.kFullCircle)
    hEffPropmts[-1].SetMarkerSize(1)

    hEffNonPrompts[-1].SetLineColor(colors[0])
    hEffNonPrompts[-1].SetLineWidth(2)
    hEffNonPrompts[-1].SetMarkerColor(colors[0])
    hEffNonPrompts[-1].SetMarkerStyle(ROOT.kFullCircle)
    hEffNonPrompts[-1].SetMarkerSize(1)

    # Prompt fraction
    hPromptFracs[-1].SetLineColor(colors[3])
    hPromptFracs[-1].SetLineWidth(2)
    hPromptFracs[-1].SetMarkerColor(colors[3])
    hPromptFracs[-1].SetMarkerStyle(ROOT.kFullCircle)
    hPromptFracs[-1].SetMarkerSize(1)

    hNonPromptFracs[-1].SetLineColor(colors[0])
    hNonPromptFracs[-1].SetLineWidth(2)
    hNonPromptFracs[-1].SetMarkerColor(colors[0])
    hNonPromptFracs[-1].SetMarkerStyle(ROOT.kFullCircle)
    hNonPromptFracs[-1].SetMarkerSize(1)

    # Raw yields
    hRawDs[-1].SetLineColor(ROOT.kBlack)
    hRawDs[-1].SetMarkerColor(ROOT.kBlack) 
    hRawDs[-1].SetMarkerStyle(ROOT.kFullCircle)
    hRawDs[-1].SetMarkerSize(1)
    hRawPrompts[-1].SetLineColor(colors[3])
    hRawPrompts[-1].SetLineWidth(2)
    hRawPrompts[-1].SetMarkerColor(colors[3])
    hRawPrompts[-1].SetMarkerStyle(ROOT.kFullCircle)
    hRawPrompts[-1].SetMarkerSize(1)
    hRawPrompts[-1].SetFillStyle(1001)
    hRawPrompts[-1].SetFillColorAlpha(colors[3],0.5)
    hRawNonPrompts[-1].SetLineColor(colors[0])
    hRawNonPrompts[-1].SetLineWidth(2)
    hRawNonPrompts[-1].SetMarkerColor(colors[0])
    hRawNonPrompts[-1].SetMarkerStyle(ROOT.kFullCircle)
    hRawNonPrompts[-1].SetMarkerSize(1)
    hRawNonPrompts[-1].SetFillStyle(1001)
    hRawNonPrompts[-1].SetFillColorAlpha(colors[0],0.5)
    hRawSums[-1].SetLineColor(colors[2])
    hRawSums[-1].SetLineWidth(2)
    hRawSums[-1].SetMarkerColor(colors[2])
    hRawSums[-1].SetMarkerStyle(ROOT.kFullDiamond)
    hRawSums[-1].SetMarkerSize(1)

infileDs.Close()


for idx, (ptMin, ptMax, hCov, hEffPropmtDs, hEffNonPromptDs, hPromptFracDs, hNonPromptFracDs, hRawDs, hRawPromptDs, hRawNonPromptDs, hRawSumDs) in enumerate(zip(ptMins, ptMaxs, hCovs, hEffPropmts, hEffNonPrompts, hPromptFracs, hNonPromptFracs, hRawDs, hRawPrompts, hRawNonPrompts, hRawSums)):
    c = ROOT.TCanvas("c","c",800,600)
    c.Divide(2,2)

    c.cd(1)

    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)

    hFrame = ROOT.gPad.DrawFrame(-0.5,-0.5,16.5,16.5,";BDT selection;BDT selection;#rho")
    hFrame.GetYaxis().SetTitleOffset(1.2)
    hFrame.GetXaxis().SetTitleOffset(1.)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetXaxis().SetLabelSize(0.045)

    #hCov.GetZaxis().SetRangeUser(0,1)
    hCov.GetZaxis().SetTitle("#rho")
    hCov.GetZaxis().SetTitleSize(0.05)
    hCov.GetZaxis().SetTitleOffset(0.6)
    hCov.Draw("colz, same")

    ROOT.gPad.RedrawAxis()

    c.cd(2)

    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.05)

    hFrame = ROOT.gPad.DrawFrame(-0.5,hEffPropmtDs.GetMinimum()/10,16.5,hEffNonPromptDs.GetMaximum()*500,";BDT selection;Acceptance #times Efficiency")
    hFrame.GetYaxis().SetTitleOffset(1.25)
    hFrame.GetXaxis().SetTitleOffset(1.)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetXaxis().SetLabelSize(0.045)
    ROOT.gPad.SetLogy()

    hEffPropmtDs.Draw("same")
    hEffNonPromptDs.Draw("same")

    legEff = ROOT.TLegend(0.15,0.15,0.45,0.3)
    legEff.SetBorderSize(0)
    legEff.SetFillStyle(0)
    legEff.SetTextSize(0.04)
    legEff.AddEntry(hEffPropmtDs,"Prompt D_{#lower[-0.2]{s}}^{+}","lp")
    legEff.AddEntry(hEffNonPromptDs,"Non-prompt D_{#lower[-0.2]{s}}^{+}","lp")
    legEff.Draw("same")

    thesisText = ROOT.TLatex(0.5, 0.865, "This Thesis")
    thesisText.SetNDC()
    thesisText.SetTextFont(42)
    thesisText.SetTextSize(0.07)
    thesisText.Draw()

    ppText = ROOT.TLatex(0.5, 0.8, "pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
    ppText.SetNDC()
    ppText.SetTextFont(42)
    ppText.SetTextSize(0.05)
    ppText.Draw()

    DecayText = ROOT.TLatex(0.5, 0.72, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
    DecayText.SetNDC()
    DecayText.SetTextFont(42)
    DecayText.SetTextSize(0.05)
    DecayText.Draw("same")

    ConjText = ROOT.TLatex(0.5, 0.67, 'and charge conjugate')
    ConjText.SetNDC()
    ConjText.SetTextFont(42)
    ConjText.SetTextSize(0.05)
    ConjText.Draw("same")
    
    ptText = ROOT.TLatex(0.5, 0.6, f'{ptMin} < #it{{p}}_{{T}} < {ptMax} Ge#kern[-0.03]{{V}}/#it{{c}}')
    ptText.SetNDC()
    ptText.SetTextFont(42)
    ptText.SetTextSize(0.05)
    ptText.Draw("same")

    ROOT.gPad.RedrawAxis()

    c.cd(3)

    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.05)

    hFrame = ROOT.gPad.DrawFrame(-0.5,0,16.5,1,";BDT selection;Fraction")
    hFrame.GetYaxis().SetTitleOffset(1.25)
    hFrame.GetXaxis().SetTitleOffset(1.)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetXaxis().SetLabelSize(0.045)

    hPromptFracDs.Draw("same")
    hNonPromptFracDs.Draw("same")

    legFrac = ROOT.TLegend(0.7,0.45,0.9,0.6)
    legFrac.SetBorderSize(0)
    legFrac.SetFillStyle(0)
    legFrac.SetTextSize(0.04)
    legFrac.AddEntry(hPromptFracDs,"Prompt D_{#lower[-0.2]{s}}^{+}","lp")
    legFrac.AddEntry(hNonPromptFracDs,"Non-prompt D_{#lower[-0.2]{s}}^{+}","lp")
    legFrac.Draw("same")

    ROOT.gPad.RedrawAxis()

    c.cd(4)

    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.05)

    hFrame = ROOT.gPad.DrawFrame(-0.5,0,16.5,hRawDs.GetMaximum()*1.2,";BDT selection;Raw yield")
    hFrame.GetYaxis().SetTitleOffset(1.25)
    hFrame.GetXaxis().SetTitleOffset(1.)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetXaxis().SetLabelSize(0.045)

    hRawDs.Draw("same")
    hRawPromptDs.Draw("hist, same")
    hRawNonPromptDs.Draw("hist, same")
    hRawSumDs.Draw("hist, same")

    chi2text = ROOT.TLatex(0.515,0.85,f"#chi^{2}/ndf = {chi2Ds[idx]}/15")
    chi2text.SetNDC()
    chi2text.SetTextFont(42)
    chi2text.SetTextSize(0.041)
    chi2text.Draw("same")

    legRY = ROOT.TLegend(0.5,0.6,0.8,0.84)
    legRY.SetBorderSize(0)
    legRY.SetFillStyle(0)
    legRY.SetTextSize(0.04)
    legRY.AddEntry(hRawDs,"Data","lp")
    legRY.AddEntry(hRawPromptDs,"Prompt D_{#lower[-0.2]{s}}^{+}","f")
    legRY.AddEntry(hRawNonPromptDs,"Non-prompt D_{#lower[-0.2]{s}}^{+}","f")
    legRY.AddEntry(hRawSumDs,"Prompt + Non-prompt D_{#lower[-0.2]{s}}^{+}","l")
    legRY.Draw("same")


    ROOT.gPad.RedrawAxis()

    c.SaveAs(f"/home/fchinu/Run3/Ds_pp_13TeV/Figures/fprompt/All_pt/Ds/DsPromptFrac{ptMin*10:.0f}_{ptMax*10:.0f}.pdf")


chi2Ds = [14.61,15.57,332.54,18.30,4.37,2.41,16.59,37.71,15.30,7.67,15.15,24.37,18.03,9.16]

# Get DPlus histos
hCovs = []
hEffPropmts = []
hEffNonPrompts = []
hPromptFracs = []
hNonPromptFracs = []
hRawDs = []
hRawPrompts = []
hRawNonPrompts = []
hRawSums = []

infileDs = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root")

for ptMin, ptMax in zip(ptMins, ptMaxs):
    # Correlation matrix
    hCovs.append(infileDs.Get("hCorrMatrixCutSets_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hCovs[-1].SetDirectory(0)

    # Efficiency
    hEffPropmts.append(infileDs.Get("hEffPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hEffPropmts[-1].SetDirectory(0)
    hEffNonPrompts.append(infileDs.Get("hEffNonPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hEffNonPrompts[-1].SetDirectory(0)


    # Prompt fraction
    hPromptFracs.append(infileDs.Get("hFracPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hPromptFracs[-1].SetDirectory(0)
    hNonPromptFracs.append(infileDs.Get("hFracNonPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hNonPromptFracs[-1].SetDirectory(0)

    # Raw yields
    hRawDs.append(infileDs.Get("hRawYieldVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hRawDs[-1].SetDirectory(0)
    hRawPrompts.append(infileDs.Get("hRawYieldPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hRawPrompts[-1].SetDirectory(0)
    hRawNonPrompts.append(infileDs.Get("hRawYieldNonPromptVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hRawNonPrompts[-1].SetDirectory(0)
    hRawSums.append(infileDs.Get("hRawYieldSumVsCut_pt{:.0f}_{:.0f}".format(ptMin*10, ptMax*10)))
    hRawSums[-1].SetDirectory(0)

    # Set style

    # Efficiency
    hEffPropmts[-1].SetLineColor(colors[3])
    hEffPropmts[-1].SetLineWidth(2)
    hEffPropmts[-1].SetMarkerColor(colors[3])
    hEffPropmts[-1].SetMarkerStyle(ROOT.kFullCircle)
    hEffPropmts[-1].SetMarkerSize(1)

    hEffNonPrompts[-1].SetLineColor(colors[0])
    hEffNonPrompts[-1].SetLineWidth(2)
    hEffNonPrompts[-1].SetMarkerColor(colors[0])
    hEffNonPrompts[-1].SetMarkerStyle(ROOT.kFullCircle)
    hEffNonPrompts[-1].SetMarkerSize(1)

    # Prompt fraction
    hPromptFracs[-1].SetLineColor(colors[3])
    hPromptFracs[-1].SetLineWidth(2)
    hPromptFracs[-1].SetMarkerColor(colors[3])
    hPromptFracs[-1].SetMarkerStyle(ROOT.kFullCircle)
    hPromptFracs[-1].SetMarkerSize(1)

    hNonPromptFracs[-1].SetLineColor(colors[0])
    hNonPromptFracs[-1].SetLineWidth(2)
    hNonPromptFracs[-1].SetMarkerColor(colors[0])
    hNonPromptFracs[-1].SetMarkerStyle(ROOT.kFullCircle)
    hNonPromptFracs[-1].SetMarkerSize(1)

    # Raw yields
    hRawDs[-1].SetLineColor(ROOT.kBlack)
    hRawDs[-1].SetMarkerColor(ROOT.kBlack) 
    hRawDs[-1].SetMarkerStyle(ROOT.kFullCircle)
    hRawDs[-1].SetMarkerSize(1)
    hRawPrompts[-1].SetLineColor(colors[3])
    hRawPrompts[-1].SetLineWidth(2)
    hRawPrompts[-1].SetMarkerColor(colors[3])
    hRawPrompts[-1].SetMarkerStyle(ROOT.kFullCircle)
    hRawPrompts[-1].SetMarkerSize(1)
    hRawPrompts[-1].SetFillStyle(1001)
    hRawPrompts[-1].SetFillColorAlpha(colors[3],0.5)
    hRawNonPrompts[-1].SetLineColor(colors[0])
    hRawNonPrompts[-1].SetLineWidth(2)
    hRawNonPrompts[-1].SetMarkerColor(colors[0])
    hRawNonPrompts[-1].SetMarkerStyle(ROOT.kFullCircle)
    hRawNonPrompts[-1].SetMarkerSize(1)
    hRawNonPrompts[-1].SetFillStyle(1001)
    hRawNonPrompts[-1].SetFillColorAlpha(colors[0],0.5)
    hRawSums[-1].SetLineColor(colors[2])
    hRawSums[-1].SetLineWidth(2)
    hRawSums[-1].SetMarkerColor(colors[2])
    hRawSums[-1].SetMarkerStyle(ROOT.kFullDiamond)
    hRawSums[-1].SetMarkerSize(1)

infileDs.Close()


for idx, (ptMin, ptMax, hCov, hEffPropmtDs, hEffNonPromptDs, hPromptFracDs, hNonPromptFracDs, hRawDs, hRawPromptDs, hRawNonPromptDs, hRawSumDs) in enumerate(zip(ptMins, ptMaxs, hCovs, hEffPropmts, hEffNonPrompts, hPromptFracs, hNonPromptFracs, hRawDs, hRawPrompts, hRawNonPrompts, hRawSums)):
    c = ROOT.TCanvas("c","c",800,600)
    c.Divide(2,2)

    c.cd(1)

    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.15)

    hFrame = ROOT.gPad.DrawFrame(-0.5,-0.5,16.5,16.5,";BDT selection;BDT selection;#rho")
    hFrame.GetYaxis().SetTitleOffset(1.2)
    hFrame.GetXaxis().SetTitleOffset(1.)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetXaxis().SetLabelSize(0.045)

    #hCov.GetZaxis().SetRangeUser(0,1)
    hCov.GetZaxis().SetTitle("#rho")
    hCov.GetZaxis().SetTitleSize(0.05)
    hCov.GetZaxis().SetTitleOffset(0.6)
    hCov.Draw("colz, same")

    ROOT.gPad.RedrawAxis()

    c.cd(2)

    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.05)

    hFrame = ROOT.gPad.DrawFrame(-0.5,hEffPropmtDs.GetMinimum()/10,16.5,hEffNonPromptDs.GetMaximum()*500,";BDT selection;Acceptance #times Efficiency")
    hFrame.GetYaxis().SetTitleOffset(1.25)
    hFrame.GetXaxis().SetTitleOffset(1.)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetXaxis().SetLabelSize(0.045)
    ROOT.gPad.SetLogy()

    hEffPropmtDs.Draw("same")
    hEffNonPromptDs.Draw("same")

    legEff = ROOT.TLegend(0.15,0.15,0.45,0.3)
    legEff.SetBorderSize(0)
    legEff.SetFillStyle(0)
    legEff.SetTextSize(0.04)
    legEff.AddEntry(hEffPropmtDs,"Prompt D^{+}","lp")
    legEff.AddEntry(hEffNonPromptDs,"Non-prompt D^{+}","lp")
    legEff.Draw("same")

    thesisText = ROOT.TLatex(0.5, 0.865, "This Thesis")
    thesisText.SetNDC()
    thesisText.SetTextFont(42)
    thesisText.SetTextSize(0.07)
    thesisText.Draw()

    ppText = ROOT.TLatex(0.5, 0.8, "pp collisions, #sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
    ppText.SetNDC()
    ppText.SetTextFont(42)
    ppText.SetTextSize(0.05)
    ppText.Draw()

    DecayText = ROOT.TLatex(0.5, 0.72, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
    DecayText.SetNDC()
    DecayText.SetTextFont(42)
    DecayText.SetTextSize(0.05)
    DecayText.Draw("same")

    ConjText = ROOT.TLatex(0.5, 0.67, 'and charge conjugate')
    ConjText.SetNDC()
    ConjText.SetTextFont(42)
    ConjText.SetTextSize(0.05)
    ConjText.Draw("same")
    
    ptText = ROOT.TLatex(0.5, 0.6, f'{ptMin} < #it{{p}}_{{T}} < {ptMax} Ge#kern[-0.03]{{V}}/#it{{c}}')
    ptText.SetNDC()
    ptText.SetTextFont(42)
    ptText.SetTextSize(0.05)
    ptText.Draw("same")

    ROOT.gPad.RedrawAxis()

    c.cd(3)

    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.05)

    hFrame = ROOT.gPad.DrawFrame(-0.5,0,16.5,1,";BDT selection;Fraction")
    hFrame.GetYaxis().SetTitleOffset(1.25)
    hFrame.GetXaxis().SetTitleOffset(1.)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetXaxis().SetLabelSize(0.045)

    hPromptFracDs.Draw("same")
    hNonPromptFracDs.Draw("same")

    legFrac = ROOT.TLegend(0.15,0.45,0.35,0.6)
    legFrac.SetBorderSize(0)
    legFrac.SetFillStyle(0)
    legFrac.SetTextSize(0.04)
    legFrac.AddEntry(hPromptFracDs,"Prompt D^{+}","lp")
    legFrac.AddEntry(hNonPromptFracDs,"Non-prompt D^{+}","lp")
    legFrac.Draw("same")

    ROOT.gPad.RedrawAxis()

    c.cd(4)

    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.05)

    hFrame = ROOT.gPad.DrawFrame(-0.5,0,16.5,hRawDs.GetMaximum()*1.2,";BDT selection;Raw yield")
    hFrame.GetYaxis().SetTitleOffset(1.25)
    hFrame.GetXaxis().SetTitleOffset(1.)
    hFrame.GetYaxis().SetTitleSize(0.05)
    hFrame.GetXaxis().SetTitleSize(0.05)
    hFrame.GetYaxis().SetLabelSize(0.045)
    hFrame.GetXaxis().SetLabelSize(0.045)

    hRawDs.Draw("same")
    hRawPromptDs.Draw("hist, same")
    hRawNonPromptDs.Draw("hist, same")
    hRawSumDs.Draw("hist, same")

    chi2text = ROOT.TLatex(0.515,0.85,f"#chi^{2}/ndf = {chi2Ds[idx]}/15")
    chi2text.SetNDC()
    chi2text.SetTextFont(42)
    chi2text.SetTextSize(0.041)
    chi2text.Draw("same")

    legRY = ROOT.TLegend(0.5,0.6,0.8,0.84)
    legRY.SetBorderSize(0)
    legRY.SetFillStyle(0)
    legRY.SetTextSize(0.04)
    legRY.AddEntry(hRawDs,"Data","lp")
    legRY.AddEntry(hRawPromptDs,"Prompt D^{+}","f")
    legRY.AddEntry(hRawNonPromptDs,"Non-prompt D^{+}","f")
    legRY.AddEntry(hRawSumDs,"Prompt + Non-prompt D^{+}","l")
    legRY.Draw("same")


    ROOT.gPad.RedrawAxis()

    c.SaveAs(f"/home/fchinu/Run3/Ds_pp_13TeV/Figures/fprompt/All_pt/Dplus/DplusPromptFrac{ptMin*10:.0f}_{ptMax*10:.0f}.pdf")
