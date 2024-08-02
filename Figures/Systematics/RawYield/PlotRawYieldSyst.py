import ROOT 
import pickle
import yaml
import numpy as np
import matplotlib as mpl

def get_discrete_matplotlib_palette(paletteName):
    cmap = mpl.colormaps[paletteName]
    colors = cmap.colors
    ROOTColorIndices = []
    ROOTColors = []
    for color in colors:
        idx = ROOT.TColor.GetFreeColorIndex()
        ROOTColors.append(ROOT.TColor(idx, color[0], color[1], color[2],"color%i" % idx))
        ROOTColorIndices.append(idx)
        
    return ROOTColorIndices, ROOTColors


def ProduceFigure(multiTrialDict, multiTrialCfg, ptMin, ptMax):
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas("canvas", "canvas", 2000, 1500)
    canvas.Divide(2, 2, 0.0001, 0.0001)

    # Colors for plotting
    colors, _ = get_discrete_matplotlib_palette('tab10')

    multiTrialCfg['bincounting']['nsigma'] = multiTrialCfg['bincounting']['nsigma'][:1] # Only 3 sigma for now

    # First pad: Raw Yields
    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.1)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.1)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(0, 2500, len(multiTrialDict["convergedTrials"])*2 + 1, 7000, "Raw Yields;Trial;Raw yield")
    hFrame.GetYaxis().SetLabelSize(0.035)
    hFrame.GetYaxis().SetTitleSize(0.04)
    hFrame.GetYaxis().SetTitleOffset(1.2)
    hFrame.GetXaxis().SetLabelSize(0.035)
    hFrame.GetXaxis().SetTitleSize(0.04)
    legend1 = ROOT.TLegend(0.3, 0.42, 0.9, 0.57)
    legend1.SetBorderSize(0)
    legend1.SetFillStyle(0)
    legend1.SetTextSize(0.04)
    legend1.SetNColumns(2)
    legend1.SetMargin(0.1)
    
    graph1 = ROOT.TGraphErrors(len(multiTrialDict["rawyieldsDs"]))
    graph2 = ROOT.TGraphErrors(len(multiTrialDict["rawyieldsDplus"]))
    for i in range(len(multiTrialDict["rawyieldsDs"])):
        graph1.SetPoint(i, i+1, multiTrialDict["rawyieldsDs"][i])
        graph1.SetPointError(i, 0, multiTrialDict["rawyieldsDs_err"][i])
        graph2.SetPoint(i, i+1, multiTrialDict["rawyieldsDplus"][i])
        graph2.SetPointError(i, 0, multiTrialDict["rawyieldsDplus_err"][i])
    
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerColor(colors[0])
    graph1.SetTitle(f'Raw Yields')
    graph1.GetXaxis().SetTitle('Trial')
    graph1.GetYaxis().SetTitle('Raw yield')
    graph1.Draw('PZ same')
    graph2.SetMarkerStyle(20)
    graph2.SetMarkerColor(colors[1])
    graph2.Draw('PZ same')
    
    legend1.AddEntry(graph1, "D_{s}^{+}, Fit", "lp")
    legend1.AddEntry(graph2, "D^{+}, Fit", "lp")

    rawyieldsBinCountDs = []
    rawyieldsBinCountDplus = []

    for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
        rawyieldsBinCountDs.append(ROOT.TGraphErrors(len(multiTrialDict["binCountDs"][i])))
        rawyieldsBinCountDplus.append(ROOT.TGraphErrors(len(multiTrialDict["binCountDplus"][i])))
        for j in range(len(multiTrialDict["binCountDs"][i])):
            rawyieldsBinCountDs[-1].SetPoint(j, len(multiTrialDict["binCountDs"][i]) * (i + 1) + j + 1, multiTrialDict["binCountDs"][i][j])
            rawyieldsBinCountDs[-1].SetPointError(j, 0, multiTrialDict["binCountDs_err"][i][j])
            rawyieldsBinCountDplus[-1].SetPoint(j, len(multiTrialDict["binCountDs"][i]) * (i + 1) + j + 1, multiTrialDict["binCountDplus"][i][j])
            rawyieldsBinCountDplus[-1].SetPointError(j, 0, multiTrialDict["binCountDs_err"][i][j])
        
        rawyieldsBinCountDs[-1].SetMarkerStyle(20)
        rawyieldsBinCountDs[-1].SetMarkerColor(colors[i + 2])
        rawyieldsBinCountDs[-1].Draw('PZ same')
        rawyieldsBinCountDplus[-1].SetMarkerStyle(20)
        rawyieldsBinCountDplus[-1].SetMarkerColor(colors[i + 4])
        rawyieldsBinCountDplus[-1].Draw('PZ same')
        
        legend1.AddEntry(rawyieldsBinCountDs[-1], f"D_{{s}}^{{+}}, Bin counting {nsigma}#sigma", "lp")
        legend1.AddEntry(rawyieldsBinCountDplus[-1], f"D^{{+}}, Bin counting {nsigma}#sigma", "lp")
    
    legend1.Draw()

    #Central values
    hLine1 = ROOT.TLine(0, multiTrialDict["hRawYieldsDsCentral"].GetBinContent(multiTrialDict["hRawYieldsDsCentral"].FindBin(ptMin+0.05)),
                        len(multiTrialDict["binCountDs"][i]) * (len(multiTrialCfg['bincounting']['nsigma']) + 1) + 1, 
                        multiTrialDict["hRawYieldsDsCentral"].GetBinContent(multiTrialDict["hRawYieldsDsCentral"].FindBin(ptMin+0.05)))
    hLine1.SetLineColor(ROOT.kRed)
    hLine1.SetLineStyle(2)
    hLine1.SetLineWidth(2)
    hLine1.Draw()

    hLine2 = ROOT.TLine(0, multiTrialDict["hRawYieldsDplusCentral"].GetBinContent(multiTrialDict["hRawYieldsDplusCentral"].FindBin(ptMin+0.05)),
                        len(multiTrialDict["binCountDs"][i]) * (len(multiTrialCfg['bincounting']['nsigma']) + 1) + 1, 
                        multiTrialDict["hRawYieldsDplusCentral"].GetBinContent(multiTrialDict["hRawYieldsDplusCentral"].FindBin(ptMin+0.05)))
    hLine2.SetLineColor(colors[2])
    hLine2.SetLineStyle(2)
    hLine2.SetLineWidth(2)
    hLine2.Draw()

    ROOT.gPad.RedrawAxis()

    # Second pad: Ratios
    canvas.cd(4)
    ROOT.gPad.SetLeftMargin(0.1)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.1)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame2 = ROOT.gPad.DrawFrame(1.6, 0., 1.8, 20, "D_{s}^{+}/D^{+} Ratio;Ratio;Entries")
    hFrame2.GetYaxis().SetLabelSize(0.035)
    hFrame2.GetYaxis().SetTitleSize(0.04)
    hFrame2.GetYaxis().SetTitleOffset(1.2)
    hFrame2.GetXaxis().SetLabelSize(0.035)
    hFrame2.GetXaxis().SetTitleSize(0.04)
    hist1 = ROOT.TH1F("hist1", "Ds+/D+ Ratio", 30, min(multiTrialDict["ratios"]), max(multiTrialDict["ratios"]))
    for value in multiTrialDict["ratios"]:
        hist1.Fill(value)
    
    hist1.SetLineColor(ROOT.kBlack)
    hist1.SetLineWidth(2)
    hist1.SetFillColorAlpha(colors[0],0.65)
    hist1.Draw("same")

    legend2 = ROOT.TLegend(0.5, 0.6, 0.7, 0.8)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0)
    legend2.SetTextSize(0.037)
    legend2.AddEntry(hist1, "Fit", "f")
    
    histosBinCount = []
    for i, nsigma in enumerate(multiTrialCfg['bincounting']['nsigma']):
        histosBinCount.append(ROOT.TH1F(f"hist2_{i}", f"Bin Counting {nsigma} sigma", 30, min(multiTrialDict["binCountRatios"][i]), max(multiTrialDict["binCountRatios"][i])))
        for value in multiTrialDict["binCountRatios"][i]:
            histosBinCount[-1].Fill(value)
        
        histosBinCount[-1].SetLineColor(colors[i + 1])
        histosBinCount[-1].SetFillColorAlpha(colors[i + 1], 0.3)
        histosBinCount[-1].Draw("same")

        legend2.AddEntry(histosBinCount[-1], f"Bin counting {nsigma}#sigma", "f")

    legend2.Draw("same")

    infofit = ROOT.TPaveText(0.7, 0.71, 0.9, 0.79, "NDC")
    infofit.SetBorderSize(0)
    infofit.SetFillColor(0)
    infofit.SetFillStyle(0)
    infofit.SetTextFont(42)
    infofit.SetTextSize(0.037)
    infofit.AddText(f'#mu = {np.mean(multiTrialDict["ratios"]):.3f}')
    infofit.AddText(f'#sigma = {np.std(multiTrialDict["ratios"]):.3f}')
    infofit.Draw("same")

    info3sig = ROOT.TPaveText(0.7, 0.61, 0.9, 0.69, "NDC")
    info3sig.SetBorderSize(0)
    info3sig.SetFillColor(0)
    info3sig.SetFillStyle(0)
    info3sig.SetTextFont(42)
    info3sig.SetTextSize(0.037)
    info3sig.AddText(f'#mu = {np.mean(multiTrialDict["binCountRatios"][0]):.3f}')
    info3sig.AddText(f'#sigma = {np.std(multiTrialDict["binCountRatios"][0]):.3f}')
    info3sig.Draw("same")

    # Central values
    central_value = multiTrialDict["hRawYieldsDsCentral"].GetBinContent(multiTrialDict["hRawYieldsDsCentral"].FindBin(ptMin+0.05)) / multiTrialDict["hRawYieldsDplusCentral"].GetBinContent(multiTrialDict["hRawYieldsDplusCentral"].FindBin(ptMin+0.05))
    vLine = ROOT.TLine(central_value, 0, central_value, 20)
    vLine.SetLineColor(ROOT.kRed)
    vLine.SetLineStyle(2)
    vLine.SetLineWidth(2)
    vLine.Draw("same")

    ROOT.gPad.RedrawAxis()
    
    # Third pad: Width
    canvas.cd(3)
    ROOT.gPad.SetLeftMargin(0.1)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.1)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame = ROOT.gPad.DrawFrame(0, 6.5, len(multiTrialDict["convergedTrials"])*2 + 1, 9.75, "Peak widths;Trial;Width (MeV/#it{c}^{2})")
    hFrame.GetYaxis().SetLabelSize(0.035)
    hFrame.GetYaxis().SetTitleSize(0.04)
    hFrame.GetYaxis().SetTitleOffset(1.2)
    hFrame.GetXaxis().SetLabelSize(0.035)
    hFrame.GetXaxis().SetTitleSize(0.04)
    legend3 = ROOT.TLegend(0.6, 0.45, 0.8, 0.55)
    legend3.SetBorderSize(0)
    legend3.SetFillStyle(0)
    legend3.SetTextSize(0.04)

    
    graph5 = ROOT.TGraphErrors(len(multiTrialDict["sigmasDs"]))
    graph6 = ROOT.TGraphErrors(len(multiTrialDict["sigmasDplus"]))
    for i in range(len(multiTrialDict["sigmasDs"])):
        graph5.SetPoint(i, i+1, multiTrialDict["sigmasDs"][i]*1000)
        graph5.SetPointError(i, 0, multiTrialDict["sigmasDs_err"][i]*1000)
        graph6.SetPoint(i, i+1, multiTrialDict["sigmasDplus"][i]*1000)
        graph6.SetPointError(i, 0, multiTrialDict["sigmasDplus_err"][i]*1000)
    
    graph5.SetMarkerStyle(20)
    graph5.SetMarkerColor(colors[0])
    graph5.Draw('PZ same')
    graph6.SetMarkerStyle(20)
    graph6.SetMarkerColor(colors[1])
    graph6.Draw('PZ same')
    
    legend3.AddEntry(graph5, "D_{s}^{+}", "lp")
    legend3.AddEntry(graph6, "D^{+}", "lp")

    # Central values
    hLine3 = ROOT.TLine(0, multiTrialDict["hSigmaDsCentral"].GetBinContent(multiTrialDict["hSigmaDsCentral"].FindBin(ptMin+0.05))*1000,
                        len(multiTrialDict["convergedTrials"])*2 + 1, 
                        multiTrialDict["hSigmaDsCentral"].GetBinContent(multiTrialDict["hSigmaDsCentral"].FindBin(ptMin+0.05))*1000)
    hLine3.SetLineColor(colors[3])
    hLine3.SetLineStyle(2)
    hLine3.SetLineWidth(2)
    hLine3.Draw("same")

    hLine4 = ROOT.TLine(0, multiTrialDict["hSigmaDplusCentral"].GetBinContent(multiTrialDict["hSigmaDplusCentral"].FindBin(ptMin+0.05))*1000,
                        len(multiTrialDict["convergedTrials"])*2 + 1, 
                        multiTrialDict["hSigmaDplusCentral"].GetBinContent(multiTrialDict["hSigmaDplusCentral"].FindBin(ptMin+0.05))*1000)
    hLine4.SetLineColor(colors[2])
    hLine4.SetLineStyle(2)
    hLine4.SetLineWidth(2)
    hLine4.Draw("same")

    legend3.Draw("same")

    ROOT.gPad.RedrawAxis()

    # Fourth pad: Chi^2/ndf
    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.1)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.1)
    ROOT.gPad.SetBottomMargin(0.12)
    hFrame3 = ROOT.gPad.DrawFrame(0, 0, len(multiTrialDict["convergedTrials"])*2 + 1, 2, "#chi^{2}/ndf;Trial;#chi^{2}/ndf")
    hFrame3.GetYaxis().SetLabelSize(0.035)
    hFrame3.GetYaxis().SetTitleSize(0.04)
    hFrame3.GetYaxis().SetTitleOffset(1.2)
    hFrame3.GetXaxis().SetLabelSize(0.035)
    hFrame3.GetXaxis().SetTitleSize(0.04)
    graph7 = ROOT.TGraph(len(multiTrialDict["chi2s"]))
    for i in range(len(multiTrialDict["chi2s"])):
        graph7.SetPoint(i, i+1, multiTrialDict["chi2s"][i])
    
    graph7.SetMarkerStyle(20)
    graph7.SetMarkerColor(colors[0])
    graph7.Draw('P same')



    thesisText = ROOT.TLatex(0.55, 0.7, "This Thesis")
    thesisText.SetNDC()
    thesisText.SetTextFont(42)
    thesisText.SetTextSize(0.07)
    thesisText.Draw()

    ppText = ROOT.TLatex(0.55, 0.65, "pp collisions")
    ppText.SetNDC()
    ppText.SetTextFont(42)
    ppText.SetTextSize(0.05)
    ppText.Draw()

    TevText = ROOT.TLatex(0.55, 0.585, "#sqrt{#it{s}} = 13.6 Te#kern[-0.03]{V}")
    TevText.SetNDC()
    TevText.SetTextFont(42)
    TevText.SetTextSize(0.05)
    TevText.Draw()

    DecayText = ROOT.TLatex(0.55, 0.5, 'D_{s}^{+}, D^{+} #rightarrow #phi#pi^{+}#rightarrow K^{+}K^{#font[122]{-}}#pi^{+}')
    DecayText.SetNDC()
    DecayText.SetTextFont(42)
    DecayText.SetTextSize(0.05)
    DecayText.Draw("same")

    ConjText = ROOT.TLatex(0.55, 0.45, 'and charge conjugate')
    ConjText.SetNDC()
    ConjText.SetTextFont(42)
    ConjText.SetTextSize(0.05)
    ConjText.Draw("same")

    ptText = ROOT.TLatex(0.55, 0.37, '1.0 < #it{p}_{T} < 1.5 Ge#kern[-0.03]{V}/#it{c}')
    ptText.SetNDC()
    ptText.SetTextFont(42)
    ptText.SetTextSize(0.05)
    ptText.Draw("same")

    ROOT.gPad.RedrawAxis()

    canvas.SaveAs(f'/home/fchinu/Run3/Ds_pp_13TeV/Figures/Systematics/RawYield/RawYieldSyst.pdf')
    canvas.SaveAs(f'/home/fchinu/Run3/Ds_pp_13TeV/Figures/Systematics/RawYield/RawYieldSyst.png')

if __name__=='__main__':
    with open('/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/results/pt10.0_15.0.pkl', 'rb') as f:
        analysisResults = pickle.load(f)
    
    with open('/home/fchinu/Run3/Ds_pp_13TeV/Systematics/RawYields/config_multi_trial_lowpt.yml', 'r') as f:
        analysisConfig = yaml.load(f, Loader=yaml.FullLoader)

    

    ProduceFigure(analysisResults, analysisConfig['multitrial'], 1.0, 1.5)


