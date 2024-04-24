import ROOT 
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from DfUtils import ConvertGraphToHist

colors = [ROOT.kRed, ROOT.kOrange-3,ROOT.kGreen+2,ROOT.kAzure+4,ROOT.kMagenta+2]
markers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kOpenCircle, ROOT.kOpenSquare]

# DD = Data Driven, TD = Theory Driven

DD_Ds_OldMCFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB.root")
DD_Ds_OldMCFrac = DD_Ds_OldMCFile.Get("hRawFracPrompt")
DD_Ds_OldMCFrac.SetDirectory(0)
DD_Ds_OldMCFile.Close()

DD_Ds_NewMCFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_newMC.root")
DD_Ds_NewMCFrac = DD_Ds_NewMCFile.Get("hRawFracPrompt")
DD_Ds_NewMCFrac.SetDirectory(0)
DD_Ds_NewMCFile.Close()

DD_Dplus_OldMCFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB.root")
DD_Dplus_OldMCFrac = DD_Dplus_OldMCFile.Get("hRawFracPrompt")
DD_Dplus_OldMCFrac.SetDirectory(0)
DD_Dplus_OldMCFile.Close()

DD_Dplus_NewMCFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB_newMC.root")
DD_Dplus_NewMCFrac = DD_Dplus_NewMCFile.Get("hRawFracPrompt")
DD_Dplus_NewMCFrac.SetDirectory(0)
DD_Dplus_NewMCFile.Close()

TD_Ds_OldMCFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/theory_driven/FD_fraction_Ds.root")
TD_Ds_OldMCFrac = TD_Ds_OldMCFile.Get("gfraction")
TD_Ds_OldMCFile.Close()
TD_Ds_OldMCFrac = ConvertGraphToHist(TD_Ds_OldMCFrac)

TD_Ds_NewMCFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/theory_driven/FD_fraction_Ds_newMC.root")
TD_Ds_NewMCFrac = TD_Ds_NewMCFile.Get("gfraction")
TD_Ds_NewMCFile.Close()
TD_Ds_NewMCFrac = ConvertGraphToHist(TD_Ds_NewMCFrac)

TD_Dplus_OldMCFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/theory_driven/FD_fraction_Dplus.root")
TD_Dplus_OldMCFrac = TD_Dplus_OldMCFile.Get("gfraction")
TD_Dplus_OldMCFile.Close()
TD_Dplus_OldMCFrac = ConvertGraphToHist(TD_Dplus_OldMCFrac)

TD_Dplus_NewMCFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/theory_driven/FD_fraction_Dplus_newMC.root")
TD_Dplus_NewMCFrac = TD_Dplus_NewMCFile.Get("gfraction")
TD_Dplus_NewMCFile.Close()
TD_Dplus_NewMCFrac = ConvertGraphToHist(TD_Dplus_NewMCFrac)

DD_Ratio_OldMC = DD_Ds_OldMCFrac.Clone("DD_Ratio_OldMC")
DD_Ratio_OldMC.Divide(DD_Dplus_OldMCFrac)

DD_Ratio_NewMC = DD_Ds_NewMCFrac.Clone("DD_Ratio_NewMC")
DD_Ratio_NewMC.Divide(DD_Dplus_NewMCFrac)

TD_Ratio_OldMC = TD_Ds_OldMCFrac.Clone("TD_Ratio_OldMC")
TD_Ratio_OldMC.Divide(TD_Dplus_OldMCFrac)

TD_Ratio_NewMC = TD_Ds_NewMCFrac.Clone("TD_Ratio_NewMC")
TD_Ratio_NewMC.Divide(TD_Dplus_NewMCFrac)

c = ROOT.TCanvas("c", "c", 800, 600)
c.DrawFrame(0, 0, 24, 1.5, ";p_{T} (GeV/c);D_{s}^{+}/D^{+} prompt fraction ratio")
DD_Ratio_OldMC.SetLineColor(colors[0])
DD_Ratio_OldMC.SetMarkerColor(colors[0])
DD_Ratio_OldMC.SetMarkerStyle(markers[0])
DD_Ratio_OldMC.Draw("same")
DD_Ratio_NewMC.SetLineColor(colors[1])
DD_Ratio_NewMC.SetMarkerColor(colors[1])
DD_Ratio_NewMC.SetMarkerStyle(markers[1])
DD_Ratio_NewMC.Draw("same")
TD_Ratio_OldMC.SetLineColor(colors[2])
TD_Ratio_OldMC.SetMarkerColor(colors[2])
TD_Ratio_OldMC.SetMarkerStyle(markers[2])
TD_Ratio_OldMC.Draw("same")
TD_Ratio_NewMC.SetLineColor(colors[3])
TD_Ratio_NewMC.SetMarkerColor(colors[3])
TD_Ratio_NewMC.SetMarkerStyle(markers[3])
TD_Ratio_NewMC.Draw("same")

leg = ROOT.TLegend(0.4, 0.2, 0.9, 0.5)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(DD_Ratio_OldMC, "Data driven, LHC22b1[a,b]", "lp")
leg.AddEntry(DD_Ratio_NewMC, "Data driven, LHC24d3", "lp")
leg.AddEntry(TD_Ratio_OldMC, "Theory driven, LHC22b1[a,b]", "lp")
leg.AddEntry(TD_Ratio_NewMC, "Theory driven, LHC24d3", "lp")
leg.Draw()

c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/CompareFractionRatios.png")
