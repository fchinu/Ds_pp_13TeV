import ROOT 

infileDs = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Ds/CutVarDs_pp13TeV_MB_LHC24d3a.root")
hRawFracDs = infileDs.Get("hRawFracPrompt")
hRawFracDs.SetDirectory(0)
infileDs.Close()

infileDplus = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/Dplus/CutVarDplus_pp13TeV_MB_LHC24d3a.root")
hRawFracDplus = infileDplus.Get("hRawFracPrompt")
hRawFracDplus.SetDirectory(0)
infileDplus.Close()

cDs = ROOT.TCanvas("cDs", "cDs", 800, 800)
cDs.SetTopMargin(0.05)
cDs.SetRightMargin(0.05)
cDs.SetLeftMargin(0.13)
cDs.SetBottomMargin(0.13)
hRawFracDs.Draw()
cDs.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/DsPromptFraction_LHC24d3a.pdf")

cDplus = ROOT.TCanvas("cDplus", "cDplus", 800, 800)
cDplus.SetTopMargin(0.05)
cDplus.SetRightMargin(0.05)
cDplus.SetLeftMargin(0.13)
cDplus.SetBottomMargin(0.13)
hRawFracDplus.Draw()
cDplus.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/FD_Fraction/data_driven/DplusPromptFraction_LHC24d3a.pdf")

