from ROOT import TFile

fileDs = TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/MC/RawYields_DsMC.root")
histoSigmaDs = fileDs.Get("hRawYieldsSigma")
histoSigmaDs.SetDirectory(0)
fileDs.Close()

fileDplus = TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/MC/RawYields_DplusMC.root")
histoSigmaDplus = fileDplus.Get("hRawYieldsSigma")
histoSigmaDplus.SetDirectory(0)
fileDplus.Close()

outputFile = TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/SigmaDsOverDplusMC.root", "RECREATE")
histoRatio = histoSigmaDs.Clone("hSigmaDsOverDplus")
histoRatio.Divide(histoSigmaDplus)
outputFile.cd()
histoRatio.Write()
outputFile.Close()