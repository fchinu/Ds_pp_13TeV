import ROOT 
import pdg

DPLUS_MASS_WINDOW = 0.015 # GeV/c2

def get_fraction_of_lost_signal():
    """
    Calculate the fraction of lost signal for Ds+ and D+ particles.

    This script loads a ROOT file containing histograms of DsAsPiKPi and DplusAsPiKPi,
    calculates the fraction of lost signal for each particle as a function of transverse momentum (pT),
    and saves the results in a new ROOT file.

    Returns:
        None
    """
    # Load the ROOT file
    f = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/LossForRemovingDplusMassForPiKPi/DsDplusToKKpi_pythia8_seed42.root")
    # Get the histos
    histoDsAsPiKPi = f.Get("histoDsAsPiKPi")
    histoDsAsPiKPi.SetDirectory(0)
    histoDplusAsPiKPi = f.Get("histoDplusAsPiKPi")
    histoDplusAsPiKPi.SetDirectory(0)
    # Close the ROOT file
    f.Close()

    histoDsLost = histoDsAsPiKPi.ProjectionX("histoDsLost")
    histoDsLost.Reset()
    histoDplusLost = histoDsLost.Clone("histoDplusLost")

    histoDsLost.SetTitle("Fraction of lost signal for D_{s}^{+};#it{p}_{T} (GeV/#it{c});Fraction of lost signal")
    histoDplusLost.SetTitle("Fraction of lost signal for D^{+};#it{p}_{T} (GeV/#it{c});Fraction of lost signal")

    dplusMin = pdg.connect().get_particle_by_name('D+').mass - DPLUS_MASS_WINDOW/2
    dplusMax = pdg.connect().get_particle_by_name('D+').mass + DPLUS_MASS_WINDOW/2

    for i in range(1, histoDsAsPiKPi.GetNbinsX()+1):
        hProjDs = histoDsAsPiKPi.ProjectionY("hProjDs", i, i)
        hProjDplus = histoDplusAsPiKPi.ProjectionY("hProjDplus", i, i)
        lostDs = hProjDs.Integral(hProjDs.FindBin(dplusMin), hProjDs.FindBin(dplusMax))
        lostDplus = hProjDplus.Integral(hProjDplus.FindBin(dplusMin), hProjDplus.FindBin(dplusMax))
        histoDsLost.SetBinContent(i, lostDs/hProjDs.Integral())
        histoDplusLost.SetBinContent(i, lostDplus/hProjDplus.Integral())

    outFile = ROOT.TFile("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/LossForRemovingDplusMassForPiKPi/LostSignalFrac.root", "RECREATE")
    histoDsLost.Write()
    histoDplusLost.Write()
    outFile.Close()

if __name__=="__main__":
    get_fraction_of_lost_signal()
