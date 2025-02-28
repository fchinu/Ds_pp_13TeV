import ROOT 
import argparse

def get_raw_prompt_fraction(corr_yields_p, corr_yields_np, effacc_p, effacc_np):
    """
    Helper function to get the raw prompt fraction given the efficiencies

    Parameters
    -----------------------------------------------------
    - effacc_p: float
        eff x acc for prompt signal
    - effacc_np: float
        eff x acc for non-prompt signal

    Returns
    -----------------------------------------------------
    - f_p, f_p_unc: (float, float)
        raw prompt fraction with its uncertainty
    """

    rawy_p = effacc_p * corr_yields_p
    rawy_np = effacc_np * corr_yields_np
    f_p = rawy_p / (rawy_p + rawy_np)

    return f_p

parser = argparse.ArgumentParser()
parser.add_argument("EffFile", type=str, help="Efficiency file path")
parser.add_argument("DsCorrYields", type=str, help="Ds corrected yields file path")
parser.add_argument("DplusCorrYields", type=str, help="Dplus corrected yields file path")
parser.add_argument("OutputFileDs", type=str, help="Ds Output file path")
parser.add_argument("OutputFileDplus", type=str, help="Ds Output file path")
args = parser.parse_args()

EffFile = ROOT.TFile(args.EffFile)
EffDsPrompt = EffFile.Get("Eff_DsPrompt")
EffDsPrompt.SetDirectory(0)
EffDsFD = EffFile.Get("Eff_DsFD")
EffDsFD.SetDirectory(0)
EffDplusPrompt = EffFile.Get("Eff_DplusPrompt")
EffDplusPrompt.SetDirectory(0)
EffDplusFD = EffFile.Get("Eff_DplusFD")
EffDplusFD.SetDirectory(0)
EffFile.Close()

DsCorrYieldsFile = ROOT.TFile(args.DsCorrYields)
hDsPrompt = DsCorrYieldsFile.Get("hCorrYieldsPrompt")
hDsPrompt.SetDirectory(0)
hDsFD = DsCorrYieldsFile.Get("hCorrYieldsNonPrompt")
hDsFD.SetDirectory(0)
DsCorrYieldsFile.Close()

DplusCorrYieldsFile = ROOT.TFile(args.DplusCorrYields)
hDplusPrompt = DplusCorrYieldsFile.Get("hCorrYieldsPrompt")
hDplusPrompt.SetDirectory(0)
hDplusFD = DplusCorrYieldsFile.Get("hCorrYieldsNonPrompt")
hDplusFD.SetDirectory(0)
DplusCorrYieldsFile.Close()

hDsFrac = hDsPrompt.Clone("hDsFrac")
hDplusFrac = hDplusPrompt.Clone("hDplusFrac")

for iPt in range(hDsFrac.GetNbinsX()):
    effacc_p = EffDsPrompt.GetBinContent(iPt+1)
    effacc_np = EffDsFD.GetBinContent(iPt+1)
    corr_yields_p = hDsPrompt.GetBinContent(iPt+1)
    corr_yields_np = hDsFD.GetBinContent(iPt+1)
    f_p = get_raw_prompt_fraction(corr_yields_p, corr_yields_np, effacc_p, effacc_np)
    hDsFrac.SetBinContent(iPt+1, f_p)
    hDsFrac.SetBinError(iPt+1, 1.e-8)

    effacc_p = EffDplusPrompt.GetBinContent(iPt+1)
    effacc_np = EffDplusFD.GetBinContent(iPt+1)
    corr_yields_p = hDplusPrompt.GetBinContent(iPt+1)
    corr_yields_np = hDplusFD.GetBinContent(iPt+1)
    f_p = get_raw_prompt_fraction(corr_yields_p, corr_yields_np, effacc_p, effacc_np)
    hDplusFrac.SetBinContent(iPt+1, f_p)
    hDplusFrac.SetBinError(iPt+1, 1.e-8)

hDsFrac.SetName("hRawFracPrompt")
hDplusFrac.SetName("hRawFracPrompt")

outputFile = ROOT.TFile(args.OutputFileDs, "RECREATE")
hDsFrac.Write()
outputFile.Close()
outputFile = ROOT.TFile(args.OutputFileDplus, "RECREATE")
hDplusFrac.Write()
outputFile.Close()


