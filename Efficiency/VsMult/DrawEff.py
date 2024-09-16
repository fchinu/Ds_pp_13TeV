import ROOT
import sys
sys.path.append("/home/fchinu/Run3/ThesisUtils")
from PlotUtils import get_discrete_matplotlib_palette

CENT_MINS = [0, 0, 1, 10, 30, 50, 70]
CENT_MAXS = [100, 1, 10, 30, 50, 70, 100]

PARTICLE_CLASSES = [("Ds", "Prompt"), ("Ds", "NonPrompt"), ("Dplus", "Prompt"), ("Dplus", "NonPrompt")]

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadLeftMargin(0.15)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadBottomMargin(0.1)

if __name__ == '__main__':
    infileVsMult = ROOT.TFile.Open('/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/VsMult/Efficiency_LHC24d3b.root')
    infileVsPt = ROOT.TFile.Open('/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/VsPt/Efficiencies_LHC24d3a.root')
    for particle, cl in PARTICLE_CLASSES:
        infileVsPt.cd()
        if cl == 'NonPrompt':
            hEffMBVsPt = infileVsPt.Get(f'Eff_{particle}FD')
        else:
            hEffMBVsPt = infileVsPt.Get(f'Eff_{particle}Prompt')
        hEffMBVsPt.SetDirectory(0)
        infileVsMult.cd()
        #hEffMB = infileVsMult.Get(f'Eff_{particle}{cl}')
        #hEffMB.SetDirectory(0)
        hEffs = {}
        hRatios = {}
        for cent_min, cent_max in zip(CENT_MINS, CENT_MAXS):
            print(f'Eff_{particle}{cl}, Cent_{cent_min}_{cent_max}')
            hEffs[f"{cent_min}_{cent_max}"] = infileVsMult.Get(f'Eff_{particle}{cl}_Cent_{cent_min}_{cent_max}')
            hEffs[f"{cent_min}_{cent_max}"].SetDirectory(0)
        colors, _ = get_discrete_matplotlib_palette("tab10")
        c = ROOT.TCanvas("c", "c", 800, 600)
        c.Divide(2, 1)
        c.cd(2).DrawFrame(0, 0.8, 24, 1.2, "Ratio to 0#font[122]{-}100%;#it{p}_{T} (GeV/#it{c});Ratio")
        c.cd(1).SetLogy()
        leg = ROOT.TLegend(0.5, 0.2, 0.9, 0.5)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.04)
        hEffMBVsPt.SetLineColor(colors[0])
        hEffMBVsPt.SetMarkerColor(colors[0])
        hEffMBVsPt.SetMarkerStyle(ROOT.kFullCircle)
        hEffMBVsPt.SetMarkerSize(0.5)
        hEffMBVsPt.Draw("hist,][")
        #hEffMB.SetLineColor(colors[1])
        #hEffMB.SetMarkerColor(colors[1])
        #hEffMB.SetMarkerStyle(ROOT.kFullCircle)
        #hEffMB.SetMarkerSize(1)
        #leg.AddEntry(hEffMB, "MB")
        #hEffMB.Draw("same")
        for i, (cent_min, cent_max) in enumerate(zip(CENT_MINS, CENT_MAXS)):
            c.cd(1)
            hEffs[f"{cent_min}_{cent_max}"].SetLineColor(colors[i+2])
            hEffs[f"{cent_min}_{cent_max}"].SetMarkerColor(colors[i+2])
            hEffs[f"{cent_min}_{cent_max}"].SetMarkerStyle(ROOT.kFullCircle)
            hEffs[f"{cent_min}_{cent_max}"].SetMarkerSize(1)
            leg.AddEntry(hEffs[f"{cent_min}_{cent_max}"], f"{cent_min}-{cent_max}%")
            hEffs[f"{cent_min}_{cent_max}"].Draw("same") 
            if i != 0:
                c.cd(2)
                hRatios[f"Ratio_{cent_min}_{cent_max}"] = hEffs[f"{cent_min}_{cent_max}"].Clone()
                hRatios[f"Ratio_{cent_min}_{cent_max}"].Divide(hEffs[f"{CENT_MINS[0]}_{CENT_MAXS[0]}"])
                hRatios[f"Ratio_{cent_min}_{cent_max}"].SetName(f"Ratio_{cent_min}_{cent_max}")
                hRatios[f"Ratio_{cent_min}_{cent_max}"].SetLineColor(colors[i+2])
                hRatios[f"Ratio_{cent_min}_{cent_max}"].SetMarkerColor(colors[i+2])
                hRatios[f"Ratio_{cent_min}_{cent_max}"].SetMarkerStyle(ROOT.kFullCircle)
                hRatios[f"Ratio_{cent_min}_{cent_max}"].SetMarkerSize(1)
                hRatios[f"Ratio_{cent_min}_{cent_max}"].Draw("same")
               
        c.cd(1)
        leg.Draw()
         

        c.SaveAs(f'/home/fchinu/Run3/Ds_pp_13TeV/Efficiency/VsMult/Eff{particle}{cl}_EffVsMult.png')

        
            