import argparse
import sys
import yaml
import ROOT
import numpy as np
from pathlib import Path
sys.path.append(f"{Path(__file__).parent.parent.parent.parent.parent}/utils")
from plot_utils import get_discrete_matplotlib_palette

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
colors, _ = get_discrete_matplotlib_palette("tab10")

def get_ratio_graphs(graph_1, graph_2, labels):
    g_ratio = graph_1.Clone(f"g_ratio_{labels[0]}_over_{labels[1]}")
    for i in range(g_ratio.GetN()):
        g_ratio.SetPointY(i, graph_1.GetY()[i] / graph_2.GetY()[i])
        # Correlated uncertainties
        g_ratio.SetPointEYlow(i, abs(graph_1.GetErrorYlow(i) / graph_1.GetY()[i] - graph_2.GetErrorYlow(i) / graph_2.GetY()[i]) * g_ratio.GetY()[i])
        g_ratio.SetPointEYhigh(i, abs(graph_1.GetErrorYhigh(i) / graph_1.GetY()[i] - graph_2.GetErrorYhigh(i) / graph_2.GetY()[i]) * g_ratio.GetY()[i])
    return g_ratio

def systematic(config_path):
    # Load configuration
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    input_files = config["inputs"]
    labels = config["labels"]
    output_dir = config["output_dir"]
    cent_mins = config["cent"]["mins"]
    cent_maxs = config["cent"]["maxs"]
    assigned_syst = config["assigned_syst"]

    with ROOT.TFile.Open(f"{output_dir}/systematic.root", "RECREATE") as _:
        pass  # Just to create/clear the file

    for i_cent, (cent_min, cent_max, files) in enumerate(zip(cent_mins, cent_maxs, input_files)):
        g_fracs_ds = []
        g_fracs_dp = []
        for file in files:
            with ROOT.TFile.Open(file, "READ") as f:
                g_fracs_ds.append(f.Get("DsPrompt_RawFraction"))
                g_fracs_dp.append(f.Get("DplusPrompt_RawFraction"))

        g_ratios = [get_ratio_graphs(g_ds, g_dp, labels) for g_ds, g_dp in zip(g_fracs_ds, g_fracs_dp)]

        c = ROOT.TCanvas(f"c_systematic_cent_{cent_min}_{cent_max}", "", 1200, 600)
        c.Divide(3)

        # Ds prompt fraction
        h_frame_ds = c.cd(1).DrawFrame(0, 0., 24, 1., f"; #it{{p}}_{{T}} (GeV/#it{{c}}); D_{{s}}^{{+}} prompt fraction")
        leg_ds = ROOT.TLegend(0.6, 0.3, 0.9, 0.4)
        leg_ds.SetBorderSize(0)
        leg_ds.SetFillStyle(0)
        leg_ds.SetTextSize(0.04)
        for i, (g, l) in enumerate(zip(g_fracs_ds, labels)):
            g.SetMarkerStyle(ROOT.kFullCircle)
            g.SetMarkerColor(colors[i])
            g.SetLineColor(colors[i])
            g.SetFillColorAlpha(colors[i], 0.7)
            g.Draw("E2 same")
            g.DrawClone("Z same")
            leg_ds.AddEntry(g, l, "lp")
        leg_ds.Draw("same")

        # Dp prompt fraction
        h_frame_dp = c.cd(2).DrawFrame(0, 0., 24, 1., f"; #it{{p}}_{{T}} (GeV/#it{{c}}); D^{{+}} prompt fraction")
        leg_dp = ROOT.TLegend(0.6, 0.3, 0.9, 0.4)
        leg_dp.SetBorderSize(0)
        leg_dp.SetFillStyle(0)
        leg_dp.SetTextSize(0.04)
        for i, (g, l) in enumerate(zip(g_fracs_dp, labels)):
            g.SetMarkerStyle(ROOT.kFullCircle)
            g.SetMarkerColor(colors[i])
            g.SetLineColor(colors[i])
            g.SetFillColorAlpha(colors[i], 0.7)
            g.Draw("E2 same")
            g.DrawClone("Z same")
            leg_dp.AddEntry(g, l, "lp")
        leg_dp.Draw("same")

        # Ratio Ds/Dp
        h_frame_ratio = c.cd(3).DrawFrame(0, 0.5, 24, 1.5, f"; #it{{p}}_{{T}} (GeV/#it{{c}}); D_{{s}}^{{+}}/D^{{+}} prompt fraction")
        leg_ratio = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
        leg_ratio.SetBorderSize(0)
        leg_ratio.SetFillStyle(0)
        leg_ratio.SetTextSize(0.04)
        for i, g in enumerate(g_ratios):
            g.SetMarkerStyle(ROOT.kFullCircle)
            g.SetMarkerColor(colors[i])
            g.SetLineColor(colors[i])
            g.SetFillColorAlpha(colors[i], 0.7)
            g.Draw("E2 same")
            g.DrawClone("Z same")
            leg_ratio.AddEntry(g, f"{labels[i]}", "lp")
        leg_ratio.Draw("same")

        c.SaveAs(f"{output_dir}/DsPrompt_Fraction_Systematic_cent_{cent_min}_{cent_max}.pdf")

        # Save results
        pt_bins = [g_ratios[0].GetX()[i] - g_ratios[0].GetErrorXlow(i) for i in range(g_ratios[0].GetN())] 
        pt_bins.append(g_ratios[0].GetX()[g_ratios[0].GetN()-1] + g_ratios[0].GetErrorXhigh(g_ratios[0].GetN()-1))
    
        g_syst_up = g_ratios[0].Clone(f"g_syst_up_{cent_min}_{cent_max}")
        g_syst_down = g_ratios[0].Clone(f"g_syst_down_{cent_min}_{cent_max}")

        h_assigned_syst_upper = ROOT.TH1F(f"assigned_syst_upper_{cent_min}_{cent_max}", ";#it{p}_{T} (GeV/#it{c}); Assigned systematic uncertainty", len(pt_bins)-1, np.asarray(pt_bins, 'd'))
        h_assigned_syst_lower = ROOT.TH1F(f"assigned_syst_lower_{cent_min}_{cent_max}", ";#it{p}_{T} (GeV/#it{c}); Assigned systematic uncertainty", len(pt_bins)-1, np.asarray(pt_bins, 'd'))
        for i_pt in range(g_ratios[0].GetN()):
            vals_upper = [g_ratios[i].GetY()[i_pt] + g_ratios[i].GetErrorYhigh(i_pt) for i in range(len(g_ratios))]
            vals_lower = [g_ratios[i].GetY()[i_pt] - g_ratios[i].GetErrorYlow(i_pt) for i in range(len(g_ratios))]
            vals = [g_ratios[i].GetY()[i_pt] for i in range(len(g_ratios))]

            max_frac = max(vals_upper)
            min_frac = min(vals_lower)
            syst_frac_up = (max_frac - vals[0]) / vals[0]
            syst_frac_down = (vals[0] - min_frac) / vals[0]

            g_syst_up.SetPointY(i_pt, syst_frac_up)
            g_syst_up.SetPointError(i_pt, g_syst_up.GetErrorXlow(i_pt), g_syst_up.GetErrorXhigh(i_pt), 0.01, 0.01)
            g_syst_down.SetPointY(i_pt, syst_frac_down)
            g_syst_down.SetPointError(i_pt, g_syst_down.GetErrorXlow(i_pt), g_syst_down.GetErrorXhigh(i_pt), 0.01, 0.01)
            h_assigned_syst_upper.SetBinContent(i_pt+1, assigned_syst[i_pt][i_cent][0])
            h_assigned_syst_lower.SetBinContent(i_pt+1, assigned_syst[i_pt][i_cent][1])

        c = ROOT.TCanvas(f"c_systematic_assigned_cent_{cent_min}_{cent_max}", "", 600, 600)
        ROOT.gPad.SetLeftMargin(0.15)
        leg = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        h_frame = c.DrawFrame(0, 0., 24, 0.5, f"; #it{{p}}_{{T}} (GeV/#it{{c}}); Systematic uncertainty on D_{{s}}^{{+}}/D^{{+}} prompt fraction")
        g_syst_up.SetLineColor(ROOT.kRed+1)
        g_syst_up.SetMarkerColor(ROOT.kRed+1)
        g_syst_up.SetMarkerSize(1.2)
        g_syst_down.SetLineColor(ROOT.kAzure-2)
        g_syst_down.SetMarkerColor(ROOT.kAzure-2)
        g_syst_down.SetMarkerSize(1.2)
        g_syst_up.Draw("PZ same")
        g_syst_down.Draw("PZ same")

        leg.AddEntry(g_syst_up, "Upper systematic", "lp")
        leg.AddEntry(g_syst_down, "Lower systematic", "lp")
        leg.Draw("same")

        c.SaveAs(f"{output_dir}/syst_{cent_min}_{cent_max}.pdf")

        with ROOT.TFile.Open(f"{output_dir}/systematic.root", "UPDATE") as f_out:
            g_syst_up.Write()
            g_syst_down.Write()
            h_assigned_syst_upper.Write()
            h_assigned_syst_lower.Write()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate systematic uncertainties for FD fraction using fc method"
    )
    parser.add_argument(
        "config",
        type=str,
        help="Path to the configuration YAML file",
    )
    args = parser.parse_args()

    systematic(args.config)

