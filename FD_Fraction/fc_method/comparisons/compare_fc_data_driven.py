import argparse
from pathlib import Path
import yaml
import ROOT
import numpy as np

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadRightMargin(0.025)
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadLeftMargin(0.13)

ROOT.gStyle.SetTitleSize(0.05, "XYZ")
ROOT.gStyle.SetLabelSize(0.05, "XYZ")
ROOT.gStyle.SetTitleOffset(1.3, "Y")
ROOT.gStyle.SetTitleOffset(1.0, "X")

def get_hist_from_graph(graph):
    """Convert a TGraphErrors to a TH1D histogram."""
    n_points = graph.GetN()
    edges = [graph.GetX()[i] - graph.GetErrorX(i) for i in range(n_points)]
    edges += [graph.GetX()[n_points - 1] + graph.GetErrorX(n_points - 1)]
    edges = np.asarray(edges, "d")
    hist = ROOT.TH1D(
        "hist_from_graph",
        "hist_from_graph",
        n_points,
        edges
    )
    for i in range(n_points):
        y = graph.GetY()[i]
        hist.SetBinContent(i + 1, y)
        hist.SetBinError(i + 1, 1.e-12)
    return hist


def compare(config_path):
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    files_fc = cfg["inputs"]["fc"]
    files_dd_ds = cfg["inputs"]["data_driven_ds"]
    files_dd_dp = cfg["inputs"]["data_driven_dplus"]
    cent_mins = cfg["inputs"]["cent"]["mins"]
    cent_maxs = cfg["inputs"]["cent"]["maxs"]
    
    
    c_ds = ROOT.TCanvas("c_ds", "c_ds", 2400, 2400)
    c_ds.Divide(3, 3, 0.01, 0.01)
    
    cent_text = ROOT.TLatex()
    cent_text.SetNDC()
    cent_text.SetTextSize(0.05)
    cent_text.SetTextFont(42)

    histos_fc = []
    legs_ds = []
    for i_cent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
        c_ds.cd(i_cent + 1).DrawFrame(0, 0, 24, 1, "; #it{p}_{T} (GeV/#it{c}); D_{s}^{+} prompt fraction")
        # cent_text.DrawLatex(0.6, 0.5, f"{cent_min}#minus{cent_max}%")
        try:
            with ROOT.TFile.Open(str(files_dd_ds[i_cent])) as f_dd:
                h_dd = f_dd.Get(f"hRawFracPrompt_cent_{cent_min}_{cent_max}")
                if (cent_min, cent_max) == (80, 90) or (cent_min, cent_max) == (70, 80):
                    # exclude last point
                    h_dd.SetBinContent(h_dd.GetNbinsX(), 1.e12)

                h_dd.SetMarkerStyle(ROOT.kOpenCircle)
                h_dd.SetMarkerColor(ROOT.kRed)
                h_dd.SetLineColor(ROOT.kRed)
            with ROOT.TFile.Open(str(files_fc[i_cent])) as f_fc:
                g_fc = f_fc.Get("DsPrompt_RawFraction")
                g_fc.SetMarkerStyle(ROOT.kFullCircle)
                g_fc.SetMarkerColor(ROOT.kAzure - 3)
                g_fc.SetLineColor(ROOT.kAzure - 3)
                g_fc.SetFillColorAlpha(ROOT.kAzure - 3, 0.3)
                g_fc.SetFillStyle(1001)
            histos_fc.append(get_hist_from_graph(g_fc))
            histos_fc[-1].Draw("E SAME")

            g_fc.Draw("E2 SAME")
            h_dd.Draw("SAME")

            legs_ds.append(ROOT.TLegend(0.4, 0.3, 0.8, 0.45))
            legs_ds[-1].SetBorderSize(0)
            legs_ds[-1].SetFillStyle(0)
            legs_ds[-1].SetTextSize(0.05)
            legs_ds[-1].SetHeader(f"{cent_min}#minus{cent_max}%")
            legs_ds[-1].AddEntry(h_dd, f"Data-driven", "pl")
            legs_ds[-1].AddEntry(g_fc, f"FC method", "lf")
            legs_ds[-1].Draw()

        except Exception as e:
            print(f"Could not process cent {cent_min}-{cent_max}%: {e}")
            continue
    
    c_dp = ROOT.TCanvas("c_dp", "c_dp", 2400, 2400)
    c_dp.Divide(3, 3, 0.01, 0.01)
    
    cent_text = ROOT.TLatex()
    cent_text.SetNDC()
    cent_text.SetTextSize(0.05)
    cent_text.SetTextFont(42)

    histos_fc = []
    legs_dp = []
    for i_cent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
        c_dp.cd(i_cent + 1).DrawFrame(0, 0, 24, 1, "; #it{p}_{T} (GeV/#it{c}); D^{+} prompt fraction")
        # cent_text.DrawLatex(0.6, 0.5, f"{cent_min}#minus{cent_max}%")
        # try:
        with ROOT.TFile.Open(str(files_dd_dp[i_cent])) as f_dd:
            h_dd = f_dd.Get(f"hRawFracPrompt_cent_{cent_min}_{cent_max}")
            if (cent_min, cent_max) == (80, 90) or (cent_min, cent_max) == (70, 80):
                # exclude last point
                h_dd.SetBinContent(h_dd.GetNbinsX(), 1.e12)
            h_dd.SetMarkerStyle(ROOT.kOpenCircle)
            h_dd.SetMarkerColor(ROOT.kRed)
            h_dd.SetLineColor(ROOT.kRed)
        with ROOT.TFile.Open(str(files_fc[i_cent])) as f_fc:
            g_fc = f_fc.Get("DplusPrompt_RawFraction")
            g_fc.SetMarkerStyle(ROOT.kFullCircle)
            g_fc.SetMarkerColor(ROOT.kAzure - 3)
            g_fc.SetLineColor(ROOT.kAzure - 3)
            g_fc.SetFillColorAlpha(ROOT.kAzure - 3, 0.3)
            g_fc.SetFillStyle(1001)
        histos_fc.append(get_hist_from_graph(g_fc))
        histos_fc[-1].Draw("E SAME")

        g_fc.Draw("E2 SAME")
        h_dd.Draw("SAME")

        legs_dp.append(ROOT.TLegend(0.4, 0.3, 0.8, 0.45))
        legs_dp[-1].SetBorderSize(0)
        legs_dp[-1].SetFillStyle(0)
        legs_dp[-1].SetTextSize(0.05)
        legs_dp[-1].SetHeader(f"{cent_min}#minus{cent_max}%")
        legs_dp[-1].AddEntry(h_dd, f"Data-driven", "pl")
        legs_dp[-1].AddEntry(g_fc, f"FC method", "lf")
        legs_dp[-1].Draw()

    c_ds.SaveAs(f"{cfg['output_dir']}/compare_fc_data_driven_ds.pdf")
    c_dp.SaveAs(f"{cfg['output_dir']}/compare_fc_data_driven_dp.pdf")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare FC data-driven results with data-driven one"
    )
    parser.add_argument(
        "config",
        type=str,
        help="Path to the configuration file",
    )
    args = parser.parse_args()

    compare(args.config)