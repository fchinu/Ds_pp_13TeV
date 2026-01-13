import argparse
import ROOT
import yaml
import itertools

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

def systematics(config_path):
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    inputs = config["inputs"]
    labels = config["labels"]
    cent_mins = config["cent"]["mins"]
    cent_maxs = config["cent"]["maxs"]
    output_dir = config["output_dir"]
    assigned_syst = config["assigned_syst"]

    particles = ["Ds", "Dplus"]
    origins = ["Prompt", "NonPrompt"]


    with ROOT.TFile.Open(f"{output_dir}/systematic.root", "RECREATE") as output_file:
        pass  # Just to create/clear the file

    for i_cent, (cent_min, cent_max) in enumerate(zip(cent_mins, cent_maxs)):
        efficiencies = {}
        weights = {}

        for (particle, origin) in itertools.product(particles, origins):
            efficiencies[(particle, origin)] = []
            weights[(particle, origin)] = []

        for input_file in inputs:
            with ROOT.TFile.Open(input_file, "READ") as file:
                for (particle, origin) in itertools.product(particles, origins):
                    efficiencies[(particle, origin)].append(file.Get(f"eff_{particle}{origin}_cent_{cent_min}_{cent_max}"))
                    weights[(particle, origin)].append(file.Get(f"Weights/PtWeight{particle}{origin}Cands_{cent_min}_{cent_max}"))

                    if origin == "Prompt":
                        efficiencies[(particle, origin)][-1].SetMarkerStyle(ROOT.kFullCircle)
                        efficiencies[(particle, origin)][-1].SetMarkerColor(ROOT.kRed)
                        efficiencies[(particle, origin)][-1].SetLineColor(ROOT.kRed)
                    else:
                        efficiencies[(particle, origin)][-1].SetMarkerStyle(ROOT.kFullSquare)
                        efficiencies[(particle, origin)][-1].SetMarkerColor(ROOT.kAzure-3)
                        efficiencies[(particle, origin)][-1].SetLineColor(ROOT.kAzure-3)

                    if particle == "Ds":
                        weights[(particle, origin)][-1].SetMarkerStyle(ROOT.kFullCircle)
                        weights[(particle, origin)][-1].SetMarkerColor(ROOT.kGreen-2)
                        weights[(particle, origin)][-1].SetLineColor(ROOT.kGreen-2)
                    else:
                        weights[(particle, origin)][-1].SetMarkerStyle(ROOT.kFullSquare)
                        weights[(particle, origin)][-1].SetMarkerColor(ROOT.kMagenta-2)
                        weights[(particle, origin)][-1].SetLineColor(ROOT.kMagenta-2)

        
        c = ROOT.TCanvas(f"c_{cent_min}_{cent_max}", f"c_{cent_min}_{cent_max}", 2400, 1800)
        c.Divide(2, 2)
        
        # Ds efficiency
        h_frame_ds = c.cd(1).DrawFrame(0, 1.e-5, 24, 1., ";#it{p}_{T} (GeV/#it{c});D_{s}^{+} efficiency")
        ROOT.gPad.SetLogy()
        for h_eff_prompt_ds in efficiencies[("Ds", "Prompt")]:
            h_eff_prompt_ds.Draw("same ][ hist")
            h_eff_prompt_ds.DrawClone("same p")
        for h_eff_nonprompt_ds in efficiencies[("Ds", "NonPrompt")]:
            h_eff_nonprompt_ds.Draw("same ][ hist")
            h_eff_nonprompt_ds.DrawClone("same p") 

        legend_ds = ROOT.TLegend(0.4, 0.3, 0.7, 0.5)
        legend_ds.SetBorderSize(0)
        legend_ds.SetFillStyle(0)
        legend_ds.SetTextSize(0.04)
        for label, h_eff in zip(labels, efficiencies[("Ds", "Prompt")]):
            legend_ds.AddEntry(h_eff, f"Prompt - {label}", "pl")
        for label, h_eff in zip(labels, efficiencies[("Ds", "NonPrompt")]):
            legend_ds.AddEntry(h_eff, f"Non-Prompt - {label}", "pl")
        legend_ds.Draw("same")

        # Dplus efficiency
        h_frame_dplus = c.cd(2).DrawFrame(0, 1.e-5, 24, 1., ";#it{p}_{T} (GeV/#it{c});D^{+} efficiency")
        ROOT.gPad.SetLogy()
        for h_eff_prompt_dplus in efficiencies[("Dplus", "Prompt")]:
            h_eff_prompt_dplus.Draw("same ][ hist")
            h_eff_prompt_dplus.DrawClone("same p")
        for h_eff_nonprompt_dplus in efficiencies[("Dplus", "NonPrompt")]:
            h_eff_nonprompt_dplus.Draw("same ][ hist")
            h_eff_nonprompt_dplus.DrawClone("same p") 

        legend_dplus = ROOT.TLegend(0.4, 0.3, 0.7, 0.5)
        legend_dplus.SetBorderSize(0)
        legend_dplus.SetFillStyle(0)
        legend_dplus.SetTextSize(0.04)
        for label, h_eff in zip(labels, efficiencies[("Dplus", "Prompt")]):
            legend_dplus.AddEntry(h_eff, f"Prompt - {label}", "pl")
        for label, h_eff in zip(labels, efficiencies[("Dplus", "NonPrompt")]):
            legend_dplus.AddEntry(h_eff, f"Non-Prompt - {label}", "pl")
        legend_dplus.Draw("same")

        # Weights ratio
        h_frame_w_ds = c.cd(3).DrawFrame(0, 0., 24, 2., ";#it{p}_{T} (GeV/#it{c});Weights ratio")
        ratios_ds_prompt = [weights[("Ds", "Prompt")][0].Clone("ratio_ds_prompt_"+str(i)) for i in range(len(weights[("Ds", "Prompt")][1:]))]
        ratios_dp_prompt = [weights[("Dplus", "Prompt")][0].Clone("ratio_dplus_prompt_"+str(i)) for i in range(len(weights[("Dplus", "Prompt")][1:]))]
        
        legend_w = ROOT.TLegend(0.4, 0.7, 0.7, 0.8)
        legend_w.SetBorderSize(0)
        legend_w.SetFillStyle(0)
        legend_w.SetTextSize(0.04)
           
        # We only draw propmt ones
        for i, (h_w_ds, h_w_dp, l) in enumerate(zip(weights[("Ds", "Prompt")][1:], weights[("Dplus", "Prompt")][1:], labels[1:])):
            ratios_ds_prompt[i].Divide(h_w_ds)
            ratios_ds_prompt[i].Draw("same ][ hist")
            ratios_ds_prompt[i].DrawClone("same p")

            ratios_dp_prompt[i].Divide(h_w_dp)
            ratios_dp_prompt[i].Draw("same ][ hist")
            ratios_dp_prompt[i].DrawClone("same p")

            legend_w.AddEntry(h_w_ds, f"Prompt D_{{s}}^{{+}} {labels[0]}/{l}", "pl")
            legend_w.AddEntry(h_w_dp, f"Prompt D^{{+}} {labels[0]}/{l}", "pl")
        legend_w.Draw("same")

        # Efficiency ratio Ds/Dplus
        h_frame_eff_ratio = c.cd(4).DrawFrame(0, 0.8, 24, 1.2, ";#it{p}_{T} (GeV/#it{c});D_{s}^{+}/D^{+} efficiency ratio")
        ratios_eff_ds_prompt = [efficiencies[("Ds", "Prompt")][0].Clone("eff_ratio_ds_prompt_"+str(i)) for i in range(len(efficiencies[("Ds", "Prompt")][1:]))]
        ratios_eff_ds_nonprompt = [efficiencies[("Ds", "NonPrompt")][0].Clone("eff_ratio_ds_nonprompt_"+str(i)) for i in range(len(efficiencies[("Ds", "NonPrompt")][1:]))]
        ratios_eff_dp_prompt = [efficiencies[("Dplus", "Prompt")][0].Clone("eff_ratio_dp_prompt_"+str(i)) for i in range(len(efficiencies[("Dplus", "Prompt")][1:]))]
        ratios_eff_dp_nonprompt = [efficiencies[("Dplus", "NonPrompt")][0].Clone("eff_ratio_dp_nonprompt_"+str(i)) for i in range(len(efficiencies[("Dplus", "NonPrompt")][1:]))]
        ratios_syst_prompt = []
        ratios_syst_nonprompt = []

        legend_eff_ratio = ROOT.TLegend(0.4, 0.7, 0.7, 0.8)
        legend_eff_ratio.SetBorderSize(0)
        legend_eff_ratio.SetFillStyle(0)
        legend_eff_ratio.SetTextSize(0.04)

        for i, (h_eff_ds_prompt, h_eff_dp_prompt, h_eff_ds_nonprompt, h_eff_dp_nonprompt, l) in enumerate(zip(
            efficiencies[("Ds", "Prompt")][1:], 
            efficiencies[("Dplus", "Prompt")][1:], 
            efficiencies[("Ds", "NonPrompt")][1:], 
            efficiencies[("Dplus", "NonPrompt")][1:], 
            labels[1:]
        )):
            ratios_eff_ds_prompt[i].Divide(h_eff_ds_prompt)
            ratios_eff_dp_prompt[i].Divide(h_eff_dp_prompt)
            ratios_eff_ds_nonprompt[i].Divide(h_eff_ds_nonprompt)
            ratios_eff_dp_nonprompt[i].Divide(h_eff_dp_nonprompt)

            ratios_eff_ds_prompt[i].Divide(ratios_eff_dp_prompt[i])
            ratios_eff_ds_nonprompt[i].Divide(ratios_eff_dp_nonprompt[i])

            ratios_syst_prompt.append(ratios_eff_ds_prompt[i].Clone("syst_eff_ratio_prompt_"+str(i)))
            ratios_syst_nonprompt.append(ratios_eff_ds_nonprompt[i].Clone("syst_eff_ratio_nonprompt_"+str(i)))

            ratios_eff_ds_prompt[i].Draw("same ][ hist")
            ratios_eff_ds_prompt[i].DrawClone("same p")
            ratios_eff_ds_nonprompt[i].Draw("same ][ hist")
            ratios_eff_ds_nonprompt[i].DrawClone("same p")

            legend_eff_ratio.AddEntry(ratios_eff_ds_prompt[i], f"Prompt {labels[0]}/{l}", "pl")
            legend_eff_ratio.AddEntry(ratios_eff_ds_nonprompt[i], f"Non-Prompt {labels[0]}/{l}", "pl")

        legend_eff_ratio.Draw("same")

        h_assigned_syst = ratios_syst_prompt[0].Clone(f"assigned_syst_{cent_min}_{cent_max}")
        for i_pt in range(ratios_syst_prompt[0].GetNbinsX()):
            vals_prompt = [ratios_syst_prompt[i].GetBinContent(i_pt+1) for i in range(len(ratios_syst_prompt))]
            vals_nonprompt = [ratios_syst_nonprompt[i].GetBinContent(i_pt+1) for i in range(len(ratios_syst_nonprompt))]

            vals_prompt.append(1.0)  # Include the reference
            vals_nonprompt.append(1.0)  # Include the reference

            max_prompt = max(vals_prompt)
            min_prompt = min(vals_prompt)
            syst_prompt = (max_prompt - min_prompt)

            max_nonprompt = max(vals_nonprompt)
            min_nonprompt = min(vals_nonprompt)
            syst_nonprompt = (max_nonprompt - min_nonprompt)

            ratios_syst_prompt[0].SetBinContent(i_pt+1, syst_prompt)
            ratios_syst_prompt[0].SetBinError(i_pt+1, 0)
            ratios_syst_nonprompt[0].SetBinContent(i_pt+1, syst_nonprompt)
            ratios_syst_nonprompt[0].SetBinError(i_pt+1, 0)
            h_assigned_syst.SetBinContent(i_pt+1, assigned_syst[i_pt][i_cent])

        with ROOT.TFile.Open(f"{output_dir}/systematic.root", "UPDATE") as output_file:
            ratios_syst_prompt[0].SetName(f"dispersion_{cent_min}_{cent_max}")
            ratios_syst_prompt[0].Write()
            h_assigned_syst.Write()


        c.SaveAs(f"{output_dir}/Efficiency_PtWeight_Cent_{cent_min}_{cent_max}.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get the systematic uncertainty due to the pt weighting."
    )
    parser.add_argument(
        "config",
        type=str,
        help="Path to the config file.",
    )
    args = parser.parse_args()

    systematics(args.config)