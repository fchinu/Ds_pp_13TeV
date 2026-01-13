import argparse
import ROOT
import yaml

ROOT.TH1.AddDirectory(False)

def merge(config_path):
    # Load configuration
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    input_files = config["inputs"]
    swap_with_fc = config["swap_with_fc"]
    output_file = config["output"]

    with ROOT.TFile.Open(output_file, "RECREATE") as f_out:
        for i_cent, cent_key in enumerate(swap_with_fc.keys()):
            with ROOT.TFile.Open(input_files["data_driven"], "READ") as f_data:
                h_dd = f_data.Get(f"assigned_syst_{cent_key}")

            with ROOT.TFile.Open(input_files["fc"], "READ") as f_fc:
                h_fc_lower = f_fc.Get(f"assigned_syst_lower_{cent_key}")
                h_fc_upper = f_fc.Get(f"assigned_syst_upper_{cent_key}")


            h_syst_lower = h_dd.Clone(f"assigned_syst_lower_{cent_key}")
            h_syst_upper = h_dd.Clone(f"assigned_syst_upper_{cent_key}")
            # Swap specified points
            for pt_interval in swap_with_fc[cent_key]:
                pt_centre = sum(pt_interval) / 2
                h_syst_lower.SetBinContent(h_syst_lower.FindBin(pt_centre), h_fc_lower.GetBinContent(h_fc_lower.FindBin(pt_centre)))
                h_syst_upper.SetBinContent(h_syst_upper.FindBin(pt_centre), h_fc_upper.GetBinContent(h_fc_upper.FindBin(pt_centre)))

            # Save merged graphs
            f_out.cd()
            h_syst_lower.Write(f"assigned_syst_lower_{cent_key}")
            h_syst_upper.Write(f"assigned_syst_upper_{cent_key}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Path to config file")
    args = parser.parse_args()

    merge(args.config)
