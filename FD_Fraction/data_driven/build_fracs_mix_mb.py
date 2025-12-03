import argparse
import os
import ROOT
import yaml


def build_fracs(config):
    with open(config, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    mb_file_ds = ROOT.TFile.Open(config["mb_files"]["ds"])
    mb_file_dplus = ROOT.TFile.Open(config["mb_files"]["dplus"])

    output_file = config["output"]["ds"] + f"_MB.root"
    os.system(f"cp {config['mb_files']['ds']} {output_file}")
    output_file = config["output"]["dplus"] + f"_MB.root"
    os.system(f"cp {config['mb_files']['dplus']} {output_file}")

    cent_mins = config["cent"]["mins"]
    cent_maxs = config["cent"]["maxs"]

    for cent_min, cent_max in zip(cent_mins, cent_maxs):
        if config["swap_with_mb"]["ds"] is None or f"{cent_min}_{cent_max}" not in config["swap_with_mb"]["ds"]:
            # just copy the file to the output
            output_file = config["output"]["ds"] + f"_{cent_min}_{cent_max}.root"
            os.system(f"cp {config['mb_files']['ds'].replace("MB.root", f"{cent_min}_{cent_max}.root")} {output_file}")
        else:
            output_filename = config["output"]["ds"] + f"_{cent_min}_{cent_max}.root"
            with ROOT.TFile.Open(output_filename, "recreate") as output_file:
                with ROOT.TFile.Open(config["mb_files"]["ds"].replace("MB.root", f"{cent_min}_{cent_max}.root")) as cent_file_ds:
                    h_prompt_ds = cent_file_ds.Get(f"hRawFracPrompt_cent_{cent_min}_{cent_max}")
                    h_prompt_ds.SetDirectory(0)
                    h_non_prompt_ds = cent_file_ds.Get(f"hRawFracNonPrompt_cent_{cent_min}_{cent_max}")
                    h_non_prompt_ds.SetDirectory(0)
                for pt_min, pt_max in config["swap_with_mb"]["ds"][f"{cent_min}_{cent_max}"]:
                    i_bin = h_prompt_ds.FindBin((pt_min+pt_max)/2)
                    h_prompt_ds_mb = mb_file_ds.Get(f"hRawFracPrompt_cent_{cent_min}_{cent_max}")
                    h_prompt_ds.SetBinContent(i_bin, h_prompt_ds_mb.GetBinContent(i_bin))
                    h_prompt_ds.SetBinError(i_bin, h_prompt_ds_mb.GetBinError(i_bin))
                    h_non_prompt_ds_mb = mb_file_ds.Get(f"hRawFracNonPrompt_cent_{cent_min}_{cent_max}")
                    h_non_prompt_ds.SetBinContent(i_bin, h_non_prompt_ds_mb.GetBinContent(i_bin))
                    h_non_prompt_ds.SetBinError(i_bin, h_non_prompt_ds_mb.GetBinError(i_bin))
                output_file.cd()
                h_prompt_ds.Write()
                h_non_prompt_ds.Write()

        if config["swap_with_mb"]["dplus"] is None or f"{cent_min}_{cent_max}" not in config["swap_with_mb"]["dplus"]:
            # just copy the file to the output
            output_file = config["output"]["dplus"] + f"_{cent_min}_{cent_max}.root"
            os.system(f"cp {config['mb_files']['dplus'].replace("MB.root", f"{cent_min}_{cent_max}.root")} {output_file}")
        else:
            output_filename = config["output"]["dplus"] + f"_{cent_min}_{cent_max}.root"
            with ROOT.TFile.Open(output_filename, "recreate") as output_file:
                with ROOT.TFile.Open(config["mb_files"]["dplus"].replace("MB.root", f"{cent_min}_{cent_max}.root")) as cent_file_dplus:
                    h_prompt_dplus = cent_file_dplus.Get(f"hRawFracPrompt_cent_{cent_min}_{cent_max}")
                    h_prompt_dplus.SetDirectory(0)
                    h_non_prompt_dplus = cent_file_dplus.Get(f"hRawFracNonPrompt_cent_{cent_min}_{cent_max}")
                    h_non_prompt_dplus.SetDirectory(0)
                for pt_min, pt_max in config["swap_with_mb"]["dplus"][f"{cent_min}_{cent_max}"]:
                    i_bin = h_prompt_dplus.FindBin((pt_min+pt_max)/2)
                    h_prompt_dplus_mb = mb_file_dplus.Get(f"hRawFracPrompt_cent_{cent_min}_{cent_max}")
                    h_prompt_dplus.SetBinContent(i_bin, h_prompt_dplus_mb.GetBinContent(i_bin))
                    h_prompt_dplus.SetBinError(i_bin, h_prompt_dplus_mb.GetBinError(i_bin))
                    h_non_prompt_dplus_mb = mb_file_dplus.Get(f"hRawFracNonPrompt_cent_{cent_min}_{cent_max}")
                    h_non_prompt_dplus.SetBinContent(i_bin, h_non_prompt_dplus_mb.GetBinContent(i_bin))
                    h_non_prompt_dplus.SetBinError(i_bin, h_non_prompt_dplus_mb.GetBinError(i_bin))
                output_file.cd()
                h_prompt_dplus.Write()
                h_non_prompt_dplus.Write()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build the prompt fractions mixing with the MB")
    parser.add_argument("config", help="Configuration file")
    args = parser.parse_args()

    build_fracs(args.config)