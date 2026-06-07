'''
Script to evaluate the efficiency starting from the output THnSparse.
'''
from importlib.metadata import files
from math import sqrt
import argparse
import itertools
from pathlib import Path
import os
from dataclasses import dataclass
from typing import Optional, Tuple, List
import numpy as np
import yaml
import ROOT
# pylint: disable=no-member
import sys # pylint: disable=wrong-import-order
sys.path.append(os.path.abspath(f"{__file__}/../../utils")) # pylint: disable=wrong-import-position
sys.path.append(os.path.abspath(f"{__file__}/../../utils/sparse_reweighter")) # pylint: disable=wrong-import-position
from utils import enforce_list, CentralityInfo, PtInfo, pt_folders_exist, get_root_files_in_directory
from sparse_reweighter import SparseReweighter

PARTICLE_CLASSES = [("Ds", "Prompt"), ("Ds", "NonPrompt"),
                    ("Dplus", "Prompt"), ("Dplus", "NonPrompt")] # Ds must be the first one
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

class EfficiencyCalculator:
    """
    Class to calculate efficiencies from THnSparse histograms.
    """
    def __init__(self, cfg):
        self.particles = ["Ds", "Dplus"]
        self.origins = ["Prompt", "NonPrompt"]
        self.particle_origin = list(itertools.product(self.particles, self.origins))

        self.cfg = self._load_config(cfg)
        self.cuts_cfg = self._load_cuts_config()
        self.cent_info = self._extract_centrality_info()
        self.pt_info = self._extract_pt_info()
        self.has_mult_dependent_cuts = False
        self.cuts, self.axes = self._extract_cuts_axes()
        self.axis_cent = self.cuts_cfg['cent']['axisnum']
        self.axis_cent_gen = self.cfg['inputs']['sparse']['axis_cent_gen']
        # Reuse the shared SparseReweighter so the npv/pt weighting (including the
        # per-candidate 10%-binned pt weights) stays consistent with preprocessing.
        self.reweighter = self._build_reweighter()
        self.histos_pv_weights = self.reweighter.histos_pv_weights if self.reweighter else None
        self.histos_pt_weights = self.reweighter.histos_pt_weights if self.reweighter else None
        self.apply_cent_sel = self.cfg['inputs']['apply_cent_sel']
        self.is_preprocessed = self.cfg['inputs']['sparse']['is_preprocessed']

        self.histos, self.histos_fine_bins = [], []
        self.histos_bdt, self.histos_cent = [], []
        self.histos_cent_fine_bins, self.histos_cent_bdt, self.canvas_cent, self.legs = [], [], [], []
        self._create_histos()


    def _load_config(self, config_file: str) -> dict:
        """Load main configuration file."""
        with open(config_file, 'r', encoding="utf-8") as f:
            return yaml.safe_load(f)

    def _load_cuts_config(self) -> dict:
        """Load cuts configuration file."""
        with open(self.cfg['inputs']['cutset'], 'r', encoding="utf-8") as f:
            return yaml.safe_load(f)

    def _extract_centrality_info(self) -> Optional[CentralityInfo]:
        """Extract centrality information from cuts config."""
        if "cent" not in self.cuts_cfg:
            return None
        cent_selections = self.cuts_cfg['cent']
        return CentralityInfo(
            mins=cent_selections['min'],
            maxs=cent_selections['max']
        )

    def _extract_pt_info(self) -> PtInfo:
        """Extract pT information from cuts config."""
        pt_mins = self.cuts_cfg['pt']['min']
        pt_maxs = self.cuts_cfg['pt']['max']
        return PtInfo(mins=pt_mins, maxs=pt_maxs)

    def _extract_cuts_axes(self) -> Tuple[List[Tuple[float, float]], List[int]]:
        """Extract variable names and cuts from the cuts config."""
        var_names = list(self.cuts_cfg.keys())
        var_names.remove('mass')
        var_names.remove('pt')

        cuts_list  = [self.cuts_cfg[v] for v in var_names]
        mins_list  = [c['min'] for c in cuts_list]
        maxs_list  = [c['max'] for c in cuts_list]

        if "cent" in var_names:
            var_names.remove('cent')
            self.has_mult_dependent_cuts = any(isinstance(mins[0], list) for mins in mins_list)
            if self.has_mult_dependent_cuts:
                cuts = {}
                for i_cent, (cent_min, cent_max) in enumerate(zip(self.cent_info.mins, self.cent_info.maxs)):
                    cuts[(cent_min, cent_max)] = []
                    for var in var_names:
                        var_min = self.cuts_cfg[var]['min']
                        var_max = self.cuts_cfg[var]['max']
                        if isinstance(var_min[0], list):
                            cuts[(cent_min, cent_max)].append(list(zip(var_min[i_cent], var_max[i_cent])))
                        else:
                            cuts[(cent_min, cent_max)].append(list(zip(var_min, var_max)))
            else:
                cuts = [list(zip(self.cuts_cfg[var]['min'],
                    self.cuts_cfg[var]['max'])) for var in var_names]
        else:
            cuts = [list(zip(self.cuts_cfg[var]['min'],
                self.cuts_cfg[var]['max'])) for var in var_names]
        axes = [self.cuts_cfg[var]['axisnum'] for var in var_names]
        return cuts, axes

    def _build_reweighter(self) -> Optional[SparseReweighter]:
        if not (self.cfg['weights']['npv']['apply'] or self.cfg['weights']['pt']['apply']):
            return None

        sparse_cfg = self.cfg['inputs']['sparse']
        axes_cfg = {
            "axis_pt_gen": sparse_cfg["axis_pt_gen"],
            "axis_ptb_gen": sparse_cfg["axis_ptb_gen"],
            "axis_npv_gen": sparse_cfg["axis_npv_gen"],
            "axis_pt_reco": sparse_cfg["axis_pt_reco"],
            "axis_ptb_reco": sparse_cfg["axis_ptb_reco"],
            "axis_npv_reco": sparse_cfg["axis_npv_reco"],
            "axis_cent_gen": sparse_cfg["axis_cent_gen"],
            "axis_cent_reco": self.axis_cent,
        }
        return SparseReweighter(self.cfg['weights'], self.cent_info.edges, axes_cfg)

    def _create_histos(self):
        """
        Create histograms and canvases for efficiency evaluation.

        Parameters:
        - cfg (dict): Configuration dictionary containing various settings.
        - cuts_cfg (dict): Configuration dictionary containing cut variables.
        Returns:
        tuple: A tuple containing the following elements:
            - histos (list):
                List of ROOT.TH1F histograms for efficiency.
            - histos_fine_bins (list):
                List of ROOT.TH1F histograms for efficiency in fine bins.
            - histos_cent (list):
                List of ROOT.TH1F histograms for efficiency with centrality selections.
            - histos_cent_fine_bins (list):
                List of ROOT.TH1F histograms for efficiency in fine bins with centrality selections.
            - canvas_cent (list):
                List of ROOT.TCanvas objects for centrality selections.
        """

        for particle, origin in self.particle_origin:

            self.histos.append(ROOT.TH1F(f"eff_{particle}{origin}",
                "Efficiency; #it{p}_{T} (GeV/c); Efficiency",
                len(self.pt_info.edges)-1, np.asarray(self.pt_info.edges, "d")))
            self.histos_fine_bins.append(ROOT.TH1F(f"eff_fine_bins_{particle}{origin}",
                "Efficiency in fine bins; #it{p}_{T} (GeV/c); Efficiency", 240, 0, 24))
            self.histos_bdt.append(ROOT.TH1F(f"bdt_eff_{particle}{origin}",
                "BDT efficiency; #it{p}_{T} (GeV/c); Efficiency", len(self.pt_info.edges)-1, np.asarray(self.pt_info.edges, "d")))

            if self.cfg["weights"]["npv"]["apply"] or self.cfg["weights"]["pt"]["apply"] or self.cfg["inputs"]["apply_cent_sel"]:
                for cent_min, cent_max in zip(self.cent_info.mins, self.cent_info.maxs):
                    self.histos_cent.append(ROOT.TH1F(f"eff_{particle}{origin}_cent_{cent_min}_{cent_max}",
                        f"Efficiency_Cent_{cent_min}_{cent_max}; #it{{p}}_{{T}} (GeV/c); Efficiency",
                        len(self.pt_info.edges)-1, np.asarray(self.pt_info.edges, "d")))
                    self.histos_cent_fine_bins.append(ROOT.TH1F(
                        f"eff_fine_bins_{particle}{origin}_cent_{cent_min}_{cent_max}",
                        f"Efficiency in fine bins_Cent_{cent_min}_{cent_max};\
                        #it{{p}}_{{T}} (GeV/c); Efficiency", 240, 0, 24))
                    self.histos_cent_bdt.append(ROOT.TH1F(
                        f"bdt_eff_{particle}{origin}_cent_{cent_min}_{cent_max}",
                        f"BDT efficiency_Cent_{cent_min}_{cent_max};\
                        #it{{p}}_{{T}} (GeV/c); Efficiency", len(self.pt_info.edges)-1, np.asarray(self.pt_info.edges, "d")))

                self.canvas_cent.append(ROOT.TCanvas(
                    f"canvas_cent_{particle}{origin}", f"canvas_cent_{particle}{origin}", 800, 600))
    
    def _load_sparses_from_filename(self,
        in_file_name: str,
        particle: str,
        origin: str,
        cent_min: int = None,
        cent_max: int = None
        ) -> Tuple[ROOT.THnSparse, ROOT.THnSparse, int]:
        """
        Load THnSpares from a single input file.
        """
        suffix = ""
        if cent_min is not None and cent_max is not None:
            suffix = f"_{cent_min}_{cent_max}"

        with ROOT.TFile.Open(in_file_name) as in_file:
            try:
                h_sparse_gen_particles = in_file.Get(f"hf-task-ds/MC/{particle}/{origin}/hSparseGen{suffix}")
                try:
                    h_sparse_gen_particles.GetEntries()
                except:  # Backward compatibility (before the renaming of the sparse)
                    h_sparse_gen_particles = in_file.Get(f"hf-task-ds/MC/{particle}/{origin}/hPtYNPvContribGen{suffix}")
                    h_sparse_gen_particles.GetEntries()
            except:
                print(f"ERROR: histogram hf-task-ds/MC/{particle}/{origin}/hSparseGen{suffix} not found in file {in_file_name}."
                      " Using the unweighted sparse instead.")
                suffix = ""
                h_sparse_gen_particles = in_file.Get(f"hf-task-ds/MC/{particle}/{origin}/hSparseGen{suffix}")
                try:
                    h_sparse_gen_particles.GetEntries()
                except:  # Backward compatibility (before the renaming of the sparse)
                    h_sparse_gen_particles = in_file.Get(f"hf-task-ds/MC/{particle}/{origin}/hPtYNPvContribGen{suffix}")
                    h_sparse_gen_particles.GetEntries()

            h_sparse_reco = in_file.Get(f"hf-task-ds/MC/{particle}/{origin}/hSparseMass{suffix}")
            if self.cfg["weights"]["npv"]["reweight_on_data_events"]:
                if self.is_preprocessed:
                    n_events_mc = in_file.Get("hf-task-ds/Data/hCounterTVXafterBCcuts").GetEntries()
                else:
                    n_events_mc = in_file.Get("lumi-task/hCounterTVXafterBCcuts").GetEntries()
            else:
                n_events_mc = None

            return h_sparse_gen_particles, h_sparse_reco, n_events_mc

    def _load_sparse_pt_cent(self, particle: str, origin: str, cent_range: Tuple[float, float], pt_range: Tuple[float, float]) -> Tuple[List[ROOT.THnSparse], List[ROOT.THnSparse], List[int]]:
        """
        Load THnSpares from input files for specific pT and centrality ranges.
        """
        # get all root files in the pt folders
        file_name = Path(self.cfg["inputs"]["sparse"]["file_names"]) / \
            f'pt_{int(pt_range[0]*10)}_{int(pt_range[1]*10)}' / \
            f'AnalysisResults_pt_{int(pt_range[0]*10)}_{int(pt_range[1]*10)}.root'

        h_sparse_gen, h_sparse_reco, n_ev_mc = self._load_sparses_from_filename(str(file_name), particle, origin, cent_range[0], cent_range[1])

        return h_sparse_gen, h_sparse_reco, n_ev_mc

    def _load_sparses(self, particle: str, origin: str) -> Tuple[List[ROOT.THnSparse], List[ROOT.THnSparse], List[int]]:
        """
        Load THnSpares from input files.
        """
        h_sparses_gen_particles = []
        h_sparses_reco = []
        n_events_mc = []

        if self.is_preprocessed:
            # check if pt folders exist
            if not pt_folders_exist(self.cfg["inputs"]["sparse"]["file_names"], (self.pt_info.mins, self.pt_info.maxs)):
                raise ValueError("Not all pt folders exist in the provided directory.")

            # get all root files in the pt folders
            folders = [
                Path(self.cfg["inputs"]["sparse"]["file_names"]) / f'pt_{int(pt_min*10)}_{int(pt_max*10)}'
                for pt_min, pt_max in zip(self.pt_info.mins, self.pt_info.maxs)
            ]
            files = []
            for f in folders:
                files += get_root_files_in_directory(f)
            for f in files:
                h_sparse_gen, h_sparse_reco, n_ev_mc = self._load_sparses_from_filename(f, particle, origin)
                h_sparses_gen_particles.append(h_sparse_gen)
                h_sparses_reco.append(h_sparse_reco)
                n_events_mc.append(n_ev_mc)
        else:
            file_names = enforce_list(self.cfg["inputs"]["sparse"]["file_names"])
            for in_file_name in file_names:
                h_sparse_gen, h_sparse_reco, n_ev_mc = self._load_sparses_from_filename(in_file_name, particle, origin)
                h_sparses_gen_particles.append(h_sparse_gen)
                h_sparses_reco.append(h_sparse_reco)
                n_events_mc.append(n_ev_mc)

        return h_sparses_gen_particles, h_sparses_reco, n_events_mc

    def _get_event_weights(self, n_events_mc: List[int]) -> Optional[List[float]]:
        event_weights = None
        if self.cfg["weights"]["npv"]["reweight_on_data_events"]:
            n_events_data = []
            for an_results in self.cfg["weights"]["npv"]["data_analysis_results"]:
                an_results = enforce_list(an_results)
                n_events_data.append(0)
                for an_result in an_results:
                    in_file = ROOT.TFile.Open(an_result)
                    n_events_data[-1] += in_file.Get("lumi-task/hCounterTVXafterBCcuts").GetEntries()
                    in_file.Close()
            event_weights = np.array(n_events_data) / np.array(n_events_mc)

        return event_weights

    def _get_reco_particles(self, h_sparse_reco_ds, i_pt, pt_min, pt_max, cent_min, cent_max, apply_bdt):
        for i, axis in enumerate(self.axes):
            if not apply_bdt and axis in self.cfg['inputs']['sparse']['axes_bdt']:  # no cut on BDT score
                continue

            if self.has_mult_dependent_cuts:
                min_val = self.cuts[(cent_min, cent_max)][i][i_pt][0]
                max_val = self.cuts[(cent_min, cent_max)][i][i_pt][1]
            else:
                min_val = self.cuts[i][i_pt][0]
                max_val = self.cuts[i][i_pt][1]
            h_sparse_reco_ds.GetAxis(axis).SetRangeUser(min_val, max_val)

        # WARNING: just like TH3::Project3D("yx") and TTree::Draw("y:x"), Projection(y,x)
        # uses the first argument to define the y-axis and the second for the x-axis!
        h_pt_npv_projected_reco_ds = h_sparse_reco_ds.Projection(
            self.cfg['inputs']['sparse']['axis_npv_reco'], self.cfg['inputs']['sparse']['axis_pt_reco'], "EO")

        return self._get_integral_of_projection(
            h_pt_npv_projected_reco_ds, "x", pt_min, pt_max)

    def _get_gen_particles(self, h_sparse_gen, pt_min, pt_max):
        # WARNING: just like TH3::Project3D("yx") and TTree::Draw("y:x"), Projection(y,x)
        # uses the first argument to define the y-axis and the second for the x-axis!
        h_pt_npv_projected_gen = h_sparse_gen.Projection(
                self.cfg['inputs']['sparse']['axis_npv_gen'], self.cfg['inputs']['sparse']['axis_pt_gen'], "EO")

        return self._get_integral_of_projection(
            h_pt_npv_projected_gen, "x", pt_min, pt_max)

    def _calculate_efficiencies_with_unc(self, sels, recos, gens, ev_weights) -> tuple:
        """
        Calculate efficiencies with uncertainty.

        Parameters:
        - sels (float or list): The number of selected events.
        - gens (float or list): The number of generated events.
        - ev_weights (float or list): The event weights.

        Returns:
        - eff (float): The efficiency.
        - unc (float): The uncertainty of the efficiency.
        """
        sels = enforce_list(sels)
        recos = enforce_list(recos)
        gens = enforce_list(gens)
        tot_sel, tot_reco, tot_gen = 0, 0, 0
        if ev_weights is None:
            ev_weights = [1] * len(sels)
        for sel, reco, gen, weight in zip(sels, recos, gens, ev_weights):
            tot_sel += sel * weight
            tot_reco += reco * weight
            tot_gen += gen * weight
        eff = tot_sel / tot_gen # pylint: disable=redefined-outer-name
        bdt_eff = tot_sel / tot_reco
        return eff, sqrt(eff * (1 - eff) / tot_gen), bdt_eff, sqrt(bdt_eff * (1 - bdt_eff) / tot_reco)

    def _get_integral_of_projection(self, hist_2d, axis, minimum, maximum):
        """
        Get the integral of a 2D histogram projection along an axis.

        Parameters:
        - hist_2d (TH2F): The 2D histogram.
        - axis (str): "x" or "y".
        - minimum (float): The minimum value of the range.
        - maximum (float): The maximum value of the range.

        Returns:
        - integral (float): The integral of the projection.
        """
        if axis == "x":
            hist_1d = hist_2d.ProjectionX("projection", 0, -1, "EO")
        elif axis == "y":
            hist_1d = hist_2d.ProjectionY("projection", 0, -1, "EO")
        else:
            raise ValueError("Invalid axis. Choose 'x' or 'y'.")

        # Right edge of the bin is not included
        return hist_1d.Integral(hist_1d.FindBin(minimum), hist_1d.FindBin(maximum-0.001))

    def _fill_eff_histos(self, idx, i_pt, i_cent, eff, unc, bdt_eff, use_centrality):
        """Set histogram bin contents and errors for efficiency results."""
        if use_centrality:
            i_hist = idx * len(self.cent_info.mins) + i_cent
            self.histos_cent[i_hist].SetBinContent(i_pt + 1, eff)
            self.histos_cent[i_hist].SetBinError(i_pt + 1, unc)
            self.histos_cent_bdt[i_hist].SetBinContent(i_pt + 1, bdt_eff)
        else:
            self.histos[idx].SetBinContent(i_pt + 1, eff)
            self.histos[idx].SetBinError(i_pt + 1, unc)
            self.histos_bdt[idx].SetBinContent(i_pt + 1, bdt_eff)

    def get_eff(self, particle_class, h_sparses_gen_particles, h_sparses_reco_ds, pt_info, cent_info = (None, None), events_weights=None): # pylint: disable=redefined-outer-name
        """
        Calculate the efficiency of particle selection.

        Parameters:
        - particle_class (tuple): Tuple containing the particle and its origin.
        - pt_info (tuple): Tuple containing the index, minimum and maximum values of pT.
        - cent_info (tuple, optional): Tuple containing the minimum and maximum values of centrality.
        - events_weights (list, optional): List of weights for each dataset.)

        Returns:
        - float: Efficiency of particle selection.
        """
        particle, origin = particle_class # pylint: disable=redefined-outer-name
        i_pt, pt_min, pt_max = pt_info
        cent_min, cent_max = cent_info
        gen_particles = []
        selected_particles = []
        reco_particles = []

        # Loop on the different input files
        for h_sparse_gen_particles, h_sparse_reco_ds in zip(h_sparses_gen_particles, h_sparses_reco_ds):
            if self.histos_pv_weights is not None or self.histos_pt_weights is not None or self.apply_cent_sel:
                h_sparse_gen_sel = h_sparse_gen_particles.Clone(
                    f"hSparseGenParticlesSel_{particle}{origin}_Cent_{cent_min}_{cent_max}_{i_pt}")
                h_sparse_reco_ds_sel = h_sparse_reco_ds.Clone(
                    f"hSparseRecoDsSel_{particle}{origin}_Cent_{cent_min}_{cent_max}_{i_pt}")
            else:
                h_sparse_gen_sel = h_sparse_gen_particles.Clone(
                    f"hSparseGenSel_{particle}{origin}")
                h_sparse_reco_ds_sel = h_sparse_reco_ds.Clone(
                    f"hSparseRecoDsSel_{particle}{origin}")

            # Apply centrality selection
            if self.apply_cent_sel:
                if cent_info[0] is None or cent_info[1] is None:
                    print("ERROR: centrality dependent requested, but no centrality bins available. Exit")
                    sys.exit(0)
                h_sparse_reco_ds_sel.GetAxis(self.axis_cent).SetRangeUser(cent_info[0], cent_info[1])
                h_sparse_gen_sel.GetAxis(self.axis_cent_gen).SetRangeUser(cent_info[0], cent_info[1])

            reco_parts = self._get_reco_particles(h_sparse_reco_ds_sel, i_pt, pt_min, pt_max, cent_info[0], cent_info[1], apply_bdt=False)
            sel_parts = self._get_reco_particles(h_sparse_reco_ds_sel, i_pt, pt_min, pt_max, cent_info[0], cent_info[1], apply_bdt=True)
            gen_parts = self._get_gen_particles(h_sparse_gen_sel, pt_min, pt_max)

            gen_particles.append(gen_parts)
            reco_particles.append(reco_parts)
            selected_particles.append(sel_parts)

        return self._calculate_efficiencies_with_unc(selected_particles, reco_particles, gen_particles, events_weights)

    def _draw_histos_with_centrality(self, colors):
        for idx, (particle, origin) in enumerate(self.particle_origin):
            leg = ROOT.TLegend(0.6, 0.2, 0.9, 0.8)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextSize(0.035)
            self.legs.append(leg)

            canvas = self.canvas_cent[idx]
            canvas.cd().SetLogy()

            # Draw min bias
            h_mb = self.histos[idx]
            h_mb.SetMarkerColor(colors[0])
            h_mb.SetLineColor(colors[0])
            h_mb.Draw()
            leg.AddEntry(h_mb, "Min. bias", "l")

            # Draw centrality bins
            for i_cent, (cent_min, cent_max) in enumerate(zip(self.cent_info.mins, self.cent_info.maxs)):
                i_hist = idx * len(self.cent_info.mins) + i_cent
                h_cent = self.histos_cent[i_hist]
                h_cent.SetMarkerColor(colors[i_cent + 1])
                h_cent.SetLineColor(colors[i_cent + 1])
                h_cent.Draw("same")
                leg.AddEntry(h_cent, f"{cent_min}#font[122]{{-}}{cent_max}%", "l")

            leg.Draw()


    def _write_histos(self, output_file, histos, fine_bins, bdt=None):
        output_file.cd()
        for i, h in enumerate(histos):
            h.Write()
            fine_bins[i].Write()
            if bdt:
                bdt[i].Write()

    def _compute_efficiency_ratios(self, use_centrality):
        hratios, hratios_fine_bins, hratios_cent, hratios_cent_fine_bins = ([] for _ in range(4))
        offset = 2

        for idx, (particle, origin) in enumerate(PARTICLE_CLASSES):
            if use_centrality:
                for i_cent, _ in enumerate(zip(self.cent_info.mins, self.cent_info.maxs)):
                    idx_cent = idx * len(self.cent_info.mins) + i_cent
                    if "Ds" in particle:
                        h_ratio = self.histos_cent[idx_cent].Clone(self.histos_cent[idx_cent].GetName().replace("eff", "effratio"))
                        h_ratio_fine = self.histos_cent_fine_bins[idx_cent].Clone(self.histos_cent_fine_bins[idx_cent].GetName().replace("eff", "effratio"))
                        h_ratio.SetTitle("Efficiency ratio")
                        h_ratio_fine.SetTitle("Efficiency ratio")
                        hratios_cent.append(h_ratio)
                        hratios_cent_fine_bins.append(h_ratio_fine)
                    else:
                        idx_ref = (idx - offset) * len(self.cent_info.mins) + i_cent
                        hratios_cent[idx_ref].Divide(self.histos_cent[idx_cent])
                        hratios_cent_fine_bins[idx_ref].Divide(self.histos_cent_fine_bins[idx_cent])
                        name_suffix = f"ompt_{particle}{origin}"
                        for h, n in zip([hratios_cent[idx_ref], hratios_cent_fine_bins[idx_ref]], ["", "_fine"]):
                            h.SetName(h.GetName().replace("ompt", name_suffix))
            else:
                if "Ds" in particle:
                    h_ratio = self.histos[idx].Clone(self.histos[idx].GetName().replace("eff", "effratio"))
                    h_ratio_fine = self.histos_fine_bins[idx].Clone(self.histos_fine_bins[idx].GetName().replace("eff", "effratio"))
                    h_ratio.SetTitle("Efficiency ratio")
                    h_ratio_fine.SetTitle("Efficiency ratio")
                    hratios.append(h_ratio)
                    hratios_fine_bins.append(h_ratio_fine)
                else:
                    hratios[idx - offset].Divide(self.histos[idx])
                    hratios_fine_bins[idx - offset].Divide(self.histos_fine_bins[idx])
                    name_suffix = f"ompt_{particle}{origin}"
                    for h, n in zip([hratios[idx - offset], hratios_fine_bins[idx - offset]], ["", "_fine"]):
                        h.SetName(h.GetName().replace("ompt", name_suffix))

        return hratios, hratios_fine_bins, hratios_cent, hratios_cent_fine_bins

    def dump_outputs(self, out_file_name, use_centrality):
        colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kOrange+7, ROOT.kGreen+2, ROOT.kAzure+4, ROOT.kBlack]

        if use_centrality:
            self._draw_histos_with_centrality(colors)

        output_file = ROOT.TFile.Open(out_file_name, "update")
        self._write_histos(output_file, self.histos, self.histos_fine_bins, self.histos_bdt)

        if use_centrality:
            self._write_histos(output_file, self.histos_cent, self.histos_cent_fine_bins, self.histos_cent_bdt)
            for canvas in self.canvas_cent:
                canvas.Write()

        # Efficiency ratios
        hratios, hratios_fine_bins, hratios_cent, hratios_cent_fine_bins = self._compute_efficiency_ratios(use_centrality)

        self._write_histos(output_file, hratios, hratios_fine_bins)
        if use_centrality:
            self._write_histos(output_file, hratios_cent, hratios_cent_fine_bins)

        output_file.Close()

    def run(self):
        """Main method to run the efficiency calculation."""
        apply_cent_sel = self.cfg['inputs']['apply_cent_sel']

        # Create output file
        out_file_name = os.path.join(
            self.cfg['output_dir'],
            f"efficiency_{self.cfg['suffix']}.root"
        )
        output_file = ROOT.TFile(out_file_name, "RECREATE")
        output_file.Close()

        if self.cfg['weights']['npv']['apply']:
            # Save the weights in the output file
            with ROOT.TFile.Open(out_file_name, "UPDATE") as output_file:
                if not output_file.GetListOfKeys().Contains("Weights"):
                    output_file.mkdir("Weights")
                output_file.cd("Weights")
                for particle, origin in self.particle_origin:
                    for i_cent, (cent_min, cent_max) in enumerate(zip(self.cent_info.mins, self.cent_info.maxs)):
                        self.histos_pv_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].Write()
        if self.cfg['weights']['pt']['apply']:
            # Save the weights in the output file
            with ROOT.TFile.Open(out_file_name, "UPDATE") as output_file:
                if not output_file.GetListOfKeys().Contains("Weights"):
                    output_file.mkdir("Weights")
                output_file.cd("Weights")
                for particle, origin in self.particle_origin:
                    for _, _, histo in self.histos_pt_weights[f"{particle}{origin}"]:
                        if histo is not None:
                            histo.Write()

        use_cent = apply_cent_sel or self.cfg["weights"]["npv"]["apply"] or self.cfg["weights"]["pt"]["apply"]
        for idx, (particle, origin) in enumerate(self.particle_origin):
            if not self.is_preprocessed:
                h_sparses_gen_particles, h_sparses_reco, n_events_mc = self._load_sparses(particle, origin)
            if self.cfg["weights"]["npv"]["reweight_on_data_events"]:
                event_weights = self._get_event_weights(n_events_mc)
            else:
                event_weights = None

            # Decide whether to iterate over centrality bins
            cent_bins = (
                list(zip(self.cent_info.mins, self.cent_info.maxs)) if use_cent else [None]
            )

            for i_cent, cent_range in enumerate(cent_bins):
                if not self.is_preprocessed:
                    h_sparses_gen_particles_for_eff = [s.Clone() for s in h_sparses_gen_particles] # After rewighting (if any)
                    h_sparses_reco_for_eff = [s.Clone() for s in h_sparses_reco] # After rewighting (if any)

                    # Apply weights if requested
                    if self.reweighter is not None:
                        for h_sparse_gen, h_sparse_reco in zip(h_sparses_gen_particles_for_eff, h_sparses_reco_for_eff):
                            self.reweighter.apply_weights(h_sparse_gen, particle, origin, cent_range[0], cent_range[1], is_gen=True)
                            self.reweighter.apply_weights(h_sparse_reco, particle, origin, cent_range[0], cent_range[1], is_gen=False)

                for i_pt, (pt_min, pt_max) in enumerate(zip(self.pt_info.mins, self.pt_info.maxs)):
                    if self.cfg["verbose"]:
                        print(f"Processing pT bin {pt_min} - {pt_max} for {particle} {origin}")
                    if self.is_preprocessed:
                        h_sparses_gen_particles_for_eff, h_sparses_reco_for_eff, n_events_mc = self._load_sparse_pt_cent(particle, origin, cent_range, (pt_min, pt_max))
                        if self.cfg["weights"]["npv"]["reweight_on_data_events"]:
                            event_weights = self._get_event_weights(n_events_mc)
                        h_sparses_gen_particles_for_eff = enforce_list(h_sparses_gen_particles_for_eff)
                        h_sparses_reco_for_eff = enforce_list(h_sparses_reco_for_eff)

                    eff, unc, bdt_eff, _ = self.get_eff(
                        (particle, origin),
                        h_sparses_gen_particles_for_eff,
                        h_sparses_reco_for_eff,
                        (i_pt, pt_min, pt_max),
                        cent_info=cent_range,
                        events_weights=event_weights,
                    )

                    self._fill_eff_histos(idx, i_pt, i_cent, eff, unc, bdt_eff, use_cent)

            if not self.is_preprocessed:
                del h_sparses_gen_particles, h_sparses_reco

        self.dump_outputs(out_file_name, use_cent)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate efficiencies for ML selection.')
    parser.add_argument('config_file', type=str, help='Path to the configuration file.')
    args = parser.parse_args()

    evaluator = EfficiencyCalculator(args.config_file)
    evaluator.run()
