'''
Script to evaluate the efficiency starting from the output THnSparse.
'''
from math import sqrt
import argparse
import itertools
import os
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict
import numpy as np
from tqdm import tqdm
import yaml
import ROOT
# pylint: disable=no-member
import sys # pylint: disable=wrong-import-order

PARTICLE_CLASSES = [("Ds", "Prompt"), ("Ds", "NonPrompt"),
                    ("Dplus", "Prompt"), ("Dplus", "NonPrompt")] # Ds must be the first one
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

@dataclass
class CentralityInfo:
    """Container for centrality bin information."""
    mins: List[float]
    maxs: List[float]


@dataclass
class PtInfo:
    """Container for pT bin information."""
    mins: List[float]
    maxs: List[float]
    edges: List[float]

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
        self.cuts, self.axes = self._extract_cuts_axes()
        self.axis_cent = self.cuts_cfg['cent']['axisnum']
        self.axis_cent_gen = self.cfg['inputs']['sparse']['axis_cent_gen']
        self.histos_pv_weights = self._load_pv_weights() if self.cfg['weights']['npv']['apply'] else None
        self.histos_pt_weights = self._load_pt_weights() if self.cfg['weights']['pt']['apply'] else None
        self.apply_cent_sel = self.cfg['inputs']['apply_cent_sel']

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
        pt_edges = pt_mins + [pt_maxs[-1]]
        return PtInfo(mins=pt_mins, maxs=pt_maxs, edges=pt_edges)

    def _extract_cuts_axes(self) -> Tuple[List[Tuple[float, float]], List[int]]:
        """Extract variable names and cuts from the cuts config."""
        var_names = list(self.cuts_cfg.keys())
        var_names.remove('mass')
        var_names.remove('pt')
        if "cent" in var_names:
            var_names.remove('cent')

        cuts = [list(zip(self.cuts_cfg[var]['min'],
            self.cuts_cfg[var]['max'])) for var in var_names]
        axes = [self.cuts_cfg[var]['axisnum'] for var in var_names]
        return cuts, axes

    def _enforce_list(self, var):
        """
        Ensure that the input variable is a list.

        Parameters:
        - var (any): The input variable.

        Returns:
        - list: The input variable as a list.
        """
        if not isinstance(var, list):
            var = [var]
        return var

    def _load_pv_weights(self) -> Dict[str, ROOT.TH1F]:
        """
        Get the weights from the input file.
        """

        histos_pv_weights = {}

        with ROOT.TFile.Open(self.cfg["weights"]["npv"]["file_name"]) as infile_weights:
            for particle, origin in self.particle_origin:
                for cent_min, cent_max in zip(self.cent_info.mins, self.cent_info.maxs):
                    histos_pv_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"] = infile_weights.Get(
                        f"centrality_{cent_min}_{cent_max}/{particle}_{origin}/hNPvContribCands_weights_{particle}_{origin}")
                    histos_pv_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].SetName(
                        f"NPvWeight{particle}{origin}Cands_{cent_min}_{cent_max}")
                    histos_pv_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].SetTitle(
                        f"NPvWeight{particle}{origin}Cands_{cent_min}_{cent_max}")
        return histos_pv_weights

    def _load_pt_weights(self) -> Dict[str, ROOT.TH1F]:
        """
        Get the weights from the input file.
        """

        histos_pt_weights = {}

        with ROOT.TFile.Open(self.cfg["weights"]["pt"]["file_name"]) as infile_weights:
            for particle, origin in self.particle_origin:
                for cent_min, cent_max in zip(self.cent_info.mins, self.cent_info.maxs):
                    histo_name = self.cfg["weights"]["pt"]["base_hist_name"][particle][origin].replace("*", f"{cent_min}_{cent_max}")
                    histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"] = infile_weights.Get(
                        f"weights/{histo_name}"
                    )
                    if isinstance(histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"], ROOT.TH1): # typically not present for 0-100%
                        histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].SetName(
                            f"PtWeight{particle}{origin}Cands_{cent_min}_{cent_max}")
                        histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].SetTitle(
                            f"PtWeight{particle}{origin}Cands_{cent_min}_{cent_max}")
                    else:
                        histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"] = None
        return histos_pt_weights

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

    def _load_sparses(self, particle: str, origin: str) -> Tuple[List[ROOT.THnSparse], List[ROOT.THnSparse], List[int]]:
        """
        Load THnSpares from input files.
        """
        h_sparses_gen_particles = []
        h_sparses_reco = []
        n_events_mc = []

        self._enforce_list(self.cfg["inputs"]["sparse"]["file_names"])
        for in_file_name in self.cfg["inputs"]["sparse"]["file_names"]:
            in_file = ROOT.TFile.Open(in_file_name)
            h_sparses_gen_particles.append(in_file.Get(f"hf-task-ds/MC/{particle}/{origin}/hSparseGen"))
            try:
                h_sparses_gen_particles[-1].GetEntries()
            except:  # Backward compatibility (before the renaming of the sparse)
                h_sparses_gen_particles.pop()
                h_sparses_gen_particles.append(in_file.Get(f"hf-task-ds/MC/{particle}/{origin}/hPtYNPvContribGen"))
                h_sparses_gen_particles[-1].GetEntries()

            h_sparses_reco.append(in_file.Get(f"hf-task-ds/MC/{particle}/{origin}/hSparseMass"))
            if self.cfg["weights"]["npv"]["reweight_on_data_events"]:
                n_events_mc.append(in_file.Get("lumi-task/hCounterTVXafterBCcuts").GetEntries())
            in_file.Close()

        return h_sparses_gen_particles, h_sparses_reco, n_events_mc

    def _get_event_weights(self, n_events_mc: List[int]) -> Optional[List[float]]:
        event_weights = None
        if self.cfg["weights"]["npv"]["reweight_on_data_events"]:
            n_events_data = []
            for an_results in self.cfg["weights"]["npv"]["data_analysis_results"]:
                an_results = self._enforce_list(an_results)
                n_events_data.append(0)
                for an_result in an_results:
                    in_file = ROOT.TFile.Open(an_result)
                    n_events_data[-1] += in_file.Get("lumi-task/hCounterTVXafterBCcuts").GetEntries()
                    in_file.Close()
            event_weights = np.array(n_events_data) / np.array(n_events_mc)

        return event_weights

    def _get_reco_particles(self, h_sparse_reco_ds, i_pt, pt_min, pt_max, apply_bdt):
        for i, axis in enumerate(self.axes):
            if not apply_bdt and axis in self.cfg['inputs']['sparse']['axes_bdt']:  # no cut on BDT score
                continue

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
        sels = self._enforce_list(sels)
        recos = self._enforce_list(recos)
        gens = self._enforce_list(gens)
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

    def _apply_weights(self, sparse, particle, origin, cent_min, cent_max, is_gen=True):

        coord =  [0] * sparse.GetNdimensions()
        coord = np.array(coord, dtype=np.int32)

        pt_weight, multi_weight = 1., 1.
        pt_axis, mult_axis = -1, -1
        histo_pt_weight, histo_pv_weight = None, None
        if self.histos_pt_weights is not None:
            histo_pt_weight = self.histos_pt_weights[
                f"{particle}{origin}_Cent_{cent_min}_{cent_max}"]
            if origin == "Prompt":
                if is_gen:
                    pt_axis = self.cfg["inputs"]["sparse"]["axis_pt_gen"]
                else:
                    pt_axis = self.cfg["inputs"]["sparse"]["axis_pt_reco"]
            else:
                if is_gen:
                    pt_axis = self.cfg["inputs"]["sparse"]["axis_ptb_gen"]
                else:
                    pt_axis = self.cfg["inputs"]["sparse"]["axis_ptb_reco"]
        if self.histos_pv_weights is not None:
            histo_pv_weight = self.histos_pv_weights[
                f"{particle}{origin}_Cent_{cent_min}_{cent_max}"]
            mult_axis = self.cfg["inputs"]["sparse"]["axis_npv_gen"] if is_gen else self.cfg["inputs"]["sparse"]["axis_npv_reco"]

        for i in tqdm(range(1, sparse.GetNbins())):
            content = sparse.GetBinContent(i, coord)
            error = np.sqrt(sparse.GetBinError2(i))

            if histo_pt_weight is not None:
                pt_weight = histo_pt_weight.GetBinContent(int(coord[pt_axis]))

            if histo_pv_weight is not None:
                multi_weight = histo_pv_weight.GetBinContent(int(coord[mult_axis]))

            weight = multi_weight * pt_weight
            sparse.SetBinContent(i, content*weight)
            sparse.SetBinError(i, error*weight)

        # reweight bin 0 (crashes if at the beginning of the loop)
        content = sparse.GetBinContent(0, coord)
        error = np.sqrt(sparse.GetBinError2(0))

        pt_weight, multi_weight = 1., 1.
        if histo_pt_weight is not None:
            pt_weight = histo_pt_weight.GetBinContent(int(coord[pt_axis]))
        if histo_pv_weight is not None:
            multi_weight = histo_pv_weight.GetBinContent(int(coord[mult_axis]))

        weight = multi_weight * pt_weight
        sparse.SetBinContent(0, content*weight)
        sparse.SetBinError(0, error*weight)


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

            reco_parts = self._get_reco_particles(h_sparse_reco_ds_sel, i_pt, pt_min, pt_max, apply_bdt=False)
            sel_parts = self._get_reco_particles(h_sparse_reco_ds_sel, i_pt, pt_min, pt_max, apply_bdt=True)
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
                    for i_cent, (cent_min, cent_max) in enumerate(zip(self.cent_info.mins, self.cent_info.maxs)):
                        if self.histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"] is not None:
                            self.histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].Write()

        use_cent = apply_cent_sel or self.cfg["weights"]["npv"]["apply"] or self.cfg["weights"]["pt"]["apply"]
        for idx, (particle, origin) in enumerate(self.particle_origin):
            h_sparses_gen_particles, h_sparses_reco, n_events_mc = self._load_sparses(particle, origin)
            event_weights = self._get_event_weights(n_events_mc)

            # Decide whether to iterate over centrality bins
            cent_bins = (
                list(zip(self.cent_info.mins, self.cent_info.maxs)) if use_cent else [None]
            )

            for i_cent, cent_range in enumerate(cent_bins):
                h_sparses_gen_particles_for_eff = [s.Clone() for s in h_sparses_gen_particles] # After rewighting (if any)
                h_sparses_reco_for_eff = [s.Clone() for s in h_sparses_reco] # After rewighting (if any)

                # Apply weights if requested
                if self.histos_pv_weights is not None or self.histos_pt_weights is not None:
                    for h_sparse_gen, h_sparse_reco in zip(h_sparses_gen_particles_for_eff, h_sparses_reco_for_eff):
                        self._apply_weights(h_sparse_gen, particle, origin, cent_range[0], cent_range[1], is_gen=True)
                        self._apply_weights(h_sparse_reco, particle, origin, cent_range[0], cent_range[1], is_gen=False)

                for i_pt, (pt_min, pt_max) in enumerate(zip(self.pt_info.mins, self.pt_info.maxs)):
                    if self.cfg["verbose"]:
                        print(f"Processing pT bin {pt_min} - {pt_max} for {particle} {origin}")

                    eff, unc, bdt_eff, _ = self.get_eff(
                        (particle, origin),
                        h_sparses_gen_particles_for_eff,
                        h_sparses_reco_for_eff,
                        (i_pt, pt_min, pt_max),
                        cent_info=cent_range,
                        events_weights=event_weights,
                    )

                    self._fill_eff_histos(idx, i_pt, i_cent, eff, unc, bdt_eff, use_cent)

            del h_sparses_gen_particles, h_sparses_reco

        self.dump_outputs(out_file_name, use_cent)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate efficiencies for ML selection.')
    parser.add_argument('config_file', type=str, help='Path to the configuration file.')
    args = parser.parse_args()

    evaluator = EfficiencyCalculator(args.config_file)
    evaluator.run()
