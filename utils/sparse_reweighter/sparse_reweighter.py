import itertools
import numpy as np
import ROOT
from tqdm import tqdm
import sys
sys.path.append(f"{__file__}/../..")

from utils import enforce_list, CentralityInfo
ROOT.TH1.AddDirectory(False)

class SparseReweighter:

    def __init__(
            self,
            cfg: dict,
            centrality_bins: list,
            axes_cfg: dict = None
        ):
        """
        Class for reweighting THnSpares.

        Args:
            sparse (dict): THnSpares to be reweighted.
            cfg (dict): Configuration dictionary. Should contain a structure like:
                npv:
                    apply: False
                    file_name: 'file.root'
                pt:
                    apply: True
                    file_name: 'ptweights.root'
                    base_hist_name: 
                        Ds: 
                            Prompt: 'hist_ds_weights_tamu_*'
                            NonPrompt: 'hist_bs_weights_tamu_*' # not exact, it should be a mixture 50%-50% of B and Bs
                        Dplus: 
                            Prompt: 'hist_d_weights_tamu_*'
                            NonPrompt: 'hist_b_weights_tamu_*'
            centrality_bins (list): List of centrality bin edges.
            axes_cfg (dict): Configuration for the axes of the THnSpares. Can be skipped if contained in cfg. Should contain:
                axis_pt_gen
                axis_ptb_gen
                axis_npv_gen
                axis_pt_reco
                axis_ptb_reco
                axis_npv_reco
                axis_cent_gen
        """
        self.particles = ["Ds", "Dplus"]
        self.origins = ["Prompt", "NonPrompt"]
        self.particle_origin = list(itertools.product(self.particles, self.origins))

        self.cfg = cfg
        self.cent_info = self._get_cent_info(centrality_bins)
        self.axes_cfg = self._get_axes_cfg(axes_cfg)
        
        self.histos_pv_weights = None
        self.histos_pt_weights = None
        self._load_weights()
        print(self.histos_pt_weights)
    
    def _get_cent_info(self, centrality_bins) -> CentralityInfo:
        """Get centrality bin info from bin edges."""
        return CentralityInfo(
            mins=centrality_bins[:-1],
            maxs=centrality_bins[1:]
        )

    def _get_axes_cfg(self, axes_cfg: dict = None) -> dict:
        """Get axes configuration from cfg or input."""
        keys = [
            "axis_pt_gen",
            "axis_ptb_gen",
            "axis_npv_gen",
            "axis_pt_reco",
            "axis_ptb_reco",
            "axis_npv_reco",
            "axis_cent_gen"
        ]
        if axes_cfg is not None:
            return {key: axes_cfg[key] for key in keys}
        else:
            return {key: self.cfg[key] for key in keys}

    def _load_pv_weights(self) -> dict[str, ROOT.TH1F]:
        """
        Get the weights from the input file.
        """

        self.histos_pv_weights = {}

        with ROOT.TFile.Open(self.cfg["npv"]["file_name"]) as infile_weights:
            for particle, origin in self.particle_origin:
                for cent_min, cent_max in zip(self.cent_info.mins, self.cent_info.maxs):
                    self.histos_pv_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"] = infile_weights.Get(
                        f"centrality_{cent_min}_{cent_max}/{particle}_{origin}/hNPvContribCands_weights_{particle}_{origin}")
                    self.histos_pv_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].SetName(
                        f"NPvWeight{particle}{origin}Cands_{cent_min}_{cent_max}")
                    self.histos_pv_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].SetTitle(
                        f"NPvWeight{particle}{origin}Cands_{cent_min}_{cent_max}")

    def _load_pt_weights(self) -> dict[str, ROOT.TH1F]:
        """
        Get the weights from the input file.
        """

        self.histos_pt_weights = {}

        with ROOT.TFile.Open(self.cfg["pt"]["file_name"]) as infile_weights:
            infile_weights.ls()
            for particle, origin in self.particle_origin:
                for cent_min, cent_max in zip(self.cent_info.mins, self.cent_info.maxs):
                    histo_name = self.cfg["pt"]["base_hist_name"][particle][origin].replace("*", f"{cent_min}_{cent_max}")
                    self.histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"] = infile_weights.Get(
                        f"weights/{histo_name}"
                    )
                    if isinstance(self.histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"], ROOT.TH1): # typically not present for 0-100%
                        self.histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].SetName(
                            f"PtWeight{particle}{origin}Cands_{cent_min}_{cent_max}")
                        self.histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].SetTitle(
                            f"PtWeight{particle}{origin}Cands_{cent_min}_{cent_max}")
                    else:
                        self.histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"] = None

    def _load_weights(self):
        """
        Load the weights histograms if reweighting is requested.
        """
        if self.cfg['npv']['apply']:
            self._load_pv_weights()

        if self.cfg['pt']['apply']:
            self._load_pt_weights()

    def _apply_weights(self, sparse, particle, origin, cent_min, cent_max, is_gen=True):

        coord =  [0] * sparse.GetNdimensions()
        coord = np.array(coord, dtype=np.int32)

        pt_weight, multi_weight = 1., 1.
        pt_axis, mult_axis = -1, -1
        histo_pt_weight, histo_pv_weight = None, None
        print(f"Applying weights for {particle} {origin} Cent {cent_min}-{cent_max}, is_gen={is_gen}")
        if self.histos_pt_weights is not None:
            print(f"Applying pt weights for {particle} {origin} Cent {cent_min}-{cent_max}, is_gen={is_gen}")
            histo_pt_weight = self.histos_pt_weights[
                f"{particle}{origin}_Cent_{cent_min}_{cent_max}"]
            if origin == "Prompt":
                if is_gen:
                    pt_axis = self.axes_cfg["axis_pt_gen"]
                else:
                    pt_axis = self.axes_cfg["axis_pt_reco"]
            else:
                if is_gen:
                    pt_axis = self.axes_cfg["axis_ptb_gen"]
                else:
                    pt_axis = self.axes_cfg["axis_ptb_reco"]

        if self.histos_pv_weights is not None:
            histo_pv_weight = self.histos_pv_weights[
                f"{particle}{origin}_Cent_{cent_min}_{cent_max}"]
            mult_axis = self.axes_cfg["axis_npv_gen"] if is_gen else self.axes_cfg["axis_npv_reco"]

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

    def reweight(
            self,
            sparse: ROOT.THnSparse,
            is_gen: bool = False,
            particle: str = "Ds",
            origin: str = "Prompt"):
        """
        Apply the reweighting to the THnSpares and save them to the output file.
        """
        sparses = {}
        for cent_range in zip(self.cent_info.mins, self.cent_info.maxs):
            h_sparse = sparse.Clone()

            # Apply weights if requested
            print(self.histos_pv_weights, self.histos_pt_weights)
            if self.histos_pv_weights is not None or self.histos_pt_weights is not None:
                self._apply_weights(h_sparse, particle, origin, cent_range[0], cent_range[1], is_gen=is_gen)
            sparses[f"{cent_range[0]}_{cent_range[1]}"] = h_sparse
        return sparses

    def dump_to_root(self, out_file_name: str, option: str = "UPDATE"):
        """
        Dump the reweighted THnSpares to a ROOT file.

        Args:
            out_file_name (str): Name of the output ROOT file.
            option (str): File open option.
        """
        if self.cfg['npv']['apply'] or self.cfg['pt']['apply']:
            # Save the weights in the output file
            with ROOT.TFile.Open(out_file_name, option) as output_file:
                if not output_file.GetListOfKeys().Contains("Weights"):
                    output_file.mkdir("Weights")
                output_file.cd("Weights")
                for particle, origin in self.particle_origin:
                    for i_cent, (cent_min, cent_max) in enumerate(zip(self.cent_info.mins, self.cent_info.maxs)):
                        if self.cfg['npv']['apply']:
                            self.histos_pv_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].Write()
                        else:
                            self.histos_pt_weights[f"{particle}{origin}_Cent_{cent_min}_{cent_max}"].Write()
