import ROOT

StdFileName = "SHM_PbPb010_5dot02TeV_DsOverD0.root"
BarEnhFileName = "SHM_baryon_enh_PbPb010_5dot02TeV_DsOverD0.root"

inFileStd = ROOT.TFile.Open(StdFileName)
hist_DsOverD0 = inFileStd.Get("Ds1968plu")
hist_DsOverD0.SetName("hist_DsOverD0")
hist_DpOverD0 = inFileStd.Get("Dc1869plu")
hist_DpOverD0.SetName("hist_DpOverD0")
hist_DsOverDp = hist_DsOverD0.Clone("hist_DsOverDp")

inFileBarEnh = ROOT.TFile.Open(BarEnhFileName)
hist_DsOverD0_BarEnh = inFileBarEnh.Get("Ds1968plu")
hist_DsOverD0_BarEnh.SetName("hist_DsOverD0_BarEnh")
hist_DpOverD0_BarEnh = inFileBarEnh.Get("Dc1869plu")
hist_DpOverD0_BarEnh.SetName("hist_DpOverD0_BarEnh")
hist_DsOverDp_BarEnh = hist_DsOverD0.Clone("hist_DsOverDp_BarEnh")

for iPt in range(1, hist_DsOverD0.GetNbinsX()+1):
    DsOverD0 = hist_DsOverD0.GetBinContent(iPt)
    DsOverD0Unc = hist_DsOverD0.GetBinError(iPt)
    DsOverD0_BarEnh = hist_DsOverD0_BarEnh.GetBinContent(iPt)
    DsOverD0Unc_BarEnh = hist_DsOverD0_BarEnh.GetBinError(iPt)
    DpOverD0 = hist_DpOverD0.GetBinContent(iPt)
    DpOverD0_BarEnh = hist_DpOverD0_BarEnh.GetBinContent(iPt)
    DsOverDp = DsOverD0 / DpOverD0
    DsOverDp_BarEnh = DsOverD0_BarEnh / DpOverD0_BarEnh
    hist_DsOverDp.SetBinContent(iPt, DsOverDp)
    hist_DsOverDp.SetBinError(iPt, DsOverDp*DsOverD0Unc/DsOverD0)
    hist_DsOverDp_BarEnh.SetBinContent(iPt, DsOverDp_BarEnh)
    hist_DsOverDp_BarEnh.SetBinError(iPt, DsOverDp_BarEnh*DsOverD0Unc_BarEnh/DsOverD0_BarEnh)

outFile = ROOT.TFile("SHM_PbPb010_5dot02TeV_DsOverDp.root", "recreate")
hist_DsOverD0.Write()
hist_DpOverD0.Write()
hist_DsOverDp.Write()
hist_DsOverD0_BarEnh.Write()
hist_DpOverD0_BarEnh.Write()
hist_DsOverDp_BarEnh.Write()
outFile.Close()
