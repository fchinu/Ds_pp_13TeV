import ROOT

# Open the file
inFile = ROOT.TFile.Open("/home/fchinu/Run3/Ds_pp_13TeV/Datasets/Ds_pp_run3_ml/Data/Train191573/Collisions.root")
hCollisions = inFile.Get("hCollisions")
hCollisions.SetDirectory(0)
inFile.Close()

hFigure  = ROOT.TH1F("hFigure", "Event selection;;Accepted collisions;", 4, 0, 4)
hFigure.SetBinContent(1, hCollisions.GetBinContent(1))
hFigure.SetBinContent(2, hCollisions.GetBinContent(3))
hFigure.SetBinContent(3, hCollisions.GetBinContent(4))
hFigure.SetBinContent(4, hCollisions.GetBinContent(7))

hFigure.SetMarkerStyle(ROOT.kFullDiamond)
hFigure.SetMarkerSize(2)
hFigure.SetMarkerColor(ROOT.kAzure-3)
hFigure.SetLineColor(ROOT.kAzure-3)
hFigure.SetLineWidth(2)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

thesisText = ROOT.TLatex(0.45, 0.75, 'This Thesis')
thesisText.SetNDC()
thesisText.SetTextFont(42)
thesisText.SetTextSize(0.052)


Ncoll = ROOT.TLatex(0.45, 0.69, 'Analysed events: 54.669 #times10^{9}')
Ncoll.SetNDC()
Ncoll.SetTextFont(42)
Ncoll.SetTextSize(0.04)

# Create the canvas
c = ROOT.TCanvas("c", "c", 800, 600)
hFrame = c.DrawFrame(0, 54e9, 4, 60e9, "Event selection;;Accepted collisions")
hFrame.GetXaxis().Set(4, 0, 4)
hFrame.GetXaxis().SetBinLabel(1, "All collisions")
hFrame.GetXaxis().SetBinLabel(2, "Trigger")
hFrame.GetXaxis().SetBinLabel(3, "TF border")
hFrame.GetXaxis().SetBinLabel(4, "#it{z}_{vtx}")
hFrame.GetXaxis().SetLabelSize(0.05)

hFigure.Draw("HIST;same")
hFigure.DrawClone("P;same")
hFigure.SetMarkerStyle(56)
hFigure.SetMarkerColor(ROOT.kBlack)
hFigure.DrawClone("P;same")
thesisText.Draw("same")
Ncoll.Draw("same")
c.SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Figures/Event_Sel/EvSel.pdf")

