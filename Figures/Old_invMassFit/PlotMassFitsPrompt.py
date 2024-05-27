'''
Script for the plots of Ds mass fits
'''
import sys
sys.path.append('/home/fchinu/DmesonAnalysis/utils')
from ROOT import TFile, TCanvas, TLegend, TLatex, TCanvas, TPad # pylint: disable=import-error,no-name-in-module
from ROOT import kBlue, kRed, kBlack # pylint: disable=import-error,no-name-in-module
from StyleFormatter import SetGlobalStyle, SetObjectStyle

inFileDs = TFile.Open('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Old_invMassFit/RawYields.root')
hRawYieldsDs = inFileDs.Get('hRawYields')
hSigmaDs = inFileDs.Get('hRawYieldsSigma')
hMeanDs = inFileDs.Get('hRawYieldsMean')
cMassDs = inFileDs.Get('cMass0')
hInvMassDs24 = cMassDs.GetPad(2).GetListOfPrimitives().FindObject('fHistoInvMass')
fMassTotDs24 = cMassDs.GetPad(2).GetListOfPrimitives().FindObject('funcmass')
fMassBkgDs24 = cMassDs.GetPad(2).GetListOfPrimitives().FindObject('funcbkgrefit')
hRawYieldsDs.SetDirectory(0)
hSigmaDs.SetDirectory(0)
hMeanDs.SetDirectory(0)
hInvMassDs24.SetDirectory(0)
SetObjectStyle(hInvMassDs24, color=kBlack, markersize=1.7, linewidth=3)
SetObjectStyle(fMassTotDs24, linewidth=6, linecolor=kBlue)
SetObjectStyle(fMassBkgDs24, linewidth=6, linecolor=kRed, linestyle=2)
fMassTotDs24.SetNpx(1000)
fMassBkgDs24.SetNpx(1000)

latALICE = TLatex()
latALICE.SetNDC()
latALICE.SetTextSize(0.06)
latALICE.SetTextFont(42)
latALICE.SetTextColor(kBlack)

latLabel = TLatex()
latLabel.SetNDC()
latLabel.SetTextSize(0.04)
latLabel.SetTextFont(42)
latLabel.SetTextColor(kBlack)
SetGlobalStyle(padbottommargin=0.14, padleftmargin=0.15,padrightmargin=0.05, padtopmargin=0.05, titleoffsety=1.4, maxdigits=3)

cMassFitDs24 = TCanvas('cMassFitDs24', '', 1000, 1000)
hFrameDs = cMassFitDs24.DrawFrame(1.75, 0, 2.1, 4000, \
    f';#it{{M}}(KK#pi) (GeV/#it{{c}}^{{2}}); Counts per {hInvMassDs24.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}')
hFrameDs.GetYaxis().SetDecimals()
hInvMassDs24.Draw('same')
fMassBkgDs24.Draw('same')
fMassTotDs24.Draw('same')
latALICE.DrawLatex(0.19, 0.88, 'This Thesis')
latLabel.DrawLatex(0.19, 0.82, 'pp, #sqrt{#it{s}} = 13.6 TeV')
latLabel.DrawLatex(0.19, 0.76, 'D_{s}^{+},D^{+}#rightarrow#phi#pi^{+}#rightarrowK^{+}K^{#font[122]{-}}#pi^{+}')
latLabel.DrawLatex(0.19, 0.70, 'and charge conjugate')
latLabel.DrawLatex(0.19, 0.64, '1.5 < #it{p}_{T} < 2 GeV/#it{c}')
#latLabel.DrawLatex(0.47, 0.3, f'#mu = ({hMeanDs.GetBinContent(1):.4f} #pm {hMeanDs.GetBinError(1):.4f}) GeV/#it{{c}}^{{2}}')
#latLabel.DrawLatex(0.47, 0.25, f'#sigma = {hSigmaDs.GetBinContent(1)*1000:.0f} MeV/#it{{c}}^{{2}} (fixed)')
#latLabel.DrawLatex(0.47, 0.2, f'S = {hRawYieldsDs.GetBinContent(1):.0f} #pm {hRawYieldsDs.GetBinError(1):.0f}')
#latLabel.DrawLatex(0.59, 0.23, f'#sigma = ({hSigmaDs.GetBinContent(1)*1000:.0f} #pm {hSigmaDs.GetBinError(1)*1000:.0f}) MeV/#it{{c}}^{{2}}')
cMassFitDs24.Update()



cMassFitDs24.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Old_invMassFit/InvMassFitDs1p5_2.pdf')
cMassFitDs24.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Old_invMassFit/InvMassFitDs1p5_2.eps')
cMassFitDs24.SaveAs('/home/fchinu/Run3/Ds_pp_13TeV/Figures/Old_invMassFit/InvMassFitDs1p5_2.png')

input('Press enter to exit')