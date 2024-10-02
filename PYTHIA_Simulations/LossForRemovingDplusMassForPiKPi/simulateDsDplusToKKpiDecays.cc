#include <TFile.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TDatabasePDG.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "Pythia8/Pythia.h"

//__________________________________________________________________________________________________
void fillHistoMassPiKPi(TH2D* histo, const float& pt, Pythia8::Pythia& pythia) {

  auto idPart = pythia.event[1].id();
  pythia.particleData.mayDecay(idPart, true);
  pythia.moreDecays();

  std::vector<ROOT::Math::PxPyPzMVector> fourVecPi{};
  ROOT::Math::PxPyPzMVector fourVecKa;
  std::vector<ROOT::Math::PxPyPzMVector> fourVecPiAsKa{};

  float massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  float massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();


  int iPi{0};
  for (int iPart{1}; iPart < pythia.event.size(); ++iPart)
  {
    if (iPi > 1) { // something went wrong with forced decays
      return;
    }
    auto absPdg = std::abs(pythia.event[iPart].id());
    if (absPdg == 211) { // pi
      fourVecPi.push_back(ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), massPi));
      fourVecPiAsKa.push_back(ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), massKa));
      iPi++;
    } else if (absPdg == 321) {
      fourVecKa = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), massKa);
    }
  }

  if (fourVecPi.size() != 1) { // something went wrong with forced decays
    return;
  }

  std::vector<ROOT::Math::PxPyPzMVector> pAsKKPi{ fourVecPi[0] + fourVecPiAsKa[1] + fourVecKa, fourVecPi[1] + fourVecPiAsKa[0] + fourVecKa };histo->Fill(pt, pAsKKPi[0].M());
  histo->Fill(pt, pAsKKPi[1].M());
}

//__________________________________________________________________________________________________
void simulateDsDplusToKKpiDecays(int nDecays, int seed) {

  gROOT->SetBatch(true);

  // init pythia object
  Pythia8::Pythia pythia;
  pythia.readString("SoftQCD:inelastic = on");
  pythia.readString("Tune:pp = 14");
  pythia.readString("411:onMode = off");
  pythia.readString("431:onMode = off");
  pythia.readString("411:onIfMatch = 321 321 211"); // D+ -> KKPi
  pythia.readString("431:onIfMatch = 321 321 211"); // Ds+ -> KKPi
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed %d", seed));
  pythia.init();

  gRandom->SetSeed(seed);

  TFile outFile(Form("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/LossForRemovingDplusMassForPiKPi/DsDplusToKKpi_pythia8_seed%d.root", seed), "recreate");
  TH2D* histoDsAsPiKPi = new TH2D("histoDsAsPiKPi", "histoDsAsPiKPi; #it{p}_{T} (GeV/#it{c}); Mass (GeV/#it{c^{2}})", 100, 0, 20, 500, 0, 5);
  TH2D* histoDplusAsPiKPi = new TH2D("histoDplusAsPiKPi", "histoDplusAsPiKPi; #it{p}_{T} (GeV/#it{c}); Mass (GeV/#it{c^{2}})", 100, 0, 20, 500, 0, 5);

  //__________________________________________________________
  // perform the simulation
  auto fDecay = new TF1("fDecay", "exp(-x/0.312)", 0., 1000.);
  float massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  float massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();

  for (auto iDecay{0}; iDecay < nDecays; ++iDecay) {
    auto ptD = gRandom->Uniform(0., 24.);
    auto yD = gRandom->Uniform(-.5, .5);
    auto phiD = gRandom->Rndm() * 2 * TMath::Pi();
    auto pxD = ptD * TMath::Cos(phiD);
    auto pyD = ptD * TMath::Sin(phiD);

    auto mtDplus = TMath::Sqrt(massDplus * massDplus + ptD * ptD);
    auto mtDs = TMath::Sqrt(massDs * massDs + ptD * ptD);
    auto pzDplus = mtDplus * TMath::SinH(yD);
    auto pzDs = mtDs * TMath::SinH(yD);
    auto pDplus = TMath::Sqrt(ptD * ptD + pzDplus * pzDplus);
    auto pDs = TMath::Sqrt(ptD * ptD + pzDs * pzDs);
    auto eDplus = TMath::Sqrt(massDplus * massDplus + pDplus * pDplus);
    auto eDs = TMath::Sqrt(massDs * massDs + pDs * pDs);

    // Dplus
    Pythia8::Particle Dplus;
    Dplus.id(411);
    Dplus.status(81);
    Dplus.m(massDplus);
    Dplus.xProd(0.);
    Dplus.yProd(0.);
    Dplus.zProd(0.);
    Dplus.tProd(0.);
    Dplus.e(eDplus);
    Dplus.px(pxD);
    Dplus.py(pyD);
    Dplus.pz(pzDplus);
    Dplus.tau(fDecay->GetRandom());

    // Ds
    Pythia8::Particle Ds;
    Ds.id(411);
    Ds.status(81);
    Ds.m(massDs);
    Ds.xProd(0.);
    Ds.yProd(0.);
    Ds.zProd(0.);
    Ds.tProd(0.);
    Ds.e(eDs);
    Ds.px(pxD);
    Ds.py(pyD);
    Ds.pz(pzDs);
    Ds.tau(fDecay->GetRandom());

    // Ds
    pythia.event.reset();
    pythia.event.append(Ds);

    fillHistoMassPiKPi(histoDsAsPiKPi, ptD, pythia);

    // Dplus
    pythia.event.reset();
    pythia.event.append(Dplus);

    fillHistoMassPiKPi(histoDplusAsPiKPi, ptD, pythia);
  }

  // save root output file
  outFile.cd();
  histoDsAsPiKPi->Write();
  histoDplusAsPiKPi->Write();
  outFile.Close();
}
