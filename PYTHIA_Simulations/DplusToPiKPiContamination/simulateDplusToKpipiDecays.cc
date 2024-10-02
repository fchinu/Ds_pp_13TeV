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
void simulateDplusToKpipiDecays(int nDecays, int seed) {

  gROOT->SetBatch(true);

  // init pythia object
  Pythia8::Pythia pythia;
  pythia.readString("SoftQCD:inelastic = on");
  pythia.readString("Tune:pp = 14");
  pythia.readString("411:onMode = off");
  pythia.readString("411:onIfMatch = 321 211 211");
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed %d", seed));
  pythia.init();

  gRandom->SetSeed(seed);

  TFile outFile(Form("/home/fchinu/Run3/Ds_pp_13TeV/Projections_RawYields/SimulateDplus/dplus_to_Kpipi_pythia8_seed%d.root", seed), "recreate");
  TNtuple* ntupleD = new TNtuple("treeD", "treeD", "ptD:yD:ptPi0:ptPi1:ptKa:etaPi0:etaPi1:etaKa:massPiPi:massKaPiPi:massKaKaPi");

  //__________________________________________________________
  // perform the simulation
  auto fDecay = new TF1("fDecay", "exp(-x/0.312)", 0., 1000.);
  float massD = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  float massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  float massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  for (auto iDecay{0}; iDecay < nDecays; ++iDecay) {
    auto ptD = gRandom->Uniform(0., 24.);
    auto yD = gRandom->Uniform(-.5, .5);
    auto phiD = gRandom->Rndm() * 2 * TMath::Pi();
    auto pxD = ptD * TMath::Cos(phiD);
    auto pyD = ptD * TMath::Sin(phiD);
    auto mt = TMath::Sqrt(massD * massD + ptD * ptD);
    auto pzD = mt * TMath::SinH(yD);
    auto pD = TMath::Sqrt(ptD * ptD + pzD * pzD);
    auto eD = TMath::Sqrt(massD * massD + pD * pD);

    // Dplus
    Pythia8::Particle Dplus;
    Dplus.id(411);
    Dplus.status(81);
    Dplus.m(massD);
    Dplus.xProd(0.);
    Dplus.yProd(0.);
    Dplus.zProd(0.);
    Dplus.tProd(0.);
    Dplus.e(eD);
    Dplus.px(pxD);
    Dplus.py(pyD);
    Dplus.pz(pzD);
    Dplus.tau(fDecay->GetRandom());

    pythia.event.reset();
    pythia.event.append(Dplus);
    auto idPart = pythia.event[1].id();
    pythia.particleData.mayDecay(idPart, true);
    pythia.moreDecays();

    std::vector<ROOT::Math::PxPyPzMVector> fourVecPi{};
    ROOT::Math::PxPyPzMVector fourVecKa;
    std::vector<ROOT::Math::PxPyPzMVector> fourVecPiAsKa{};
    int iPi{0};
    for (int iPart{1}; iPart < pythia.event.size(); ++iPart)
    {
      if (iPi > 2) { // something went wrong with forced decays
        continue;
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
    if (fourVecPi.size() != 2) { // something went wrong with forced decays
      continue;
    }

    ROOT::Math::PxPyPzMVector fourVecPiPi = fourVecPi[0] + fourVecPi[1];
    ROOT::Math::PxPyPzMVector fourVecKaPiPi = fourVecKa + fourVecPi[0] + fourVecPi[1];
    ROOT::Math::PxPyPzMVector fourVecKaKaPi = fourVecKa + fourVecPiAsKa[0] + fourVecPi[1];
    ntupleD->Fill(ptD, yD, fourVecPi[0].Pt(), fourVecPi[1].Pt(), fourVecKa.Pt(), fourVecPi[0].Eta(), fourVecPi[1].Eta(),  fourVecKa.Eta(), fourVecPiPi.M(), fourVecKaPiPi.M(), fourVecKaKaPi.M());
  }

  // save root output file
  ntupleD->Write();
  outFile.Close();
}
