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
void simulateLcTopKpiDecays(int nDecays, int seed) {

  gROOT->SetBatch(true);

  // init pythia object
  Pythia8::Pythia pythia;
  pythia.readString("SoftQCD:inelastic = on");
  pythia.readString("Tune:pp = 14");
  pythia.readString("4122:onMode = off");
  pythia.readString("4122:addChannel = 1 1 0 2212 -321 211");
  pythia.readString("4122:onIfMatch = 2212 321 211");
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed %d", seed));
  pythia.init();

  gRandom->SetSeed(seed);

  TFile outFile(Form("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/LcToPKPiContamination/lc_to_pkpi_pythia8_seed%d.root", seed), "recreate");
  TNtuple* ntupleD = new TNtuple("treeLc", "treeLc", "ptLc:yLc:ptPr:ptKa:ptPi:massPrKaPi:massKaKaPi");

  //__________________________________________________________
  // perform the simulation
  auto fDecay = new TF1("fDecay", "exp(-x/0.312)", 0., 1000.);
  float massLc = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  float massPr = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  float massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  float massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  for (auto iDecay{0}; iDecay < nDecays; ++iDecay) {
    auto ptLc = gRandom->Uniform(0., 24.);
    auto yLc = gRandom->Uniform(-.5, .5);
    auto phiLc = gRandom->Rndm() * 2 * TMath::Pi();
    auto pxLc = ptLc * TMath::Cos(phiLc);
    auto pyLc = ptLc * TMath::Sin(phiLc);
    auto mt = TMath::Sqrt(massLc * massLc + ptLc * ptLc);
    auto pzLc = mt * TMath::SinH(yLc);
    auto pLc = TMath::Sqrt(ptLc * ptLc + pzLc * pzLc);
    auto eLc = TMath::Sqrt(massLc * massLc + pLc * pLc);

    // Lc
    Pythia8::Particle Lc;
    Lc.id(4122);
    Lc.status(81);
    Lc.m(massLc);
    Lc.xProd(0.);
    Lc.yProd(0.);
    Lc.zProd(0.);
    Lc.tProd(0.);
    Lc.e(eLc);
    Lc.px(pxLc);
    Lc.py(pyLc);
    Lc.pz(pzLc);
    Lc.tau(fDecay->GetRandom());

    pythia.event.reset();
    pythia.event.append(Lc);
    auto idPart = pythia.event[1].id();
    pythia.particleData.mayDecay(idPart, true);
    pythia.moreDecays();

    ROOT::Math::PxPyPzMVector fourVecPr;
    ROOT::Math::PxPyPzMVector fourVecPi;
    ROOT::Math::PxPyPzMVector fourVecKa;
    ROOT::Math::PxPyPzMVector fourVecPrAsKa;
    int iPr{0};
    int iKa{0};
    int iPi{0};
    for (int iPart{1}; iPart < pythia.event.size(); ++iPart)
    {
      auto absPdg = std::abs(pythia.event[iPart].id());
      if (absPdg == 2212) { // pr
        fourVecPr = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), massPr);
        fourVecPrAsKa = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), massKa);
        iPr++;
      } else if (absPdg == 321) {
        fourVecKa = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), massKa);
        iKa++;
      } else if (absPdg == 211) {
        fourVecPi = ROOT::Math::PxPyPzMVector(pythia.event[iPart].px(), pythia.event[iPart].py(), pythia.event[iPart].pz(), massPi);
        iPi++;
      }
    }

    if (iPi != 1 || iKa != 1 || iPr != 1) {
      continue;
    }
    ROOT::Math::PxPyPzMVector fourVecPrKaPi = fourVecPr + fourVecKa + fourVecPi;
    ROOT::Math::PxPyPzMVector fourVecKaKaPi = fourVecPrAsKa + fourVecKa + fourVecPi;
    ntupleD->Fill(ptLc, yLc, fourVecPr.Pt(), fourVecKa.Pt(), fourVecPi.Pt(), fourVecPrKaPi.M(), fourVecKaKaPi.M());
  }

  // save root output file
  ntupleD->Write();
  outFile.Close();
}
