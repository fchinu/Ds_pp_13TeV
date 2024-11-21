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
void testTofTime(int nDecays, int seed) {

  gROOT->SetBatch(true);

  // init pythia object
  Pythia8::Pythia pythia;
  pythia.readString("SoftQCD:inelastic = on");
  pythia.readString("Tune:pp = 14");
  pythia.readString("3122:onMode = off");
  pythia.readString("3122:onIfMatch = 2212 -211");
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed %d", seed));
  pythia.init();

  gRandom->SetSeed(seed);
  const double R = 3700; // mm

  TFile outFile(Form("/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/TOF_time_direct_LF/testTofTime_seed%d.root", seed), "recreate");
  TNtuple* ntupleLam = new TNtuple("treeFromLambda", "treeFromLambda", "ptLam:ptPi:ptPr:phiLam:phiPi:phiPr:timeLam:xLam:yLam:zLam:timePi:xPr:yPr:timePr");
  TNtuple* ntupleDirectPi = new TNtuple("treeDirectPi", "treeDirect", "ptPi:phiPi:timePi:xPi:yPi");
  TNtuple* ntupleDirectPr = new TNtuple("treeDirectPr", "treeDirect", "ptPr:phiPr:timePr:xPr:yPr");

  //__________________________________________________________
  // perform the simulation
  auto fDecay = new TF1("fDecay", "exp(-x/78)", 0., 1000.);
  float massLam = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  float massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  float massPr = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  for (auto iDecay{0}; iDecay < nDecays; ++iDecay) {
    auto ptLam = gRandom->Uniform(0., 24.);
    auto rapLam = gRandom->Uniform(-.5, .5);
    auto phiLambda = gRandom->Rndm() * 2 * TMath::Pi();
    auto pxLam = ptLam * TMath::Cos(phiLambda);
    auto pyLam = ptLam * TMath::Sin(phiLambda);
    auto mt = TMath::Sqrt(massLam * massLam + ptLam * ptLam);
    auto pzLam = mt * TMath::SinH(rapLam);
    auto pLam = TMath::Sqrt(ptLam * ptLam + pzLam * pzLam);
    auto eLam = TMath::Sqrt(massLam * massLam + pLam * pLam);

    // Lambda
    Pythia8::Particle Lambda;
    Lambda.id(3122);
    Lambda.status(81);
    Lambda.m(massLam);
    Lambda.xProd(0.);
    Lambda.yProd(0.);
    Lambda.zProd(0.);
    Lambda.tProd(0.);
    Lambda.e(eLam);
    Lambda.px(pxLam);
    Lambda.py(pyLam);
    Lambda.pz(pzLam);
    Lambda.tau(fDecay->GetRandom());

    pythia.event.reset();
    pythia.event.append(Lambda);
    auto idPart = pythia.event[1].id();
    pythia.particleData.mayDecay(idPart, true);
    pythia.moreDecays();

    float timeLam{0.}, xLam{0.}, yLam{0.}, zLam{0.}, phiLam{0.};
    float timePi{0.}, phiPi{0.}, ptPi{0.}, xPi{0.}, yPi{0.};
    float timePr{0.}, phiPr{0.}, ptPr{0.}, xPr{0.}, yPr{0.};
    int iLam{0};
    int iPi{0};
    int iPr{0};
    for (int iPart{1}; iPart < pythia.event.size(); ++iPart)
    {
      auto absPdg = std::abs(pythia.event[iPart].id());
      if (absPdg == 3122) { // lambda

        timeLam = pythia.event[iPart].tDec();
        xLam = pythia.event[iPart].xDec();
        yLam = pythia.event[iPart].yDec();
        zLam = pythia.event[iPart].zDec();
        phiLam = pythia.event[iPart].phi();
        iLam++;
      } else if (absPdg == 211) { // pi

        double x0 = pythia.event[iPart].xProd(); // Production vertex x
        double y0 = pythia.event[iPart].yProd(); // Production vertex y
        double z0 = pythia.event[iPart].zProd(); // Production vertex z
        // Calculate coefficients of quadratic equation
        double vx = pythia.event[iPart].px() / pythia.event[iPart].m();
        double vy = pythia.event[iPart].py() / pythia.event[iPart].m();
        double vz = pythia.event[iPart].pz() / pythia.event[iPart].m();
        double a = vx * vx + vy * vy;
        double b = 2 * (x0 * vx + y0 * vy);
        double c = x0 * x0 + y0 * y0 - R * R;
        double discriminant = b * b - 4 * a * c;
        double t = (-b + std::sqrt(discriminant)) / (2 * a);

        timePi = t; //pythia.event[iPart].tDec();
        xPi = vx * t + x0;
        yPi = vy * t + y0;
        phiPi = pythia.event[iPart].phi();
        ptPi = pythia.event[iPart].pT();
        iPi++;      
      } else if (absPdg == 2212) { // pr

        double x0 = pythia.event[iPart].xProd(); // Production vertex x
        double y0 = pythia.event[iPart].yProd(); // Production vertex y
        double z0 = pythia.event[iPart].zProd(); // Production vertex z
        // Calculate coefficients of quadratic equation
        double vx = pythia.event[iPart].px() / pythia.event[iPart].m();
        double vy = pythia.event[iPart].py() / pythia.event[iPart].m();
        double vz = pythia.event[iPart].pz() / pythia.event[iPart].m();
        double a = vx * vx + vy * vy;
        double b = 2 * (x0 * vx + y0 * vy);
        double c = x0 * x0 + y0 * y0 - R * R;
        double discriminant = b * b - 4 * a * c;
        double t = (-b + std::sqrt(discriminant)) / (2 * a);

        timePr = t; //pythia.event[iPart].tDec();
        xPr = vx * t + x0;
        yPr = vy * t + y0;
        phiPr = pythia.event[iPart].phi();
        ptPr = pythia.event[iPart].pT();
        iPr++;
      }
    }
    if (iLam != 1 || iPi != 1 || iPr != 1) { // something went wrong with forced decays
      std::cout << "iLam = " << iLam << ", iPi = " << iPi << ", iPr = " << iPr << std::endl;
      continue;
    }
    timePi = timePi + timeLam;
    timePr = timePr + timeLam;
    ntupleLam->Fill(ptLam, ptPi, ptPr, phiLam, phiPi, phiPr, timeLam, xLam, yLam, zLam, timePi, xPr, yPr, timePr);
  }
  std::cout << "Done with the simulation" << std::endl;

  // save root output file
  ntupleLam->Write();

  for (auto iDecay{0}; iDecay < nDecays; ++iDecay) {
    auto ptParticle = gRandom->Uniform(0., 24.);
    auto rapParticle = gRandom->Uniform(-.5, .5);
    auto phiParticle = gRandom->Rndm() * 2 * TMath::Pi();
    auto pxParticle = ptParticle * TMath::Cos(phiParticle);
    auto pyParticle = ptParticle * TMath::Sin(phiParticle);

    // Pi
    auto mt = TMath::Sqrt(massPi * massPi + ptParticle * ptParticle);
    auto pzParticle = mt * TMath::SinH(rapParticle);
    auto pParticle = TMath::Sqrt(ptParticle * ptParticle + pzParticle * pzParticle);
    auto eParticle = TMath::Sqrt(massPi * massPi + pParticle * pParticle);

    Pythia8::Particle Pion;
    Pion.id(211);
    Pion.status(81);
    Pion.m(massPi);
    Pion.xProd(0.);
    Pion.yProd(0.);
    Pion.zProd(0.);
    Pion.tProd(0.);
    Pion.e(eParticle);
    Pion.px(pxParticle);
    Pion.py(pyParticle);
    Pion.pz(pzParticle);

    pythia.event.reset();
    pythia.event.append(Pion);
    auto idPart = pythia.event[1].id();
    float timePi{0.}, phiPi{0.}, ptPi{0.}, xPi{0.}, yPi{0.};
    int iPi{0};
    for (int iPart{1}; iPart < pythia.event.size(); ++iPart)
    {
      auto absPdg = std::abs(pythia.event[iPart].id());
      if (absPdg == 211) { // pi

        double x0 = pythia.event[iPart].xProd(); // Production vertex x
        double y0 = pythia.event[iPart].yProd(); // Production vertex y
        double z0 = pythia.event[iPart].zProd(); // Production vertex z
        // Calculate coefficients of quadratic equation
        double vx = pythia.event[iPart].px() / pythia.event[iPart].m();
        double vy = pythia.event[iPart].py() / pythia.event[iPart].m();
        double vz = pythia.event[iPart].pz() / pythia.event[iPart].m();
        double a = vx * vx + vy * vy;
        double b = 2 * (x0 * vx + y0 * vy);
        double c = x0 * x0 + y0 * y0 - R * R;
        double discriminant = b * b - 4 * a * c;
        double t = (-b + std::sqrt(discriminant)) / (2 * a);

        timePi = t; //pythia.event[iPart].tDec();
        xPi = vx * t + x0;
        yPi = vy * t + y0;
        phiPi = pythia.event[iPart].phi();
        ptPi = pythia.event[iPart].pT();
        iPi++;      
      }
    }
    timePi = timePi;
    ntupleDirectPi->Fill(ptPi, phiPi, timePi, xPi, yPi);

    // Pr
    mt = TMath::Sqrt(massPr * massPr + ptParticle * ptParticle);
    pzParticle = mt * TMath::SinH(rapParticle);
    pParticle = TMath::Sqrt(ptParticle * ptParticle + pzParticle * pzParticle);
    eParticle = TMath::Sqrt(massPr * massPr + pParticle * pParticle);

    Pythia8::Particle Proton;
    Proton.id(2212);
    Proton.status(81);
    Proton.m(massPr);
    Proton.xProd(0.);
    Proton.yProd(0.);
    Proton.zProd(0.);
    Proton.tProd(0.);
    Proton.e(eParticle);
    Proton.px(pxParticle);
    Proton.py(pyParticle);
    Proton.pz(pzParticle);

    pythia.event.reset();
    pythia.event.append(Proton);
    idPart = pythia.event[1].id();
    float timePr{0.}, phiPr{0.}, ptPr{0.}, xPr{0.}, yPr{0.};
    int iPr{0};
    for (int iPart{1}; iPart < pythia.event.size(); ++iPart)
    {
      auto absPdg = std::abs(pythia.event[iPart].id());
      if (absPdg == 2212) { // pr

        double x0 = pythia.event[iPart].xProd(); // Production vertex x
        double y0 = pythia.event[iPart].yProd(); // Production vertex y
        double z0 = pythia.event[iPart].zProd(); // Production vertex z
        // Calculate coefficients of quadratic equation
        double vx = pythia.event[iPart].px() / pythia.event[iPart].m();
        double vy = pythia.event[iPart].py() / pythia.event[iPart].m();
        double vz = pythia.event[iPart].pz() / pythia.event[iPart].m();
        double a = vx * vx + vy * vy;
        double b = 2 * (x0 * vx + y0 * vy);
        double c = x0 * x0 + y0 * y0 - R * R;
        double discriminant = b * b - 4 * a * c;
        double t = (-b + std::sqrt(discriminant)) / (2 * a);

        timePr = t; //pythia.event[iPart].tDec();
        xPr = vx * t + x0;
        yPr = vy * t + y0;
        phiPr = pythia.event[iPart].phi();
        ptPr = pythia.event[iPart].pT();
        iPr++;      
      }
    }
    timePr = timePr;
    ntupleDirectPr->Fill(ptPr, phiPr, timePr, xPr, yPr);
  }

  // save root output file
  ntupleDirectPi->Write();
  ntupleDirectPr->Write();
  outFile.Close();

}
