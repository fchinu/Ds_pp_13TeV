#include <iostream>
#include <vector>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TStyle.h"

std::vector<bool> GeneratePiK(int nTracks, double fracPi) {
  std::vector<bool> isPi(nTracks, false);
  for (int i = 0; i < nTracks; ++i) {
    if (gRandom->Rndm() < fracPi) isPi[i] = true;
  }
  return isPi;
}

std::vector<bool> DoCombinatorial(std::vector<bool> Prong0, std::vector<bool> Prong1, std::vector<bool> Prong2) {
    std::vector<bool> isPiKK;
    for (std::vector<bool>::size_type i = 0; i < Prong0.size(); ++i)
        for (std::vector<bool>::size_type j = 0; j < Prong1.size(); ++j)
            for (std::vector<bool>::size_type k = i; k < Prong2.size(); ++k) {
                if (Prong0[i] && !Prong1[j] && !Prong2[k]) {
                    isPiKK.push_back(true);
                } else if (!Prong0[i] && !Prong1[j] && Prong2[k]) {
                    isPiKK.push_back(false);
                }
            }
    return isPiKK;
}


std::vector<long> Event(int nTracks, double fracPi) {
    std::vector<bool> Prong0 = GeneratePiK(nTracks, fracPi);
    std::vector<bool> Prong1 = GeneratePiK(nTracks, fracPi);
    std::vector<bool> isPiKK = DoCombinatorial(Prong0, Prong1, Prong0);
    //PiKK, KKPi
    return {std::count(isPiKK.begin(), isPiKK.end(), true), std::count(isPiKK.begin(), isPiKK.end(), false)};
}

void ToyMC(int nEvents, double fracPi) {
    TRandom3 rndm;
    rndm.SetSeed(42);
    gRandom = &rndm;
    TH2D *h = new TH2D("h", "h", 100, -5.e4, 5.e4, 100, 0, 400);
    long totnPiKK = 0;
    long totnKKPi = 0;
    for (int i = 0; i < nEvents; ++i) {
        int nTracks = gRandom->Uniform(10, 400);
        std::vector<long> nPiKKEvent = Event(nTracks, fracPi);
        long nPiKK = nPiKKEvent[0];
        long nKKPi = nPiKKEvent[1];
        totnPiKK += nPiKK;
        totnKKPi += nKKPi;
        h->Fill(nPiKK - nKKPi, nTracks);
    }
    std::cout << "nPiKK = " << totnPiKK << ", nKKPi = " << totnKKPi << std::endl;
    TCanvas *cDifference = new TCanvas("cDifference", "cDifference", 800, 600);
    gStyle->SetOptStat(0);
    h->Draw("colz");
    cDifference->SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Tests/CompareKKPi_PiKK/PiKK_KKPi_Difference_Vs_NTracks.png"); 
    TCanvas *cProjection = new TCanvas("cProjection", "cProjection", 800, 600);
    TH1F *hProjection = (TH1F*)h->ProjectionX();
    hProjection->Draw();
    cProjection->SaveAs("/home/fchinu/Run3/Ds_pp_13TeV/Tests/CompareKKPi_PiKK/PiKK_KKPi_Difference_Projection.png");
}

