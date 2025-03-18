#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#endif
#include <iostream>

void CreateTree(TString filnam="charmHadr.txt", int ispp=kFALSE, float ymin=-9999, float ymax=9999){
  //std::cout << "CreateTree(" << filnam << ", " << ispp << ", " << ymin << ", " << ymax << ")" << std::endl;
  const int maxFiles=1090;
  TString nameFils[maxFiles];
  int firstJob=1;
  int lastJob=1090;
  int nFiles=0;
  for(int jf=firstJob; jf<=lastJob; jf++){
    TString command=Form("ls JobPbPb%d/%s > /dev/null 2>&1",jf,filnam.Data());
    if(gSystem->Exec(command.Data())==0){
      nameFils[nFiles++]=Form("JobPbPb%d/%s",jf,filnam.Data());
    }
  }
  printf("Number of files to be read = %d\n",nFiles);

  char line[200];
  const int maxHad=200;
  int id,btag;
  float px,py,pz,en;
  int arrayId[maxHad],arrayBtag[maxHad];
  float arrayPx[maxHad],arrayPy[maxHad],arrayPz[maxHad],arrayEn[maxHad];
  float pt,y;
  int jEv,nChTot,nChMid,nChFwd,nPhysPrimMid,nPhysPrimFwd;
  int totEvents=0;
  float avexSec=0;
  TFile* outfil=new TFile("CharmHadrons_merged.root","recreate");
  int nMulBins=101;
  double minMul=-0.5;
  double maxMul=100.5;
  int nDiffBins=21;
  double minDiff=-10.5;
  double maxDiff=10.5;
  if(!ispp){
    nMulBins=100;
    minMul=0;
    maxMul=2500;
    nDiffBins=20;
    minDiff=-100;
    maxDiff=100;
  }
  TH2F* hMultCorrel = new TH2F("hMultCorrel"," ; n_{ch}^{|#eta|<0.8} ; n_{ch}^{fwd}",nMulBins,minMul,maxMul,nMulBins,minMul,maxMul);
  TH1F* hMultDiffMid = new TH1F("hMultDiffMid", " ; #Deltan_{ch} (meth1-meth2)",nDiffBins,minDiff,maxDiff);
  TH1F* hMultDiffFwd = new TH1F("hMultDiffFwd", " ; #Deltan_{ch} (meth1-meth2)",nDiffBins,minDiff,maxDiff);
  TTree* treeHadr=new TTree("treeHadr","hadrons");
  treeHadr->Branch("pdgCode",&id,"pdgCode/I");
  treeHadr->Branch("promptTag",&btag,"promptTag/I");
  //  treeHadr->Branch("px",&px,"px/F");
  //  treeHadr->Branch("py",&py,"py/F");
  //  treeHadr->Branch("pz",&pz,"pz/F");
  treeHadr->Branch("pt",&pt,"pt/F");
  treeHadr->Branch("y",&y,"y/F");
  treeHadr->Branch("nChMid",&nChMid,"nChMid/I");
  treeHadr->Branch("nChFwd",&nChFwd,"nChFwd/I");
  //  treeHadr->Branch("energy",&en,"energy/F");

  for(int jf=0; jf<nFiles; jf++){
    FILE* infilpPb=fopen(nameFils[jf].Data(),"r");
    cout<<infilpPb<<endl;
    int nD0=0;
    int nChHad=0;
    while(1){
      fgets(line,200,infilpPb);
      if(strstr(line,"-----")!=0) break;
      if(strstr(line,"Event")!=0){
       	sscanf(line,"Event %d",&jEv);
	nChHad=0;
       	printf("---> Event %d\n",jEv);
      }else if(strstr(line,"nCharged")!=0){
       	sscanf(line,"nCharged %d %d %d %d %d",&nChTot,&nChMid,&nChFwd,&nPhysPrimMid,&nPhysPrimFwd);
	hMultCorrel->Fill(nChMid,nChFwd);
	hMultDiffMid->Fill(nChMid-nPhysPrimMid);
	hMultDiffFwd->Fill(nChFwd-nPhysPrimFwd);
	//printf("Charm hadrons = %d Mults mid: %d %d fwd: %d %d\n",nChHad,nChMid,nPhysPrimMid,nChFwd,nPhysPrimFwd);
	for(int ih=0; ih<nChHad; ih++){
	  id=arrayId[ih];
	  btag=arrayBtag[ih];
	  pt=TMath::Sqrt(arrayPx[ih]*arrayPx[ih]+arrayPy[ih]*arrayPy[ih]);
	  y=0.5*TMath::Log((arrayEn[ih]+arrayPz[ih])/(arrayEn[ih]-arrayPz[ih]));
	  if(y>ymin && y<ymax) treeHadr->Fill();
	}
      }else{
	sscanf(line,"%d %d %f %f %f %f",&id,&btag,&px,&py,&pz,&en);
	arrayId[nChHad]=id;
	arrayBtag[nChHad]=btag;
	arrayPx[nChHad]=px;
	arrayPy[nChHad]=py;
	arrayPz[nChHad]=pz;
	arrayEn[nChHad]=en;
	//printf("Ch Hadr %d %d %d %f\n",nChHad,id,btag,px);
	++nChHad;
	if(TMath::Abs(id)==421) nD0++;
      }
    }
    int nEvents;
    float xSec;
    fgets(line,200,infilpPb);
    if(strstr(line,"sectionHI")!=0){
      float accEv,xSecDum,nAccDum;
      sscanf(line,"N events = %d  Cross sectionHI = %f nAcceptedHI = %f Cross section = %f nAccepted = %f\n",&nEvents,&xSecDum,&nAccDum,&xSec,&accEv);
      if(TMath::Abs(accEv-nEvents)>0.5){ 
	printf("MISMATCH in event count: %d vs %f\n",nEvents,accEv);
	nEvents=(int)(accEv+0.001);
      }      
    }else{
      float accEv;
      sscanf(line,"N events = %d  Cross section = %f nAccepted = %f",&nEvents,&xSec,&accEv);
      if(TMath::Abs(accEv-nEvents)>0.5){ 
	printf("MISMATCH in event count: %d vs %f\n",nEvents,accEv);
	nEvents=(int)(accEv+0.001);
      }
    }
    printf("File %d %s Number of events: %d xSec=%f Number of D0=%d\n",jf,nameFils[jf].Data(),nEvents,xSec,nD0);
    totEvents+=nEvents;
    avexSec+=xSec*nEvents;
    fclose(infilpPb);
  }
  std::cout << "here";
  avexSec/=(float)totEvents;
  printf("Total number of events=%d   <xSec>=%f\n",totEvents,avexSec);
  TH1F* hNorm=new TH1F("hNorm","",2,-0.5,1.5);
  hNorm->GetXaxis()->SetBinLabel(1,"N Events");
  hNorm->GetXaxis()->SetBinLabel(2,"Cross sec");
  hNorm->SetBinContent(1,totEvents);
  hNorm->SetBinContent(2,avexSec);


  outfil->cd();
  treeHadr->Write();
  hNorm->Write();
  hMultCorrel->Write();
  hMultDiffMid->Write();
  hMultDiffFwd->Write();
  outfil->Close();
}
