void NormalizeToXsec(TH1F* histo, double normFac){
    for(int ib=1; ib<=histo->GetNbinsX(); ib++){
      double c0=histo->GetBinContent(ib);
      double ec0=TMath::Sqrt(c0)/c0;
      c0*=(normFac/histo->GetBinWidth(ib));
      histo->SetBinContent(ib,c0);
      histo->SetBinError(ib,ec0*c0);
    }
  }
  
  void FillParticleRatio(TH1F* hnum, TH1F* hden, TH1F* hrat){
    for(int ib=1; ib<=hnum->GetNbinsX(); ib++){
      double num=hnum->GetBinContent(ib);
      double ernum=0;
      if(num>0) ernum=hnum->GetBinError(ib)/num;
      double den=hden->GetBinContent(ib);
      double erden=0;
      if(den>0) erden=hden->GetBinError(ib)/den;
      double r=0;
      if(den>0) r=num/den;
      double er=r*TMath::Sqrt(ernum*ernum+erden*erden);
      hrat->SetBinContent(ib,r);
      hrat->SetBinError(ib,er);
    }
    hrat->SetStats(0);
  }
  
  void ComputeRpPb(TH1F* hpPb, TH1F* hpp, TH1F* hRpPb){
    for(int ib=1; ib<=hpPb->GetNbinsX(); ib++){
      double num=hpPb->GetBinContent(ib);
      double den=hpp->GetBinContent(ib);
      double ernum=hpPb->GetBinError(ib)/num;
      double erden=hpp->GetBinError(ib)/den;
      if(den>0){
        double r=num/208./den;
        double er=r*TMath::Sqrt(ernum*ernum+erden*erden);
        hRpPb->SetBinContent(ib,r);
        hRpPb->SetBinError(ib,er);
      }else{
        hRpPb->SetBinContent(ib,0.);
        hRpPb->SetBinError(ib,0.);
      }
    }
  }
  
  
  void PlotRatiosVsMult(TString inpfilnam="CharmHadrons_merged.root",
                TString optPart="charm",
                bool isPbPb=true,
                bool usePercentiles=true){
  
    TFile* infil=new TFile(inpfilnam.Data());
    TH2F* hMultCorrel=(TH2F*)infil->Get("hMultCorrel");
    TH1D* hnChMid=hMultCorrel->ProjectionX();
    TH1D* hnChFwd=hMultCorrel->ProjectionY();
    double aveMid=hnChMid->GetMean();
    double aveFwd=hnChFwd->GetMean();
    hnChMid->GetXaxis()->SetTitle("n_{ch}");
    hnChMid->SetName("hnChMid");
    hnChFwd->SetName("hnChFwd");
    const int percBins=6;
    Double_t percF[percBins]={0};
    Double_t percM[percBins]={0};
    Double_t prob[percBins]={0.,0.8,0.9,0.95,0.98,1.};
    if(isPbPb){
      prob[1]=0.3;
      prob[2]=0.5;
      prob[3]=0.7;
      prob[4]=0.9;
      prob[5]=0.99;
    }
    hnChFwd->GetQuantiles(percBins,percF,prob);
    hnChMid->GetQuantiles(percBins,percM,prob);
    printf("Percentile Limits Mid = ");
    for(int i=0; i<percBins; i++) printf("%f ",percM[i]);
    printf("\n");
    printf("Percentile Limits Fwd = ");
    for(int i=0; i<percBins; i++) printf("%f ",percF[i]);
    printf("\n");
    
    TH1F* hmdm=(TH1F*)infil->Get("hMultDiffMid");
    TH1F* hmdf=(TH1F*)infil->Get("hMultDiffFwd");
  
    TCanvas* cmd = new TCanvas("cmd","",800,800);
    cmd->SetLogy();
    cmd->SetTickx();
    cmd->SetTicky();
    hmdm->SetLineWidth(2);
    hmdm->Draw();
    hmdf->SetLineWidth(2);
    hmdf->SetLineColor(4);
    hmdf->Draw("sames");
    gPad->Update();
    TPaveStats* tpsm=(TPaveStats*)hmdm->GetListOfFunctions()->FindObject("stats");
    tpsm->SetTextColor(hmdm->GetLineColor());
    tpsm->SetY1NDC(0.72);
    tpsm->SetY2NDC(0.92);
    gPad->Modified();
    TPaveStats* tpsf=(TPaveStats*)hmdf->GetListOfFunctions()->FindObject("stats");
    tpsf->SetTextColor(hmdf->GetLineColor());
    tpsf->SetY1NDC(0.51);
    tpsf->SetY2NDC(0.71);
    gPad->Modified();
    cmd->SaveAs("MultCountsDiff.png");
  
    TCanvas* cm = new TCanvas("cm","",1200,600);
    cm->Divide(2,1);
    cm->cd(1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    hMultCorrel->SetStats(0);
    hMultCorrel->Draw("colz");
    TLatex* t1=new TLatex(0.15,0.82,"PYTHIA8, CR mode 2 + ropes, pp, #sqrt{s} = 13.6 TeV");
    t1->SetNDC();
    t1->Draw();
    cm->cd(2);
    gPad->SetLogy();
    hnChMid->SetLineWidth(2);
    hnChMid->Draw();
    hnChFwd->SetLineWidth(2);
    hnChFwd->SetLineColor(4);
    hnChFwd->Draw("sames");
    gPad->Update();
    TPaveStats* tpsm2=(TPaveStats*)hnChMid->GetListOfFunctions()->FindObject("stats");
    tpsm2->SetTextColor(hnChMid->GetLineColor());
    tpsm2->SetY1NDC(0.72);
    tpsm2->SetY2NDC(0.92);
    gPad->Modified();
    TPaveStats* tpsf2=(TPaveStats*)hnChFwd->GetListOfFunctions()->FindObject("stats");
    tpsf2->SetTextColor(hnChFwd->GetLineColor());
    tpsf2->SetY1NDC(0.51);
    tpsf2->SetY2NDC(0.71);
    gPad->Modified();
    double maxheigh;
    for(int i=0; i<percBins; i++){
      maxheigh=1.5*hnChFwd->GetBinContent(hnChFwd->FindBin(percF[i]));
      TLine* lf=new TLine(percF[i],0.,percF[i],maxheigh);
      lf->SetLineStyle(2);
      lf->SetLineColor(hnChFwd->GetLineColor());
      lf->Draw();
      if(prob[i]>0 && prob[i]<1){
        TLatex *tf=new TLatex(percF[i],maxheigh,Form("%.0f",100.-prob[i]*100.));
        tf->SetTextColor(hnChFwd->GetLineColor());
        tf->SetTextFont(43);
        tf->SetTextSize(18);
        tf->Draw();
      }
      maxheigh=1.8*hnChMid->GetBinContent(hnChMid->FindBin(percM[i]));
      TLine* lm=new TLine(percM[i],0.,percM[i],maxheigh);
      lm->SetLineStyle(2);
      lm->SetLineColor(hnChMid->GetLineColor());
      lm->Draw();
      if(prob[i]>0 && prob[i]<1){
        TLatex *tm=new TLatex(percM[i],maxheigh,Form("%.0f",100.-prob[i]*100.));
        tm->SetTextColor(hnChMid->GetLineColor());
        tm->SetTextFont(43);
        tm->SetTextSize(18);
        tm->Draw();
      }
    }
    cm->SaveAs("MultDistr-Percentiles.png");
  
    
    int id,prtag;
    float pt,y;
    int nChMid,nChFwd;
    TTree* treeHadr=(TTree*)infil->Get("treeHadr");
    treeHadr->SetBranchAddress("pdgCode",&id);
    treeHadr->SetBranchAddress("promptTag",&prtag);
    treeHadr->SetBranchAddress("pt",&pt);
    treeHadr->SetBranchAddress("y",&y);
    treeHadr->SetBranchAddress("nChMid",&nChMid);
    treeHadr->SetBranchAddress("nChFwd",&nChFwd);
    TH1F* hNorm=(TH1F*)infil->Get("hNorm");
  
    int nMultBins=6;
    double multLims[7]={-0.5,10.5,20.5,30.5,40.5,50.5,75.5};
    const int nPtBins=6;
    double ptLims[7]={1.,2.,4.,6.,8.,12.,24.};
    double multLimsF[100];
    double multLimsM[100];
    if(usePercentiles){
      nMultBins=percBins-1;
      for(int j=0; j<percBins; j++){
        multLimsF[j]=percF[j];
        multLimsM[j]=percM[j];
      }
    }else{
      for(int j=0; j<nMultBins+1; j++){
        multLimsF[j]=multLims[j];
        multLimsM[j]=multLims[j];
      }
    }
  
    const int maxPart=10;
    int pdgArray[maxPart];
    TString partNames[maxPart];
    TString partTitles[maxPart];
    double ratioMin[maxPart],ratioMax[maxPart];
    int nHadr=4;
    if(optPart=="charm"){
      nHadr=4;
      pdgArray[0]=421;
      pdgArray[1]=411;
      pdgArray[2]=431;
      pdgArray[3]=4122;
      partNames[0]="D0";
      partNames[1]="Dp";
      partNames[2]="Ds";
      partNames[3]="Lc";
      partTitles[0]="D^{0}";
      partTitles[1]="D^{+}";
      partTitles[2]="D_{s}^{+}";
      partTitles[3]="#Lambda_{c}^{+}";
      ratioMin[0]=0.1;
      ratioMax[0]=0.85;
      ratioMin[1]=0.3;
      ratioMax[1]=0.65;
      ratioMin[2]=0.05;
      ratioMax[2]=0.41;
      ratioMin[3]=0.02;
      ratioMax[3]=0.85;
      if(isPbPb) ratioMax[3]=1.25;
    }else{
      nHadr=5;
      pdgArray[0]=211;
      pdgArray[1]=321;
      pdgArray[2]=3122;
      pdgArray[3]=3312;
      pdgArray[4]=3334;
      partNames[0]="pi";
      partNames[1]="K";
      partNames[2]="L";
      partNames[3]="Xi";
      partNames[4]="Om";
      partTitles[0]="#pi^{+}";
      partTitles[1]="K^{+}";
      partTitles[2]="#Lambda";
      partTitles[3]="#Xi^{-}";
      partTitles[4]="#Omega^{-}";
      ratioMin[1]=0.;
      ratioMax[1]=0.3;
      ratioMin[2]=0.;
      ratioMax[2]=0.09;
      ratioMin[3]=0.;
      ratioMax[3]=0.01;
      ratioMin[4]=0.;
      ratioMax[4]=0.001;
    }
      
  
    TH1F* hYieldMid[maxPart][nPtBins];
    TH1F* hYieldFwd[maxPart][nPtBins];
    for(int ih=0; ih<nHadr; ih++){
      for(int ip=0; ip<nPtBins; ip++){
        hYieldMid[ih][ip] = new TH1F(Form("hYieldMidPrompt%spt%d-%d",partNames[ih].Data(),(int)ptLims[ip],(int)ptLims[ip+1])," ; n_{ch}^{mid}",nMultBins,multLimsM);
        hYieldFwd[ih][ip] = new TH1F(Form("hYieldFwdPrompt%spt%d-%d",partNames[ih].Data(),(int)ptLims[ip],(int)ptLims[ip+1])," ; n_{ch}^{fwd}",nMultBins,multLimsF);
      }
    }
  
    for(int j=0; j<treeHadr->GetEntries(); j++){
      treeHadr->GetEntry(j);
      if(TMath::Abs(y)<1. && prtag==0){
        int thePtBin=TMath::BinarySearch(nPtBins,ptLims,(double)pt);
        if(thePtBin>=0 && thePtBin<nPtBins){
      for(int ih=0; ih<nHadr; ih++){
        if(TMath::Abs(id)==pdgArray[ih]){
          hYieldMid[ih][thePtBin]->Fill(nChMid);
          hYieldFwd[ih][thePtBin]->Fill(nChFwd);
        }
      }
        }
      }
    }
    printf("Entries for Ds in pt bins:");
    for(int ip=0; ip<nPtBins; ip++) printf(" %d ", (int)hYieldFwd[2][ip]->GetEntries());
    printf("\n");
    
  
    
    TH1F* hRatioMid[maxPart][nPtBins];
    TH1F* hRatioFwd[maxPart][nPtBins];
    for(int ih=1; ih<nHadr; ih++){
      for(int ip=0; ip<nPtBins; ip++){   
        hRatioMid[ih][ip]=(TH1F*)hYieldMid[ih][ip]->Clone(Form("h%s%sRatioMidMultPt%d%d",partNames[ih].Data(),partNames[0].Data(),(int)ptLims[ip],(int)ptLims[ip+1]));
        FillParticleRatio(hYieldMid[ih][ip],hYieldMid[0][ip],hRatioMid[ih][ip]);
        hRatioMid[ih][ip]->GetYaxis()->SetTitle(Form("%s / %s",partTitles[ih].Data(),partTitles[0].Data()));
        hRatioMid[ih][ip]->SetMinimum(ratioMin[ih]);
        hRatioMid[ih][ip]->SetMaximum(ratioMax[ih]);
        hRatioFwd[ih][ip]=(TH1F*)hYieldFwd[ih][ip]->Clone(Form("h%s%sRatioFwdMultPt%d%d",partNames[ih].Data(),partNames[0].Data(),(int)ptLims[ip],(int)ptLims[ip+1]));
        FillParticleRatio(hYieldFwd[ih][ip],hYieldFwd[0][ip],hRatioFwd[ih][ip]);
        hRatioFwd[ih][ip]->GetYaxis()->SetTitle(Form("%s / %s",partTitles[ih].Data(),partTitles[0].Data()));
        hRatioFwd[ih][ip]->SetMinimum(ratioMin[ih]);
        hRatioFwd[ih][ip]->SetMaximum(ratioMax[ih]);
      }
    }
    if(optPart=="charm"){
      for(int ip=0; ip<nPtBins; ip++){   
        hRatioMid[0][ip]=(TH1F*)hYieldMid[2][ip]->Clone(Form("h%s%sRatioMidMultPt%d%d",partNames[2].Data(),partNames[1].Data(),(int)ptLims[ip],(int)ptLims[ip+1]));
        FillParticleRatio(hYieldMid[2][ip],hYieldMid[1][ip],hRatioMid[0][ip]);
        hRatioMid[0][ip]->GetYaxis()->SetTitle(Form("%s / %s",partTitles[2].Data(),partTitles[1].Data()));
        hRatioMid[0][ip]->SetMinimum(ratioMin[0]);
        hRatioMid[0][ip]->SetMaximum(ratioMax[0]);
        hRatioFwd[0][ip]=(TH1F*)hYieldFwd[2][ip]->Clone(Form("h%s%sRatioFwdMultPt%d%d",partNames[2].Data(),partNames[1].Data(),(int)ptLims[ip],(int)ptLims[ip+1]));
        FillParticleRatio(hYieldFwd[2][ip],hYieldFwd[1][ip],hRatioFwd[0][ip]);
        hRatioFwd[0][ip]->GetYaxis()->SetTitle(Form("%s / %s",partTitles[2].Data(),partTitles[1].Data()));
        hRatioFwd[0][ip]->SetMinimum(ratioMin[0]);
        hRatioFwd[0][ip]->SetMaximum(ratioMax[0]);
      }
    }
  
    int ptBinToDraw=2;
    double ptMin=ptLims[ptBinToDraw];
    double ptMax=ptLims[ptBinToDraw+1];
    
    TCanvas* cy=new TCanvas("cy","",1500,800);
    cy->Divide(3,2);
    for(int ih=0; ih<nHadr; ih++){
      cy->cd(ih+1);
      hYieldMid[ih][ptBinToDraw]->Draw();
    }
    
    TLatex* t2=new TLatex(0.15,0.72,Form("%.1f <p_{T} < %.1f GeV/c, |y|<1",ptMin,ptMax));
    t2->SetNDC();
  
    TCanvas* cratf=new TCanvas("cratf","Fwd",1500,800);
    cratf->Divide(2,2);
    for(int ih=1; ih<nHadr; ih++){
      cratf->cd(ih);
      hRatioFwd[ih][ptBinToDraw]->Draw();
      if(ih==1){
        t1->Draw();
        t2->Draw();
      }
    }
    if(optPart=="charm"){
      cratf->cd(nHadr);
      hRatioFwd[0][ptBinToDraw]->Draw();
    }
    cratf->SaveAs("CharmHadRatiosVsFwdMult.png");
    
    TCanvas* cratm=new TCanvas("cratm","Mid",1500,800);
    cratm->Divide(2,2);
    for(int ih=1; ih<nHadr; ih++){
      cratm->cd(ih);
      hRatioMid[ih][ptBinToDraw]->Draw();
      if(ih==1){
        t1->Draw();
        t2->Draw();
      }
    }
    if(optPart=="charm"){
      cratm->cd(nHadr);
      hRatioMid[0][ptBinToDraw]->Draw();
    }
    cratm->SaveAs("CharmHadRatiosVsMidMult.png");
  
    TFile* outroot = new TFile("CharmHadRatiosVsMult_merged.root","recreate");
    for(int ip=0; ip<nPtBins; ip++){
      for(int ih=1; ih<nHadr; ih++){
        hRatioMid[ih][ip]->Write();
        hRatioFwd[ih][ip]->Write();
      }
      if(optPart=="charm"){
        hRatioMid[0][ip]->Write();
        hRatioFwd[0][ip]->Write();
      }
    }
  }
  