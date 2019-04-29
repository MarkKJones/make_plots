#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void make_hist_shms_ngcer_cosmic(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_shms_ngcer_cosmic_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
   Double_t adcamp[4];
   tsimc->SetBranchAddress("P.ngcer.goodAdcPulseAmp",&adcamp) ;
   Double_t adcint[4];
   tsimc->SetBranchAddress("P.ngcer.goodAdcPulseInt",&adcint) ;
   Double_t adctime[4];
   tsimc->SetBranchAddress("P.ngcer.goodAdcPulseTime",&adctime) ;
   Double_t adctdcdifftime[4];
   tsimc->SetBranchAddress("P.ngcer.goodAdcTdcDiffTime",&adctdcdifftime) ;

 Double_t npeSum;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&npeSum) ;
 Double_t aeroSum;
   tsimc->SetBranchAddress("P.aero.npeSum",&aeroSum) ;
 Double_t xAtCer;
   tsimc->SetBranchAddress("P.ngcer.xAtCer",&xAtCer) ;
 Double_t yAtCer;
   tsimc->SetBranchAddress("P.ngcer.yAtCer",&yAtCer) ;
 Double_t ntrack;
   tsimc->SetBranchAddress("P.dc.ntrack",&ntrack) ;
 Double_t etottracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etottracknorm) ;
   //
   Double_t ycut[4]={9,-10,10,-8};
   Double_t xcut[4]={40,20,-20,-20};
   //
   TH1F *h_adcamp[4];
   TH1F *h_adctimediff[4];
   TH1F *h_adcint[4];
   TH1F *h_adcamp_cut[4];
   TH1F *h_adcint_cut[4];
   TH1F *h_adcamp_cutlow[4];
   TH1F *h_adcint_cutlow[4];
   TH2F *h_cerXY_adcamp[4];
   TH2F *h_cerXY_adcamp_now[4];
   TH2F *h_cerXY_adcint[4];
   TH2F *h_cerXY_adcint_now[4];
   //   TString temp;
   for (Int_t n=0; n<4;n++) {
   temp=Form("Run %d ; PMT %d  Pulse Amp ; Counts",nrun,n);
   h_adcamp[n] = new TH1F(Form("h_adcamp_%d",n),temp,50,0.,50.);
   HList.Add(h_adcamp[n]);
   temp=Form("Run %d ; PMT %d  Pulse Time diff ; Counts",nrun,n);
   h_adctimediff[n] = new TH1F(Form("h_adctimediff_%d",n),temp,400,-200.,200.);
   HList.Add(h_adctimediff[n]);
   temp=Form("Run %d ; PMT %d  Pulse Int ; Counts",nrun,n);
   h_adcint[n] = new TH1F(Form("h_adcint_%d",n),temp,50,0.,50.);
   HList.Add(h_adcint[n]);
   temp=Form("Run %d Cut ; PMT %d  Pulse Amp ; Counts",nrun,n);
   h_adcamp_cut[n] = new TH1F(Form("h_adcamp_cut_%d",n),temp,50,0.,50.);
   HList.Add(h_adcamp_cut[n]);
   temp=Form("Run %d Cut; PMT %d  Pulse Int ; Counts",nrun,n);
   h_adcint_cut[n] = new TH1F(Form("h_adcint_cut_%d",n),temp,50,0.,50.);
   HList.Add(h_adcint_cut[n]);
   temp=Form("Run %d Cutlow ; PMT %d  Pulse Amp ; Counts",nrun,n);
   h_adcamp_cutlow[n] = new TH1F(Form("h_adcamp_cutlow_%d",n),temp,50,0.,50.);
   HList.Add(h_adcamp_cutlow[n]);
   temp=Form("Run %d Cutlow; PMT %d  Pulse Int ; Counts",nrun,n);
   h_adcint_cutlow[n] = new TH1F(Form("h_adcint_cutlow_%d",n),temp,50,0.,50.);
   HList.Add(h_adcint_cutlow[n]);
    temp=Form("Run %d Pulse amp ; PMT %d  Y ; X ",nrun,n);
   h_cerXY_adcamp[n] = new TH2F(Form("hcerXY_adcamp_%d",n),temp,160,-40,40.,100,-60,60.0);
   h_cerXY_adcamp_now[n] = new TH2F(Form("hcerXY_adcamp_now_%d",n),temp,160,-40,40.,100,-60,60.0);
   HList.Add(h_cerXY_adcamp[n]);
   HList.Add(h_cerXY_adcamp_now[n]);
    temp=Form("Run %d Pulse int ; PMT %d  Y ; X ",nrun,n);
   h_cerXY_adcint[n] = new TH2F(Form("hcerXY_adcint_%d",n),temp,160,-40,40.,100,-60,60.0);
   h_cerXY_adcint_now[n] = new TH2F(Form("hcerXY_adcint_now_%d",n),temp,160,-40,40.,100,-60,60.0);
   HList.Add(h_cerXY_adcint[n]);
   HList.Add(h_cerXY_adcint_now[n]);
   }
   TString temp=Form("Run %d ; Etottracknorm  ; Counts",nrun);
   TH1F *hetotnorm = new TH1F("hetotnorm",temp,100,0,2.0);
   HList.Add(hetotnorm);
   temp=Form("Run %d ; NpeSUm  ; Counts",nrun);
   TH1F *hcernpeSum = new TH1F("hcernpeSum",temp,160,0,40.);
   HList.Add(hcernpeSum);
  temp=Form("Run %d ; Ctime_eK  ; Counts",nrun);
   TH1F *hCtime_eK = new TH1F("hCtime_eK",temp,320,0,80.);
   HList.Add(hCtime_eK);
   temp=Form("Run %d ; NGCer NpeSUm  ; Counts",nrun);
   TH1F *hcernpeSum_aerocut = new TH1F("hcernpeSum_aerocut",temp,160,0,40.);
   HList.Add(hcernpeSum_aerocut);
   temp=Form("Run %d ; NGCer NpeSUm (mmiss cut)  ; Counts",nrun);
   TH1F *hcernpeSum_mmcut = new TH1F("hcernpeSum_mmcut",temp,160,0,40.);
   HList.Add(hcernpeSum_mmcut);
   temp=Form("Run %d ; Aero NpeSUm  ; Counts",nrun);
   TH1F *haeronpeSum = new TH1F("haeronpeSum",temp,160,0,40.);
   HList.Add(haeronpeSum);
    temp=Form("Run %d ; Mmiss  ; Counts",nrun);
   TH1F *hmmiss = new TH1F("hmmiss",temp,150,0.5,2.0);
   HList.Add(hmmiss);
   temp=Form("Run %d ; NpeSUm  ; Etottracknorm",nrun);
   TH2F *hcernpeSum_cal = new TH2F("hcernpeSum_cal",temp,160,0,40.,100,0,2.0);
   HList.Add(hcernpeSum_cal);
   temp=Form("Run %d ; NpeSUm  ; X at Cer",nrun);
   TH2F *hcernpeSum_cerX = new TH2F("hcernpeSum_cerX",temp,160,0,40.,12,-60,60);
   HList.Add(hcernpeSum_cerX);
   temp=Form("Run %d ; Y at CEr ; X at Cer",nrun);
   TH2F *hcerXY = new TH2F("hcerXY",temp,160,-40,40.,100,-60,60.0);
   HList.Add(hcerXY);
   TH2F *hcerXYnpe = new TH2F("hcerXYnpe",temp,160,-40,40.,100,-60,60.0);
   HList.Add(hcerXYnpe);

//
   Int_t hit_index=0;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		hetotnorm->Fill(etottracknorm);
		if (1==1) {
		  hit_index=-1;
                 for (Int_t n=0; n<4;n++) {
		   h_adcamp[n]->Fill(adcamp[n]);
		   h_adctimediff[n]->Fill(adctdcdifftime[n]);
		   h_adcint[n]->Fill(adcint[n]);
		   if (adcamp[n]>5 && xAtCer>40 && xAtCer<80) {
		   h_adcamp_cut[n]->Fill(adcamp[n]);
		   h_adcint_cut[n]->Fill(adcint[n]);
		   hit_index=n;
		   h_adcamp_cutlow[n]->Fill(adcamp[n]);
		   h_adcint_cutlow[n]->Fill(adcint[n]);
		       }
		   h_cerXY_adcamp_now[n]->Fill(yAtCer,xAtCer);
		   h_cerXY_adcamp[n]->Fill(yAtCer,xAtCer,adcamp[n]);
		   h_cerXY_adcint_now[n]->Fill(yAtCer,xAtCer);
		   h_cerXY_adcint[n]->Fill(yAtCer,xAtCer,adcint[n]);
		 }
                 }
		 /*		 if (hit_index >-1) {
		   if (hit_index==0) h_adcamp_cutlow[3]->Fill(adcamp[3]);
		   if (hit_index==1) h_adcamp_cutlow[2]->Fill(adcamp[2]);
		   if (hit_index==2) h_adcamp_cutlow[1]->Fill(adcamp[1]);
		   if (hit_index==3) h_adcamp_cutlow[0]->Fill(adcamp[0]);
		 }
		 */
		hcernpeSum->Fill(npeSum);
	}
 TFile hsimc(outputhist,"recreate");
        HList.Write();
//
}
