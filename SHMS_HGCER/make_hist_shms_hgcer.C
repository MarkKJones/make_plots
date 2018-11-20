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

void make_hist_shms_hgcer(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_hgcer_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t npeSum;
   tsimc->SetBranchAddress("P.hgcer.npeSum",&npeSum) ;
 Double_t aeroSum;
   tsimc->SetBranchAddress("P.aero.npeSum",&aeroSum) ;
 Double_t xAtCer;
   tsimc->SetBranchAddress("P.hgcer.xAtCer",&xAtCer) ;
 Double_t yAtCer;
   tsimc->SetBranchAddress("P.hgcer.yAtCer",&yAtCer) ;
 Double_t ntrack;
   tsimc->SetBranchAddress("P.dc.ntrack",&ntrack) ;
 Double_t emiss;
   tsimc->SetBranchAddress("P.kin.secondary.emiss",&emiss) ;
 Double_t pmiss;
   tsimc->SetBranchAddress("P.kin.secondary.pmiss",&pmiss) ;
 Double_t ctime_eK;
   tsimc->SetBranchAddress("CTime.eKCoinTime_ROC1",&ctime_eK) ;
 Double_t etottracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etottracknorm) ;
   //
   TString temp=Form("Run %d ; Etottracknorm  ; Counts",nrun);
   TH1F *hetotnorm = new TH1F("hetotnorm",temp,100,0,2.0);
   HList.Add(hetotnorm);
   temp=Form("Run %d ; NpeSUm  ; Counts",nrun);
   TH1F *hcernpeSum = new TH1F("hcernpeSum",temp,160,0,40.);
   HList.Add(hcernpeSum);
  temp=Form("Run %d ; Ctime_eK  ; Counts",nrun);
   TH1F *hCtime_eK = new TH1F("hCtime_eK",temp,320,0,80.);
   HList.Add(hCtime_eK);
   temp=Form("Run %d ; HGCer NpeSUm  ; Counts",nrun);
   TH1F *hcernpeSum_aerocut = new TH1F("hcernpeSum_aerocut",temp,160,0,40.);
   HList.Add(hcernpeSum_aerocut);
   temp=Form("Run %d ; HGCer NpeSUm (mmiss cut)  ; Counts",nrun);
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

Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		if (1==1) {
		hetotnorm->Fill(etottracknorm);
		if (etottracknorm > 0.05) {
		hcernpeSum->Fill(npeSum);
		if (aeroSum>2.) hcernpeSum_aerocut->Fill(npeSum);
		haeronpeSum->Fill(aeroSum);
		hcernpeSum_cal->Fill(npeSum,etottracknorm);
		hcerXY->Fill(yAtCer,xAtCer);
		hcerXYnpe->Fill(yAtCer,xAtCer,npeSum);
                hCtime_eK->Fill(ctime_eK);
		Double_t mmiss = TMath::Sqrt(emiss*emiss-pmiss*pmiss);
                 if (TMath::Abs(ctime_eK-44.5)<2.) {
		  hmmiss->Fill(mmiss);
                  if (TMath::Abs(mmiss-0.9)<0.1) {
                   hcernpeSum_mmcut->Fill(npeSum);
                   hcernpeSum_cerX->Fill(npeSum,xAtCer);
		  }
		 }
		}
		}
}
//
}
