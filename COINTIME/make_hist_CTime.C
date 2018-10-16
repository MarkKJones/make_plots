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
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
using namespace std;

void make_hist_CTime(TString basename="",Int_t nrun=1267){
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
   outputhist= "hist/"+basename+"_CTime_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  pTRIG6_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTimeRaw",&pTRIG6_ROC1_tdcTimeRaw);
 Double_t  pTRIG6_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTime",&pTRIG6_ROC1_tdcTime);
 Double_t  pTRIG6_ROC2_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTime",&pTRIG6_ROC2_tdcTime);
 Double_t  pTRIG4_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTime",&pTRIG4_ROC1_tdcTime);
 Double_t  pTRIG1_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTime",&pTRIG1_ROC1_tdcTime);
 Double_t  CoinTime_RAW_ROC1;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC1",&CoinTime_RAW_ROC1);
 Double_t  CoinTime_RAW_ROC2;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC2",&CoinTime_RAW_ROC2);
 Double_t  epCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC1",&epCoinTime_ROC1);
 Double_t  epCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC2",&epCoinTime_ROC2);
 Double_t  betanotrack;
   tsimc->SetBranchAddress("P.hod.betanotrack",&betanotrack);
 Double_t  betatrack;
   tsimc->SetBranchAddress("P.hod.beta",&betatrack);
 Double_t  delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
 Double_t  pntr;
   tsimc->SetBranchAddress("P.dc.ntrack",&pntr);
 Double_t  hntr;
   tsimc->SetBranchAddress("H.dc.ntrack",&hntr);
 Double_t  starttime;
   tsimc->SetBranchAddress("P.hod.starttime",&starttime);
   // Define histograms
    TH1F *hbetanotrack = new TH1F("hbetanotrack",Form("Run %d ; Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetanotrack);
    TH2F *hbetanotrack_delta = new TH2F("hbetanotrack_delta",Form("Run %d ; Beta notrack;Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hbetanotrack_delta);
    TH1F *hbetatrack = new TH1F("hbetatrack",Form("Run %d ; Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetatrack);
    TH1F *hstarttime = new TH1F("hstarttime",Form("Run %d ; Starttime;Counts",nrun), 280, -10.,60.0);
    HList.Add(hstarttime);
    TH2F *hstarttime_pdelta = new TH2F("hstarttime_pdelta",Form("Run %d ; Starttime; SHMS delta",nrun), 280, -10.,60.0,100,-10,20);
    HList.Add(hstarttime_pdelta);
    TH1F *hCoinTime_RAW_ROC1 = new TH1F("hcoinTime_RAW_ROC1",Form("Run %d ; CoinTime_RAW_ROC1  ;Counts",nrun),600, -50.,100.0);
    HList.Add(hCoinTime_RAW_ROC1);
    TH2F *hCoinTime_RAW_ROC1_pdelta = new TH2F("hcoinTime_RAW_ROC1_pdelta",Form("Run %d ; CoinTime_RAW_ROC1  ; SHMS delta",nrun),600, -50.,100.0,100,-10,20);
    HList.Add(hCoinTime_RAW_ROC1_pdelta);
    TH1F *hCoinTime_RAW_ROC2 = new TH1F("hcoinTime_RAW_ROC2",Form("Run %d ; CoinTime_RAW_ROC2  ;Counts",nrun),600, -50.,100.0);
    HList.Add(hCoinTime_RAW_ROC2);
    TH1F *hepCoinTime_ROC1 = new TH1F("hepcoinTime_ROC1",Form("Run %d ; epCoinTimeW_ROC1  ;Counts",nrun),800, -200.,200.0);
    HList.Add(hepCoinTime_ROC1);
    TH1F *hepCoinTime_ROC2 = new TH1F("hepcoinTime_ROC2",Form("Run %d ; epCoinTimeW_ROC2  ;Counts",nrun),800, -200.,200.0);
    HList.Add(hepCoinTime_ROC2);
    TH1F *hptrig1_ROC1 = new TH1F("hptrig1_ROC1",Form("Run %d ; ptrig1_ROC1 (SHMS trig in SHMS ROC) ;Counts",nrun),800, 0.,200.0);
    HList.Add(hptrig1_ROC1);
    TH1F *hptrig4_ROC1 = new TH1F("hptrig4_ROC1",Form("Run %d ; ptrig4_ROC1 (HMS trig in SHMS ROC)  ;Counts",nrun),800, -50.,300.0);
    HList.Add(hptrig4_ROC1);
    TH2F *hptrig4_ptrig1_ROC1 = new TH2F("hptrig4_ptrig1_ROC1",Form("Run %d ; ptrig4_ROC1 (HMS trig in SHMS ROC)  ; ptrig1_ROC1 (SHMS trig in SHMS ROC) ",nrun),200, -50.,300.0,200, 0.,200.0);
    HList.Add(hptrig4_ptrig1_ROC1);
    TH2F *hptrig4_ptrig1_ROC1_cut = new TH2F("hptrig4_ptrig1_ROC1_cut",Form("Run %d SHMS set time ; ptrig4_ROC1 (HMS trig in SHMS ROC)  ; ptrig1_ROC1 (SHMS trig in SHMS ROC) ",nrun),200, -50.,300.0,200, 0.,200.0);
    HList.Add(hptrig4_ptrig1_ROC1_cut);
    TH2F *hptrig6shms_ptrig6hms = new TH2F("hptrig6shms_ptrig6hms",Form("Run %d ; coin trig in SHMS ROC  ; coin trig in HMS ROC ",nrun),200, -50.,300.0,200, 0.,200.0);
    HList.Add(hptrig6shms_ptrig6hms);
 // loop over entries
    Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (pTRIG6_ROC1_tdcTimeRaw && pntr>0 && hntr >0) {
		hbetanotrack->Fill(betanotrack);
		hbetatrack->Fill(betatrack);
		hstarttime->Fill(starttime);
		hCoinTime_RAW_ROC1->Fill(CoinTime_RAW_ROC1);
		hCoinTime_RAW_ROC1_pdelta->Fill(CoinTime_RAW_ROC1,delta);
 		hptrig1_ROC1->Fill(pTRIG1_ROC1_tdcTime);
 		hptrig4_ROC1->Fill(pTRIG4_ROC1_tdcTime);
 		hptrig4_ptrig1_ROC1->Fill(pTRIG4_ROC1_tdcTime,pTRIG1_ROC1_tdcTime);
 		if (TMath::Abs(CoinTime_RAW_ROC1+2.5)<5) hptrig4_ptrig1_ROC1_cut->Fill(pTRIG4_ROC1_tdcTime,pTRIG1_ROC1_tdcTime);
 		hCoinTime_RAW_ROC2->Fill(CoinTime_RAW_ROC2);
		hepCoinTime_ROC1->Fill(epCoinTime_ROC1);
 		hepCoinTime_ROC2->Fill(epCoinTime_ROC2);
               hbetanotrack_delta->Fill(betanotrack,delta);
	       hptrig6shms_ptrig6hms->Fill(pTRIG6_ROC2_tdcTime,pTRIG6_ROC1_tdcTime);
	       hstarttime_pdelta->Fill(starttime,delta);
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
