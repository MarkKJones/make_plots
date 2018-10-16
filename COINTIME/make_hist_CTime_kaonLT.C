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

void make_hist_CTime_kaonLT(TString basename="",Int_t nrun=1267){
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
   outputhist= "hist/"+basename+"_CTime_kaonLT_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Int_t  nhits_RF;
   tsimc->SetBranchAddress("Ndata.D.pRF_ROC2",&nhits_RF);
 Double_t  RFtime[10];
   tsimc->SetBranchAddress("D.pRF_ROC2",&RFtime);
 Double_t  CoinTime_RAW_ROC1;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC1",&CoinTime_RAW_ROC1);
 Double_t  CoinTime_RAW_ROC2;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC2",&CoinTime_RAW_ROC2);
 Double_t  epCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC1",&epCoinTime_ROC1);
 Double_t  epCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC2",&epCoinTime_ROC2);
 Double_t  eKCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.eKCoinTime_ROC1",&eKCoinTime_ROC1);
 Double_t  eKCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.eKCoinTime_ROC2",&eKCoinTime_ROC2);
 Double_t  epiCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC1",&epiCoinTime_ROC1);
 Double_t  epiCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC2",&epiCoinTime_ROC2);
 Double_t  pbetanotrack;
   tsimc->SetBranchAddress("P.hod.betanotrack",&pbetanotrack);
 Double_t  pbetatrack;
   tsimc->SetBranchAddress("P.hod.beta",&pbetatrack);
 Double_t  hbetanotrack;
   tsimc->SetBranchAddress("H.hod.betanotrack",&hbetanotrack);
 Double_t  hbetatrack;
   tsimc->SetBranchAddress("H.hod.beta",&hbetatrack);
 Double_t  pdelta;
   tsimc->SetBranchAddress("P.gtr.dp",&pdelta);
 Double_t  hdelta;
   tsimc->SetBranchAddress("H.gtr.dp",&hdelta);
 Double_t  pntr;
   tsimc->SetBranchAddress("P.dc.ntrack",&pntr);
 Double_t  hntr;
   tsimc->SetBranchAddress("H.dc.ntrack",&hntr);
 Double_t  hstarttime;
   tsimc->SetBranchAddress("H.hod.starttime",&hstarttime);
 Double_t  pstarttime;
   tsimc->SetBranchAddress("P.hod.starttime",&pstarttime);
   // Define histograms
    TH1F *hRFhits = new TH1F("hRFhits",Form("Run %d ; Number of RF hits;Counts",nrun), 10, 0.,10.0);
    HList.Add(hRFhits);
    TH1F *hRFtime[10];
    for (Int_t h=0;h<10;h++) {
      hRFtime[h] = new TH1F(Form("hRFtime_hit_%d",h),Form("Run %d ; Raw Time RF hit %d;Counts",nrun,h), 600., 0.,1200.0);
    HList.Add(hRFtime[h]);
    }
      hRF_tdiff_2hits_12 = new TH1F("hRF_tdiff_2hits_12",Form("Run %d RF mult=2; Time Diff hits 1-2;Counts",nrun), 800.,500,520.0);
      Double_t tdcchanperns = 0.09766;
    TH1F *hhbetanotrack = new TH1F("hhbetanotrack",Form("Run %d ; HMS Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hhbetanotrack);
    TH2F *hhbetanotrack_delta = new TH2F("hhbetanotrack_delta",Form("Run %d ; HMS Beta notrack;Delta",nrun), 200, .5,1.5,100,-10,10);
    HList.Add(hhbetanotrack_delta);
    TH1F *hhbetatrack = new TH1F("hhbetatrack",Form("Run %d ; HMS Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hhbetatrack);
    TH1F *hpbetanotrack = new TH1F("hpbetanotrack",Form("Run %d ; SHMS Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetanotrack);
    TH2F *hpbetanotrack_delta = new TH2F("hpbetanotrack_delta",Form("Run %d ; SHMS Beta notrack;Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hpbetanotrack_delta);
    TH1F *hpbetatrack = new TH1F("hpbetatrack",Form("Run %d ; SHMS Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetatrack);
    TH1F *hhstarttime = new TH1F("hhstarttime",Form("Run %d ; HMS Starttime;Counts",nrun), 280, -10.,80.0);
    HList.Add(hhstarttime);
    TH1F *hpstarttime = new TH1F("hpstarttime",Form("Run %d ; SHMS Starttime;Counts",nrun), 280, -10.,80.0);
    HList.Add(hpstarttime);
    TH2F *hpstarttime_pdelta = new TH2F("hpstarttime_pdelta",Form("Run %d ; SHMS Starttime; SHMS delta",nrun), 280, -10.,80.0,100,-10,20);
    HList.Add(hpstarttime_pdelta);
    TH2F *hhstarttime_hdelta = new TH2F("hhstarttime_pdelta",Form("Run %d ; HMS Starttime; HMS delta",nrun), 280, -10.,80.0,100,-10,10);
    HList.Add(hhstarttime_hdelta);
    TH1F *hCoinTime_RAW_ROC1 = new TH1F("hcoinTime_RAW_ROC1",Form("Run %d ; CoinTime_RAW_ROC1  ;Counts",nrun),600, -50.,100.0);
    HList.Add(hCoinTime_RAW_ROC1);
    TH2F *hCoinTime_RAW_ROC1_pdelta = new TH2F("hcoinTime_RAW_ROC1_pdelta",Form("Run %d ; CoinTime_RAW_ROC1  ; SHMS delta",nrun),600, -50.,100.0,100,-10,20);
    HList.Add(hCoinTime_RAW_ROC1_pdelta);
    TH1F *hCoinTime_RAW_ROC2 = new TH1F("hcoinTime_RAW_ROC2",Form("Run %d ; CoinTime_RAW_ROC2  ;Counts",nrun),600, -50.,100.0);
    HList.Add(hCoinTime_RAW_ROC2);
    TH1F *hepCoinTime_ROC1 = new TH1F("hepcoinTime_ROC1",Form("Run %d ; epCoinTime_ROC1  ;Counts",nrun),800,0.,100.0);
    HList.Add(hepCoinTime_ROC1);
    TH1F *hepCoinTime_ROC2 = new TH1F("hepcoinTime_ROC2",Form("Run %d ; epCoinTime_ROC2  ;Counts",nrun),800,0.,100.0);
    HList.Add(hepCoinTime_ROC2);
    TH2F *hepCoinTime_ROC2_pdelta = new TH2F("hepcoinTime_ROC2_pdelta",Form("Run %d ; epCoinTime_ROC2  ; SHMS delta",nrun),100,40.,50.0,60,-10,20);
    HList.Add(hepCoinTime_ROC2_pdelta);
    TH2F *hepCoinTime_ROC2_hdelta = new TH2F("hepcoinTime_ROC2_hdelta",Form("Run %d ; epCoinTime_ROC2  ; HMS delta",nrun),100,40.,50.0,40,-10,10);
    HList.Add(hepCoinTime_ROC2_pdelta);
    TH1F *hepiCoinTime_ROC2 = new TH1F("hepicoinTime_ROC2",Form("Run %d ; epiCoinTime_ROC2  ;Counts",nrun),800,0.,100.0);
    HList.Add(hepiCoinTime_ROC2);
    TH1F *heKCoinTime_ROC2 = new TH1F("heKcoinTime_ROC2",Form("Run %d ; eKCoinTime_ROC2  ;Counts",nrun),800,0.,100.0);
    HList.Add(heKCoinTime_ROC2);
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
		if (pntr>0 && hntr >0) {
		  hRFhits->Fill(nhits_RF);
		    for (Int_t h=0;h<nhits_RF;h++) {
		      hRFtime[h]->Fill(RFtime[h]*tdcchanperns);
		    } 
		    if (nhits_RF==2) hRF_tdiff_2hits_12->Fill((RFtime[1]-RFtime[0])*tdcchanperns); 
		hhbetanotrack->Fill(hbetanotrack);
		hhbetatrack->Fill(hbetatrack);
		hhstarttime->Fill(hstarttime);
		hpbetanotrack->Fill(pbetanotrack);
		hpbetatrack->Fill(pbetatrack);
		hpstarttime->Fill(pstarttime);
                hhbetanotrack_delta->Fill(hbetanotrack,hdelta);
                hpbetanotrack_delta->Fill(pbetanotrack,pdelta);
		hCoinTime_RAW_ROC1->Fill(CoinTime_RAW_ROC1);
		hCoinTime_RAW_ROC1_pdelta->Fill(CoinTime_RAW_ROC1,pdelta);
  		hCoinTime_RAW_ROC2->Fill(CoinTime_RAW_ROC2);
		hepCoinTime_ROC1->Fill(epCoinTime_ROC1);
 		hepCoinTime_ROC2->Fill(epCoinTime_ROC2);
 		hepCoinTime_ROC2_pdelta->Fill(epCoinTime_ROC2,pdelta);
 		hepCoinTime_ROC2_hdelta->Fill(epCoinTime_ROC2,hdelta);
 		hepiCoinTime_ROC2->Fill(epiCoinTime_ROC2);
 		heKCoinTime_ROC2->Fill(eKCoinTime_ROC2);
	       hpstarttime_pdelta->Fill(pstarttime,pdelta);
	       hhstarttime_hdelta->Fill(hstarttime,hdelta);
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
