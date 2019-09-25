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

void make_hist_lumi(TString basename="",Int_t nrun=2043,Double_t mean_current=2.0){
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
   outputhist= "hist/"+basename+"_lumi_shms_jpsi_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
TTree *tscal = (TTree*) fsimc->Get("TSP");
//
 Double_t  Scal_evNumber;
   tscal->SetBranchAddress("evNumber",&Scal_evNumber);
 Double_t  Scal_BCM4B_charge;
   tscal->SetBranchAddress("P.BCM4B.scalerCharge",&Scal_BCM4B_charge);
 Double_t  Scal_BCM4B_current;
   tscal->SetBranchAddress("P.BCM4B.scalerCurrent",&Scal_BCM4B_current);
 Double_t  Scal_BCM1_charge;
   tscal->SetBranchAddress("P.BCM1.scalerCharge",&Scal_BCM1_charge);
 Double_t  Scal_BCM1_current;
   tscal->SetBranchAddress("P.BCM1.scalerCurrent",&Scal_BCM1_current);
 Double_t  Scal_time;
   tscal->SetBranchAddress("P.1MHz.scalerTime",&Scal_time);
 Double_t  Scal_EDTM;
   tscal->SetBranchAddress("P.EDTM.scaler",&Scal_EDTM);
Double_t  Scal_hEL_CLEAN;
   tscal->SetBranchAddress("P.hEL_CLEAN.scaler",&Scal_hEL_CLEAN);
Double_t  Scal_hEL_REAL;
   tscal->SetBranchAddress("P.hEL_REAL.scaler",&Scal_hEL_REAL);
Double_t  Scal_TRIG2;
   tscal->SetBranchAddress("P.pTRIG2.scaler",&Scal_TRIG2);
 Double_t  Scal_TRIG3;
   tscal->SetBranchAddress("P.pTRIG3.scaler",&Scal_TRIG3);
 Double_t  Scal_TRIG1;
   tscal->SetBranchAddress("P.pTRIG1.scaler",&Scal_TRIG1);
 Double_t  Scal_TRIG4;
   tscal->SetBranchAddress("P.pTRIG4.scaler",&Scal_TRIG4);
 Double_t  Scal_TRIG5;
   tscal->SetBranchAddress("P.pTRIG5.scaler",&Scal_TRIG5);
 Double_t  Scal_Splane[4];
   tscal->SetBranchAddress("P.S1X.scaler",&Scal_Splane[0]);
   tscal->SetBranchAddress("P.S1Y.scaler",&Scal_Splane[1]);
   tscal->SetBranchAddress("P.S2X.scaler",&Scal_Splane[2]);
   tscal->SetBranchAddress("P.S2Y.scaler",&Scal_Splane[3]);
//loop through scalers
     Int_t nscal_reads=0;
     Int_t nscal_reads_cut=0;
     Double_t prev_read=-1;
     Double_t ave_current=0;
     Double_t ave_current_cut=0;
     Double_t charge_sum=0;
     Double_t charge_sum_cut=0;
     Double_t prev_charge=0;
     Double_t event_flag[10000];
     Double_t scal_event_number[10000];
     Double_t tot_scal_EDTM=0;
     Double_t tot_scal_cut_EDTM=0;
     Double_t prev_EDTM=0;
     Double_t tot_scal_hEL_CLEAN=0;
     Double_t tot_scal_cut_hEL_CLEAN=0;
     Double_t prev_hEL_CLEAN=0;
     Double_t tot_scal_hEL_REAL=0;
     Double_t tot_scal_cut_hEL_REAL=0;
     Double_t prev_hEL_REAL=0;
     Double_t tot_scal_TRIG2=0;
     Double_t tot_scal_TRIG3=0;
     Double_t prev_TRIG2=0;
     Double_t prev_TRIG3=0;
     Double_t tot_scal_cut_TRIG2=0;
     Double_t tot_scal_cut_TRIG3=0;
     Double_t tot_scal_TRIG1=0;
     Double_t tot_scal_TRIG4=0;
     Double_t prev_TRIG1=0;
     Double_t prev_TRIG4=0;
     Double_t tot_scal_cut_TRIG1=0;
     Double_t tot_scal_cut_TRIG4=0;
     Double_t tot_scal_cut_time=0;
     Double_t tot_scal_TRIG5=0;
     Double_t prev_TRIG5=0;
    Double_t tot_scal_cut_TRIG5=0;
    Double_t threshold_cut=3.;
			   //
    Double_t tot_scal_Splane[4]={0,0,0,0};
     Double_t prev_Splane[4]={0,0,0,0};
    Double_t tot_scal_cut_Splane[4]={0,0,0,0};
			   //
      Double_t tot_scal_time=0;
     Double_t prev_time=0;
     //
Long64_t scal_entries = tscal->GetEntries();
 cout << " scal ent = " << scal_entries << endl;
 Double_t nlast=float(scal_entries);
 TH1F *h_cur_entry = new TH1F("h_cur_entry","; ENtry;current",nlast,0,nlast);
Long64_t data_entries = tsimc->GetEntries();
	for (int i = 0; i < scal_entries; i++) {
      		tscal->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
          event_flag[nscal_reads] = 0;
             scal_event_number[nscal_reads] = Scal_evNumber;
          ave_current+=Scal_BCM1_current;
	  h_cur_entry->Fill(float(i),Scal_BCM1_current);
	  Bool_t   special_skip=kFALSE;
	  if (i ==178 && nrun==7572) special_skip=kTRUE;
	  if (TMath::Abs(Scal_BCM1_current-mean_current) < threshold_cut && !special_skip) {
             event_flag[nscal_reads] = 1;
             ave_current_cut+=Scal_BCM4B_current;
 	     tot_scal_cut_time+=(Scal_time-prev_time);
 	     tot_scal_cut_EDTM+=(Scal_EDTM-prev_EDTM);
 	     tot_scal_cut_hEL_CLEAN+=(Scal_hEL_CLEAN-prev_hEL_CLEAN);
	     tot_scal_cut_hEL_REAL+=(Scal_hEL_REAL-prev_hEL_REAL);
 	     tot_scal_cut_TRIG2+=(Scal_TRIG2-prev_TRIG2);
 	     tot_scal_cut_TRIG3+=(Scal_TRIG3-prev_TRIG3);
 	     tot_scal_cut_TRIG1+=(Scal_TRIG1-prev_TRIG1);
 	     tot_scal_cut_TRIG4+=(Scal_TRIG4-prev_TRIG4);
 	     tot_scal_cut_TRIG5+=(Scal_TRIG5-prev_TRIG5);
	     for (Int_t s=0;s<4;s++) tot_scal_cut_Splane[s]+=(Scal_Splane[s]-prev_Splane[s]);
	     //	     cout << i << " " << tot_scal_cut_Splane[0] << " " << Scal_Splane[0] << " " << prev_Splane[0] << endl;

             charge_sum_cut+=(Scal_BCM1_charge-prev_charge);
             nscal_reads_cut++;
	  }
	  prev_charge = Scal_BCM1_charge;
	  prev_time = Scal_time;
	  prev_EDTM = Scal_EDTM;
	  prev_hEL_CLEAN = Scal_hEL_CLEAN;
	  prev_hEL_REAL = Scal_hEL_REAL;
	  prev_TRIG2 = Scal_TRIG2;
	  prev_TRIG3 = Scal_TRIG3;
	  prev_TRIG1 = Scal_TRIG1;
	  prev_TRIG4 = Scal_TRIG4;
	  prev_TRIG5 = Scal_TRIG5;
	     for (Int_t s=0;s<4;s++) prev_Splane[s]=Scal_Splane[s];
	  // cout <<  nscal_reads <<  " " << Scal_BCM4B_current << " " << event_flag[nscal_reads] << " " << Scal_TRIG1 << endl;
          nscal_reads++;
          charge_sum=Scal_BCM1_charge;
	  tot_scal_EDTM=Scal_EDTM;
	  tot_scal_hEL_CLEAN=Scal_hEL_CLEAN;
	  tot_scal_hEL_REAL=Scal_hEL_REAL;
	  tot_scal_TRIG2=Scal_TRIG2;
	  tot_scal_TRIG3=Scal_TRIG3;
	  tot_scal_TRIG1=Scal_TRIG1;
	  tot_scal_TRIG4=Scal_TRIG4;
	  tot_scal_TRIG5=Scal_TRIG5;
          tot_scal_time=Scal_time;
	  for (Int_t s=0;s<4;s++) tot_scal_Splane[s]=Scal_Splane[s];
	}
	//
	Double_t EDTM_timeraw;
	tsimc->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw",&EDTM_timeraw);
	Double_t gevnum;
	tsimc->SetBranchAddress("g.evnum",&gevnum);
	Double_t pgoodscinhit;
	tsimc->SetBranchAddress("P.hod.goodscinhit",&pgoodscinhit);
	Double_t hgoodscinhit;
	tsimc->SetBranchAddress("H.hod.goodscinhit",&hgoodscinhit);
	Double_t gevtyp;
	tsimc->SetBranchAddress("g.evtyp",&gevtyp);
	Double_t pntrack;
	tsimc->SetBranchAddress("P.dc.ntrack",&pntrack);
	Double_t hntrack;
	tsimc->SetBranchAddress("H.dc.ntrack",&hntrack);
 	Double_t hbetanotrack;
	tsimc->SetBranchAddress("H.hod.betanotrack",&hbetanotrack);
 	Double_t petotnorm;
	tsimc->SetBranchAddress("P.cal.etotnorm",&petotnorm);
	//
	TH1F* h_EDTM_CUT = new TH1F("h_EDTM"," ; EDTM (beam cut)",1000,0.,10000.);
	TH1F* h_hScinHit = new TH1F("h_hScinhit"," ; HMS Scinhit",10,0.,10.);
	TH1F* h_hScinHit_track = new TH1F("h_hScinhit_track"," ; HMS Scinhit_track",10,0.,10.);
	TH1F* h_pScinHit = new TH1F("h_pScinhit"," ; SHMS Scinhit",10,0.,10.);
	TH1F* h_pScinHit_track = new TH1F("h_pScinhit_track"," ; SHMS Scinhit_track",10,0.,10.);
	TH1F* h_ev1 = new TH1F("h_ev1"," ; hEv1",10,0.,10.);
	TH1F* h_ev4 = new TH1F("h_ev4"," ; hEv4",10,0.,10.);
	//
  Int_t nscal_reads_2=0;
   prev_read=-1;
	for (int i = 0; i < data_entries; i++) {
      		tsimc->GetEntry(i);
                if (i%100000==0) cout << " Entry = " << i << endl;
	        if (event_flag[nscal_reads_2]==1&&gevtyp==1) h_ev1->Fill(gevtyp);
	        if (event_flag[nscal_reads_2]==1&&gevtyp==4) h_ev4->Fill(gevtyp);
                if (event_flag[nscal_reads_2]==1&&gevtyp==4 && EDTM_timeraw>0) h_EDTM_CUT->Fill(EDTM_timeraw);
		if (event_flag[nscal_reads_2]==1&&gevtyp==4) {
		  if ( hbetanotrack < 0.9) h_hScinHit->Fill(hgoodscinhit);
		  if ( hbetanotrack < 0.9 && hntrack>0) h_hScinHit_track->Fill(hgoodscinhit);
		  if ( petotnorm > .7) h_pScinHit->Fill(pgoodscinhit);
		  if ( petotnorm > .7 && pntrack>0) h_pScinHit_track->Fill(pgoodscinhit);
		}
		//
		if (gevnum>scal_event_number[nscal_reads_2]) {
                   nscal_reads_2++;
		}
	}
	//
	Int_t nEDTM = h_EDTM_CUT->Integral();
	Double_t edtm_lt=  nEDTM/tot_scal_cut_EDTM;
	Int_t good_ev4 = h_ev4->Integral();
	Int_t good_ev1 = h_ev1->Integral();
	cout << " Select events with current  = " << mean_current << " +/- " << threshold_cut << " charge  " << charge_sum_cut<< endl;
	cout << " good_ev4 = " << good_ev4 << " " << " good_ev1 = " << good_ev1 << " " << endl;
	cout << "EV4  Computer Livetime = " << (good_ev4-nEDTM)/(tot_scal_cut_TRIG5-tot_scal_cut_EDTM) << endl;
         cout << " EDTM LT = "  << edtm_lt << endl;
	//
}
