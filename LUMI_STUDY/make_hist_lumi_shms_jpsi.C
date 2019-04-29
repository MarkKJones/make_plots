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

void make_hist_lumi_shms_jpsi(TString basename="",Int_t nrun=2043,Double_t mean_current=2.0){
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
 Double_t  Scal_time;
   tscal->SetBranchAddress("P.1MHz.scalerTime",&Scal_time);
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
          ave_current+=Scal_BCM4B_current;
	  h_cur_entry->Fill(float(i),Scal_BCM4B_current);
	  Bool_t   special_skip=kFALSE;
	  if (i ==178 && nrun==7572) special_skip=kTRUE;
	  if (TMath::Abs(Scal_BCM4B_current-mean_current) < threshold_cut && !special_skip) {
             event_flag[nscal_reads] = 1;
             ave_current_cut+=Scal_BCM4B_current;
 	     tot_scal_cut_time+=(Scal_time-prev_time);
 	     tot_scal_cut_TRIG2+=(Scal_TRIG2-prev_TRIG2);
 	     tot_scal_cut_TRIG3+=(Scal_TRIG3-prev_TRIG3);
 	     tot_scal_cut_TRIG1+=(Scal_TRIG1-prev_TRIG1);
 	     tot_scal_cut_TRIG4+=(Scal_TRIG4-prev_TRIG4);
 	     tot_scal_cut_TRIG5+=(Scal_TRIG5-prev_TRIG5);
	     for (Int_t s=0;s<4;s++) tot_scal_cut_Splane[s]+=(Scal_Splane[s]-prev_Splane[s]);
	     //	     cout << i << " " << tot_scal_cut_Splane[0] << " " << Scal_Splane[0] << " " << prev_Splane[0] << endl;

             charge_sum_cut+=(Scal_BCM4B_charge-prev_charge);
             nscal_reads_cut++;
	  }
	  prev_charge = Scal_BCM4B_charge;
	  prev_time = Scal_time;
	  prev_TRIG2 = Scal_TRIG2;
	  prev_TRIG3 = Scal_TRIG3;
	  prev_TRIG1 = Scal_TRIG1;
	  prev_TRIG4 = Scal_TRIG4;
	  prev_TRIG5 = Scal_TRIG5;
	     for (Int_t s=0;s<4;s++) prev_Splane[s]=Scal_Splane[s];
	  // cout <<  nscal_reads <<  " " << Scal_BCM4B_current << " " << event_flag[nscal_reads] << " " << Scal_TRIG1 << endl;
          nscal_reads++;
          charge_sum=Scal_BCM4B_charge;
	  tot_scal_TRIG2=Scal_TRIG2;
	  tot_scal_TRIG3=Scal_TRIG3;
	  tot_scal_TRIG1=Scal_TRIG1;
	  tot_scal_TRIG4=Scal_TRIG4;
	  tot_scal_TRIG5=Scal_TRIG5;
          tot_scal_time=Scal_time;
	  for (Int_t s=0;s<4;s++) tot_scal_Splane[s]=Scal_Splane[s];
	}
   //
	Double_t gevtyp;
	tsimc->SetBranchAddress("g.evtyp",&gevtyp);
	Double_t gevnum;
	tsimc->SetBranchAddress("g.evnum",&gevnum);
	Double_t ntrack;
	tsimc->SetBranchAddress("P.dc.ntrack",&ntrack);
	Double_t goodscinhit;
	tsimc->SetBranchAddress("P.hod.goodscinhit",&goodscinhit);
	Double_t starttime;
	tsimc->SetBranchAddress("P.hod.starttime",&starttime);
 	Double_t betanotrack;
	tsimc->SetBranchAddress("P.hod.betanotrack",&betanotrack);
  Double_t ngcer_npeSum;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&ngcer_npeSum);
   Double_t delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
   Double_t etotnorm;
   tsimc->SetBranchAddress("P.cal.etotnorm",&etotnorm);
   Double_t etottracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etottracknorm);
	//
	TH1F* h_starttime = new TH1F("h_starttime"," ; Startime (elec good scinhit)",60,0.,120.);
	TH1F* h_starttime_ntrack = new TH1F("h_starttime_ntrack"," ; Startime (elec good scinhit ntrack)",60,0.,120.);
	TH1F* h_etotnorm_all = new TH1F("h_etotnorm_all"," All events; etot  norm",100,0.,3.);
	TH1F* h_etottracknorm_trig1 = new TH1F("h_etottracknorm_trig1"," Trig 1 ; etot track norm",100,0.,3.);
	TH1F* h_etottracknorm_trig4 = new TH1F("h_etottracknorm_trig4"," Trig 4 ; etot track norm",100,0.,3.);
	TH1F* h_etottracknorm_all = new TH1F("h_etottracknorm_all"," All events; etot track norm",100,0.,3.);
	TH1F* h_ngcernpeSum_all = new TH1F("h_ngcernpeSum_all"," All events; NG cer npe sum",120,0.,30.);
	TH1F* h_ngcernpeSum_all_etotcut = new TH1F("h_ngcernpeSum_all_etotcut"," All events; NG cer npe sum",120,0.,30.);
	TH1F* h_etottracknorm_curcut = new TH1F("h_etottracknorm_curcut"," Current cut; etot track norm",100,0.,3.);
	//
   Int_t nscal_reads_2=0;
   prev_read=-1;
   Int_t ps2;
   Double_t trackeff;
   if (nrun==7570) ps2=9;
   if (nrun==7570) trackeff=0.9915;
   if (nrun==7571) ps2=3;
   if (nrun==7571) trackeff=0.9857;
   if (nrun==7572) ps2=2;
   if (nrun==7572) trackeff=0.9807;
   if (nrun==7574) ps2=9;
   if (nrun==7574) trackeff=0.9691;
   if (nrun==7499) ps2=33;
   if (nrun==7499) trackeff=0.9616;
   if (nrun==7500) ps2=33;
   if (nrun==7500) trackeff=0.9756;
   if (nrun==7501) ps2=33;
   if (nrun==7501) trackeff=0.9756;
   if (nrun==7419) ps2=129;
   if (nrun==7419) trackeff=1.0;
   if (nrun==7420) ps2=129;
   if (nrun==7420) trackeff=1.0;
	for (int i = 0; i < data_entries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		h_ngcernpeSum_all->Fill(ngcer_npeSum);
		if (etotnorm > 0.6) h_ngcernpeSum_all_etotcut->Fill(ngcer_npeSum);
		if (etotnorm>0.6&&ngcer_npeSum!=0&& delta>-10&& delta<22) h_etottracknorm_all->Fill(etottracknorm);
 		if (etotnorm > 0.6&&ngcer_npeSum!=0) h_etotnorm_all->Fill(etotnorm);
		  if (event_flag[nscal_reads_2]==1&&etotnorm>0.6&&ngcer_npeSum!=0&&goodscinhit==1) {
		  h_starttime->Fill(starttime);
		  if (ntrack >0) h_starttime_ntrack->Fill(starttime);
		  }
                if (event_flag[nscal_reads_2]==1&&etotnorm>0.6&&ngcer_npeSum!=0&& delta>-10&& delta<22) {
		  h_etottracknorm_curcut->Fill(etottracknorm);
		  if (gevtyp==1) {
		    h_etottracknorm_trig1->Fill(etottracknorm);
		      }
		  if (gevtyp==4) {
		    h_etottracknorm_trig4->Fill(etottracknorm);
		      }
 		}
	       	//  cout << nscal_reads_2 << " " << gevtyp << " " << gevnum << " " << scal_event_number[nscal_reads_2] << endl;
		
		if (gevnum>scal_event_number[nscal_reads_2]) {
                   nscal_reads_2++;
		}
	}
	//
	Double_t calc_treff= float(h_starttime_ntrack->Integral())/float(h_starttime->Integral()) ;
	Int_t good_ev1 = h_etottracknorm_trig1->Integral();
	Double_t good_ev1_err = TMath::Sqrt(good_ev1);
	Int_t good_ev4 = h_etottracknorm_trig4->Integral();
	Double_t good_ev4_err = TMath::Sqrt(good_ev4);
	Double_t good_ev = (good_ev1*ps2+good_ev4)/calc_treff;
	Double_t good_ev_err = TMath::Sqrt( (good_ev1_err*ps2)*(good_ev1_err*ps2) + (good_ev4_err)*(good_ev4_err))/calc_treff;
	//
        Double_t err_trig2 = 1./TMath::Sqrt(tot_scal_cut_TRIG2);
        Double_t err_trig3 = 1./TMath::Sqrt(tot_scal_cut_TRIG3);
        Double_t err_trig1 = 1./TMath::Sqrt(tot_scal_cut_TRIG1);
        Double_t err_trig4 = 1./TMath::Sqrt(tot_scal_cut_TRIG4);
        Double_t err_trig5 = 1./TMath::Sqrt(tot_scal_cut_TRIG5);
			   cout << nrun << " " << charge_sum_cut << " " << mean_current << " "<< charge_sum_cut/tot_scal_cut_time << " " << tot_scal_cut_time/tot_scal_time << " " << tot_scal_cut_TRIG2/tot_scal_cut_time << " " << tot_scal_cut_TRIG4/tot_scal_cut_time << " " << good_ev/tot_scal_cut_time  << " " << tot_scal_cut_TRIG2/charge_sum_cut << " " << err_trig2*tot_scal_cut_TRIG2/charge_sum_cut << " " << tot_scal_cut_TRIG4/charge_sum_cut << " " << err_trig4*tot_scal_cut_TRIG4/charge_sum_cut<< " " << good_ev/charge_sum_cut<< " " << good_ev_err/charge_sum_cut << " " << calc_treff<< endl;
cout << nrun << " " << mean_current << " "<< charge_sum_cut/tot_scal_cut_time << " "<< tot_scal_cut_Splane[0]/tot_scal_cut_time<< " "<< tot_scal_cut_Splane[1]/tot_scal_cut_time<< " "<< tot_scal_cut_Splane[2]/tot_scal_cut_time<< " "<< tot_scal_cut_Splane[3]/tot_scal_cut_time<< " " << calc_treff<< endl;

// loop through data
Long64_t nentries = tsimc->GetEntries();
//
}
