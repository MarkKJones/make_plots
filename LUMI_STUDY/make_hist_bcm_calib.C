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

void make_hist_bcm_calib(TString basename="",Int_t nrun=20){
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
   inputroot="ROOTfiles_vcs/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_bcm_calib_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 

TTree *tscal = (TTree*) fsimc->Get("TSP");
//
 Double_t  Scal_evNumber;
   tscal->SetBranchAddress("evNumber",&Scal_evNumber);
 Double_t  Scal_BCM4B_rate;
   tscal->SetBranchAddress("P.BCM4B.scalerRate",&Scal_BCM4B_rate);
 Double_t  Scal_BCM4B_charge;
   tscal->SetBranchAddress("P.BCM4B.scalerCharge",&Scal_BCM4B_charge);
 Double_t  Scal_BCM4B_current;
   tscal->SetBranchAddress("P.BCM4B.scalerCurrent",&Scal_BCM4B_current);
 Double_t  Scal_BCM4C_rate;
   tscal->SetBranchAddress("P.BCM4C.scalerRate",&Scal_BCM4C_rate);
 Double_t  Scal_BCM4C_charge;
   tscal->SetBranchAddress("P.BCM4C.scalerCharge",&Scal_BCM4C_charge);
 Double_t  Scal_BCM4C_current;
   tscal->SetBranchAddress("P.BCM4C.scalerCurrent",&Scal_BCM4C_current);
 Double_t  Scal_BCM4A_rate;
   tscal->SetBranchAddress("P.BCM4A.scalerRate",&Scal_BCM4A_rate);
 Double_t  Scal_BCM4A_charge;
   tscal->SetBranchAddress("P.BCM4A.scalerCharge",&Scal_BCM4A_charge);
 Double_t  Scal_BCM4A_current;
   tscal->SetBranchAddress("P.BCM4A.scalerCurrent",&Scal_BCM4A_current);
 Double_t  Scal_Unser_rate;
   tscal->SetBranchAddress("P.Unser.scalerRate",&Scal_Unser_rate);
 Double_t  Scal_Unser_charge;
   tscal->SetBranchAddress("P.Unser.scalerCharge",&Scal_Unser_charge);
 Double_t  Scal_Unser_current;
   tscal->SetBranchAddress("P.Unser.scalerCurrent",&Scal_Unser_current);
 Double_t  Scal_time;
   tscal->SetBranchAddress("P.1MHz.scalerTime",&Scal_time);
   //
Long64_t scal_entries = tscal->GetEntries();
 cout << " scal ent = " << scal_entries << endl;
 Double_t nlast=float(scal_entries);
 TH1F *h_cur4A_entry = new TH1F("h_cur4A_entry","; ENtry; 4A current",nlast,0,nlast);
 TH1F *h_cur4B_entry = new TH1F("h_cur4B_entry","; ENtry; 4B current",nlast,0,nlast);
 TH1F *h_cur4C_entry = new TH1F("h_cur4C_entry","; ENtry; 4C current",nlast,0,nlast);
 TH1F *h_curUnserRate_entry = new TH1F("h_curUnserRate_entry","; ENtry; Unser Rate",nlast,0,nlast);
 TH1F *h_curUnserCurrent_entry = new TH1F("h_curUnserCurrent_entry","; ENtry; Unser Current",nlast,0,nlast);
 TH1F *h_cur4A_time = new TH1F("h_cur4A_time","; Time; 4A current",2200,0,22000);
 TH1F *h_cur4B_time = new TH1F("h_cur4B_time","; Time; 4B current",2200,0,22000);
 TH1F *h_cur4C_time = new TH1F("h_cur4C_time","; Time; 4C current",2200,0,22000);
 TH1F *h_curUnserRate_time = new TH1F("h_curUnserRate_time","; Time; Unser Rate",2200,0,22000);
 TH1F *h_curUnserCurrent_time = new TH1F("h_curUnserCurrent_time","; Time; Unser Current",2200,0,22000);
 //
 Double_t UnserRateLowEntry[1000];
 Double_t UnserRateHighEntry[1000];
 Double_t UnserGoodLowEntry[1000];
 Double_t UnserGoodHighEntry[1000];
 Int_t NZero_regions;
 Int_t NBeam_regions;
 Double_t Unser_zero = 323000;
 Double_t Unser_zero_range = 10000;
 Int_t BeamRange=8; // number of entries

		Double_t save_unser1=0;
		Double_t save_unser2=0;
		Double_t unser_thres=0;
		Int_t nreg=0;
	for (int i = 0; i < scal_entries; i++) {
                if (i%50000==0) cout << " Entry = " << i << endl;
		UnserRateLowEntry[nreg]=i;
    		tscal->GetEntry(i++);
	        save_unser1=Scal_Unser_rate;
     		tscal->GetEntry(i);
		save_unser2=Scal_Unser_rate;
		//cout << "entry = " << i << " abs  = " << abs((save_unser1-save_unser2)/save_unser1) << " " <<  save_unser1 << " " <<  save_unser2 << endl;
		if (Scal_Unser_rate-323000 < 3000) unser_thres=.04;
		if (Scal_Unser_rate-323000 >= 3000) unser_thres=.02;
		while (abs((save_unser1-save_unser2)/save_unser1) < unser_thres && i<scal_entries) {
     		tscal->GetEntry(i++);
	        save_unser1=Scal_Unser_rate;
     		tscal->GetEntry(i);
		save_unser2=Scal_Unser_rate;
		//		cout << nreg << " entry  " << i << " " << abs((save_unser1-save_unser2)/save_unser1) << " " << save_unser1 << " " <<  save_unser2 << endl;
		
		  }
		UnserRateHighEntry[nreg++]=i-1;
		//cout << nreg-1 << " " << UnserRateLowEntry[nreg-1] << " " <<  UnserRateHighEntry[nreg-1] << endl;
		i=i-1;
		if (nreg>1000) nreg=1000;
		//		cin>> nreg;
	}
	//
	Int_t ngood=0;
	Int_t nrange=6;
	for (int i = 0; i < nreg; i++) {
	  if  ( (UnserRateHighEntry[i]-UnserRateLowEntry[i]) > nrange) {
	    UnserGoodHighEntry[ngood]=UnserRateHighEntry[i];
	    UnserGoodLowEntry[ngood++]=UnserRateLowEntry[i];
	    cout << " good = " << ngood-1 << " " << UnserGoodLowEntry[ngood-1]<< " " << UnserGoodHighEntry[ngood-1] << endl;
	  }
	}
	//
	Int_t ng=0;
	Double_t UnserNumGood[1000];
	Double_t UnserAveGood[1000];
	for (Int_t j=0;j<1000;j++) {
	  UnserNumGood[j]=0;
	  UnserAveGood[j]=0.0;
	}
	for (int i = 0; i < scal_entries; i++) {
      		tscal->GetEntry(i);
		if (i>=UnserGoodLowEntry[ng] && i<=UnserGoodHighEntry[ng]) {
		  UnserNumGood[ng]++;
		  if (i==UnserGoodLowEntry[ng]) {
		    UnserAveGood[ng]=Scal_Unser_rate;
		  } else {
                      UnserAveGood[ng]+=Scal_Unser_rate;
		  }
		  if (i==UnserGoodHighEntry[ng]) {
		    UnserAveGood[ng]=UnserAveGood[ng]/UnserNumGood[ng];
		    cout << " ng = " << ng << " Ave = " << UnserAveGood[ng] << " numgood = " << UnserNumGood[ng] << endl;
		    ng++;
		    UnserAveGood[ng]=0;
		  }		    
		}
	  h_cur4C_entry->Fill(float(i),Scal_BCM4C_rate);
	  h_cur4B_entry->Fill(float(i),Scal_BCM4B_rate);
	  h_cur4A_entry->Fill(float(i),Scal_BCM4A_rate);
	  h_curUnserRate_entry->Fill(float(i),Scal_Unser_rate);
	  h_curUnserCurrent_entry->Fill(float(i),Scal_Unser_current);
	  h_cur4C_time->Fill(Scal_time,Scal_BCM4C_rate);
	  h_cur4B_time->Fill(Scal_time,Scal_BCM4B_rate);
	  h_cur4A_time->Fill(Scal_time,Scal_BCM4A_rate);
	  h_curUnserRate_time->Fill(Scal_time,Scal_Unser_rate);
	  h_curUnserCurrent_time->Fill(Scal_time,Scal_Unser_current);
	}
	//
	cout << " number of good regions = " << ng-1 << endl;
        Int_t ng_cur_regions=0;
        Int_t ng_zero_regions=0;
	Int_t UnserCurNumGood[1000];
	Int_t UnserZeroNumGood[1000];
	for (int i = 0; i < ng; i++) {
	  if (UnserAveGood[i]>Unser_zero*1.03) {
	    UnserCurNumGood[ng_cur_regions]=i;
	    cout << " Cur ng = " << i << " " << ng_cur_regions << " " << UnserAveGood[i] << endl;
	    ng_cur_regions++;
	  } else {
	    UnserZeroNumGood[ng_zero_regions]=i;
	    cout << " zero ng = " << i << " " << ng_zero_regions << " " << UnserAveGood[i] << endl;
	    ng_zero_regions++;
	  }
	}
	// Calculate global average for zero Unser
	Double_t global_zero_unser=0;
	for (Int_t i = 0; i < ng_zero_regions; i++) {
	  global_zero_unser+=UnserAveGood[UnserZeroNumGood[i]];
	  //cout << UnserZeroNumGood[i] << " " << UnserAveGood[UnserZeroNumGood[i]] << endl;
	}      
	global_zero_unser=global_zero_unser/(ng_zero_regions);
	cout << " Global unser zero = " <<  global_zero_unser << endl;
	// match zero regions with current regions
	Int_t cur_test;
	Int_t lowgood;
	Int_t highgood;
	Int_t zero_test;
	Int_t lowzero;
	Int_t highzero;
	Double_t UnserCurGoodZero[1000];
	Double_t UnserCurGoodNumZero[1000];
	Int_t Zero_region_range=8;
	for (int i = 0; i < 1000; i++) {
	  UnserCurGoodZero[i]=0.;
	  UnserCurGoodNumZero[i]=0.;
	}
	for (int i = 0; i < ng_cur_regions; i++) {
	  cur_test=UnserCurNumGood[i];
	  lowgood=UnserGoodLowEntry[cur_test];
	  highgood=UnserGoodHighEntry[cur_test];
	     for (int j = 0; j < ng_zero_regions; j++) {
	       zero_test=UnserZeroNumGood[j];
	       lowzero=UnserGoodLowEntry[zero_test];
	       highzero=UnserGoodHighEntry[zero_test];
	       if ( lowgood > highzero && abs(lowgood-highzero) < Zero_region_range) {
		 cout << " Match low good " << " ng_cur = " << i <<" " << zero_test << " " << highzero << " " << lowgood << endl;
                 UnserCurGoodZero[i]+=UnserAveGood[zero_test];
		 UnserCurGoodNumZero[i]++;
	       }
	       if ( lowzero > highgood && abs(lowzero-highgood) < Zero_region_range) {
		 cout << " Match high good " << " ng_cur = " << i  <<" "<< zero_test << " " << lowzero << " " << highgood << endl;
                 UnserCurGoodZero[i]+=UnserAveGood[zero_test];
		 UnserCurGoodNumZero[i]++;
	       }
	     }	  
	}
	//
	//
	for (int i = 0; i < ng_cur_regions; i++) {
	  Double_t unser_zero;
	  if (UnserCurGoodNumZero[i]==0) {
	    unser_zero=global_zero_unser;
	      } else {
	    unser_zero=UnserCurGoodZero[i]/UnserCurGoodNumZero[i];
	  }
	  cout << " Cur reg = " << i << " " << unser_zero << " " << UnserCurNumGood[i] << " " << UnserAveGood[UnserCurNumGood[i]] << " " << (UnserAveGood[UnserCurNumGood[i]]-unser_zero)/5020.<< endl;
	}
	//
  TCanvas *cUnser = new TCanvas("cUnser","Unser",700,700);
  TLine *sline;
  cUnser->Divide(1,1);
  cUnser->cd(1);
  h_curUnserRate_entry->Draw("HIST ");
  for (Int_t n=0;n<ng;n++) {
    sline = new TLine(UnserGoodLowEntry[n]+.5,UnserAveGood[n],UnserGoodHighEntry[n]+.5,UnserAveGood[n]);
    sline->SetLineColor(2);
    sline->SetLineWidth(2);
    sline->Draw();
  }
   //
}

