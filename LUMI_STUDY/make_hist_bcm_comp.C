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

void make_hist_bcm_comp(TString basename="",Int_t nrun=2043,Double_t mean_current=2.0){
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
   outputhist= "hist/"+basename+"_bcm_comp_hist.root";
 TObjArray HList(0);
  TString outputpdf;
  outputpdf= Form("plots/bcm_comp_run_%d.pdf",nrun);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tscal = (TTree*) fsimc->Get("TSP");
//
 Double_t  Scal_evNumber;
   tscal->SetBranchAddress("evNumber",&Scal_evNumber);
 Double_t  Scal_time;
   tscal->SetBranchAddress("P.1MHz.scalerTime",&Scal_time);
 Double_t  Scal_BCM4B_charge;
   tscal->SetBranchAddress("P.BCM4B.scalerCharge",&Scal_BCM4B_charge);
 Double_t  Scal_BCM4B_current;
   tscal->SetBranchAddress("P.BCM4B.scalerCurrent",&Scal_BCM4B_current);
 Double_t  Scal_BCM4A_charge;
   tscal->SetBranchAddress("P.BCM4A.scalerCharge",&Scal_BCM4A_charge);
 Double_t  Scal_BCM4A_current;
   tscal->SetBranchAddress("P.BCM4A.scalerCurrent",&Scal_BCM4A_current);
 Double_t  Scal_BCM1_charge;
   tscal->SetBranchAddress("P.BCM1.scalerCharge",&Scal_BCM1_charge);
 Double_t  Scal_BCM1_current;
   tscal->SetBranchAddress("P.BCM1.scalerCurrent",&Scal_BCM1_current);
 Double_t  Scal_BCM2_charge;
   tscal->SetBranchAddress("P.BCM2.scalerCharge",&Scal_BCM2_charge);
 Double_t  Scal_BCM2_current;
   tscal->SetBranchAddress("P.BCM2.scalerCurrent",&Scal_BCM2_current);
   //
   vector<Double_t> rat_bcm1_bcm2;
   vector<Double_t> rat_bcm1_bcm4a;
   vector<Double_t> rat_bcm1_bcm4b;
   vector<Double_t> rat_bcm4b_bcm4a;
   vector<Double_t> err;
   vector<Double_t> nent;
   
    Long64_t scal_entries = tscal->GetEntries();
	for (int i = 0; i < scal_entries; i++) {
      		tscal->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		//		cout << " entry = " << i << " " << Scal_BCM4A_current << " " << Scal_BCM4B_current/Scal_BCM4A_current  <<" " << Scal_BCM1_current/Scal_BCM4A_current <<" " << Scal_BCM2_current/Scal_BCM4A_current << " "  << Scal_BCM4B_current/Scal_BCM1_current  <<" " << Scal_BCM2_current/Scal_BCM1_current << endl; 
		if (Scal_BCM4A_current > 12.) {
		rat_bcm1_bcm2.push_back(Scal_BCM1_current/Scal_BCM2_current);
		rat_bcm1_bcm4a.push_back(Scal_BCM1_current/Scal_BCM4A_current);
		rat_bcm1_bcm4b.push_back(Scal_BCM1_current/Scal_BCM4B_current);
		rat_bcm4b_bcm4a.push_back(Scal_BCM4B_current/Scal_BCM4A_current);
	        err.push_back(.00001);
		nent.push_back(Scal_time/60);
		}
	}
   TMultiGraph *mgt2 = new TMultiGraph();
   TLegend* legt2 = new TLegend(0.2,0.75,0.48,0.89);
   TGraphErrors *gr_ph1;
   gr_ph1 = new TGraphErrors(rat_bcm1_bcm2.size(),&(nent[0]),&(rat_bcm1_bcm2[0]),0,&(err[0]));
      gr_ph1->SetMarkerStyle(20);
     gr_ph1->SetMarkerColor(2);
   mgt2->Add(gr_ph1);
   TGraphErrors *gr_ph2;
   gr_ph2 = new TGraphErrors(rat_bcm1_bcm4a.size(),&(nent[0]),&(rat_bcm1_bcm4a[0]),0,&(err[0]));
     gr_ph2->SetMarkerStyle(21);
     gr_ph2->SetMarkerColor(3);
    mgt2->Add(gr_ph2);
   TGraphErrors *gr_ph3;
   gr_ph3 = new TGraphErrors(rat_bcm4b_bcm4a.size(),&(nent[0]),&(rat_bcm4b_bcm4a[0]),0,&(err[0]));
     gr_ph3->SetMarkerStyle(22);
     gr_ph3->SetMarkerColor(4);
    mgt2->Add(gr_ph3);
   TGraphErrors *gr_ph4;
   gr_ph4 = new TGraphErrors(rat_bcm1_bcm4b.size(),&(nent[0]),&(rat_bcm1_bcm4b[0]),0,&(err[0]));
     gr_ph4->SetMarkerStyle(23);
     gr_ph4->SetMarkerColor(6);
    mgt2->Add(gr_ph4);
  legt2->SetBorderSize(1);
  legt2->SetFillColor(0);
   legt2->AddEntry(gr_ph1,"BCM1/BCM2 ","p");
   legt2->AddEntry(gr_ph2,"BCM1/BCM4A ","p");
   legt2->AddEntry(gr_ph4,"BCM1/BCM4B ","p");
   legt2->AddEntry(gr_ph3,"BCM4B/BCM4A ","p");
    TCanvas *ctreff = new TCanvas("ctreff","ctreff",700,700);
   ctreff->Divide(1,1);
   ctreff->cd(1);
    mgt2->Draw("AP");
    mgt2->SetTitle(Form("Run %d; Time (min); Current ratio ",nrun));
   mgt2->SetMaximum(1.02);
   mgt2->SetMinimum(.98);
     legt2->Draw();
  ctreff->Print(outputpdf);
    //
}
