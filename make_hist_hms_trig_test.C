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

void make_hist_hms_trig_test(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_trig_test_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 static const Int_t trignum=11;
 const char* trigname[trignum]={
 "D.hdcref1",
 "D.hdcref2",
 "D.hdcref3",
 "D.hdcref4",
 "D.hhodoref1",
 "D.htrigref1",
 "D.htrigref2",
 "D.htrigref3",
 "D.htrig1",
 "D.htrig2",
 "D.htrig3",
};
 Double_t data[trignum][50];
 Int_t ndata[trignum];
 for (Int_t ip=0;ip<trignum;ip++) {
   tsimc->SetBranchAddress(Form("%s",trigname[ip]),&data[ip]) ;
   tsimc->SetBranchAddress(Form("Ndata.%s",trigname[ip]),&ndata[ip]) ;
 }

   // Define histograms
 TH1F *hist[trignum][5][5];
 TH1F *nhist[trignum];
 Int_t nbins=400.;
 Int_t xhi=2000.;
 for (Int_t ip=0;ip<trignum;ip++) {
   xhi=2000;
   nbins=400;
   if ( ip >4) xhi=600;
   if ( ip >4) nbins=300;
   nhist[ip]= new TH1F(Form("hist_ndata_%s",trigname[ip]),Form("; %s Mult  ; Counts ",trigname[ip]),10,0,10);
   HList.Add( nhist[ip]);
   for (Int_t ic=0;ic<5;ic++) {
     for (Int_t ih=0;ih<ic+1;ih++) {
       hist[ip][ic][ih]= new TH1F(Form("hist_%s_m%d_h%d",trigname[ip],ic+1,ih+1),Form("; %s (ns) hit = %d for Mult = %d ; Counts ",trigname[ip],ih+1,ic+1),nbins,0,xhi);
     HList.Add( hist[ip][ic][ih]);
   }
   }
 }
 TH1F *dhist_mult2;
 TH1F *dhist_mult3_12;
 TH1F *dhist_mult3_23;
 TH1F *dhist_mult4_12;
 TH1F *dhist_mult4_23;
 TH1F *dhist_mult4_34;
 nbins=125;
 xhi=250;
 Int_t xlo=   0;
 dhist_mult2= new TH1F(Form("hist_%s_mult_2_diff_1_2",trigname[5] ) ,Form("; %s Hit 2- Hit 1 Time( ns) Mult = 2 ; Counts ",trigname[5]),nbins,xlo,xhi); 
 HList.Add(dhist_mult2);
 dhist_mult3_12= new TH1F(Form("hist_%s_mult_3_diff_1_2",trigname[5]),Form("; %s Hit 2- Hit 1  time (ns) Mult = 3 ; Counts ",trigname[5]),nbins,xlo,xhi); 
 HList.Add(dhist_mult3_12);
 dhist_mult3_23= new TH1F(Form("hist_%s_mult_3_diff_2_3",trigname[5]),Form("; %s  Hit 3- Hit 2 time (ns) Mult = 3  ; Counts ",trigname[5]),nbins,xlo,xhi); 
 HList.Add(dhist_mult3_23);
 dhist_mult4_12= new TH1F(Form("hist_%s_mult_4_diff_1_2",trigname[5]),Form("; %s  Hit 2- Hit 1 time (ns) Mult = 4 ; Counts ",trigname[5]),nbins,xlo,xhi); 
 HList.Add(dhist_mult4_12);
 dhist_mult4_23= new TH1F(Form("hist_%s_mult_4_diff_2_3",trigname[5]),Form("; %s Hit 3- Hit 2 time (ns) Mult = 4; Counts ",trigname[5]),nbins,xlo,xhi); 
 HList.Add(dhist_mult4_23);
 dhist_mult4_34= new TH1F(Form("hist_%s_mult_4_diff_3_4",trigname[5]),Form("; %s Hit 4- Hit 3 time (ns) Mult = 4 ; Counts ",trigname[5]),nbins,xlo,xhi); 
 HList.Add(dhist_mult4_34);
 TH2F *mult_trig1_2;
 mult_trig1_2 =new TH2F(Form("hist_mult_%s_%s",trigname[5],trigname[6]),Form("; %s Mult ; %s Mult ",trigname[5],trigname[6]),10,0,10,10,0,10); 
 HList.Add(mult_trig1_2);
 TH2F *mult_trig1_3;
 mult_trig1_3 =new TH2F(Form("hist_mult_%s_%s",trigname[5],trigname[7]),Form("; %s Mult ; %s Mult ",trigname[5],trigname[7]),10,0,10,10,0,10); 
 HList.Add(mult_trig1_3);
 TH2F *mult_trig1_dc1;
 mult_trig1_dc1 =new TH2F(Form("hist_mult_%s_%s",trigname[5],trigname[0]),Form("; %s Mult ; %s Mult ",trigname[5],trigname[0]),10,0,10,10,0,10); 
 HList.Add(mult_trig1_dc1);
 TH2F *mult_trig2_dc1;
 mult_trig2_dc1 =new TH2F(Form("hist_mult_%s_%s",trigname[6],trigname[0]),Form("; %s Mult ; %s Mult ",trigname[6],trigname[0]),10,0,10,10,0,10); 
 HList.Add(mult_trig2_dc1);
 TH1F *hist_tref2_mult;
 hist_tref2_mult= new TH1F(Form("hist_%s_mult_cut",trigname[6]),Form("; %s Mult for Trigref1 mult==1 and time<100; Counts ",trigname[6]),10,0,10); 
 HList.Add(hist_tref2_mult);
 TH1F *hist_tref2_mult3_diff12;
 hist_tref2_mult3_diff12= new TH1F(Form("hist_%s_mult3_tdiff12_cut",trigname[6]),Form("; %s Mult=3 Hit 2-1 time for T1 mult==1 and time<100; Counts ",trigname[6]),nbins,xlo,xhi); 
 HList.Add(hist_tref2_mult3_diff12);
 TH1F *hist_tref2_mult3_diff23;
 hist_tref2_mult3_diff23= new TH1F(Form("hist_%s_mult3_tdiff23_cut",trigname[6]),Form("; %s Mult=3 Hit 3-2 time for T1 mult==1 and time<100; Counts ",trigname[6]),nbins,xlo,xhi); 
 HList.Add(hist_tref2_mult3_diff23);
 TH1F *hist_tref2_mult3_time[3];
		       for (Int_t ih=0;ih<3;ih++) {
			 hist_tref2_mult3_time[ih]= new TH1F(Form("hist_%s_mult3_hit%d_time_cut",trigname[6],ih),Form("T1 mult==1 and time<100; %s Mult=3 Hit %d Time ; Counts ",trigname[6],ih),nbins,xlo,500); 
 HList.Add(hist_tref2_mult3_time[ih]);
		       }
		       xhi=550.;
 TH2F *tref2time_mult3_hit12;
 tref2time_mult3_hit12 =new TH2F(Form("hist_time_%s_hit12",trigname[6]),Form("%s Mult =3 ; Hit 1 time ; Hit 2 - Hit 1 time ",trigname[6]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref2time_mult3_hit12);
 TH2F *tref2time_mult3_hit23;
 tref2time_mult3_hit23 =new TH2F(Form("hist_time_%s_hit23",trigname[6]),Form("%s Mult =3 ; Hit 2 time ; Hit 3 - Hit 1 time ",trigname[6]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref2time_mult3_hit23);
 TH2F *tref2time_mult3_hit13;
 tref2time_mult3_hit13 =new TH2F(Form("hist_time_%s_hit13",trigname[6]),Form("%s Mult =3 ; Hit 1 time ; Hit 3 - Hit 1 time ",trigname[6]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref2time_mult3_hit13);
 TH2F *tref2time_mult3_hit31;
 tref2time_mult3_hit31 =new TH2F(Form("hist_time_%s_hit31",trigname[6]),Form("%s Mult =3 ; Hit 3 time ; Hit 3 - Hit 1 time ",trigname[6]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref2time_mult3_hit31);
 TH2F *tref1time_mult3_hit12;
 tref1time_mult3_hit12 =new TH2F(Form("hist_time_%s_hit12",trigname[5]),Form("%s Mult =3 ; Hit 1 time ; Hit 2 - Hit 1 time ",trigname[5]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit12);
 TH2F *tref1time_mult3_hit23;
 tref1time_mult3_hit23 =new TH2F(Form("hist_time_%s_hit23",trigname[5]),Form("%s Mult =3 ; Hit 2 time ; Hit 3 - Hit 2 time ",trigname[5]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit23);
 TH2F *tref1time_mult3_hit13;
 tref1time_mult3_hit13 =new TH2F(Form("hist_time_%s_hit13",trigname[5]),Form("%s Mult =3 ; Hit 1 time ; Hit 3 - Hit 1 time ",trigname[5]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit13);
 TH2F *tref1time_mult3_hit31;
 tref1time_mult3_hit31 =new TH2F(Form("hist_time_%s_hit31",trigname[5]),Form("%s Mult =3 ; Hit 3 time ; Hit 3 - Hit 1 time ",trigname[5]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit31);

 // loop over entries
Long64_t nentries = tsimc->GetEntries();
          for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		mult_trig1_2->Fill(ndata[5],ndata[6]);
		mult_trig1_dc1->Fill(ndata[5],ndata[0]);
		mult_trig2_dc1->Fill(ndata[6],ndata[0]);
		mult_trig1_3->Fill(ndata[5],ndata[7]);
                   for (Int_t ip=0;ip<trignum;ip++) {
                       nhist[ip]->Fill(ndata[ip]);
		       if (ndata[ip] <6) {
		       for (Int_t ih=0;ih<ndata[ip];ih++) {
                              hist[ip][ndata[ip]-1][ih]->Fill(.1*data[ip][ih]);
			  }
		       }
		   }
		   if (ndata[5]==3) {
		     tref1time_mult3_hit12->Fill(.1*data[5][0],.1*data[5][1]-.1*data[5][0]);
		     tref1time_mult3_hit23->Fill(.1*data[5][1],.1*data[5][2]-.1*data[5][1]);
		     tref1time_mult3_hit13->Fill(.1*data[5][0],.1*data[5][2]-.1*data[5][0]);
		     tref1time_mult3_hit31->Fill(.1*data[5][2],.1*data[5][2]-.1*data[5][0]);
		   }
		   if (ndata[6]==3) {
		     tref2time_mult3_hit12->Fill(.1*data[6][0],.1*data[6][1]-.1*data[6][0]);
		     tref2time_mult3_hit23->Fill(.1*data[6][1],.1*data[6][2]-.1*data[6][1]);
		     tref2time_mult3_hit13->Fill(.1*data[6][0],.1*data[6][2]-.1*data[6][0]);
		     tref2time_mult3_hit31->Fill(.1*data[6][2],.1*data[6][2]-.1*data[6][0]);
		   }
		   if (ndata[5]==1 && .1*data[5][0]<100) {
		     hist_tref2_mult->Fill(ndata[6]);
		     if (ndata[6]==3) {
		       for (Int_t ih=0;ih<3;ih++) {
			 hist_tref2_mult3_time[ih]->Fill(.1*data[6][ih]);
		       }
		     hist_tref2_mult3_diff12->Fill(.1*data[6][1]-.1*data[6][0]);
		     hist_tref2_mult3_diff23->Fill(.1*data[6][2]-.1*data[6][1]);
		     }
		   }
		   if (ndata[5]==2) {
		     dhist_mult2->Fill(.1*data[5][1]-.1*data[5][0]);
                   }
		   if (ndata[5]==3) {
		     dhist_mult3_12->Fill(.1*data[5][1]-.1*data[5][0]);
		     dhist_mult3_23->Fill(.1*data[5][2]-.1*data[5][1]);
                   }
		   if (ndata[5]==4) {
		     dhist_mult4_12->Fill(.1*data[5][1]-.1*data[5][0]);
		     dhist_mult4_23->Fill(.1*data[5][2]-.1*data[5][1]);
		     dhist_mult4_34->Fill(.1*data[5][3]-.1*data[5][2]);
                   }
	
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
