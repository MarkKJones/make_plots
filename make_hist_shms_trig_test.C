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

void make_hist_shms_trig_test(TString basename="",Int_t nrun=2043){
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
 static const Int_t trignum=15;
 const char* trigname[trignum]={
 "D.pdcref1",
 "D.pdcref2",
 "D.pdcref3",
 "D.pdcref4",
 "D.pdcref5",
 "D.pdcref6",
 "D.pdcref7",
 "D.pdcref8",
 "D.pdcref9",
 "D.pdcref10",
 "D.phodoref1",
"D.phodoref2",
 "D.ptrigref1",
 "D.ptrigref2",
"D.pedtm"
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
   if ( ip >9) xhi=1000;
   if ( ip >9) nbins=500;
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
 nbins=250;
 xhi=500;
 Int_t xlo=   0;
 dhist_mult2= new TH1F(Form("hist_%s_mult_2_diff_1_2",trigname[12]),Form("; %s Hit 2- Hit 1 Time( ns) Mult = 2 ; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_mult2);
 dhist_mult3_12= new TH1F(Form("hist_%s_mult_3_diff_1_2",trigname[12]),Form("; %s Hit 2- Hit 1  time (ns) Mult = 3 ; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_mult3_12);
 dhist_mult3_23= new TH1F(Form("hist_%s_mult_3_diff_2_3",trigname[12]),Form("; %s  Hit 3- Hit 2 time (ns) Mult = 3  ; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_mult3_23);
 dhist_mult4_12= new TH1F(Form("hist_%s_mult_4_diff_1_2",trigname[12]),Form("; %s  Hit 2- Hit 1 time (ns) Mult = 4 ; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_mult4_12);
 dhist_mult4_23= new TH1F(Form("hist_%s_mult_4_diff_2_3",trigname[12]),Form("; %s Hit 3- Hit 2 time (ns) Mult = 4; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_mult4_23);
 dhist_mult4_34= new TH1F(Form("hist_%s_mult_4_diff_3_4",trigname[12]),Form("; %s Hit 4- Hit 3 time (ns) Mult = 4 ; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_mult4_34);
 TH1F *dhist_cut_mult2;
 TH1F *dhist_cut_mult3_12;
 TH1F *dhist_cut_mult3_23;
 TH1F *dhist_cut_mult4_12;
 TH1F *dhist_cut_mult4_23;
 TH1F *dhist_cut_mult4_34;
 nbins=250;
 xhi=500;
 xlo=   0;
 dhist_cut_mult2= new TH1F(Form("hist_cut_%s_mult_2_diff_1_2",trigname[12]),Form("; %s Hit 2- Hit 1 Time( ns) Mult = 2 NDC1==0; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_cut_mult2);
 dhist_cut_mult3_12= new TH1F(Form("hist_cut_%s_mult_3_diff_1_2",trigname[12]),Form("; %s Hit 2- Hit 1  time (ns) Mult = 3 NDC1==0; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_cut_mult3_12);
 dhist_cut_mult3_23= new TH1F(Form("hist_cut_%s_mult_3_diff_2_3",trigname[12]),Form("; %s  Hit 3- Hit 2 time (ns) Mult = 3 NDC1==0 ; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_cut_mult3_23);
 dhist_cut_mult4_12= new TH1F(Form("hist_cut_%s_mult_4_diff_1_2",trigname[12]),Form("; %s  Hit 2- Hit 1 time (ns) Mult = 4 NDC1==0; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_cut_mult4_12);
 dhist_cut_mult4_23= new TH1F(Form("hist_cut_%s_mult_4_diff_2_3",trigname[12]),Form("; %s Hit 3- Hit 2 time (ns) Mult = 4 NDC1==0 ; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_cut_mult4_23);
 dhist_cut_mult4_34= new TH1F(Form("hist_cut_%s_mult_4_diff_3_4",trigname[12]),Form("; %s Hit 4- Hit 3 time (ns) Mult = 4 NDC1==0; Counts ",trigname[12]),nbins,xlo,xhi); 
 HList.Add(dhist_cut_mult4_34);
 TH2F *dhist_cut_mult3_12_23;
 dhist_cut_mult3_12_23= new TH2F(Form("hist_cut_%s_mult_3_diff_1_2_diff_2_3",trigname[12]),Form("; %s Hit 2- Hit 1  time (ns) Mult = 3 NDC1==0; Hit 3- Hit 2  time (ns) ",trigname[12]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(dhist_cut_mult3_12_23);
 xhi=600.;
 TH2F *tref2time_mult3_hit12;
 tref2time_mult3_hit12 =new TH2F(Form("hist_time_%s_hit12",trigname[13]),Form("%s Mult =3 ; Hit 1 time ; Hit 2 - Hit 1 time ",trigname[13]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref2time_mult3_hit12);
 TH2F *tref2time_mult3_hit23;
 tref2time_mult3_hit23 =new TH2F(Form("hist_time_%s_hit23",trigname[13]),Form("%s Mult =3 ; Hit 2 time ; Hit 3 - Hit 1 time ",trigname[13]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref2time_mult3_hit23);
 TH2F *tref2time_mult3_hit13;
 tref2time_mult3_hit13 =new TH2F(Form("hist_time_%s_hit13",trigname[13]),Form("%s Mult =3 ; Hit 1 time ; Hit 3 - Hit 1 time ",trigname[13]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref2time_mult3_hit13);
 xhi=2000.;
 TH2F *tdc1time_mult3_hit12;
 tdc1time_mult3_hit12 =new TH2F(Form("hist_time_%s_hit12",trigname[0]),Form("%s Mult =3 ; Hit 1 time ; Hit 2 - Hit 1 time ",trigname[0]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tdc1time_mult3_hit12);
 TH2F *tdc1time_mult3_hit23;
 tdc1time_mult3_hit23 =new TH2F(Form("hist_time_%s_hit23",trigname[0]),Form("%s Mult =3 ; Hit 2 time ; Hit 3 - Hit 1 time ",trigname[0]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tdc1time_mult3_hit23);
 TH2F *tdc1time_mult3_hit13;
 tdc1time_mult3_hit13 =new TH2F(Form("hist_time_%s_hit13",trigname[0]),Form("%s Mult =3 ; Hit 1 time ; Hit 3 - Hit 1 time ",trigname[0]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tdc1time_mult3_hit13);
 xhi=700.;
 TH2F *tref1time_mult3_hit12;
 tref1time_mult3_hit12 =new TH2F(Form("hist_time_%s_hit12",trigname[12]),Form("%s Mult =3 ; Hit 1 time ; Hit 2 - Hit 1 time ",trigname[12]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit12);
 TH2F *tref1time_mult3_hit23;
 tref1time_mult3_hit23 =new TH2F(Form("hist_time_%s_hit23",trigname[12]),Form("%s Mult =3 ; Hit 2 time ; Hit 3 - Hit 2 time ",trigname[12]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit23);
 TH2F *tref1time_mult3_hit23_cut1;
 tref1time_mult3_hit23_cut1 =new TH2F(Form("hist_time_%s_hit23_cut1",trigname[12]),Form("%s Mult =3 Hit 1 time < 300; Hit 2 time ; Hit 3 - Hit 2 time ",trigname[12]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit23_cut1);
 TH2F *tref1time_mult3_hit13;
 tref1time_mult3_hit13 =new TH2F(Form("hist_time_%s_hit13",trigname[12]),Form("%s Mult =3 ; Hit 1 time ; Hit 3 - Hit 1 time ",trigname[12]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit13);
 TH2F *tref1time_mult3_hit13_cut1;
 tref1time_mult3_hit13_cut1=new TH2F(Form("hist_time_%s_hit13_cut1",trigname[12]),Form("%s Mult =3 Hit 1 < 300 ; Hit 1 time ; Hit 3 - Hit 1 time ",trigname[12]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit13_cut1);
 TH2F *tref1time_mult3_hit31;
 tref1time_mult3_hit31 =new TH2F(Form("hist_time_%s_hit31",trigname[12]),Form("%s Mult =3 ; Hit 3 time ; Hit 3 - Hit 1 time ",trigname[12]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit31);
 TH2F *tref1time_mult3_hit31_cut;
 tref1time_mult3_hit31_cut =new TH2F(Form("hist_time_%s_hit31_cut",trigname[12]),Form("%s Mult =3 Hit 1 time < 300; Hit 3 time ; Hit 3 - Hit 1 time ",trigname[12]),nbins,xlo,xhi,nbins,xlo,xhi); 
 HList.Add(tref1time_mult3_hit31_cut);
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
                for (Int_t ip=0;ip<trignum;ip++) {
                       nhist[ip]->Fill(ndata[ip]);
		       if (ndata[ip] <6) {
                          for (Int_t ih=0;ih<ndata[ip];ih++) {
                              hist[ip][ndata[ip]-1][ih]->Fill(.1*data[ip][ih]);
			  }
		       }
		}
		   if (ndata[12]==3) {
		     tref1time_mult3_hit12->Fill(.1*data[12][0],.1*data[12][1]-.1*data[12][0]);
		     tref1time_mult3_hit23->Fill(.1*data[12][1],.1*data[12][2]-.1*data[12][1]); 
		     if (.1*data[12][0]<300) tref1time_mult3_hit23_cut1->Fill(.1*data[12][1],.1*data[12][2]-.1*data[12][1]);
		     tref1time_mult3_hit13->Fill(.1*data[12][0],.1*data[12][2]-.1*data[12][0]);
		     if (.1*data[12][0]<300) tref1time_mult3_hit13_cut1->Fill(.1*data[12][0],.1*data[12][2]-.1*data[12][0]);
		     tref1time_mult3_hit31->Fill(.1*data[12][2],.1*data[12][2]-.1*data[12][0]);
		     if (.1*data[12][0]<300) tref1time_mult3_hit31_cut->Fill(.1*data[12][2],.1*data[12][2]-.1*data[12][0]);
		   }
		   if (ndata[13]==3) {
		     tref2time_mult3_hit12->Fill(.1*data[13][0],.1*data[13][1]-.1*data[13][0]);
		     tref2time_mult3_hit23->Fill(.1*data[13][1],.1*data[13][2]-.1*data[13][1]);
		     tref2time_mult3_hit13->Fill(.1*data[13][0],.1*data[13][2]-.1*data[13][0]);
		   }
		   if (ndata[0]==3) {
		     tdc1time_mult3_hit12->Fill(.1*data[0][0],.1*data[0][1]-.1*data[0][0]);
		     tdc1time_mult3_hit23->Fill(.1*data[0][1],.1*data[0][2]-.1*data[0][1]);
		     tdc1time_mult3_hit13->Fill(.1*data[0][0],.1*data[0][2]-.1*data[0][0]);
		   }
		   if (ndata[12]==2) {
		     dhist_mult2->Fill(.1*data[12][1]-.1*data[12][0]);
		     if(ndata[0]==0) dhist_cut_mult2->Fill(.1*data[12][1]-.1*data[12][0]);
                   }
		   if (ndata[12]==3) {
		     dhist_mult3_12->Fill(.1*data[12][1]-.1*data[12][0]);
		     dhist_mult3_23->Fill(.1*data[12][2]-.1*data[12][1]);
		     if(ndata[0]==0) dhist_cut_mult3_12->Fill(.1*data[12][1]-.1*data[12][0]);
		     if(ndata[0]==0) dhist_cut_mult3_23->Fill(.1*data[12][2]-.1*data[12][1]);
		     if(ndata[0]==0) dhist_cut_mult3_12_23->Fill(.1*data[12][1]-.1*data[12][0],.1*data[12][2]-.1*data[12][1]);
                   }
		   if (ndata[12]==4) {
		     dhist_mult4_12->Fill(.1*data[12][1]-.1*data[12][0]);
		     dhist_mult4_23->Fill(.1*data[12][2]-.1*data[12][1]);
		     dhist_mult4_34->Fill(.1*data[12][3]-.1*data[12][2]);
		     if(ndata[0]==0) dhist_cut_mult4_12->Fill(.1*data[12][1]-.1*data[12][0]);
		     if(ndata[0]==0) dhist_cut_mult4_23->Fill(.1*data[12][2]-.1*data[12][1]);
		     if(ndata[0]==0) dhist_cut_mult4_34->Fill(.1*data[12][3]-.1*data[12][2]);
                   }
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
