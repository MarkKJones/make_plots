#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void comp_hist_CTime(Int_t nr1,Int_t nr2,Int_t nev=100000) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
 //
     TString outputpdf;
     outputpdf=Form("plots/comp_%d_%d_CTime.pdf",nr1,nr2);
  TString inputroot[2];
   TFile *fhistroot[2];
   inputroot[0]=Form("hist/coin_replay_production_%d_%d_CTime_hist.root",nr1,nev);
   inputroot[1]=Form("hist/coin_replay_production_%d_%d_CTime_hist.root",nr2,nev);
     cout << " infile root = " << inputroot[0] << endl;
     cout << " infile root = " << inputroot[1] << endl;
   fhistroot[0] =  new TFile(inputroot[0]);
   fhistroot[1] =  new TFile(inputroot[1]);
   static const Int_t nhist=1;
   //  TString hname[nhist]={"hepicoinTime_ROC2"};
   TString hname[nhist]={"hCTcalc_all"};
    TH1F *hist[2][nhist];
       Double_t fac[2]={1.,1.}; 
       if (nr1==5434) {
	 fac[0]=1./8.07/.93/.98;
	 fac[1]=1./5.38/.92/.99;
       }
       if (nr1==5435) {
	 fac[0]=1./5.38/.92/.99;
	 fac[1]=1./5.38/.92/.99;
       }
       if (nr1==5369) {
         if (hname[0]=="hCTcalc_all") {
	 fac[0]=1./20.27/.99/.90;
	 fac[1]=1./6.14/.98/.78/.95;
         } else {
	 fac[0]=1./20.27/.99/.90;
	 fac[1]=1./6.14/.98/.78/.95;
	 }
       }
       if (nr1==5371) {
	 fac[0]=1./6.14/.99/.84;
	 fac[1]=1./6.14/.99/.84;
       }
  for (Int_t nf=0;nf<2;nf++) {
  for (Int_t ip=0;ip<nhist;ip++) {
    // cout <<  hname[ip] << endl;
       hist[nf][ip] = (TH1F*)fhistroot[nf]->Get(hname[ip]);
       if (!hist[nf][ip]) cout << " no hist = " << hname[ip] << endl;
  }
  }
  //
  TCanvas *can1 = new TCanvas("can1","can1",700,700);
  can1->Divide(1,1);
  can1->cd(1);
  hist[0][0]->Draw();
  Double_t pbin=43.1;
  Double_t bbin=51.3;
   pbin=37.8;
  bbin=45.6;
  if (nr1==5369&&nev==100001&&hname[0]=="hCTcalc_all") {
  pbin=24.73; // ct_calc_all 5371
  bbin=32.88;
  }
  if (nr1==5369&&nev==100002&&hname[0]=="hCTcalc_all") {
  pbin=24.73; // ct_calc_all 5371
  bbin=32.88;
  }
  if (nr1==5369&&nev==100003&&hname[0]=="hCTcalc_all") {
  pbin=24.73; // ct_calc_all 5371
  bbin=32.88;
  }
  if (nr1==5369&&nev==100002&&hname[0]=="hepicoinTime_ROC2") {
  pbin=37.5; 
  bbin=41.48;
  }
  if (nr1==5369&&nev==100000&&hname[0]=="hCTcalc_all") {
    pbin=30.44;
  bbin=38.53;
 }
  Double_t wid=2.0;
  Double_t sum = hist[0][0]->Integral(hist[0][0]->FindBin(pbin-wid),hist[0][0]->FindBin(pbin+wid));
  Double_t back = hist[0][0]->Integral(hist[0][0]->FindBin(bbin-wid),hist[0][0]->FindBin(bbin+wid));
  cout << nr1 << " Integral = " << sum << " Yield = " << sum*fac[0] << " back = " << back*fac[0] << " sum-back = " << (sum-back)*fac[0]<< endl ;
  //
  TCanvas *can2 = new TCanvas("can2","can2",700,700);
  can2->Divide(1,1);
  can2->cd(1);
  hist[1][0]->Draw();
  sum = hist[1][0]->Integral(hist[1][0]->FindBin(pbin-wid),hist[1][0]->FindBin(pbin+wid));
  back = hist[1][0]->Integral(hist[1][0]->FindBin(bbin-wid),hist[1][0]->FindBin(bbin+wid));
  cout << nr2 << "  Integral = " << sum << " Yield = " << sum*fac[1]<< " back = " << back*fac[1] << " sum-back = " << (sum-back)*fac[1]<< endl ;

   
   //
}
