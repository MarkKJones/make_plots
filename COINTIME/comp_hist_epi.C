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
   inputroot[0]=Form("hist/coin_replay_coin_hElec_pProt_%d_%d_CTime_hist.root",nr1,nev);
   inputroot[1]=Form("hist/coin_replay_coin_hElec_pProt_%d_%d_CTime_hist.root",nr2,nev);
     cout << " infile root = " << inputroot[0] << endl;
     cout << " infile root = " << inputroot[1] << endl;
   fhistroot[0] =  new TFile(inputroot[0]);
   fhistroot[1] =  new TFile(inputroot[1]);
   static const Int_t nhist=1;
   //  TString hname[nhist]={"hepicoinTime_ROC2"};
   TString hname[nhist]={"hCTcalc_all"};
    TH1F *hist[2][nhist];
  for (Int_t nf=0;nf<2;nf++) {
  for (Int_t ip=0;ip<nhist;ip++) {
    // cout <<  hname[ip] << endl;
       hist[nf][ip] = (TH1F*)fhistroot[nf]->Get(hname[ip]);
       if (!hist[nf][ip]) cout << " no hist = " << hname[ip] << endl;
  }
  }
  //
}
