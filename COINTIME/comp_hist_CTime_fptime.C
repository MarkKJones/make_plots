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

void comp_hist_CTime_fptime(Int_t nr1,Int_t nr2,Int_t nev=100000) {
 gROOT->Reset();
 gStyle->SetOptStat(1000011);
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
   inputroot[0]=Form("hist/coin_replay_coin_hElec_pProt_%d_200000_CTime_hist.root",nr1);
   inputroot[1]=Form("hist/coin_replay_production_%d_500000_CTime_hist.root",nr2);
     cout << " infile root = " << inputroot[0] << endl;
     cout << " infile root = " << inputroot[1] << endl;
   fhistroot[0] =  new TFile(inputroot[0]);
   fhistroot[1] =  new TFile(inputroot[1]);
   static const Int_t nhist=2;
   static const Int_t nhist1d=4;
   //  TString hname[nhist]={"hepicoinTime_ROC2"};
   TString hname[nhist]={"hps1xfptime_trig14diffraw","hs1xdiff_ptrigdiff_ROC1_cut"};
   TString hname1d[nhist1d]={"hptrig1_ROC2","hptrig4_ROC2","hptrigd14iffroc2","hptrigd14iff"};
    TH2F *hist[2][nhist];
    TH2F *hist1d[2][nhist1d];
  for (Int_t nf=0;nf<2;nf++) {
  for (Int_t ip=0;ip<nhist;ip++) {
    // cout <<  hname[ip] << endl;
       hist[nf][ip] = (TH2F*)fhistroot[nf]->Get(hname[ip]);
       if (!hist[nf][ip]) cout << " no hist = " << hname[ip] << endl;
  }
  for (Int_t ip=0;ip<nhist1d;ip++) {
    // cout <<  hname[ip] << endl;
       hist1d[nf][ip] = (TH2F*)fhistroot[nf]->Get(hname1d[ip]);
       if (!hist1d[nf][ip]) cout << " no hist = " << hname1d[ip] << endl;
  }
  }
  //
  TCanvas *can2d[nhist];
  for (Int_t ip=0;ip<nhist;ip++) {
  can2d[ip] = new TCanvas(Form("can2d_%d",ip),Form("can2d_%d",ip),700,700);
  can2d[ip]->Divide(1,2);
  can2d[ip]->cd(1);
  hist[0][ip]->Draw("colz");
  can2d[ip]->cd(2);
  hist[1][ip]->Draw("colz");
  if (ip==0) {
    can2d[ip]->Print(outputpdf+"(");
  } else {
    can2d[ip]->Print(outputpdf);
  }
  }
  //
  TCanvas *can1d[nhist1d];
  for (Int_t ip=0;ip<nhist1d;ip++) {
    can1d[ip] = new TCanvas(Form("can1d_%d",ip),Form("can1d_%d",ip),700,700);
  can1d[ip]->Divide(1,2);
  can1d[ip]->cd(1);
  hist1d[0][ip]->Draw();
  can1d[ip]->cd(2);
  hist1d[1][ip]->Draw();
   if (ip==nhist1d-1) {
    can1d[ip]->Print(outputpdf+")");
  } else {
    can1d[ip]->Print(outputpdf);
  }
  } 
   
   //
}
