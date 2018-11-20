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

void comp_peak_5371() {
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
     outputpdf="plots/comp_cointime_5371.pdf";
    static const Int_t nftot=4;
 TString inputroot[nftot];
   TFile *fhistroot[nftot];
   TString nfName[nftot]={"5371_200007","5371_200008","5369_200007","5369_200008"};
   TString nfTitle[nftot]={"5371  new","5371 old","5369  new","5369 old"};
   for (Int_t n=0;n<nftot;n++) {
   inputroot[n]="hist/coin_replay_production_"+nfName[n]+"_CTime_hist.root";
     cout << " infile root = " << inputroot[n] << endl;
   fhistroot[n] =  new TFile(inputroot[n]);
   }
   static const Int_t nhist2d=4;
   static const Int_t nhist1d=3;
   TString hname2d[nhist2d]={"hepicoinTime_pdelta","hepicoinTime_hdelta","hepicoinTime_hxptar","hepicoinTime_pxptar"};
   TString hname1d[nhist1d]={"hepicoinTime_ROC2","hpbetatrack","hhbetatrack"};
    TH2F *hist2d[nftot][nhist2d];
    TH2F *hist1d[nftot][nhist1d];
  for (Int_t nf=0;nf<nftot;nf++) {
  for (Int_t ip=0;ip<nhist2d;ip++) {
    // cout <<  hname[ip] << endl;
       hist2d[nf][ip] = (TH2F*)fhistroot[nf]->Get(hname2d[ip]);
       hist2d[nf][ip]->SetTitle(nfTitle[nf]);
       if (!hist2d[nf][ip]) cout << " no hist = " << hname2d[ip] << endl;
  }
  for (Int_t ip=0;ip<nhist1d;ip++) {
    // cout <<  hname[ip] << endl;
       hist1d[nf][ip] = (TH2F*)fhistroot[nf]->Get(hname1d[ip]);
        hist1d[nf][ip]->SetTitle(nfTitle[nf]);
      if (!hist1d[nf][ip]) cout << " no hist = " << hname1d[ip] << endl;
  }
  }
  //
  TCanvas *can2d[nhist2d];
   for (Int_t ip=0;ip<nhist2d;ip++) {
     can2d[ip]=new TCanvas(Form("can2d_%d",ip),Form("can2d_%d",ip),700,700);
     can2d[ip]->Divide(2,2);
   for (Int_t nf=0;nf<nftot;nf++) {
     can2d[ip]->cd(nf+1);
     hist2d[nf][ip]->Draw("colz");
   }    
   if (ip==0) can2d[ip]->Print(outputpdf+"(");
   if (ip!=0) can2d[ip]->Print(outputpdf);
   } 
  //
  TCanvas *can1d[nhist1d];
  Double_t ct_cent[nftot][nhist1d];
  Double_t ct_error[nftot][nhist1d];
  Double_t lim[nhist1d]={1.,.1,.1};
   for (Int_t ip=0;ip<nhist1d;ip++) {
     can1d[ip]=new TCanvas(Form("can1d_%d",ip),Form("can1d_%d",ip),700,700);
     can1d[ip]->Divide(2,2);
   for (Int_t nf=0;nf<nftot;nf++) {
     can1d[ip]->cd(nf+1);
     if (ip==0) hist1d[nf][ip]->GetXaxis()->SetRangeUser(30,60);
     if (ip==1) hist1d[nf][ip]->GetXaxis()->SetRangeUser(.7,1.4);
     if (ip==2) hist1d[nf][ip]->GetXaxis()->SetRangeUser(.7,1.4);
     hist1d[nf][ip]->Draw();
     Double_t xlo=hist1d[nf][ip]->GetBinCenter(hist1d[nf][ip]->GetMaximumBin())-lim[ip];
     Double_t xhi=hist1d[nf][ip]->GetBinCenter(hist1d[nf][ip]->GetMaximumBin())+lim[ip];
 	  hist1d[nf][ip]->Fit("gaus","","",xlo,xhi);
	  TF1 *fit = hist1d[nf][ip]->GetFunction("gaus");
	  ct_cent[nf][ip] = fit->GetParameter(1);
	  ct_error[nf][ip] = fit->GetParError(1);
  }    
   if (ip!=nhist1d-1) can1d[ip]->Print(outputpdf); 
   if (ip==nhist1d-1) can1d[ip]->Print(outputpdf+")"); 
   }
  //
}
