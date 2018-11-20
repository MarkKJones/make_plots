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

void comp_peak_runs() {
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
     outputpdf="plots/comp_cointime_runs.pdf";
    static const Int_t nftot=4;
 TString inputroot[nftot];
   TFile *fhistroot[nftot];
   TString nfName[nftot]={"3420_200007","3543_200007","3923_200007","3968_200007"};
   TString nfTitle[nftot]={"run 3420","run 3543","run 3923","run 3968"};
   for (Int_t n=0;n<nftot;n++) {
   inputroot[n]="hist/coin_replay_coin_hElec_pProt_"+nfName[n]+"_CTime_hist.root";
     cout << " infile root = " << inputroot[n] << endl;
   fhistroot[n] =  new TFile(inputroot[n]);
   }
   static const Int_t nhist2d=4;
   static const Int_t nhist1d=5;
   TString hname2d[nhist2d]={"hepicoinTime_pdelta","hepicoinTime_hdelta","hepicoinTime_hxptar","hepicoinTime_pxptar"};
   TString hname1d[nhist1d]={"hepicoinTime_ROC2","hps1xfptime","hps1yfptime","hps2xfptime","hps2yfptime"};
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
   } 
  //
  TCanvas *can1d[nhist1d];
  Double_t ct_cent[nftot];
  Double_t ct_error[nftot];
   for (Int_t ip=0;ip<nhist1d;ip++) {
     can1d[ip]=new TCanvas(Form("can1d_%d",ip),Form("can1d_%d",ip),700,700);
     can1d[ip]->Divide(2,2);
   for (Int_t nf=0;nf<nftot;nf++) {
     can1d[ip]->cd(nf+1);
     if (ip==0) {
     hist1d[nf][ip]->GetXaxis()->SetRangeUser(0,20);
     hist1d[nf][ip]->Draw();
     Double_t xlo=hist1d[nf][ip]->GetBinCenter(hist1d[nf][ip]->GetMaximumBin())-.6;
     Double_t xhi=hist1d[nf][ip]->GetBinCenter(hist1d[nf][ip]->GetMaximumBin())+.6;
 	  hist1d[nf][ip]->Fit("gaus","","",xlo,xhi);
	  TF1 *fit = hist1d[nf][ip]->GetFunction("gaus");
	  ct_cent[nf] = fit->GetParameter(1);
	  ct_error[nf] = fit->GetParError(1);
     } else {
       gStyle->SetOptStat(1000111);
    hist1d[nf][ip]->Draw();
     }
  }    
   }
  //
}
