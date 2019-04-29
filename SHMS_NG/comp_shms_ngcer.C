#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
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
void comp_shms_ngcer( TString fn1 ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/comp_shms_ngcer_"+fn1+".pdf";
    //
   TFile *fhistroot;
 //
 //
 TString inputroot;
 inputroot = "hist/"+fn1+"_shms_ngcer_hist.root";
 cout << " root file = " << inputroot << endl;
  //
 //
   TH1F *h_adcamp[4];
   TH1F *h_adcint[4];
   TH1F *h_adcamp_cut[4];
   TH1F *h_adcint_cut[4];
   TH1F *h_adcamp_cutlow[4];
   TH1F *h_adcint_cutlow[4];
   TH1F *h_adcamp_cutlow2[4];
   TH1F *h_adcint_cutlow2[4];
   TH2F *h_cerXY_adcamp[4];
   TH2F *h_cerXY_adcamp_now[4];
   TH2F *h_cerXY_adcamp_rat[4];
   TH2F *h_cerXY_adcint[4];
   TH2F *h_cerXY_adcint_now[4];
   TH2F *h_cerXY_adcint_rat[4];
  //
   fhistroot =  new TFile(inputroot);
   for (Int_t n=0; n<4;n++) {
     h_adcamp[n] = (TH1F*)fhistroot->Get(Form("h_adcamp_%d",n));
     h_adcint[n] = (TH1F*)fhistroot->Get(Form("h_adcint_%d",n));
      h_adcamp_cut[n] = (TH1F*)fhistroot->Get(Form("h_adcamp_cut_%d",n));
     h_adcint_cut[n] = (TH1F*)fhistroot->Get(Form("h_adcint_cut_%d",n));
      h_adcamp_cutlow[n] = (TH1F*)fhistroot->Get(Form("h_adcamp_cutlow_%d",n));
     h_adcint_cutlow[n] = (TH1F*)fhistroot->Get(Form("h_adcint_cutlow_%d",n));
      h_adcamp_cutlow2[n] = (TH1F*)fhistroot->Get(Form("h_adcamp_cutlow2_%d",n));
     h_adcint_cutlow2[n] = (TH1F*)fhistroot->Get(Form("h_adcint_cutlow2_%d",n));
    h_cerXY_adcamp[n] = (TH2F*)fhistroot->Get(Form("hcerXY_adcamp_%d",n));
     h_cerXY_adcamp_now[n] = (TH2F*)fhistroot->Get(Form("hcerXY_adcamp_now_%d",n));
     h_cerXY_adcamp_rat[n] = (TH2F*)h_cerXY_adcamp_now[n]->Clone();
     h_cerXY_adcamp_rat[n]->Divide(h_cerXY_adcamp[n],h_cerXY_adcamp_now[n],1,1,"B");
   } 
   //
   TCanvas *c2d;
      c2d= new TCanvas("c2d","2damp", 1000,700);
     c2d->Divide(2,2);

    for (Int_t n=0; n<4;n++) {
      c2d->cd(n+1);
      h_cerXY_adcamp_rat[n]->Draw("colz");
    }  
    c2d->Print(outputpdf+"(");
   //
    TF1* fadcamplow[4];
    Double_t adcamplow_peak[4];
    Double_t adcamplow_sigma[4];
   TCanvas *campcutlow;
      campcutlow= new TCanvas("campcutlow","cutlo", 1000,700);
     campcutlow->Divide(2,2);

    for (Int_t n=0; n<4;n++) {
      Double_t WMax = h_adcamp_cutlow[n]->GetBinCenter(h_adcamp_cutlow[n]->GetMaximumBin());
      fadcamplow[n] = new TF1(Form("fadcamplow_%d",n),"gaus",h_adcamp_cutlow[n]->GetMean()-h_adcamp_cutlow[n]->GetRMS(),h_adcamp_cutlow[n]->GetMean()+h_adcamp_cutlow[n]->GetRMS());
      campcutlow->cd(n+1);
      //h_adcamplow_cut[n]->Draw("");
      h_adcamp_cutlow[n]->Fit(Form("fadcamplow_%d",n),"QR");
      adcamplow_peak[n]=fadcamplow[n]->GetParameter(1);
      adcamplow_sigma[n]=fadcamplow[n]->GetParameter(2);
      cout << " Pmt = " << n << " peak = " << adcamplow_peak[n] << " sigma = " << adcamplow_sigma[n] 
       << endl;
     }  
    campcutlow->Print(outputpdf);
   //
   //
    TF1* fadcintlow[4];
    Double_t adcintlow_peak[4];
    Double_t adcintlow_sigma[4];
   TCanvas *cintcutlow;
      cintcutlow= new TCanvas("cintcutlow","cutintlo", 1000,700);
     cintcutlow->Divide(2,2);

    for (Int_t n=0; n<4;n++) {
      Double_t WMax = h_adcint_cutlow[n]->GetBinCenter(h_adcint_cutlow[n]->GetMaximumBin());
      fadcintlow[n] = new TF1(Form("fadcintlow_%d",n),"gaus",h_adcint_cutlow[n]->GetMean()-h_adcint_cutlow[n]->GetRMS(),h_adcint_cutlow[n]->GetMean()+h_adcint_cutlow[n]->GetRMS());
      cintcutlow->cd(n+1);
      //h_adcintlow_cut[n]->Draw("");
      h_adcint_cutlow[n]->Fit(Form("fadcintlow_%d",n),"QR");
      adcintlow_peak[n]=fadcintlow[n]->GetParameter(1);
      adcintlow_sigma[n]=fadcintlow[n]->GetParameter(2);
      cout << " Pmt = " << n << " peak = " << adcintlow_peak[n] << " sigma = " << adcintlow_sigma[n] 
       << endl;
     }  
    cintcutlow->Print(outputpdf);
   //
   //
    TF1* fadcamplow2[4];
    Double_t adcamplow2_peak[4];
    Double_t adcamplow2_sigma[4];
   TCanvas *cintcutlow2;
      cintcutlow2= new TCanvas("cintcutlow2","cutintlow2", 1000,700);
     cintcutlow2->Divide(2,2);

    for (Int_t n=0; n<4;n++) {
      Double_t WMax = h_adcamp_cutlow2[n]->GetBinCenter(h_adcamp_cutlow2[n]->GetMaximumBin());
      fadcamplow2[n] = new TF1(Form("fadcamplow2_%d",n),"gaus",h_adcamp_cutlow2[n]->GetMean()-h_adcamp_cutlow2[n]->GetRMS(),h_adcamp_cutlow2[n]->GetMean()+h_adcamp_cutlow2[n]->GetRMS());
      cintcutlow2->cd(n+1);
      //h_adcamplow2_cut[n]->Draw("");
      h_adcamp_cutlow2[n]->Fit(Form("fadcamplow2_%d",n),"QR");
      adcamplow2_peak[n]=fadcamplow2[n]->GetParameter(1);
      adcamplow2_sigma[n]=fadcamplow2[n]->GetParameter(2);
      cout << " Pmt = " << n << " peak = " << adcamplow2_peak[n] << " sigma = " << adcamplow2_sigma[n] 
       << endl;
     }  
    cintcutlow2->Print(outputpdf);
   //
   //
    TF1* fadcint[4];
    Double_t adcint_peak[4];
    Double_t adcint_sigma[4];
   TCanvas *cintcut;
      cintcut= new TCanvas("cintcut","1dint", 1000,700);
     cintcut->Divide(2,2);

    for (Int_t n=0; n<4;n++) {
      Double_t WMax = h_adcint_cut[n]->GetBinCenter(h_adcint_cut[n]->GetMaximumBin());
      fadcint[n] = new TF1(Form("fadcint_%d",n),"gaus",h_adcint_cut[n]->GetMean()-h_adcint_cut[n]->GetRMS(),h_adcint_cut[n]->GetMean()+h_adcint_cut[n]->GetRMS());
      cintcut->cd(n+1);
      //h_adcint_cut[n]->Draw("");
      h_adcint_cut[n]->Fit(Form("fadcint_%d",n),"QR");
      adcint_peak[n]=fadcint[n]->GetParameter(1);
      adcint_sigma[n]=fadcint[n]->GetParameter(2);
      Double_t npe = adcint_peak[n]/adcint_sigma[n]*adcint_peak[n]/adcint_sigma[n];
      cout << " Pmt = " << n << " peak = " << adcint_peak[n] << " sigma = " << adcint_sigma[n] <<
	" npe = " << npe<< " conversion = " << npe/adcint_peak[n] << endl;
     }  
    cintcut->Print(outputpdf);
   //
   //
    TF1* fadcamp[4];
    Double_t adcamp_peak[4];
    Double_t adcamp_sigma[4];
   TCanvas *campcut;
      campcut= new TCanvas("campcut","1damp", 1000,700);
     campcut->Divide(2,2);

    for (Int_t n=0; n<4;n++) {
      Double_t WMax = h_adcamp_cut[n]->GetBinCenter(h_adcamp_cut[n]->GetMaximumBin());
      fadcamp[n] = new TF1(Form("fadcamp_%d",n),"gaus",h_adcamp_cut[n]->GetMean()-h_adcamp_cut[n]->GetRMS(),h_adcamp_cut[n]->GetMean()+h_adcamp_cut[n]->GetRMS());
      campcut->cd(n+1);
      //h_adcamp_cut[n]->Draw("");
      h_adcamp_cut[n]->Fit(Form("fadcamp_%d",n),"QR");
      adcamp_peak[n]=fadcamp[n]->GetParameter(1);
      adcamp_sigma[n]=fadcamp[n]->GetParameter(2);
      Double_t npe = adcamp_peak[n]/adcamp_sigma[n]*adcamp_peak[n]/adcamp_sigma[n];
      cout << " Pmt = " << n << " peak = " << adcamp_peak[n] << " sigma = " << adcamp_sigma[n] <<
	" npe = " << npe<< " conversion = " << npe/adcamp_peak[n] << endl;
     }  
    campcut->Print(outputpdf+")");
   //

}
