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
void fit_shms_preshower( TString fn1 ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/comp_shms_preshower_fit_"+fn1+".pdf";
    //
   TFile *fhistroot;
 //
 //
 TString inputroot;
 inputroot = "hist/"+fn1+"_shms_ngcer_hist.root";
 cout << " root file = " << inputroot << endl;
  //
 //
  //
   fhistroot =  new TFile(inputroot);
    TH1F *h_preshow_energy_pion[2][14];
   for (Int_t ns=0;ns<2;ns++) {
   for (Int_t nr=0;nr<14;nr++) {
    h_preshow_energy_pion[ns][nr] = (TH1F*)fhistroot->Get(Form("h_preshow_energy_pion_%d_%d",ns,nr));
   }}
   //
    TF1* fene[2][14];
    Double_t ene_peak[2][14];
    Double_t ene_sigma[2][14];
   TCanvas *cNegEnergy;
   Int_t ns=0;
      cNegEnergy= new TCanvas("cNegEnergy","cNegEnergy", 1000,700);
     cNegEnergy->Divide(2,7);
    for (Int_t nr=0; nr<14;nr++) {
      //    h_preshow_energy_pion[ns][nr]->Draw();
      fene[ns][nr] = new TF1(Form("fene_%d_%d",ns,nr),"gaus",.05,.15);
      cNegEnergy->cd(nr+1);
      
      h_preshow_energy_pion[ns][nr]->Fit(Form("fene_%d_%d",ns,nr),"QR");
      ene_peak[ns][nr]=fene[ns][nr]->GetParameter(1);
      ene_sigma[ns][nr]=fene[ns][nr]->GetParameter(2);
      cout << "Neg  Pmt = " << nr << " peak = " << ene_peak[ns][nr] << " sigma = " <<  ene_sigma[ns][nr]
       << endl;
    }
   //

}
