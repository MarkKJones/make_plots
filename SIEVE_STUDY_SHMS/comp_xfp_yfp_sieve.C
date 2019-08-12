#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
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

void comp_xfp_yfp_sieve(TString basename="") {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
 //
 outputpdf="plots/"+basename+"_comp_shms_xfp_yfp_sieve_hist.pdf";
 //
 static const Int_t nftot=5;
 TString flab[nftot]={"4787_-1","6620_-1","6868_-1","7947_-1","8376_-1"};
 TString tlab[nftot]={"Run 4787 Eb=10.6 p=6.3 th=9.5","Run 6620 Eb=3.8 p=3.0 th=9.5","Run 6868 Eb=4.9 p=3.5 th=9.5","Run 7947 Eb=8.2 p=5.745 th=9.5","Run 8376 Eb=2.7 p=1.5 th=9.5"};
 TString inputroot[nftot];
   TFile *fhistroot[nftot];
    TH2F *fhist[nftot];
 for (Int_t nf=0;nf<nftot;nf++) {
     inputroot[nf]="hist/shms_replay_matrixopt_"+flab[nf]+"_shms_ztar_sieve_hist.root";
     cout << " infile root = " << inputroot[nf] << endl;
   fhistroot[nf] =  new TFile(inputroot[nf]);
   fhist[nf] = (TH2F*)fhistroot[nf]->Get("hyfp_yxfp_cent_foil");
   fhist[nf]->SetTitle(tlab[nf]);
   fhist[nf]->GetXaxis()->SetTitle("SHMS Xfp (cm)");
   fhist[nf]->GetYaxis()->SetTitle("SHMS Yfp (cm)");
 } //
 //
       TCanvas *c1 = new TCanvas("c1","xs_ys", 1000,1200);
      c1->Divide(2,3);
  for (Int_t nf=0;nf<nftot;nf++) {
    c1->cd(nf+1);
gPad->SetGridx();
gPad->SetGridy();
    fhist[nf]->DrawClone("colz");
  }     
  c1->Print(outputpdf+"(");
 //
       TCanvas *c2 = new TCanvas("c2","xs_ys_2", 1000,1200);
      c2->Divide(2,3);
  for (Int_t nf=0;nf<nftot;nf++) {
    c2->cd(nf+1);
gPad->SetGridx();
gPad->SetGridy();
 fhist[nf]->SetAxisRange(-5,5,"Y");
 fhist[nf]->SetAxisRange(-10,10,"X");
    fhist[nf]->Draw("colz");
  }     
   c2->Print(outputpdf+")");
//
}
