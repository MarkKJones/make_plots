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

void comp_hms_xfp_yfp_sieve(TString basename="") {
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
 outputpdf="plots/comp_hms_xfp_yfp_sieve_hist.pdf";
 //
 static const Int_t nftot=5;
 /* 
 TString flab[nftot]={"4787_-1","6868_-1","7947_-1","8376_-1","8388_-1"};
 TString tlab[nftot]={"Run 4787 Eb=10.6 p=5.322 th=13","Run 6868  Eb=4.9 p=3.0 th=13","Run 7947 Eb=8.2 p=5.745 th=13","Run 8376 Eb=2.7 p=1.5 th=11.5","Run 8388 Eb=2.7 p=0.4 th=11.5"};
 */
 TString flab[nftot]={"6868_-1","8376_-1","8397_-1","8395_-1","8388_-1"};
 TString tlab[nftot]={"Run 6868  Eb=4.9 p=3.0 th=13","Run 8376 Eb=2.7 p=1.5 th=11.5","Run 8397 Eb=2.7 p=1.0 th=11.5","Run 8395 Eb=2.7 p=0.8 th=11.5","Run 8388 Eb=2.7 p=0.4 th=11.5"};
 Int_t colind[nftot]={1,2,4,6,7};
TString inputroot[nftot];
   TFile *fhistroot[nftot];
    TH2F *fhist[nftot];
    TH2F *fsieve[nftot];
    TH1F *fxs[nftot];
    TH1F *fys[nftot];
 for (Int_t nf=0;nf<nftot;nf++) {
     inputroot[nf]="hist/hms_replay_matrixopt_"+flab[nf]+"_hms_ztar_sieve_hist.root";
     cout << " infile root = " << inputroot[nf] << endl;
   fhistroot[nf] =  new TFile(inputroot[nf]);
   fhist[nf] = (TH2F*)fhistroot[nf]->Get("hyfp_yxfp_foil_1");
   fhist[nf]->SetTitle(tlab[nf]);
   fhist[nf]->GetXaxis()->SetTitle("HMS Xfp (cm)");
   fhist[nf]->GetYaxis()->SetTitle("HMS Yfp (cm)");
    fsieve[nf] = (TH2F*)fhistroot[nf]->Get("hys_xs_foil_1");
   fsieve[nf]->SetTitle(tlab[nf]);
    fxs[nf] = (TH1F*)fhistroot[nf]->Get("hxs_foil_1");
   fxs[nf]->SetTitle(tlab[nf]);
    fys[nf] = (TH1F*)fhistroot[nf]->Get("hys_foil_1");
   fys[nf]->SetTitle(tlab[nf]);
} //
 //
       TCanvas *c1 = new TCanvas("c1","xfp_yfp", 1000,1000);
       c1->Divide(2,3);
  for (Int_t nf=0;nf<nftot;nf++) {
    c1->cd(nf+1);
gPad->SetGridx();
gPad->SetGridy();
    fhist[nf]->DrawClone("colz");
    
  }     
  c1->Print(outputpdf+"(");
 //
 //
       TCanvas *c2 = new TCanvas("c2","xfp_yfp_2", 1500,1000);
      c2->Divide(3,2);
  for (Int_t nf=0;nf<nftot;nf++) {
    c2->cd(nf+1);
gPad->SetGridx();
gPad->SetGridy();
 fhist[nf]->SetAxisRange(-5,5,"Y");
 fhist[nf]->SetAxisRange(-10,10,"X");
    fhist[nf]->Draw("colz");
  }     
  c2->Print(outputpdf);
  //
 //
       TCanvas *c4 = new TCanvas("c4","xs", 1500,1000);
       TLegend *lxs = new TLegend(.15,.75,.45,.95,"");
       c4->Divide(1,1);
    c4->cd(1);
  for (Int_t nf=0;nf<nftot;nf++) {
gPad->SetGridx();
gPad->SetGridy();
 fxs[nf]->SetLineColor(colind[nf]);
 if (nf==0)  fxs[nf]->DrawNormalized();
 if (nf>0)  fxs[nf]->DrawNormalized("same");
 lxs->AddEntry(fxs[nf],tlab[nf]);
  } 
  lxs->Draw();
  c4->Print(outputpdf);
  //
 //
       TCanvas *c5 = new TCanvas("c5","ys", 1500,1000);
        TLegend *lys = new TLegend(.15,.75,.45,.95,"");
      c5->Divide(1,1);
    c5->cd(1);
  for (Int_t nf=0;nf<nftot;nf++) {
gPad->SetGridx();
gPad->SetGridy();
 fys[nf]->SetLineColor(colind[nf]);
 if (nf==0)  fys[nf]->DrawNormalized();
 if (nf>0)  fys[nf]->DrawNormalized("same");
 lys->AddEntry(fys[nf],tlab[nf]);
  }     
  lys->Draw();
  c5->Print(outputpdf);
  //
       TCanvas *c3 = new TCanvas("c3","xs_ys", 1500,1000);
       c3->Divide(3,2);
  for (Int_t nf=0;nf<nftot;nf++) {
    c3->cd(nf+1);
gPad->SetGridx();
gPad->SetGridy();
    fsieve[nf]->DrawClone("colz");
    
  }     
  c3->Print(outputpdf+")");
 //
 //
}
