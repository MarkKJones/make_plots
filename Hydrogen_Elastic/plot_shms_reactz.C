#include <TString.h>
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
void plot_shms_reactz( TString fn1 ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/plot_shms_reactz.pdf";
    //
const UInt_t nftot=1;
 Int_t colind[nftot]={1};
   TFile *fhistroot[nftot];
 //
const UInt_t nEplots=3;
 TString hEname[nEplots]={"h_ereactz","h_ereactx","h_ereacty"};
const UInt_t nHplots=3;
 TString hHname[nHplots]={"h_hreactz","h_hreactx","h_hreacty"};
const UInt_t nRplots=2;
 TString hRname[nRplots]={"h_eRasterx","h_eRastery"};
const UInt_t nEpicsplots=2;
 TString hEpicsname[nEpicsplots]={"h_eEpicsx","h_eEpicsy"};
 //
 TString inputroot[2];
 inputroot[0] = "hist/"+fn1+"_shms_ep_elastic_hist.root";
  TH1F *fEhist[nEplots];
  TH1F *fRhist[nRplots];
  TH1F *fEpicshist[nEpicsplots];
 TH1F *fHhist[nHplots];
 TH1F *h_diffreactz;
  //
   for (UInt_t nf=0;nf<nftot;nf++) {
     cout << " infile root = " << inputroot[nf] << endl;
   fhistroot[nf] =  new TFile(inputroot[nf]);
    h_diffreactz = (TH1F*)fhistroot[nf]->Get("h_diffreactz");
  for (UInt_t nh=0;nh<nEpicsplots;nh++) {
     fEpicshist[nh] = (TH1F*)fhistroot[nf]->Get(hEpicsname[nh]);
   }
   for (UInt_t nh=0;nh<nRplots;nh++) {
     fRhist[nh] = (TH1F*)fhistroot[nf]->Get(hRname[nh]);
   }
   for (UInt_t nh=0;nh<nEplots;nh++) {
     fEhist[nh] = (TH1F*)fhistroot[nf]->Get(hEname[nh]);
   }
   for (UInt_t nh=0;nh<nEplots;nh++) {
     fEhist[nh] = (TH1F*)fhistroot[nf]->Get(hEname[nh]);
   }
  for (UInt_t nh=0;nh<nHplots;nh++) {
     fHhist[nh] = (TH1F*)fhistroot[nf]->Get(hHname[nh]);
   }
   }
   // W histograms
    TCanvas *creact;
     creact = new TCanvas("creact","creact", 1000,700);
     creact->Divide(1,3);
     Int_t ncd=1;
     TLegend *lW[3];
     TString laxis[3]={"; React_z (cm) ","; React_x (cm) ","; React_y (cm) "};
     TString tt;
    for (UInt_t nh=0;nh<nEplots;nh++) {
      lW[nh]= new TLegend(.8,.75,.99,.95,"");
      creact->cd(ncd++);
      fEhist[nh]->Draw();
      fEhist[nh]->SetTitle(laxis[nh]);
       fHhist[nh]->Draw("same");
       fHhist[nh]->SetLineColor(2);
       Double_t mean=fEhist[nh]->GetMean();
        tt=Form("SHMS Mean = %f",mean);
  	lW[nh]->AddEntry(fEhist[nh],tt);
        mean=fHhist[nh]->GetMean();
        tt=Form("HMS Mean = %f",mean);
  	lW[nh]->AddEntry(fHhist[nh],tt);
        lW[nh]->Draw();
   }
      creact->Print(outputpdf+"(");
      //
    TCanvas *cplot;
 gStyle->SetOptStat(11111);
     cplot = new TCanvas("cplot","cplot", 1000,700);
     cplot->Divide(2,3);
     ncd=1;
     cplot->cd(ncd++);
     h_diffreactz->Draw();
     ncd++;
     for (UInt_t nh=0;nh<nRplots;nh++) {
    cplot->cd(ncd++);
      fRhist[nh]->Draw();
    }
    for (UInt_t nh=0;nh<nEpicsplots;nh++) {
    cplot->cd(ncd++);
      fEpicshist[nh]->Draw();
    }
      cplot->Print(outputpdf+")");
}
	   //
