#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
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

void comp_coin_cosy_model() {
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
const UInt_t nruntot=4;
 TString frun[nruntot]={"3288","3371","3374","3377"};
 TString fang[nruntot]={"12.2 deg  P_e = 8.44","14 deg P_e = 7.94","10 deg P_e = 9.05","8.5 deg P_e =9.43"};
const UInt_t nver=4;
 TString vername[nver]={"orig","p15","p18","p20"};
 Int_t colind[nver]={1,2,3,5};
 TString lname[nver]={"orig","p15","p18","p20"};
const UInt_t nplots=4;
 //
  TString inputroot;
   TFile *fhistroot[nruntot][nver];
   TH2F *fhist[nruntot][nver][nplots];
   TH1F *fW[nruntot][nver];
   const char* plname[nplots]={"hDeltaDiffXfp","hDeltaDiffXpfp","hDeltaDiffYfp","hDeltaDiffYpfp"};
//
   for (UInt_t nf=0;nf<nruntot;nf++) {
   for (UInt_t nv=0;nv<nver;nv++) {
     inputroot="hist/coin_replay_coin_pElec_hProt_"+frun[nf]+"_quads_"+vername[nv]+"_coin_hist.root";
     cout << " infile root = " << inputroot << endl;
     fhistroot[nf][nv] =  new TFile(inputroot);
     for (UInt_t nh=0;nh<nplots;nh++) {
     fhist[nf][nv][nh] = (TH2F*)fhistroot[nf][nv]->Get(plname[nh]);
     }
     fW[nf][nv] = (TH1F*)fhistroot[nf][nv]->Get("hWcalc");
     if (!fW[nf][nv]) cout << " No W_calc " << endl;
   }
   }
   //
    outputpdf="plots/comp_coin_cosy_model_W.pdf";
    TCanvas *ctime[nruntot];
    TLegend *ltime[nruntot];
   for (UInt_t nr=0;nr<nruntot;nr++) {
     ctime[nr] = new TCanvas(Form("ctime_%d",nr),frun[nr], 1000,700);
     ctime[nr]->Divide(1,1);
     ctime[nr]->cd(1);
     //gPad->SetLogy();
     Double_t max=0;
      for (UInt_t nv=0;nv<nver;nv++) {
	if (fW[nr][nv]->GetMaximum()>max) max=fW[nr][nv]->GetMaximum();
      }
      ltime[nr] = new TLegend(.5,.6,.8,.8,"");
      for (UInt_t nv=0;nv<nver;nv++) {
	if (nv==0) fW[nr][nv]->Draw();
	if (nv!=0) fW[nr][nv]->Draw("same");
	fW[nr][nv]->SetLineColor(colind[nv]);
	fW[nr][nv]->SetMaximum(max*1.2);
	TString title="Run "+frun[nr]+"  Ang "+fang[nr]+" Fix Pcent = 8.7";
	fW[nr][nv]->SetTitle(title);
        TString tt=vername[nv];
	ltime[nr]->AddEntry(fW[nr][nv],tt);
      }
      ltime[nr]->Draw();
            if (nr==0) ctime[nr]->Print(outputpdf+"(");
            if (nr==(nruntot-1)) ctime[nr]->Print(outputpdf+")");
	    if (nr!=0&&nr!=(nruntot-1)) ctime[nr]->Print(outputpdf);
   }
   //
   //
    outputpdf="plots/comp_coin_cosy_model_xptar.pdf";
    TCanvas *ctime2[nver];
      for (UInt_t nv=0;nv<nver;nv++) {
     ctime2[nv] = new TCanvas(Form("ctime2_%d",nv),vername[nv], 1000,700);
     ctime2[nv]->Divide(2,2);
      for (UInt_t nr=0;nr<nruntot;nr++) {
         ctime2[nv]->cd(nr+1);
	fhist[nr][nv][1]->Draw("colz");
	TString title="Run "+frun[nr]+"  Optics version = "+vername[nv];
	fhist[nr][nv][1]->SetTitle(title);
      }
      if (nv==0) ctime2[nv]->Print(outputpdf+"(");
      if (nv!=0&&nv!=(nver-1)) ctime2[nv]->Print(outputpdf);
      if (nv==(nver-1)) ctime2[nv]->Print(outputpdf+")");
   }
}
