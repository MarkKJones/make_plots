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

void comp_coin_orig_newfit() {
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
 TString frun[nruntot]={"1711","1713","1718","1719"};
 TString fang[nruntot]={"25.01 deg P = 1.7851","30. deg P = 1.7851","30. deg P = 1.6573","30. deg P = 1.5285"};
const UInt_t nver=2;
 TString vername[nver]={"p18","newfit2"};
 Int_t colind[nver]={1,2};
 TString lname[nver]={"p18","newfit"};
const UInt_t nplots=4;
 //
  TString inputroot;
   TFile *fhistroot[nruntot][nver];
   TH2F *fhist[nruntot][nver][nplots];
   TH2F *fEmfp[nruntot][nver][nplots];
   TH2F *fDeltaDiff_fp[nruntot][nver][nplots];
   TProfile *fprof[nruntot][nver][nplots];
   TH1F *fW[nruntot][nver];
   TH1F *fEm[nruntot][nver];
   const char* plname[nplots]={"hWXfp","hWXpfp","hWYfp","hWYpfp"};
   const char* pl2name[nplots]={"hEmissXfp","hEmissXpfp","hEmissYfp","hEmissYpfp"};
   const char* pl3name[nplots]={"hDeltaDiffXfp","hDeltaDiffXpfp","hDeltaDiffYfp","hDeltaDiffYpfp"};
//
   for (UInt_t nf=0;nf<nruntot;nf++) {
   for (UInt_t nv=0;nv<nver;nv++) {
     inputroot="hist/coin_replay_coin_pElec_hProt_"+frun[nf]+"_500000_"+vername[nv]+"_coin_hist.root";
     cout << " infile root = " << inputroot << endl;
     fhistroot[nf][nv] =  new TFile(inputroot);
     for (UInt_t nh=0;nh<nplots;nh++) {
     fhist[nf][nv][nh] = (TH2F*)fhistroot[nf][nv]->Get(plname[nh]);
     fEmfp[nf][nv][nh] = (TH2F*)fhistroot[nf][nv]->Get(pl2name[nh]);
     fEmfp[nf][nv][nh]->SetName(Form("%s_%d_%d_%d",plname[nh],nf,nv,nh));
     fDeltaDiff_fp[nf][nv][nh] = (TH2F*)fhistroot[nf][nv]->Get(pl3name[nh]);
     fprof[nf][nv][nh] = fhist[nf][nv][nh]->ProfileY(Form("%s_%d_%d_%d_py",plname[nh],nf,nv,nh),1,-1,"s");
     if (!fhist[nf][nv][nh]) cout<< nf  << " No 2d " << endl;
     if (!fEmfp[nf][nv][nh]) cout<< nf  << " No Emfp " << endl;
     if (!fDeltaDiff_fp[nf][nv][nh]) cout<< nf  << " No DeltaDiff_fp " << endl;
     }
     fW[nf][nv] = (TH1F*)fhistroot[nf][nv]->Get("hW");
     if (!fW[nf][nv]) cout << " No W" << endl;
     fEm[nf][nv] = (TH1F*)fhistroot[nf][nv]->Get("hEmiss");
     if (!fEm[nf][nv]) cout << " No Emiss " << endl;
   }
   }
   //
    outputpdf="plots/comp_coin_orig_newfit_W.pdf";
    TCanvas *ctime[nruntot];
    TLegend *ltime[nruntot];
   for (UInt_t nr=0;nr<nruntot;nr++) {
     ctime[nr] = new TCanvas(Form("ctime_%d",nr),frun[nr], 1000,700);
     ctime[nr]->Divide(1,2);
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
	fW[nr][nv]->SetAxisRange(.9,1.05,"X");
	TString title="Run "+frun[nr]+"  Ang "+fang[nr];
	fW[nr][nv]->SetTitle(title);
        TString tt=vername[nv];
	ltime[nr]->AddEntry(fW[nr][nv],tt);
      }
      ltime[nr]->Draw();
     ctime[nr]->cd(1);
     //gPad->SetLogy();
     ctime[nr]->cd(2);
     max=0;
      for (UInt_t nv=0;nv<nver;nv++) {
	if (fEm[nr][nv]->GetMaximum()>max) max=fEm[nr][nv]->GetMaximum();
      }
      for (UInt_t nv=0;nv<nver;nv++) {
	if (nv==0) fEm[nr][nv]->Draw();
	if (nv!=0) fEm[nr][nv]->Draw("same");
	fEm[nr][nv]->SetLineColor(colind[nv]);
	fEm[nr][nv]->SetMaximum(max*1.2);
	fEm[nr][nv]->SetAxisRange(-.05,0.15,"X");
	TString title="Run "+frun[nr]+"  Ang "+fang[nr];
	fEm[nr][nv]->SetTitle(title);
      }
            if (nr==0) ctime[nr]->Print(outputpdf+"(");
            if (nr==(nruntot-1)) ctime[nr]->Print(outputpdf+")");
	    if (nr!=0&&nr!=(nruntot-1)) ctime[nr]->Print(outputpdf);
   }
   //
   //
    outputpdf="plots/comp_coin_orig_newfit_2d.pdf";
    TCanvas *ctime2[nver][nplots];
      for (UInt_t nv=0;nv<nver;nv++) {
      for (UInt_t np=0;np<nplots;np++) {
	ctime2[nv][np] = new TCanvas(Form("ctime2_%d_%d",nv,np),vername[nv]+plname[np], 1000,700);
     ctime2[nv][np]->Divide(2,2);
      for (UInt_t nr=0;nr<nruntot;nr++) {
         ctime2[nv][np]->cd(nr+1);
	fhist[nr][nv][np]->Draw("colz");
	TString title="Run "+frun[nr]+"  Replay version = "+vername[nv];
	fhist[nr][nv][np]->SetTitle(title);
      }
      
          if (nv==0&&np==0) { 
          ctime2[nv][np]->Print(outputpdf+"(");
	  } else {
            if (nv==(nver-1)&&np==(nplots-1)) { 
              ctime2[nv][np]->Print(outputpdf+")");
            }else {
              ctime2[nv][np]->Print(outputpdf);
            }
	  }
      }
      }
   //
   //
    outputpdf="plots/comp_coin_orig_newfit_emiss2d.pdf";
    TCanvas *ctime4[nver][nplots];
      for (UInt_t nv=0;nv<nver;nv++) {
      for (UInt_t np=0;np<nplots;np++) {
	ctime4[nv][np] = new TCanvas(Form("ctime4_%d_%d",nv,np),vername[nv]+pl2name[np], 1000,700);
     ctime4[nv][np]->Divide(2,2);
      for (UInt_t nr=0;nr<nruntot;nr++) {
         ctime4[nv][np]->cd(nr+1);
	fEmfp[nr][nv][np]->Draw("colz");
	TString title="Run "+frun[nr]+"  Replay version = "+vername[nv];
	fEmfp[nr][nv][np]->SetTitle(title);
      }
      
          if (nv==0&&np==0) { 
          ctime4[nv][np]->Print(outputpdf+"(");
	  } else {
            if (nv==(nver-1)&&np==(nplots-1)) { 
              ctime4[nv][np]->Print(outputpdf+")");
            }else {
              ctime4[nv][np]->Print(outputpdf);
            }
	  }
      }
      }
   //
   //
    outputpdf="plots/comp_coin_orig_newfit_deltadiff.pdf";
    TCanvas *ctime5[nver][nplots];
      for (UInt_t nv=0;nv<nver;nv++) {
      for (UInt_t np=0;np<nplots;np++) {
	ctime5[nv][np] = new TCanvas(Form("ctime5_%d_%d",nv,np),vername[nv]+pl3name[np], 1000,700);
     ctime5[nv][np]->Divide(2,2);
      for (UInt_t nr=0;nr<nruntot;nr++) {
         ctime5[nv][np]->cd(nr+1);
	fDeltaDiff_fp[nr][nv][np]->Draw("colz");
	TString title="Run "+frun[nr]+"  Replay version = "+vername[nv];
	fDeltaDiff_fp[nr][nv][np]->SetTitle(title);
      }
      
          if (nv==0&&np==0) { 
          ctime5[nv][np]->Print(outputpdf+"(");
	  } else {
            if (nv==(nver-1)&&np==(nplots-1)) { 
              ctime5[nv][np]->Print(outputpdf+")");
            }else {
              ctime5[nv][np]->Print(outputpdf);
            }
	  }
      }
      }
   //
    outputpdf="plots/comp_coin_orig_newfit_prof.pdf";
    TCanvas *ctime3[nver][nplots];
      for (UInt_t nv=0;nv<nver;nv++) {
      for (UInt_t np=0;np<nplots;np++) {
	ctime3[nv][np] = new TCanvas(Form("ctime3_%d_%d",nv,np),vername[nv]+plname[np], 1000,700);
     ctime3[nv][np]->Divide(2,2);
      for (UInt_t nr=0;nr<nruntot;nr++) {
         ctime3[nv][np]->cd(nr+1);
	 fEmfp[nr][nv][np]->FitSlicesX(0,0,-1,20);
         TH1D *hfit = (TH1D*)gDirectory->Get(Form("%s_%d_%d_%d_1",plname[np],nr,nv,np));
	TString title="Run "+frun[nr]+"  Replay version = "+vername[nv];
	hfit->SetTitle(title);
	hfit->SetMinimum(-.02);
	hfit->SetMaximum(+.02);
	hfit->Draw();
      }
      
          if (nv==0&&np==0) { 
          ctime3[nv][np]->Print(outputpdf+"(");
	  } else {
            if (nv==(nver-1)&&np==(nplots-1)) { 
              ctime3[nv][np]->Print(outputpdf+")");
            }else {
              ctime3[nv][np]->Print(outputpdf);
            }
	  }
      }
      }
   //
}
