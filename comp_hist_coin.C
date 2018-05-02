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

void comp_hist_coin(TString basename, TString basename1) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
 outputpdf="plots/"+basename+"_coin.pdf";
  TString inputroot;
 static const Int_t nftot=2;
   TFile *fhistroot[nftot];
     inputroot="hist/"+basename+"_coin_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[0] =  new TFile(inputroot);
     inputroot="hist/"+basename1+"_coin_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[1] =  new TFile(inputroot);
   //
  static const Int_t hnum=1;
  TString hname[hnum]={"hW"};
   TH1F *hist[2][1];
  for (Int_t ifn=0;ifn<2;ifn++) {
       hist[ifn][0]=(TH1F*)fhistroot[ifn]->Get(hname[0]);
  }
   // 3244,3247
  Double_t scal[nftot]={26.2*0.6344*0.86,36.85*0.831*.86};
       for (Int_t nft=0;nft<2;nft++) {
	 hist[nft][0]->Scale(1./scal[nft]);
       }
      TCanvas *cplot[hnum];
     for (Int_t nh=0;nh<hnum;nh++) {
       cplot[nh] = new TCanvas(Form("cplot_%d",nh),hname[nh], 700,700);
       cplot[nh]->Divide(1,2);
       cplot[nh]->cd(1);
       for (Int_t nft=0;nft<nftot;nft++) {
	 if( nft==0) hist[nft][nh]->Draw();
	 if( nft==1) hist[nft][nh]->Draw("same");
	 if( nft==1) hist[nft][nh]->SetLineColor(2);
       }   
       cplot[nh]->cd(2);
       TH1F * hd = (TH1F*)hist[0][0]->Clone();
       hd->Divide(hist[0][0],hist[1][0],1,1);
       hd->Draw();
     }
}
