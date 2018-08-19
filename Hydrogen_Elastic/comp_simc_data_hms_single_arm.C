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
void comp_simc_data_hms_single_arm( TString fn1 , TString fn2 ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/comp_simc_hms_single_arm.pdf";
    //
const UInt_t nftot=2;
 Int_t colind[nftot]={1,2};
 TString lname[nftot]={"data","simc"};
   TFile *fhistroot[nftot];
 //
const UInt_t nplots=1;
 TString h1name[nplots]={"hW"};
 //
 TString inputroot[2];
 inputroot[0] = "hist/"+fn1+".root";
 inputroot[1] = "hist/"+fn2+".root";
  TH1F *fhist[nftot][nplots];
  //
   for (UInt_t nf=0;nf<nftot;nf++) {
     cout << " infile root = " << inputroot[nf] << endl;
   fhistroot[nf] =  new TFile(inputroot[nf]);
   for (UInt_t nh=0;nh<nplots;nh++) {
     fhist[nf][nh] = (TH1F*)fhistroot[nf]->Get(h1name[nh]);
   }
   }
   // W histograms
    TCanvas *cW;
    TLegend *lW;
    TF1* fW[nftot];
    Double_t Wsum[nftot];
     cW = new TCanvas("cW","W", 1000,700);
     cW->Divide(1,1);
     cW->cd(1);
     Int_t nh=0;
      lW = new TLegend(.59,.75,.99,.95,"");
      for (UInt_t nf=0;nf<nftot;nf++) {
	if (nf==0) fhist[nf][nh]->Draw();
	if (nf!=0) fhist[nf][nh]->Draw("same");
	Double_t WMax = fhist[nf][nh]->GetBinCenter(fhist[nf][nh]->GetMaximumBin());
	fW[nf] = new TF1(Form("fW_%d",nf),"gaus",WMax-.02,WMax+.02);
	fhist[nf][nh]->SetLineColor(colind[nf]);
	fhist[nf][nh]->Fit(Form("fW_%d",nf),"QR");
	Wsum[nf]=fhist[nf][nh]->Integral(fhist[nf][nh]->FindBin(WMax-.02),fhist[nf][nh]->FindBin(WMax+.1));
        TString tt=lname[nf]+Form(" Int= %f",Wsum[nf]);
	lW->AddEntry(fhist[nf][nh],tt);
	cout << nf << " " << Wsum[nf] << WMax << " " << fW[nf]->GetParameter(1) << " " << fW[nf]->GetParameter(2) << endl;
       }
      lW->Draw();
      cW->Print(outputpdf+"(");
      cout << " ratio = " << Wsum[0]/Wsum[1] << endl;
      //
}
	   //
