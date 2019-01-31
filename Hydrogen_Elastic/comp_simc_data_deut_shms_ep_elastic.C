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
void comp_simc_data_deut_shms_ep_elastic( TString fn1 , TString fn2 ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/comp_simc_deut_shms_ep_elastic_"+fn2+".pdf";
    //
const UInt_t nftot=2;
 Int_t colind[nftot]={1,2};
 TString lname[nftot]={"data","simc"};
   TFile *fhistroot[nftot];
 //
const UInt_t nplots=5;
 TString h1name[nplots]={"hW","h_edelta","h_pdelta","h_exptar","h_eyptar"};
 //
const UInt_t n2plots=3;
 TString h2name[n2plots]={"hxptar","hyptar","hdelta"};
const UInt_t n3plots=4;
 TString h3name[n3plots]={"pxptar","pyptar","pdelta","hprot_mom_calc"};
const UInt_t n4plots=4;
 TString h4name[n4plots]={"h_exfp","h_eyfp","h_expfp","h_eypfp"};
const UInt_t n5plots=4;
 TString h5name[n5plots]={"h_pxfp","h_pyfp","h_pxpfp","h_pypfp"};
 //
 TString inputroot[2];
 inputroot[0] = "hist/"+fn1+"_deut_shms_ep_elastic_hist.root";
 inputroot[1] = "hist/"+fn2+"_simc_deut_shms_ep_elastic_hist.root";
  TH1F *fhist[nftot][nplots];
  TH1F *fhist2[nftot][n2plots];
  TH1F *fhist3[nftot][n3plots];
  TH1F *fhist4[nftot][n4plots];
  TH1F *fhist5[nftot][n5plots];
  //
   for (UInt_t nf=0;nf<nftot;nf++) {
     cout << " infile root = " << inputroot[nf] << endl;
   fhistroot[nf] =  new TFile(inputroot[nf]);
   for (UInt_t nh=0;nh<nplots;nh++) {
     fhist[nf][nh] = (TH1F*)fhistroot[nf]->Get(h1name[nh]);
   }
   for (UInt_t nh=0;nh<n4plots;nh++) {
     fhist4[nf][nh] = (TH1F*)fhistroot[nf]->Get(h4name[nh]);
   }
   for (UInt_t nh=0;nh<n5plots;nh++) {
     fhist5[nf][nh] = (TH1F*)fhistroot[nf]->Get(h5name[nh]);
   }
    }
   // W histograms
    TCanvas *cW;
    TLegend *lW;
    Double_t ymax;
    TF1* fW[nftot][3];
    Double_t Wsum[2][3];
    Double_t Wpeak[2][3];
    Double_t Wpeakerr[nftot][3];
     cW = new TCanvas("cW","W", 1000,700);
     cW->Divide(2,3);
      lW = new TLegend(.65,.75,.99,.95,"");
      for (UInt_t nh=0;nh<nplots;nh++) {
      cW->cd(nh+1);
      ymax=0;
     for (UInt_t nf=0;nf<nftot;nf++) {
       if (fhist[nf][nh]->GetMaximum()> ymax) ymax = fhist[nf][nh]->GetMaximum();
     }
     for (UInt_t nf=0;nf<nftot;nf++) {
       //	if (nf==0) fhist[nf][nh]->DrawNormalized();
       //	if (nf!=0) fhist[nf][nh]->DrawNormalized("same");
       	if (nf==0) fhist[nf][nh]->SetMaximum(ymax);
       	if (nf==0) fhist[nf][nh]->Draw();
       	if (nf!=0) fhist[nf][nh]->Draw("same");
      	if (nf==0) fhist[nf][nh]->SetLineColor(2);
	if (nh<1) {
 	Double_t WMax = fhist[nf][nh]->GetBinCenter(fhist[nf][nh]->GetMaximumBin());
	fW[nf][nh] = new TF1(Form("fW_%d_%d",nf,nh),"gaus",WMax-.02,WMax+.02);
	fhist[nf][nh]->Fit(Form("fW_%d_%d",nf,nh),"QR");
        Wpeak[nf][nh]=fW[nf][nh]->GetParameter(1);
        Wpeakerr[nf][nh]=fW[nf][nh]->GetParError(1);
	Wsum[nf][nh]=fhist[nf][nh]->Integral(fhist[nf][nh]->FindBin(WMax-.1),fhist[nf][nh]->FindBin(WMax+.1));
	     }
	if (nh==0 )lW->AddEntry(fhist[nf][nh],lname[nf]);
       }
     if (nh==0 ) lW->Draw();
     if (nh<1) {
       cout << " Data sum = " << Wsum[0][nh] << " SIMC sum = " << Wsum[1][nh] << endl ;
cout << " ratio = " << Wsum[0][nh]/Wsum[1][nh] << " =/- " << Wsum[0][nh]/Wsum[1][nh]*TMath::Sqrt(1./Wsum[0][nh])<< endl;
 	cout << " Data W = " << Wpeak[0][nh] << " Simc W = " << Wpeak[1][nh] << " Data W - SIMC W = " << Wpeak[0][nh]-Wpeak[1][nh] << endl;
	cout << Wpeak[0][nh] << " " << Wpeakerr[0][nh] << " " << Wpeak[1][nh]  << " " << Wpeakerr[1][nh] << " " << (Wpeak[0][nh]-Wpeak[1][nh]) << " " << TMath::Sqrt(Wpeakerr[0][nh]*Wpeakerr[0][nh]+Wpeakerr[1][nh]*Wpeakerr[1][nh])<< endl;
     }
      }
      cW->Print(outputpdf+"(");
   // W histograms
    TCanvas *cpfp;
    TLegend *lpfp;
     cpfp= new TCanvas("cpfp","Hadron fp", 1000,700);
     cpfp->Divide(2,2);
      lpfp = new TLegend(.65,.75,.99,.95,"");
      for (UInt_t nh=0;nh<n5plots;nh++) {
      cpfp->cd(nh+1);
      ymax=0;
     for (UInt_t nf=0;nf<nftot;nf++) {
       if (fhist5[nf][nh]->GetMaximum()> ymax) ymax = fhist5[nf][nh]->GetMaximum();
     }
     for (UInt_t nf=0;nf<nftot;nf++) {
       	if (nf==0) fhist5[nf][nh]->SetMaximum(ymax*1.1);
       	if (nf==0) fhist5[nf][nh]->Draw();
       	if (nf!=0) fhist5[nf][nh]->Draw("same");
      	if (nf==0) fhist5[nf][nh]->SetLineColor(2);
	if (nh==0 )lpfp->AddEntry(fhist5[nf][nh],lname[nf]);
       }
     if (nh==0 ) lpfp->Draw();
      }
      cpfp->Print(outputpdf);
   // W histograms
    TCanvas *cefp;
    TLegend *lefp;
     cefp= new TCanvas("cefp","Electron fp", 1000,700);
     cefp->Divide(2,2);
      lefp = new TLegend(.65,.75,.99,.95,"");
      for (UInt_t nh=0;nh<n4plots;nh++) {
      cefp->cd(nh+1);
      ymax=0;
     for (UInt_t nf=0;nf<nftot;nf++) {
       if (fhist4[nf][nh]->GetMaximum()> ymax) ymax = fhist4[nf][nh]->GetMaximum();
     }
     for (UInt_t nf=0;nf<nftot;nf++) {
       	if (nf==0) fhist4[nf][nh]->SetMaximum(ymax*1.1);
       	if (nf==0) fhist4[nf][nh]->Draw();
       	if (nf!=0) fhist4[nf][nh]->Draw("same");
      	if (nf==0) fhist4[nf][nh]->SetLineColor(2);
	if (nh==0 )lefp->AddEntry(fhist4[nf][nh],lname[nf]);
       }
     if (nh==0 ) lefp->Draw();
      }
      cefp->Print(outputpdf+")");

}

