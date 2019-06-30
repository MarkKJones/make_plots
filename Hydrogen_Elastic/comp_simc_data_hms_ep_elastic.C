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
void comp_simc_data_hms_ep_elastic( TString fn1 , TString fn2 ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/comp_simc_hms_ep_elastic_"+fn2+".pdf";
    //
const UInt_t nftot=2;
 Int_t colind[nftot]={1,2};
 TString lname[nftot]={"data","simc"};
   TFile *fhistroot[nftot];
 //
const UInt_t nplots=6;
 TString h1name[nplots]={"hW","hEmiss","hPmiss","hPmissx","hPmissy","hPmissz"};
 //
const UInt_t n2plots=3;
 TString h2name[n2plots]={"hxptar","hyptar","hdelta"};
const UInt_t n3plots=4;
 TString h3name[n3plots]={"pxptar","pyptar","pdelta","hprot_mom_calc"};
const UInt_t n4plots=4;
 TString h4name[n4plots]={"hxfp","hyfp","hxpfp","hypfp"};
const UInt_t n5plots=4;
 TString h5name[n5plots]={"pxfp","pyfp","pxpfp","pypfp"};
 //
 TString inputroot[2];
 inputroot[0] = "hist/"+fn1+"_hms_ep_elastic_hist.root";
 inputroot[1] = "hist/"+fn2+"_simc_hms_ep_elastic_hist.root";
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
   for (UInt_t nh=0;nh<n2plots;nh++) {
     fhist2[nf][nh] = (TH1F*)fhistroot[nf]->Get(h2name[nh]);
   }
   for (UInt_t nh=0;nh<n3plots;nh++) {
     fhist3[nf][nh] = (TH1F*)fhistroot[nf]->Get(h3name[nh]);
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
    TF1* fW[nftot];
    Double_t Wsum[2];
    Double_t Wpeak[2];
    Double_t Wpeakerr[nftot];
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
	if (nh==0) {
 	Double_t WMax = fhist[nf][nh]->GetBinCenter(fhist[nf][nh]->GetMaximumBin());
	fW[nf] = new TF1(Form("fW_%d",nf),"gaus",WMax-.02,WMax+.02);
	fhist[nf][nh]->Fit(Form("fW_%d",nf),"QR");
        Wpeak[nf]=fW[nf]->GetParameter(1);
        Wpeakerr[nf]=fW[nf]->GetParError(1);
	Wsum[nf]=fhist[nf][nh]->Integral(fhist[nf][nh]->FindBin(WMax-.02),fhist[nf][nh]->FindBin(WMax+.1));
	     }
	if (nh==0 )lW->AddEntry(fhist[nf][nh],lname[nf]);
       }
     if (nh==0 ) lW->Draw();
      }
      cW->Print(outputpdf+"(");
      cout << " ratio = " << Wsum[0]/Wsum[1] << " =/- " << Wsum[0]/Wsum[1]*TMath::Sqrt(1./Wsum[0])<< endl;
 	cout << " Data W = " << Wpeak[0] << " Simc W = " << Wpeak[1] << " Data W - SIMC W = " << Wpeak[0]-Wpeak[1] << endl;
	cout << Wpeak[0] << " " << Wpeakerr[0] << " " << Wpeak[1]  << " " << Wpeakerr[1] << " " << (Wpeak[0]-Wpeak[1]) << " " << TMath::Sqrt(Wpeakerr[0]*Wpeakerr[0]+Wpeakerr[1]*Wpeakerr[1])<< endl;
     //
  //  histograms
    TCanvas *cSHMS;
    TLegend *lSHMS;
     cSHMS = new TCanvas("cSHMS","SHMS", 1000,700);
     cSHMS->Divide(2,2);
      lSHMS = new TLegend(.65,.75,.99,.95,"");
      for (UInt_t nh=0;nh<n3plots;nh++) {
      cSHMS->cd(nh+1);
     for (UInt_t nf=0;nf<nftot;nf++) {
       	if (nf==0) fhist3[nf][nh]->Draw();
       	if (nf!=0) fhist3[nf][nh]->Draw("same");
      	if (nf==0) fhist3[nf][nh]->SetLineColor(2);
	if (nh==0 )lSHMS->AddEntry(fhist3[nf][nh],lname[nf]);
       }
      if (nh==0 )lSHMS->Draw();
      }
      cSHMS->Print(outputpdf);
      //
 //  histograms
    TCanvas *cHMSfp;
    TLegend *lHMSfp;
     cHMSfp = new TCanvas("cHMSfp","HMSfp", 1000,700);
     cHMSfp->Divide(2,2);
      lHMSfp = new TLegend(.65,.75,.99,.95,"");
      for (UInt_t nh=0;nh<n4plots;nh++) {
      cHMSfp->cd(nh+1);
     for (UInt_t nf=0;nf<nftot;nf++) {
       	if (nf==0) fhist4[nf][nh]->Draw();
       	if (nf!=0) fhist4[nf][nh]->Draw("same");
      	if (nf==0) fhist4[nf][nh]->SetLineColor(2);
	if (nh==0 )lHMSfp->AddEntry(fhist4[nf][nh],lname[nf]);
       }
      if (nh==0 )lHMSfp->Draw();
      }
      cHMSfp->Print(outputpdf);
      //
//  histograms
    TCanvas *cSHMSfp;
    TLegend *lSHMSfp;
     cSHMSfp = new TCanvas("cSHMSfp","SHMSfp", 1000,700);
     cSHMSfp->Divide(2,2);
      lSHMSfp = new TLegend(.65,.75,.99,.95,"");
      for (UInt_t nh=0;nh<n5plots;nh++) {
      cSHMSfp->cd(nh+1);
     for (UInt_t nf=0;nf<nftot;nf++) {
       	if (nf==0) fhist5[nf][nh]->Draw();
       	if (nf!=0) fhist5[nf][nh]->Draw("same");
      	if (nf==0) fhist5[nf][nh]->SetLineColor(2);
	if (nh==0 )lSHMSfp->AddEntry(fhist5[nf][nh],lname[nf]);
       }
      if (nh==0 )lSHMSfp->Draw();
      }
      cSHMSfp->Print(outputpdf);
      //
  //  histograms
    TCanvas *cHMS;
    TLegend *lHMS;
     cHMS = new TCanvas("cHMS","HMS", 1000,700);
     cHMS->Divide(2,2);
      lHMS = new TLegend(.65,.75,.99,.95,"");
      for (UInt_t nh=0;nh<n2plots;nh++) {
      cHMS->cd(nh+1);
     for (UInt_t nf=0;nf<nftot;nf++) {
       	if (nf==0) fhist2[nf][nh]->Draw();
       	if (nf!=0) fhist2[nf][nh]->Draw("same");
      	if (nf==0) fhist2[nf][nh]->SetLineColor(2);
	if (nh==0 )lHMS->AddEntry(fhist2[nf][nh],lname[nf]);
       }
      if (nh==0 )lHMS->Draw();
      }
      cHMS->Print(outputpdf+")");
      //
}
	   //
