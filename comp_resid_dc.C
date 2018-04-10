#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
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

void comp_resid_dc() {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
    TString basename = "shms_replay_dc_calib_1791_test";
     TString outputpdf;
const UInt_t nftot=2;
 TString setname[nftot]={"3","4"};
 Int_t colind[nftot]={1,2};
 // TString lname[nftot]={"SA=0","SA=1"};
 TString lname[nftot]={"SA=1 Zero","SA=1 Zero NoTC"};
 //
 outputpdf="plots/"+basename+"_comp";
   for (UInt_t nf=0;nf<nftot;nf++) {
     outputpdf=outputpdf+"_"+setname[nf];
   }
    outputpdf= outputpdf+"_resid.pdf";
  TString inputroot;
   TFile *fhistroot[nftot];
   TH2F *fhistresid[nftot][12];
 static const Int_t plnum=12;
 const char* plname[plnum]={"u1","u2","x1","x2","v1","v2","v2","v1","x2","x1","u2","u1"};
 const char* chname[plnum]={"1","1","1","1","1","1","2","2","2","2","2","2"};
 TString rname[plnum];
   for (UInt_t nf=0;nf<nftot;nf++) {
     inputroot="hist/"+basename+setname[nf]+"_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[nf] =  new TFile(inputroot);
   for (UInt_t ip=0;ip<plnum;ip++) {
     rname[ip]=Form("exclresid_%s%s",chname[ip],plname[ip]);
     fhistresid[nf][ip] = (TH2F*)fhistroot[nf]->Get(Form("exclresid_%s%s",chname[ip],plname[ip]));
   }
   }
   //
    //
    TCanvas *cresid[plnum];
    TLegend *lresid[plnum];
    TF1 *gfit[nftot][plnum];
    Double_t mean[nftot][plnum];
    Double_t sig[nftot][plnum];
   for (UInt_t nh=0;nh<plnum;nh++) {
     cresid[nh] = new TCanvas(Form("cresid_%d",nh),rname[nh], 1000,700);
     cresid[nh]->Divide(1,1);
     cresid[nh]->cd(1);
     gStyle->SetOptFit(0);
      lresid[nh] = new TLegend(.59,.75,.99,.95,"");
      for (UInt_t nf=0;nf<nftot;nf++) {
	if (nf==0) fhistresid[nf][nh]->Draw();
	if (nf!=0) fhistresid[nf][nh]->Draw("same");
        gfit[nf][nh]= new TF1(Form("gaus_%d",nh),"gaus",-.05,.05);
	gfit[nf][nh]->SetLineColor(colind[nf]);
	fhistresid[nf][nh]->Fit(Form("gaus_%d",nh),"Q","",-.05,.05);
	fhistresid[nf][nh]->SetLineColor(colind[nf]);
	  Double_t temp=gfit[nf][nh]->GetParameter(1);
	mean[nf][nh]=gfit[nf][nh]->GetParameter(1);
        TString tt=lname[nf]+Form(" Mean=%5.4f",mean[nf][nh]);
	  temp=gfit[nf][nh]->GetParameter(2);
	sig[nf][nh]=temp;
	lresid[nh]->AddEntry(fhistresid[nf][nh],tt);
      }
      lresid[nh]->Draw();
      if (nh==plnum-1) cresid[nh]->Print(outputpdf+")");
      if (nh!=plnum-1) cresid[nh]->Print(outputpdf);
   }
   for (UInt_t nh=0;nh<plnum;nh++) {
     cout << chname[nh] << plname[nh] ;
     for (UInt_t nf=0;nf<nftot;nf++) {
       cout << " "<< mean[nf][nh] << " " << sig[nf][nh];
     }
     cout << endl;
   }
   //
}
