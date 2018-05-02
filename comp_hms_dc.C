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

void comp_hms_dc() {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
    TString basename = "hms_replay_dc_calib_3454_test";
     TString outputpdf;
    outputpdf="plots/"+basename+"_comp_dc.pdf";
    /*
const UInt_t nftot=4;
const UInt_t nplots=4;
 TString setname[nftot]={"15","17","16","18b"};
 Int_t colind[nftot]={1,2,3,6};
 TString lname[nftot]={"xt<100,yt<100,xpt<1. xpdiff<1000","xt<20,yt<5,xpt<.4 xpdiff<1000","xt<100,yt<100,xpt<1. xpdiff<.2","xt<20,yt<5,xpt<.4 xpdiff<0.2"};
    */
    //
const UInt_t nftot=2;
const UInt_t nplots=4;
 TString setname[nftot]={"2","3"};
 Int_t colind[nftot]={1,2};
 TString lname[nftot]={"test2","test3"};
 //
 TString h1name[nplots]={"hstubx_diff","hstuby_diff","hstubxp_diff","hstub1_xp_diff"};
  TString inputroot;
   TFile *fhistroot[nftot];
   TH2F *fhist[nftot][nplots];
   TH2F *fhistresid[nftot][12];
 static const Int_t plnum=12;
 const char* plname[plnum]={"u1","u2","x1","x2","v1","v2","v2","v1","x2","x1","u2","u1"};
 const char* chname[plnum]={"1","1","1","1","1","1","2","2","2","2","2","2"};
 TString rname[plnum];
   for (UInt_t nf=0;nf<nftot;nf++) {
     inputroot="hist/"+basename+setname[nf]+"_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[nf] =  new TFile(inputroot);
   for (UInt_t nh=0;nh<nplots;nh++) {
     fhist[nf][nh] = (TH2F*)fhistroot[nf]->Get(h1name[nh]);
   }
   for (UInt_t ip=0;ip<plnum;ip++) {
     rname[ip]=Form("resid_%s%s",chname[ip],plname[ip]);
     fhistresid[nf][ip] = (TH2F*)fhistroot[nf]->Get(Form("resid_%s%s",chname[ip],plname[ip]));
   }
   }
   //
    TCanvas *ctime[nplots];
    TLegend *ltime[nplots];
   for (UInt_t nh=0;nh<nplots;nh++) {
     ctime[nh] = new TCanvas(Form("ctime_%d",nh),h1name[nh], 1000,700);
     ctime[nh]->Divide(1,1);
     ctime[nh]->cd(1);
     //gPad->SetLogy();
      ltime[nh] = new TLegend(.59,.75,.99,.95,"");
      for (UInt_t nf=0;nf<nftot;nf++) {
	if (nf==0) fhist[nf][nh]->Draw();
	if (nf!=0) fhist[nf][nh]->Draw("same");
	fhist[nf][nh]->SetLineColor(colind[nf]);
	Int_t temp=fhist[nf][nh]->Integral();
        TString tt=lname[nf]+Form(" Int=%d",temp);
	ltime[nh]->AddEntry(fhist[nf][nh],tt);
      }
      ltime[nh]->Draw();
      if (nh==0) ctime[nh]->Print(outputpdf+"(");
      if (nh!=0) ctime[nh]->Print(outputpdf);
   }
   //
    TCanvas *cresid[plnum];
    TLegend *lresid[plnum];
   for (UInt_t nh=0;nh<plnum;nh++) {
     cresid[nh] = new TCanvas(Form("cresid_%d",nh),rname[nh], 1000,700);
     cresid[nh]->Divide(1,1);
     cresid[nh]->cd(1);
     //gPad->SetLogy();
      lresid[nh] = new TLegend(.59,.75,.99,.95,"");
      for (UInt_t nf=0;nf<nftot;nf++) {
	if (nf==0) fhistresid[nf][nh]->Draw();
	if (nf!=0) fhistresid[nf][nh]->Draw("same");
	fhistresid[nf][nh]->SetLineColor(colind[nf]);
	Int_t temp=fhistresid[nf][nh]->Integral();
        TString tt=lname[nf]+Form(" Int=%d",temp);
	lresid[nh]->AddEntry(fhistresid[nf][nh],tt);
      }
      lresid[nh]->Draw();
      if (nh==plnum-1) cresid[nh]->Print(outputpdf+")");
      if (nh!=plnum-1) cresid[nh]->Print(outputpdf);
   }

   //
}
