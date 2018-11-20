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

void plot_aero(TString basename) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
 //
     TString outputpdf;
 outputpdf="plots/"+basename+"_aero.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_shms_aero_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
 static const Int_t npad=7;
 static const Int_t iside=2;
 const char* sidename[iside]={"Neg","Pos"};
 const char* sname[iside]={"neg","pos"};
    // cout <<  hname[ip] << endl;
 TH1F *hAero_npe[iside][npad];
 TH1F *hAero_mult[iside][npad];
  TH1F *hAero_adctdcdiff[iside][npad];
TH2F *hAero_adctdcdiff_x[iside][npad];
 TH2F *hAero_adctdcdiff_y[iside][npad];
TH2F *hAero_npe_x[iside][npad];
 TH2F *hAero_npe_y[iside][npad];
 for (Int_t is=0;is<iside;is++) {
 for (Int_t ipad=0;ipad<npad;ipad++) {
       hAero_npe[is][ipad] = (TH1F*)fhistroot->Get(Form("hAero_npe_%s_pad_%d",sidename[is],ipad));
       hAero_mult[is][ipad] = (TH1F*)fhistroot->Get(Form("hAero_mult_%s_pad_%d",sidename[is],ipad));
       if (!hAero_npe[is][ipad]) cout << " no hist = " << Form("hAero_npe_%s_pad_%d",sidename[is],ipad) << endl;
       hAero_adctdcdiff[is][ipad] = (TH1F*)fhistroot->Get(Form("hAero_tdiff_%s_pad_%d",sidename[is],ipad));
       hAero_adctdcdiff_x[is][ipad] = (TH2F*)fhistroot->Get(Form("hAeroX_%s_pad_%d",sidename[is],ipad));
       hAero_adctdcdiff_y[is][ipad] = (TH2F*)fhistroot->Get(Form("hAeroY_%s_pad_%d",sidename[is],ipad));
 }}
   //
       //
       Int_t n=1;
       TCanvas *cAero_npe[2] ;
       for (Int_t is=0;is<iside;is++) {
	 cAero_npe[is] = new TCanvas(Form("cAero_npe_%d",is),Form("cAero_npe_%d",is),1000,700);
       cAero_npe[is]->Divide(4,2);
       n=1;
       for (Int_t ipad=0;ipad<npad;ipad++) {
	 cAero_npe[is]->cd(n++);
	 hAero_npe[is][ipad]->Draw();
       }
       if (is==0) cAero_npe[is]->Print(outputpdf+"(");
       if (is!=0) cAero_npe[is]->Print(outputpdf);
       }
       
       TCanvas *cAero_mult[2] ;
       for (Int_t is=0;is<iside;is++) {
	 cAero_mult[is] = new TCanvas(Form("cAero_mult_%d",is),Form("cAero_mult_%d",is),1000,700);
       cAero_mult[is]->Divide(4,2);
       n=1;
       for (Int_t ipad=0;ipad<npad;ipad++) {
	 cAero_mult[is]->cd(n++);
          gPad->SetLogy();
	 hAero_mult[is][ipad]->Draw();
       }
       cAero_mult[is]->Print(outputpdf);
       }
      
       TCanvas *cAero_adctdcdiff[2] ;
       for (Int_t is=0;is<iside;is++) {
	 cAero_adctdcdiff[is] = new TCanvas(Form("cAero_adctdcdiff_%d",is),Form("cAero_adctdcdiff_%d",is),1000,700);
       cAero_adctdcdiff[is]->Divide(4,2);
       n=1;
       for (Int_t ipad=0;ipad<npad;ipad++) {
	 cAero_adctdcdiff[is]->cd(n++);
	 hAero_adctdcdiff[is][ipad]->Draw();
       }
       cAero_adctdcdiff[is]->Print(outputpdf);
       }
       
       TCanvas *cAero_adctdcdiff_x[2] ;
       for (Int_t is=0;is<iside;is++) {
	 cAero_adctdcdiff_x[is] = new TCanvas(Form("cAero_adctdcdiff_x_%d",is),Form("cAero_adctdcdiff_x_%d",is),1000,700);
       cAero_adctdcdiff_x[is]->Divide(4,2);
       n=1;
       for (Int_t ipad=0;ipad<npad;ipad++) {
	 cAero_adctdcdiff_x[is]->cd(n++);
	 hAero_adctdcdiff_x[is][ipad]->Draw("colz");
       }
       cAero_adctdcdiff_x[is]->Print(outputpdf);
}

       TCanvas *cAero_adctdcdiff_y[2] ;
       for (Int_t is=0;is<iside;is++) {
	 cAero_adctdcdiff_y[is] = new TCanvas(Form("cAero_adctdcdiff_y_%d",is),Form("cAero_adctdcdiff_y_%d",is),1000,700);
       cAero_adctdcdiff_y[is]->Divide(4,2);
       n=1;
       for (Int_t ipad=0;ipad<npad;ipad++) {
	 cAero_adctdcdiff_y[is]->cd(n++);
	 hAero_adctdcdiff_y[is][ipad]->Draw("colz");
       }
       if (is==0) cAero_adctdcdiff_x[is]->Print(outputpdf);
       if (is==1) cAero_adctdcdiff_x[is]->Print(outputpdf+")");
       }
       
}
