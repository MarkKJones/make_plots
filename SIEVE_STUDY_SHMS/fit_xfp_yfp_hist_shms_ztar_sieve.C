#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
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

void fit_xfp_yfp_hist_shms_ztar_sieve(TString basename,TString label) {
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
 outputpdf="plots/"+basename+"fit_xfp_yfp_shms_ztar_sieve_hist.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_shms_ztar_sieve_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   TH2F* hyfp_yxfp_cent_foil = (TH2F*)fhistroot->Get("hyfp_yxfp_cent_foil");
  TH2F* hyfp_yxfp_cent_foil_ypcut[5];
  	for (Int_t iz = 0; iz < 5; iz++) {
	  hyfp_yxfp_cent_foil_ypcut[iz] = (TH2F*)fhistroot->Get(Form("hyfp_yxfp_cent_foil_ypcut_%d",iz));
	}
	    //
 TF1 *linfit[5];
	    TCanvas *cfit = new TCanvas("cfit"," Ind fits",900,700);
	  cfit->Divide(4,2);
	  for (Int_t n=0;n<5;n++) {
            linfit[n] = new TF1(Form("linfit_%d",n),"pol1",-10,10);
            linfit[n]->SetLineColor(2);
            linfit[n]->SetLineWidth(2);
	    cfit->cd(n+1);
	    hyfp_yxfp_cent_foil_ypcut[n]->Draw("colz");
            hyfp_yxfp_cent_foil_ypcut[n]->Fit(Form("linfit_%d",n),"Q");
  cout << linfit[n]->GetParameter(0) << " " << linfit[n]->GetParError(0) << " " <<linfit[n]->GetParameter(1) << " " <<linfit[n]->GetParError(1) <<endl;
	  }
	  //
	    TCanvas *call = new TCanvas("call","yfp v xfp",700,700);
	    call->Divide(1,1);
	    call->cd(1);
	    hyfp_yxfp_cent_foil->Draw("colz");
	  for (Int_t n=0;n<5;n++) {
            linfit[n]->Draw("same");
	  }
	//
}
