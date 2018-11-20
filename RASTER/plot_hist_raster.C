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
#include <TEllipse.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_hist_raster(TString basename,Double_t ebeam=8.518,Double_t rsize=2.) {
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
 outputpdf="plots/"+basename+"_raster.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_raster_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   static const Int_t nh1=8;
   TString h1name[nh1]={"hrastxaRawAdc","hrastxbRawAdc","hrastyaRawAdc","hrastybRawAdc","hrastxaVolts","hrastxbVolts","hrastyaVolts","hrastybVolts"};
   const char* rname[4]={"gfrxa","gfrxb","gfrya","gfryb"};
   TH1F* hist1d[nh1];
  for (Int_t ip=0;ip<nh1;ip++) {
    // cout <<  hname[ip] << endl;
       hist1d[ip] = (TH1F*)fhistroot->Get(h1name[ip]);
       if (!hist1d[ip]) cout << " no hist = " << h1name[ip] << endl;
  }
  Double_t near_edge[4]={0.};
	Double_t far_edge[4]={0.};
       Double_t cbin=5;
  for (Int_t nh=0;nh<4;nh++) {
	Int_t nbin = hist1d[nh]->GetNbinsX();
	Bool_t looking_for_near_edge = kTRUE;
        if (cbin<hist1d[nh]->GetMaximum()*.1) cbin = hist1d[nh]->GetMaximum()*.1;
	Bool_t looking_for_far_edge = kFALSE;
 	for (Int_t i=0;i<nbin;i++) {
	  if(looking_for_near_edge &&hist1d[nh]->GetBinContent(i)>cbin ) {
	    looking_for_near_edge= kFALSE;
	    looking_for_far_edge= kTRUE;
	    near_edge[nh] = hist1d[nh]->GetBinCenter(i);
	  } 
	  if(looking_for_far_edge &&hist1d[nh]->GetBinContent(i)<cbin ) {
	    looking_for_far_edge= kFALSE;
	    far_edge[nh] = hist1d[nh]->GetBinCenter(i);
	  } 
	}
  }
  cout << "gfr_cal_mom  = " << ebeam << endl;
  for (Int_t nh=0;nh<4;nh++) {
    cout << Form("%s_adc_zero_offset = ",rname[nh]) <<  (far_edge[nh]+near_edge[nh])/2. << endl;
  }
  for (Int_t nh=0;nh<4;nh++) {
    cout << Form("%s_adcpercm = ",rname[nh]) <<  (far_edge[nh]-near_edge[nh])/(rsize/10) << endl;
  }
  for (Int_t nh=0;nh<4;nh++) {
    cout << Form("%s_adcpercm divide by cal_mom = ",rname[nh]) <<  (far_edge[nh]-near_edge[nh])/(rsize/10)/ebeam << endl;
  }
  //
  TCanvas *craster = new TCanvas("craster","Raster",800,700);
  craster->Divide(2,2);
  for (Int_t nh=0;nh<4;nh++) {
     craster->cd(nh+1);
     hist1d[nh]->Draw();
     TLine *ll = new TLine(near_edge[nh],0,near_edge[nh],hist1d[nh]->GetMaximum()/2.);
     ll->Draw();
     ll->SetLineColor(2);
     ll->SetLineWidth(2);
     TLine *hl = new TLine(far_edge[nh],0,far_edge[nh],hist1d[nh]->GetMaximum()/2.);
     hl->Draw();
     hl->SetLineColor(2);
     hl->SetLineWidth(2);
  }  
  craster->Print(outputpdf+"(");
  //
  TCanvas *c2raster = new TCanvas("c2raster","Raster",800,700);
  c2raster->Divide(2,2);
  for (Int_t nh=0;nh<4;nh++) {
     c2raster->cd(nh+1);
     hist1d[nh+4]->Draw();
  }  
  c2raster->Print(outputpdf+")");
  //
}
