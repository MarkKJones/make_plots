#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_jpsi_hms() {
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.2,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.16);
   TString outputpdf;
   outputpdf= "plots/lumi_jpsi_phase4_hms.pdf";
   TTree *T = new TTree("ntuple","data");
   Double_t nrun;
   Double_t curthres,ch,avecur,ptime;
   Double_t trig1rate,trig2rate,trig3rate,trig4rate,trig1yield,trig1yielderr,trig2yield,trig2yielderr,trig3yield,trig3yielderr,trig4yield,trig4yielderr,datayield,datayielderr;
      Long64_t nlines = T->ReadFile("LUMI_STUDY/lumi_jpsi_phase4_hms.dat","nrun/D:ch/D:curthres:avecur:ptime:trig2rate:trig4rate:datarate:trig2yield:trig2yielderr:trig4yield:trig4yielderr:datayield:datayielderr");
   printf(" found %lld points\n",nlines);
    T->SetBranchAddress("nrun",&nrun);
    T->SetBranchAddress("avecur",&avecur);
    T->SetBranchAddress("curthres",&curthres);
   T->SetBranchAddress("trig2rate",&trig2rate);
   T->SetBranchAddress("trig4rate",&trig4rate);
   T->SetBranchAddress("trig2yield",&trig2yield);
   T->SetBranchAddress("trig4yield",&trig4yield);
   T->SetBranchAddress("datayield",&datayield);
   T->SetBranchAddress("trig2yielderr",&trig2yielderr);
   T->SetBranchAddress("trig4yielderr",&trig4yielderr);
   T->SetBranchAddress("datayielderr",&datayielderr);
   vector<Double_t> v_avecur_low;
   vector<Double_t> v_avecur_low_err;
   vector<Double_t> v_avecur;
   vector<Double_t> v_avecur_err;
   vector<Double_t> v_trig2yield;
   vector<Double_t> v_trig2yield_err;
   vector<Double_t> v_trig4yield;
   vector<Double_t> v_trig4yield_err;
   vector<Double_t> v_datayield;
   vector<Double_t> v_datayield_err;
   Double_t max_trig2=0;
   Double_t max_trig4=0;
   Double_t max_data=0;
   for (Int_t n=0;n<nlines;n++) {
     T->GetEntry(n);
     cout << nrun << " " << curthres << " " << avecur  << " " << trig2yield<< " " << trig2yielderr<< " " << trig3yield<< " " << trig3yielderr << endl;
      v_avecur.push_back(avecur);
     v_avecur_err.push_back(.001);
     v_trig2yield.push_back(trig2yield);
     v_trig2yield_err.push_back(trig2yielderr);
     v_trig4yield.push_back(trig4yield);
     v_trig4yield_err.push_back(trig4yielderr);
     v_datayield.push_back(datayield);
     v_datayield_err.push_back(datayielderr);
     if (trig2yield>max_trig2) max_trig2=trig2yield;
     if (trig4yield>max_trig4) max_trig4=trig4yield;
     if (datayield>max_data) max_data=datayield;
    }
   // plot yields
  TMultiGraph *mgt2 = new TMultiGraph();
  TLegend* legt2 = new TLegend(0.2,0.75,0.48,0.89);
   TGraphErrors *gr_trig2_low;
   TGraphErrors *gr_trig2;
   gr_trig2 = new TGraphErrors(v_avecur.size(),&(v_avecur[0]),&(v_trig2yield[0]),0,&(v_trig2yield_err[0]));
     gr_trig2->SetMarkerStyle(20);
     gr_trig2->SetMarkerColor(2);
     mgt2->Add(gr_trig2);
  legt2->SetBorderSize(1);
  legt2->SetFillColor(0);
   legt2->AddEntry(gr_trig2,"HMS ELCLEAN scal I > 5","p");
  TCanvas *ctrig2 = new TCanvas("ctrig2"," trig2",700,700);
   ctrig2->Divide(1,1);
   ctrig2->cd(1);
   mgt2->Draw("AP");
   mgt2->SetTitle("Jpsi ; Current; TRIG3 HMS scal Yield ");
   mgt2->SetMaximum(max_trig2*1.1);
   mgt2->SetMinimum(max_trig2*0.9);
     legt2->Draw();
  ctrig2->Print(outputpdf+"(");
    // plot yields
  TMultiGraph *mgt3 = new TMultiGraph();
  TLegend* legt3 = new TLegend(0.2,0.75,0.48,0.89);
   TGraphErrors *gr_trig4;
   gr_trig4 = new TGraphErrors(v_avecur.size(),&(v_avecur[0]),&(v_trig4yield[0]),0,&(v_trig4yield_err[0]));
     gr_trig4->SetMarkerStyle(22);
     gr_trig4->SetMarkerColor(2);
     mgt3->Add(gr_trig4);
  legt3->SetBorderSize(1);
  legt3->SetFillColor(0);
   legt3->AddEntry(gr_trig4," HMS trig 4 coin I >5","p");
  TCanvas *ctrig4 = new TCanvas("ctrig4"," trig4",700,700);
   ctrig4->Divide(1,1);
   ctrig4->cd(1);
   mgt3->Draw("AP");
   mgt3->SetTitle("Jpsi; Current; TRIG4 Coin scal Yield ");
   mgt3->SetMaximum(max_trig4*1.1);
   mgt3->SetMinimum(max_trig4*0.9);
    legt3->Draw();
  ctrig4->Print(outputpdf);
    // plot yields
  TMultiGraph *mgdata = new TMultiGraph();
  TLegend* legdata = new TLegend(0.2,0.75,0.48,0.89);
   TGraphErrors *gr_data;
   gr_data = new TGraphErrors(v_avecur.size(),&(v_avecur[0]),&(v_datayield[0]),0,&(v_datayield_err[0]));
     gr_data->SetMarkerStyle(22);
     gr_data->SetMarkerColor(2);
     mgdata->Add(gr_data);
  legdata->SetBorderSize(1);
  legdata->SetFillColor(0);
   legdata->AddEntry(gr_data,"Data I > 1.5","p");
  TCanvas *cdata = new TCanvas("cdata"," data",700,700);
   cdata->Divide(1,1);
   cdata->cd(1);
   mgdata->Draw("AP");
   mgdata->SetTitle("Jpsi; Current; DATA HMS Yield ");
   mgdata->SetMaximum(max_data*1.1);
   mgdata->SetMinimum(max_data*0.8);
    legdata->Draw();
  cdata->Print(outputpdf+")");
  //
}
