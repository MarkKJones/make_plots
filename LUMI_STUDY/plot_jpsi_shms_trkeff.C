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

void plot_jpsi_shms_trkeff() {
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.2,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.16);
   TString outputpdf;
   outputpdf= "plots/lumi_jpsi_shms_treff.pdf";
   TTree *T = new TTree("ntuple","data");
   Double_t nrun;
   Double_t datayield,datayielderr,avecur;
   Double_t trig1rate,s1xrate,s1yrate,s2xrate,s2yrate,treff,dt;
      Long64_t nlines = T->ReadFile("LUMI_STUDY/lumi_jpsi_shms_rates.dat","nrun/D:datayield:datayielderr:avecur:trig1rate:s1xrate:s1yrate:s2xrate:s2yrate:treff:dt");
   printf(" found %lld points\n",nlines);
    T->SetBranchAddress("nrun",&nrun);
    T->SetBranchAddress("avecur",&avecur);
    T->SetBranchAddress("datayield",&datayield);
    T->SetBranchAddress("datayielderr",&datayielderr);
    T->SetBranchAddress("trig1rate",&trig1rate);
    T->SetBranchAddress("s1xrate",&s1xrate);
   T->SetBranchAddress("s2xrate",&s2xrate);
   T->SetBranchAddress("s1yrate",&s1yrate);
   T->SetBranchAddress("s2yrate",&s2yrate);
   T->SetBranchAddress("treff",&treff);
   vector<Double_t> v_avecur_ph2;
   vector<Double_t> v_avecur_err_ph2;
    vector<Double_t> v_datayield_ph2;
   vector<Double_t> v_datayield_err_ph2;
    vector<Double_t> v_datayield_ph1;
   vector<Double_t> v_datayield_err_ph1;
    vector<Double_t> v_datayield_ph4;
   vector<Double_t> v_datayield_err_ph4;
  vector<Double_t> v_trig1rate_ph2;
   vector<Double_t> v_trig1rate_err_ph2;
   vector<Double_t> v_s1xrate_ph2;
   vector<Double_t> v_s2xrate_ph2;
   vector<Double_t> v_s1yrate_ph2;
   vector<Double_t> v_s2yrate_ph2;
   vector<Double_t> v_treff_ph2;
   vector<Double_t> v_s1xrate_err_ph2;
   vector<Double_t> v_s2xrate_err_ph2;
   vector<Double_t> v_s1yrate_err_ph2;
   vector<Double_t> v_s2yrate_err_ph2;
   vector<Double_t> v_treff_err_ph2;
   //
   vector<Double_t> v_avecur_ph1;
   vector<Double_t> v_avecur_err_ph1;
   vector<Double_t> v_trig1rate_ph1;
   vector<Double_t> v_trig1rate_err_ph1;
   vector<Double_t> v_s1xrate_ph1;
   vector<Double_t> v_s2xrate_ph1;
   vector<Double_t> v_s1yrate_ph1;
   vector<Double_t> v_s2yrate_ph1;
   vector<Double_t> v_treff_ph1;
   vector<Double_t> v_s1xrate_err_ph1;
   vector<Double_t> v_s2xrate_err_ph1;
   vector<Double_t> v_s1yrate_err_ph1;
   vector<Double_t> v_s2yrate_err_ph1;
   vector<Double_t> v_treff_err_ph1;
   //
   vector<Double_t> v_avecur_ph4;
   vector<Double_t> v_avecur_err_ph4;
    vector<Double_t> v_trig1rate_ph4;
    vector<Double_t> v_trig1rate_err_ph4;
  vector<Double_t> v_s1xrate_ph4;
   vector<Double_t> v_s2xrate_ph4;
   vector<Double_t> v_s1yrate_ph4;
   vector<Double_t> v_s2yrate_ph4;
   vector<Double_t> v_treff_ph4;
   vector<Double_t> v_s1xrate_err_ph4;
   vector<Double_t> v_s2xrate_err_ph4;
   vector<Double_t> v_s1yrate_err_ph4;
   vector<Double_t> v_s2yrate_err_ph4;
   vector<Double_t> v_treff_err_ph4;
   Double_t yield_norm_ph4 = 2.075;
   Double_t yield_norm_ph2 = 25.03;
   Double_t yield_norm_ph1 = 71.0;
   for (Int_t n=0;n<nlines;n++) {
     T->GetEntry(n);
     //phase 4
     if (nrun==7570||nrun==7571||nrun==7572||nrun==7574) { 
      v_avecur_ph4.push_back(avecur);
     v_avecur_err_ph4.push_back(.001);
      v_datayield_ph4.push_back(datayield/yield_norm_ph4);
     v_datayield_err_ph4.push_back(datayielderr/yield_norm_ph4);
     v_trig1rate_ph4.push_back(trig1rate);
     v_trig1rate_err_ph4.push_back(.001);
     v_s1xrate_ph4.push_back(s1xrate);
     v_s2xrate_ph4.push_back(s2xrate);
     v_s1yrate_ph4.push_back(s1yrate);
     v_s2yrate_ph4.push_back(s2yrate);
     v_treff_ph4.push_back(treff);
     v_s1xrate_err_ph4.push_back(.001);
     v_s2xrate_err_ph4.push_back(.001);
     v_s1yrate_err_ph4.push_back(.001);
     v_s2yrate_err_ph4.push_back(.001);
     v_treff_err_ph4.push_back(.001);
     }
     //phase 2
     if (nrun==7500||nrun==7499||nrun==7501) {
      v_avecur_ph2.push_back(avecur);
       v_datayield_ph2.push_back(datayield/yield_norm_ph2);
       v_datayield_err_ph2.push_back(datayielderr/yield_norm_ph2);
    v_avecur_err_ph2.push_back(.001);
      v_trig1rate_ph2.push_back(trig1rate);
      v_trig1rate_err_ph2.push_back(.001);
   v_s1xrate_ph2.push_back(s1xrate);
     v_s2xrate_ph2.push_back(s2xrate);
     v_s1yrate_ph2.push_back(s1yrate);
     v_s2yrate_ph2.push_back(s2yrate);
     v_treff_ph2.push_back(treff);
     v_s1xrate_err_ph2.push_back(.001);
     v_s2xrate_err_ph2.push_back(.001);
     v_s1yrate_err_ph2.push_back(.001);
     v_s2yrate_err_ph2.push_back(.001);
     v_treff_err_ph2.push_back(.001);
    }
     // phase 1
     if (nrun==7420||nrun==7419) {
      v_avecur_ph1.push_back(avecur);
     v_avecur_err_ph1.push_back(.001);
      v_datayield_ph1.push_back(datayield/yield_norm_ph1);
     v_datayield_err_ph1.push_back(datayielderr/yield_norm_ph1);
     v_trig1rate_ph1.push_back(trig1rate);
      v_trig1rate_err_ph1.push_back(.001);
    v_s1xrate_ph1.push_back(s1xrate);
     v_s2xrate_ph1.push_back(s2xrate);
     v_s1yrate_ph1.push_back(s1yrate);
     v_s2yrate_ph1.push_back(s2yrate);
     v_treff_ph1.push_back(treff);
     v_s1xrate_err_ph1.push_back(.001);
     v_s2xrate_err_ph1.push_back(.001);
     v_s1yrate_err_ph1.push_back(.001);
     v_s2yrate_err_ph1.push_back(.001);
     v_treff_err_ph1.push_back(.001);
    }
   }
   // plot yields
  TMultiGraph *mgt2 = new TMultiGraph();
  TLegend* legt2 = new TLegend(0.2,0.75,0.48,0.89);
   TGraphErrors *gr_ph1;
   gr_ph1 = new TGraphErrors(v_s1xrate_ph1.size(),&(v_s1xrate_ph1[0]),&(v_treff_ph1[0]),0,&(v_treff_err_ph1[0]));
     gr_ph1->SetMarkerStyle(20);
     gr_ph1->SetMarkerColor(2);
     mgt2->Add(gr_ph1);
   TGraphErrors *gr_ph2;
   gr_ph2 = new TGraphErrors(v_s1xrate_ph2.size(),&(v_s1xrate_ph2[0]),&(v_treff_ph2[0]),0,&(v_treff_err_ph2[0]));
     gr_ph2->SetMarkerStyle(21);
     gr_ph2->SetMarkerColor(3);
     mgt2->Add(gr_ph2);
   TGraphErrors *gr_ph4;
   gr_ph4 = new TGraphErrors(v_s1xrate_ph4.size(),&(v_s1xrate_ph4[0]),&(v_treff_ph4[0]),0,&(v_treff_err_ph4[0]));
     gr_ph4->SetMarkerStyle(22);
     gr_ph4->SetMarkerColor(4);
     mgt2->Add(gr_ph4);
  legt2->SetBorderSize(1);
  legt2->SetFillColor(0);
   legt2->AddEntry(gr_ph1," Phase 1 ","p");
   legt2->AddEntry(gr_ph2," Phase 2 ","p");
   legt2->AddEntry(gr_ph4," Phase 4 ","p");
  TCanvas *ctreff = new TCanvas("ctreff","ctreff",700,700);
   ctreff->Divide(1,1);
   ctreff->cd(1);
   mgt2->Draw("AP");
   mgt2->SetTitle("Jpsi ; S1X rate; SHMS Track eff ");
   mgt2->SetMaximum(1.2);
   mgt2->SetMinimum(.8);
     legt2->Draw();
  ctreff->Print(outputpdf+"(");
   //
   // plot yields
  TMultiGraph *mgt3 = new TMultiGraph();
  TLegend* legt3 = new TLegend(0.2,0.75,0.48,0.89);
   TGraphErrors *gr_y_ph1;
   gr_y_ph1 = new TGraphErrors(v_s1xrate_ph1.size(),&(v_s1xrate_ph1[0]),&(v_datayield_ph1[0]),0,&(v_datayield_err_ph1[0]));
     gr_y_ph1->SetMarkerStyle(20);
     gr_y_ph1->SetMarkerColor(2);
     mgt3->Add(gr_y_ph1);
   TGraphErrors *gr_y_ph2;
   gr_y_ph2 = new TGraphErrors(v_s1xrate_ph2.size(),&(v_s1xrate_ph2[0]),&(v_datayield_ph2[0]),0,&(v_datayield_err_ph2[0]));
     gr_y_ph2->SetMarkerStyle(21);
     gr_y_ph2->SetMarkerColor(3);
     mgt3->Add(gr_y_ph2);
   TGraphErrors *gr_y_ph4;
   gr_y_ph4 = new TGraphErrors(v_s1xrate_ph4.size(),&(v_s1xrate_ph4[0]),&(v_datayield_ph4[0]),0,&(v_datayield_err_ph4[0]));
     gr_y_ph4->SetMarkerStyle(22);
     gr_y_ph4->SetMarkerColor(4);
     mgt3->Add(gr_y_ph4);
  legt3->SetBorderSize(1);
  legt3->SetFillColor(0);
   legt3->AddEntry(gr_y_ph1," Phase 1 ","p");
   legt3->AddEntry(gr_y_ph2," Phase 2 ","p");
   legt3->AddEntry(gr_y_ph4," Phase 4 ","p");
  TCanvas *cdatayield = new TCanvas("cdatayield","cdatayield",700,700);
   cdatayield->Divide(1,1);
   cdatayield->cd(1);
   mgt3->Draw("AP");
   mgt3->SetTitle("Jpsi ; S1X rate; Normalized yield ");
   mgt3->SetMaximum(1.2);
   mgt3->SetMinimum(.8);
     legt3->Draw();
  cdatayield->Print(outputpdf);
   //
     // plot yields
  TMultiGraph *mgt1 = new TMultiGraph();
  TLegend* legt1 = new TLegend(0.2,0.75,0.48,0.89);
   TGraphErrors *gr_trig1_ph1;
   gr_trig1_ph1 = new TGraphErrors(v_avecur_ph1.size(),&(v_avecur_ph1[0]),&(v_trig1rate_ph1[0]),0,&(v_trig1rate_err_ph1[0]));
     gr_trig1_ph1->SetMarkerStyle(20);
     gr_trig1_ph1->SetMarkerColor(2);
     mgt1->Add(gr_trig1_ph1);
   TGraphErrors *gr_trig1_ph2;
   gr_trig1_ph2 = new TGraphErrors(v_avecur_ph2.size(),&(v_avecur_ph2[0]),&(v_trig1rate_ph2[0]),0,&(v_trig1rate_err_ph2[0]));
     gr_trig1_ph2->SetMarkerStyle(21);
     gr_trig1_ph2->SetMarkerColor(3);
     mgt1->Add(gr_trig1_ph2);
   TGraphErrors *gr_trig1_ph4;
   gr_trig1_ph4 = new TGraphErrors(v_avecur_ph4.size(),&(v_avecur_ph4[0]),&(v_trig1rate_ph4[0]),0,&(v_trig1rate_err_ph4[0]));
     gr_trig1_ph4->SetMarkerStyle(22);
     gr_trig1_ph4->SetMarkerColor(4);
     mgt1->Add(gr_trig1_ph4);
  legt1->SetBorderSize(1);
  legt1->SetFillColor(0);
   legt1->AddEntry(gr_trig1_ph1," Phase 1 ","p");
   legt1->AddEntry(gr_trig1_ph2," Phase 2 ","p");
   legt1->AddEntry(gr_trig1_ph4," Phase 4 ","p");
  TCanvas *ctrig1 = new TCanvas("ctrig1","ctrig1",700,700);
   ctrig1->Divide(1,1);
   ctrig1->cd(1);
   mgt1->Draw("AP");
   mgt1->SetTitle("Jpsi ; Current; SHMS Trig1 rate ");
   // mgt1->SetMaximum(1.2);
   // mgt1->SetMinimum(.8);
     legt1->Draw();
  ctrig1->Print(outputpdf+")");
   //
   }
