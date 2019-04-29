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

void plot_comp_data() {
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.2,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.16);
   TString outputpdf;
   outputpdf= "plots/lumi_jpsi_phase4.pdf";
   TTree *T = new TTree("ntuple","data");
   Double_t nrun;
   Double_t curthres,ch,avecur,ptime;
   Double_t trig1rate,trig2rate,trig3rate,trig4rate,trig1yield,trig1yielderr,trig2yield,trig2yielderr,trig3yield,trig3yielderr,trig4yield,trig4yielderr;
      Long64_t nlines = T->ReadFile("LUMI_STUDY/lumi_jpsi_phase4.dat","nrun/D:ch/D:curthres:avecur:ptime:trig1rate:trig2rate:trig3rate:trig4rate:trig1yield:trig1yielderr:trig2yield:trig2yielderr:trig3yield:trig3yielderr:trig4yield:trig4yielderr");
   printf(" found %lld points\n",nlines);
    T->SetBranchAddress("nrun",&nrun);
    T->SetBranchAddress("avecur",&avecur);
    T->SetBranchAddress("curthres",&curthres);
   T->SetBranchAddress("trig2rate",&trig2rate);
   T->SetBranchAddress("trig3rate",&trig3rate);
   T->SetBranchAddress("trig2yield",&trig2yield);
   T->SetBranchAddress("trig3yield",&trig3yield);
   T->SetBranchAddress("trig2yielderr",&trig2yielderr);
   T->SetBranchAddress("trig3yielderr",&trig3yielderr);
   vector<Double_t> v_avecur_low;
   vector<Double_t> v_avecur_low_err;
   vector<Double_t> v_avecur;
   vector<Double_t> v_avecur_err;
   vector<Double_t> v_trig2yield;
   vector<Double_t> v_trig2yield_err;
   vector<Double_t> v_trig3yield;
   vector<Double_t> v_trig3yield_err;
   vector<Double_t> v_trig2yield_low;
   vector<Double_t> v_trig2yield_low_err;
   vector<Double_t> v_trig3yield_low;
   vector<Double_t> v_trig3yield_low_err;
   Double_t max_trig2=0;
   Double_t max_trig3=0;
   for (Int_t n=0;n<nlines;n++) {
     T->GetEntry(n);
     cout << nrun << " " << curthres << " " << avecur  << " " << trig2yield<< " " << trig2yielderr<< " " << trig3yield<< " " << trig3yielderr << endl;
    if (curthres <= 1.6) {
     v_avecur_low.push_back(avecur);
     v_avecur_low_err.push_back(.001);
     v_trig2yield_low.push_back(trig2yield);
     v_trig2yield_low_err.push_back(trig2yielderr);
     v_trig3yield_low.push_back(trig3yield);
     v_trig3yield_low_err.push_back(trig3yielderr);
     } else {
      v_avecur.push_back(avecur);
     v_avecur_err.push_back(.001);
     v_trig2yield.push_back(trig2yield);
     v_trig2yield_err.push_back(trig2yielderr);
     v_trig3yield.push_back(trig3yield);
     v_trig3yield_err.push_back(trig3yielderr);
     if (trig2yield>max_trig2) max_trig2=trig2yield;
     if (trig3yield>max_trig3) max_trig3=trig3yield;
       }
    }
   // plot yields
  TMultiGraph *mgt2 = new TMultiGraph();
  TLegend* legt2 = new TLegend(0.2,0.75,0.48,0.89);
   TGraphErrors *gr_trig2_low;
   gr_trig2_low = new TGraphErrors(v_avecur_low.size(),&(v_avecur_low[0]),&(v_trig2yield_low[0]),0,&(v_trig2yield_low_err[0]));
     gr_trig2_low->SetMarkerStyle(24);
     gr_trig2_low->SetMarkerColor(1);
     mgt2->Add(gr_trig2_low);
   TGraphErrors *gr_trig2;
   gr_trig2 = new TGraphErrors(v_avecur.size(),&(v_avecur[0]),&(v_trig2yield[0]),0,&(v_trig2yield_err[0]));
     gr_trig2->SetMarkerStyle(20);
     gr_trig2->SetMarkerColor(2);
     mgt2->Add(gr_trig2);
  legt2->SetBorderSize(1);
  legt2->SetFillColor(0);
   legt2->AddEntry(gr_trig2_low," SHMS ELCLEAN I > 1.5","p");
  legt2->AddEntry(gr_trig2," SHMS ELCLEAN I > near max","p");
  TCanvas *ctrig2 = new TCanvas("ctrig2"," trig2",700,700);
   ctrig2->Divide(1,1);
   ctrig2->cd(1);
   mgt2->Draw("AP");
   mgt2->SetTitle("Jpsi ; Current; TRIG2 SHMS Yield ");
   mgt2->SetMaximum(max_trig2*1.05);
   mgt2->SetMinimum(max_trig2*0.8);
     legt2->Draw();
  ctrig2->Print(outputpdf+"(");
    // plot yields
  TMultiGraph *mgt3 = new TMultiGraph();
  TLegend* legt3 = new TLegend(0.2,0.75,0.48,0.89);
   TGraphErrors *gr_trig3_low;
   gr_trig3_low = new TGraphErrors(v_avecur_low.size(),&(v_avecur_low[0]),&(v_trig3yield_low[0]),0,&(v_trig3yield_low_err[0]));
     gr_trig3_low->SetMarkerStyle(26);
     gr_trig3_low->SetMarkerColor(1);
     mgt3->Add(gr_trig3_low);
   TGraphErrors *gr_trig3;
   gr_trig3 = new TGraphErrors(v_avecur.size(),&(v_avecur[0]),&(v_trig3yield[0]),0,&(v_trig3yield_err[0]));
     gr_trig3->SetMarkerStyle(22);
     gr_trig3->SetMarkerColor(2);
     mgt3->Add(gr_trig3);
  legt3->SetBorderSize(1);
  legt3->SetFillColor(0);
   legt3->AddEntry(gr_trig3_low," HMS ELCLEAN I > 1.5","p");
  legt3->AddEntry(gr_trig3," HMS ELCLEAN I > near max","p");
  TCanvas *ctrig3 = new TCanvas("ctrig3"," trig3",700,700);
   ctrig3->Divide(1,1);
   ctrig3->cd(1);
   mgt3->Draw("AP");
   mgt3->SetTitle("Jpsi; Current; TRIG3 HMS Yield ");
   mgt3->SetMaximum(max_trig3*1.05);
   mgt3->SetMinimum(max_trig3*0.8);
    legt3->Draw();
  ctrig3->Print(outputpdf);
  //
   cout << v_avecur_low.size() << " " << v_avecur.size() << endl;
   if (v_avecur_low.size() >0) {
   Double_t norm2_low=v_trig2yield_low[0];
   Double_t norm2_low_err=v_trig2yield_low_err[0]/v_trig2yield_low[0];
   Double_t norm3_low_err=v_trig3yield_low_err[0]/v_trig3yield_low[0];
    Double_t norm3_low=v_trig3yield_low[0];
   for (int j=0;j<v_avecur_low.size();j++) {
     v_trig2yield_low_err[j]=v_trig2yield_low[j]/norm2_low*TMath::Sqrt((v_trig2yield_low_err[j]/v_trig2yield_low[j])*(v_trig2yield_low_err[j]/v_trig2yield_low[j])+norm2_low_err*norm2_low_err);
     if (j==0) v_trig2yield_low_err[j]=0.000001;
     v_trig2yield_low[j]=v_trig2yield_low[j]/norm2_low;
     v_trig3yield_low_err[j]=v_trig3yield_low[j]/norm3_low*TMath::Sqrt((v_trig3yield_low_err[j]/v_trig3yield_low[j])*(v_trig3yield_low_err[j]/v_trig3yield_low[j])+norm3_low_err*norm3_low_err);
    v_trig3yield_low[j]=v_trig3yield_low[j]/norm3_low;
     if (j==0) v_trig3yield_low_err[j]=0.000001;
     cout << v_trig2yield_low[j]<< " " << v_trig2yield_low_err[j]<< " " << v_trig3yield_low[j]<< " " << v_trig3yield_low_err[j] << endl;
   }
   }
   Double_t norm2=v_trig2yield[0];
    Double_t norm3=v_trig3yield[0];
   Double_t norm2_err=v_trig2yield_err[0]/v_trig2yield[0];
   Double_t norm3_err=v_trig3yield_err[0]/v_trig3yield[0];
  for (int j=0;j<v_avecur.size();j++) {
      v_trig2yield_err[j]=v_trig2yield[j]/norm2*TMath::Sqrt((v_trig2yield_err[j]/v_trig2yield[j])*(v_trig2yield_err[j]/v_trig2yield[j])+norm2_err*norm2_err);
     if (j==0) v_trig2yield_err[j]=0.000001;
      v_trig3yield_err[j]=v_trig3yield[j]/norm3*TMath::Sqrt((v_trig3yield_err[j]/v_trig3yield[j])*(v_trig3yield_err[j]/v_trig3yield[j])+norm3_err*norm3_err);
     if (j==0) v_trig3yield_err[j]=0.000001;
   v_trig2yield[j]=v_trig2yield[j]/norm2;
     v_trig3yield[j]=v_trig3yield[j]/norm3;
     cout << v_trig2yield[j]<< " " << v_trig2yield_err[j]<< " " << v_trig3yield[j]<< " " << v_trig3yield_err[j] << endl;
   }
  //
 
  TMultiGraph *mg1 = new TMultiGraph();
  TLegend* leg1 = new TLegend(0.2,0.75,0.48,0.89);
  if (v_avecur_low.size()>0) {
     TGraphErrors *gr_rat_trig2_low;
   gr_rat_trig2_low = new TGraphErrors(v_avecur_low.size(),&(v_avecur_low[0]),&(v_trig2yield_low[0]),0,&(v_trig2yield_low_err[0]));
     gr_rat_trig2_low->SetMarkerStyle(24);
     gr_rat_trig2_low->SetMarkerColor(1);
     mg1->Add(gr_rat_trig2_low);
  leg1->AddEntry(gr_rat_trig2_low," SHMS ELCLEAN I > 1.5","p");
     
   TGraphErrors *gr_rat_trig3_low;
   gr_rat_trig3_low = new TGraphErrors(v_avecur_low.size(),&(v_avecur_low[0]),&(v_trig3yield_low[0]),0,&(v_trig3yield_low_err[0]));
     gr_rat_trig3_low->SetMarkerStyle(26);
     gr_rat_trig3_low->SetMarkerColor(1);
     mg1->Add(gr_rat_trig3_low);
  leg1->AddEntry(gr_rat_trig3_low,"  HMS ELCLEAN I > 1.5","p");
  }
     TGraphErrors *gr_rat_trig2;
   gr_rat_trig2 = new TGraphErrors(v_avecur.size(),&(v_avecur[0]),&(v_trig2yield[0]),0,&(v_trig2yield_err[0]));
     gr_rat_trig2->SetMarkerStyle(20);
     gr_rat_trig2->SetMarkerColor(2);
     mg1->Add(gr_rat_trig2);
   TGraphErrors *gr_rat_trig3;
   gr_rat_trig3 = new TGraphErrors(v_avecur.size(),&(v_avecur[0]),&(v_trig3yield[0]),0,&(v_trig3yield_err[0]));
     gr_rat_trig3->SetMarkerStyle(22);
     gr_rat_trig3->SetMarkerColor(2);
     mg1->Add(gr_rat_trig3);
  leg1->SetBorderSize(1);
  leg1->SetFillColor(0);
  leg1->AddEntry(gr_rat_trig2," SHMS ELCLEAN I > near max","p");
  leg1->AddEntry(gr_rat_trig3,"  HMS ELCLEAN I > near max","p");
  TCanvas *cplot = new TCanvas("cplot"," ratio",700,700);
   cplot->Divide(1,1);
   cplot->cd(1);
   mg1->Draw("AP");
   mg1->SetTitle("Jpsi ; Current; Relative Yield ");
   mg1->SetMaximum(1.025);
   mg1->SetMinimum(0.95);
    leg1->Draw();
  cplot->Print(outputpdf+")");
}
