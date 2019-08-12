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
void comp_shms_cal( TString fn1 ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/comp_shms_cal_"+fn1+".pdf";
    //
   TFile *fhistroot;
 //
 //
 TString inputroot;
 inputroot = "hist/"+fn1+"_shms_cal_hist.root";
 cout << " root file = " << inputroot << endl;
  //
 Double_t pcal_arr_gain_cor[224]={0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                     0.00,  0.00,  0.00, 34.51, 50.03, 66.05, 31.82,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                     0.00,  0.00, 35.50, 54.70, 38.62, 27.51, 78.94, 60.79,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                     0.00,  0.00, 31.04, 33.74, 65.11, 51.95, 43.52, 33.58, 38.10, 78.75,  0.00, 36.19, 38.67, 69.22,  0.00,  0.00,
                     0.00,  0.00, 31.22, 82.58, 19.73, 39.89, 40.95, 56.92, 82.86, 41.83, 63.32, 52.60, 74.48, 15.85, 17.64,  0.00,
                     0.00,  0.00, 55.28, 54.67, 26.83, 44.92, 56.34, 48.71, 45.56, 41.57, 37.31, 23.68, 38.52, 35.82,  6.31,  0.00,
                     0.00,  0.00, 21.84, 17.69, 19.08, 14.41, 21.67, 22.10, 11.43, 23.99, 26.04, 21.50, 29.16, 16.52, 10.21,  0.00,
                     0.00,  0.00,  0.00, 40.24, 55.30, 47.48, 27.91, 46.77, 61.59, 70.41, 24.73, 65.84, 27.59, 25.26, 23.23,  0.00,
                     0.00,  0.00, 18.22, 48.76, 26.97, 62.92, 49.71, 59.54, 32.70,  0.00, 38.02, 38.81, 54.45, 32.43, 14.69,  0.00,
                     0.00,  0.00,  0.00, 47.47, 35.80, 44.88, 23.16, 31.24,  0.00,  0.00,  0.00,  0.00, 46.67, 25.56,  0.00,  0.00,
                     0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                     0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
				  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00};
 //
  //
   fhistroot =  new TFile(inputroot);
   TH1F *hetotnorm_npecut; 
   hetotnorm_npecut = (TH1F*)fhistroot->Get("hetotnorm_npecut");
   TH2F *h_calXY ; 
   TH2F *h_calXY_now; 
   TH2F *h_calXY_rat; 
     h_calXY = (TH2F*)fhistroot->Get("h_calXY");
     h_calXY_now = (TH2F*)fhistroot->Get("h_calXY_now");
      h_calXY_rat = (TH2F*)h_calXY->Clone();
     h_calXY_rat->Divide(h_calXY,h_calXY_now,1,1,"B");
   TH2F *h_calXY_le ; 
   TH2F *h_calXY_le_now; 
   TH2F *h_calXY_le_rat; 
     h_calXY_le = (TH2F*)fhistroot->Get("h_calXY_le");
     h_calXY_le_now = (TH2F*)fhistroot->Get("h_calXY_le_now");
      h_calXY_le_rat = (TH2F*)h_calXY_le->Clone();
     h_calXY_le_rat->Divide(h_calXY_le,h_calXY_le_now,1,1,"B");
   TH2F *h_calXY_he ; 
   TH2F *h_calXY_he_now; 
   TH2F *h_calXY_he_rat; 
     h_calXY_he = (TH2F*)fhistroot->Get("h_calXY_he");
     h_calXY_he_now = (TH2F*)fhistroot->Get("h_calXY_he_now");
      h_calXY_he_rat = (TH2F*)h_calXY_he->Clone();
     h_calXY_he_rat->Divide(h_calXY_he,h_calXY_he_now,1,1,"B");
   //
    TCanvas *ccalXY= new TCanvas("ccalXY","2damp", 1000,700);
     ccalXY->Divide(1,1);
     ccalXY->cd();
      h_calXY_rat->Draw("colz");
      h_calXY_rat->SetMaximum(2);  
      h_calXY_rat->SetMinimum(0);  
   Double_t Xpos,Ypos,XStep=9.,YStep=9.;
  TText* coltext= new TText(-75,77,"Col");
  coltext->Draw();
    coltext->SetTextAlign(22);
  TText* rowtext= new TText(+70,-75,"Row");
  rowtext->Draw();
    rowtext->SetTextAlign(22);
  for (Int_t nc=0; nc < 14;nc++) {
    Ypos=(14-1)*YStep/2.-YStep*nc;
    TText* coltext= new TText(Ypos,77,Form("%d",nc+1));
    coltext->Draw();
    coltext->SetTextAlign(22);
    for (Int_t nr=0; nr < 16;nr++) {
    Xpos=-(16-1)*XStep/2.+XStep*nr;
     TText* rowtext= new TText(+70,Xpos,Form("%d",nr+1));
    rowtext->Draw();
    rowtext->SetTextAlign(22);
   TBox* mod = new TBox(Ypos-YStep/2.,Xpos-XStep/2.,Ypos+YStep/2.,Xpos+XStep/2.);
    mod->Draw();
    mod->SetLineWidth(2);
    mod->SetFillStyle(0);
  }
  }
   ccalXY->Print(outputpdf+"(");
   //
  //
    TCanvas *ccalXY_he= new TCanvas("ccalXY_he",".9 < E/p < 1.1", 1000,700);
     ccalXY_he->Divide(1,1);
     ccalXY_he->cd();
      h_calXY_he_rat->Draw("colz");
      h_calXY_he_rat->SetMaximum(1.2);  
      h_calXY_he_rat->SetMinimum(0.8);  
      //  Double_t Xpos,Ypos,XStep=9.,YStep=9.;
  coltext= new TText(-75,77,"Col");
  coltext->Draw();
    coltext->SetTextAlign(22);
  rowtext= new TText(+70,-75,"Row");
  rowtext->Draw();
    rowtext->SetTextAlign(22);
  for (Int_t nc=0; nc < 14;nc++) {
    Ypos=(14-1)*YStep/2.-YStep*nc;
    TText* coltext= new TText(Ypos,77,Form("%d",nc+1));
    coltext->Draw();
    coltext->SetTextAlign(22);
    for (Int_t nr=0; nr < 16;nr++) {
    Xpos=-(16-1)*XStep/2.+XStep*nr;
     TText* rowtext= new TText(+70,Xpos,Form("%d",nr+1));
    rowtext->Draw();
    rowtext->SetTextAlign(22);
   TBox* mod = new TBox(Ypos-YStep/2.,Xpos-XStep/2.,Ypos+YStep/2.,Xpos+XStep/2.);
    mod->Draw();
    mod->SetLineWidth(2);
    mod->SetFillStyle(0);
  }
  }
   ccalXY_he->Print(outputpdf);
  //
     TCanvas *ccaletot= new TCanvas("ccaletot","etot", 1000,700);
     ccaletot->Divide(1,1);
     ccaletot->cd();
     hetotnorm_npecut->Draw();
     ccaletot->Print(outputpdf);
  //
    TCanvas *ccalXY_le= new TCanvas("ccalXY_le",".5 < E/p < .9", 1000,700);
     ccalXY_le->Divide(1,1);
     ccalXY_le->cd();
      h_calXY_le_rat->Draw("colz");
      h_calXY_le_rat->SetMaximum(1.);  
      h_calXY_le_rat->SetMinimum(0.4);  
      //  Double_t Xpos,Ypos,XStep=9.,YStep=9.;
  coltext= new TText(-75,77,"Col");
  coltext->Draw();
    coltext->SetTextAlign(22);
  rowtext= new TText(+70,-75,"Row");
  rowtext->Draw();
    rowtext->SetTextAlign(22);
  for (Int_t nc=0; nc < 14;nc++) {
    Ypos=(14-1)*YStep/2.-YStep*nc;
    TText* coltext= new TText(Ypos,77,Form("%d",nc+1));
    coltext->Draw();
    coltext->SetTextAlign(22);
    for (Int_t nr=0; nr < 16;nr++) {
    Xpos=-(16-1)*XStep/2.+XStep*nr;
     TText* rowtext= new TText(+70,Xpos,Form("%d",nr+1));
    rowtext->Draw();
    rowtext->SetTextAlign(22);
   TBox* mod = new TBox(Ypos-YStep/2.,Xpos-XStep/2.,Ypos+YStep/2.,Xpos+XStep/2.);
    mod->Draw();
    mod->SetLineWidth(2);
    mod->SetFillStyle(0);
  }
  }
   ccalXY_le->Print(outputpdf+")");
 
   //

}
