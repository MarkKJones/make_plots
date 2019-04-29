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
 inputroot = "hist/"+fn1+"_shms_ngcer_hist.root";
 cout << " root file = " << inputroot << endl;
  //
 //
  //
   fhistroot =  new TFile(inputroot);
   TH2F *h_calXY ; 
   TH2F *h_calXY_now; 
   TH2F *h_calXY_rat; 
     h_calXY = (TH2F*)fhistroot->Get("h_calXY");
     h_calXY_now = (TH2F*)fhistroot->Get("h_calXY_now");
      h_calXY_rat = (TH2F*)h_calXY->Clone();
     h_calXY_rat->Divide(h_calXY,h_calXY_now,1,1,"B");
   TH2F *h_cerXY ; 
   TH2F *h_cerXY_now; 
   TH2F *h_cerXY_rat; 
     h_cerXY = (TH2F*)fhistroot->Get("h_cerXY");
     h_cerXY_now = (TH2F*)fhistroot->Get("h_cerXY_now");
      h_cerXY_rat = (TH2F*)h_cerXY->Clone();
     h_cerXY_rat->Divide(h_cerXY,h_cerXY_now,1,1,"B");
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
    TCanvas *ccerXY= new TCanvas("ccerXY","2dcer", 1000,700);
     ccerXY->Divide(1,1);
     ccerXY->cd();
      h_cerXY_rat->Draw("colz");
      h_cerXY_rat->SetMaximum(30);  
      h_cerXY_rat->SetMinimum(0);  
   ccerXY->Print(outputpdf+")");

   //

}
