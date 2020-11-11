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
#include <TText.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
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

Double_t dt_extended(Double_t *x, Double_t *par){
  Double_t dtts = 0;
  Float_t xx=x[0];
  Double_t muts;
  Int_t istart=int(par[0]);
  muts=par[1]*xx*1000;
  for (Int_t i = istart; i < 20; i++) {
    dtts = dtts +  ( TMath::Power(muts,i) * TMath::Exp(-1.0*muts) )/ TMath::Factorial(i);   }
  //  cout << xx << " " << dtts << endl;
  //  dtts=100*dtts;
  dtts=1-dtts;
  return dtts;
}

void plot_f2_lt(Int_t ntar=4) {
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(0);
 gStyle->SetTitleOffset(1.2,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.16);
   TString outputpdf;
   outputpdf= "plots/F2_shms_lt.pdf";
   Bool_t DumSub= kTRUE;
   //
   //const Int_t ntar=3;
   TString DataFile[ntar];
   TString TarLabel[ntar];
   TString TarNames[4]={"1","2","5","17"};
   Double_t PreScale[4] = {1,2,5,17};
   for (Int_t i=0;i<ntar;i++) {
     TarLabel[i]=" Prescal = "+TarNames[i];
      DataFile[i] = "FALL18_LUMI/shms_LT_test_preScale"+TarNames[i]+".data";
   }
   //
   vector<vector<Double_t> >   CompLT;
   vector<vector<Double_t> >   EDTMLT;
   vector<vector<Double_t> >   TrigRate;
   vector<vector<Double_t> >   TrigRateErr;
   //
   CompLT.resize(ntar);
   EDTMLT.resize(ntar);
   TrigRate.resize(ntar);
   TrigRateErr.resize(ntar);
  // nrun curr currcorr clt lt treff ev2rate yield_pid err yield_track err Yield_EL_CLEAN err Yield_EL_REAL err TrigRate
   Int_t Ndata_get[3]={1,2,9};
   Int_t nrun=0;
   for (Int_t nt=0;nt<ntar;nt++) {
     ifstream file_Data(DataFile[nt].Data());
     nrun=0;
      if (file_Data.is_open()) {
	Int_t nr;
	file_Data >> nr;
	
	for (Int_t ny=0;ny<nr;ny++)  {
		  Double_t temp;
		  for (Int_t i=0;i<19;i++) {
		    file_Data >> temp;
		    if (i==1) CompLT[nt].push_back(temp);
		    if (i==2) EDTMLT[nt].push_back(temp);
		    if (i==9) TrigRate[nt].push_back(temp/1000./PreScale[nt]);
		  }
		  nrun++;
	}
	file_Data.close();
      } else {
	cout << DataFile[nt].Data() << endl;
	return;
      }
   }
      //
        TGraph *gCompLT[ntar]; 
       TGraph *gEDTMLT[ntar]; 
     for (Int_t nt=0;nt<ntar;nt++) {
	cout << TarLabel[nt] << " Nrun  = " << TrigRate[nt].size() << endl;
	//
	gCompLT[nt]= new TGraph(TrigRate[nt].size(), &TrigRate[nt][0], &CompLT[nt][0]);
	gCompLT[nt]->SetMarkerStyle(23);
	gCompLT[nt]->SetTitle(TarLabel[nt]);
	gCompLT[nt]->GetXaxis()->SetTitle("Trig Rate (kHz)");
	gCompLT[nt]->GetYaxis()->SetTitle("Computer LT");
	//
	gEDTMLT[nt]= new TGraph(TrigRate[nt].size(), &TrigRate[nt][0], &EDTMLT[nt][0]);
	gEDTMLT[nt]->SetMarkerStyle(23);
	gEDTMLT[nt]->SetTitle(TarLabel[nt]);
	gEDTMLT[nt]->GetXaxis()->SetTitle("Trig Rate (kHz)");
	gEDTMLT[nt]->GetYaxis()->SetTitle("EDTM LT");
	//
	//
	      }
      //
  TF1 *deadt3 = new TF1("deadt3",dt_extended,0,10.,2);
  deadt3->SetParameters(2,152.e-6);
  //  deadt3->Draw("same");
  deadt3->SetLineColor(2);
  TF1 *deadt5 = new TF1("deadt5",dt_extended,0,10.,2);
  deadt5->SetParameters(5,180.e-6);
  //  deadt3->Draw("same");
  deadt5->SetLineColor(2);
       //
       TF1* fLT = new TF1("fLT","1./(1+[0]*(x-[1]))");
       TF1* fLT2 = new TF1("fLT2","1./(1+[0]*(x-[1]))");
       TF1* fLT5 = new TF1("fLT5","1./(1+[0]*(x-[1]))");
       vector<Double_t > CLT_AllRuns;
       vector<Double_t > EDTMLT_AllRuns;
       vector<Double_t > TrigRate_AllRuns;
      TCanvas *cTarRate[2];
	 cTarRate[0] = new TCanvas(Form("CTarRate_%d",0),Form("CTarRatel_%d",0), 700,700);
	 cTarRate[0]->Divide(2,2);
       for (Int_t nt=0;nt<ntar;nt++) {
	 for (Int_t nr=0;nr<TrigRate[nt].size();nr++) {
	   CLT_AllRuns.push_back(CompLT[nt][nr]);
	   EDTMLT_AllRuns.push_back(EDTMLT[nt][nr]);
	   TrigRate_AllRuns.push_back(TrigRate[nt][nr]/PreScale[nt]);
	 }
	 cTarRate[0]->cd(nt+1);
	 gCompLT[nt]->Draw("AP");
	 if (nt==0) gCompLT[nt]->Fit("fLT","","");
	 if (nt==1) deadt3->Draw("same");
	 if (nt==1) gCompLT[nt]->Fit("fLT2","","",0,10);
	 if (nt==2) gCompLT[nt]->Fit("fLT5","","");
	 if (nt==2) deadt5->Draw("same");
       }
       cTarRate[0]->Print(outputpdf);
       /*
       TF1* fLT = new TF1("fLT","1./(1+[0]*(x-[1]))");
       TF1* fEDTMLT = new TF1("fEDTMLT","1./(1+[0]*(x-[1]))");
       fEDTMLT->SetLineColor(3);
	TGraph *gCompLT_AllRuns= new TGraph(TrigRate_AllRuns.size(), &TrigRate_AllRuns[0], &CLT_AllRuns[0]);
	gCompLT_AllRuns->SetMarkerStyle(23);
	gCompLT_AllRuns->SetTitle(Form("All runs"));
	gCompLT_AllRuns->GetXaxis()->SetTitle("Trig Rate (kHz)");
	gCompLT_AllRuns->GetYaxis()->SetTitle("Computer LT");
	TGraph *gEDTMLT_AllRuns= new TGraph(TrigRate_AllRuns.size(), &TrigRate_AllRuns[0], &EDTMLT_AllRuns[0]);
	gEDTMLT_AllRuns->SetMarkerStyle(22);
	gEDTMLT_AllRuns->SetMarkerColor(3);
	gEDTMLT_AllRuns->SetTitle(Form("All runs"));
	gEDTMLT_AllRuns->GetXaxis()->SetTitle("Trig Rate (kHz)");
	gEDTMLT_AllRuns->GetYaxis()->SetTitle("EDTM LT");
      TCanvas *cTarRate_AllRuns;
      cTarRate_AllRuns = new TCanvas("cTarRate_AllRuns","cTarRate_AllRuns",700,700);
      cTarRate_AllRuns->Divide(1,1);
      cTarRate_AllRuns->cd(1);
      gCompLT_AllRuns->Draw("AP");
      gEDTMLT_AllRuns->Draw("P same");
      //      gCompLT_AllRuns->Fit("fLT","","",5,20);
      */
}
