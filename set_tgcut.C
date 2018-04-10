#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMarker.h>
#include <TCutG.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
using namespace std;
void set_tgcut(TString basename,Int_t flag=0) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
 //
     TString inputroot;
    inputroot="hist/"+basename+"_hist.root";
      TFile *fhistroot;
     fhistroot =  new TFile(inputroot);
     TString outputcut="hist/"+basename+"_cut.root";
     TFile *fcut= new TFile(outputcut,"UPDATE");
       TH2F *fhist;
      fhist = (TH2F*)fhistroot->Get("hytar_yptar");
      //
   TCanvas *ccompcut; 
   ccompcut= new TCanvas("ccompcut","cut",900,500);
   ccompcut->Divide(1,1);
   ccompcut->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
    gPad->SetLogz();    
      fhist->SetTitle(basename);
      fhist->SetMinimum(10);
       fhist->Draw("colz");
       //	 	  gPad->WaitPrimitive();
        
		  //       TCutG *c1;
		  //		  c1 = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
    
   //
       Int_t nfoils=2;
       //             cout <<"enter number of foils = " << endl;
       //      cin >> nfoils ;
  
       TCutG *cutg[nfoils];
       Int_t nf=0;
       fcut->cd();
       while (nf<nfoils) {
	 	  TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
	 			  TString cutname=Form("ytar_yptar_cut_%d",nf);
		  cutg[nf]=(TCutG*)(tempg->Clone(cutname));
		  cutg[nf]->SetName(cutname);
		  cutg[nf]->Print();
		  cutg[nf]->Write();
                  ccompcut->Update();
		  nf++;
		  cout << " nf = " << nf << endl;
                  //tempg->Delete();
	     }
       //

}
