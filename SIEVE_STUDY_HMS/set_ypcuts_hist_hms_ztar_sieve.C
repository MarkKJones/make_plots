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

void set_ypcuts_hist_hms_ztar_sieve(TString basename) {
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
 outputpdf="plots/"+basename+"_set_ypcuts_hms_ztar_sieve_hist.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_hms_ztar_sieve_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   TH1F *hyptar_cent_foil= (TH1F*)fhistroot->Get("hyptar_foil_1");
   if (!hyptar_cent_foil) return;
   //
   TCanvas *cyp = new TCanvas("cyp","Yp center foil",700,700);
   cyp->Divide(1,1);
   cyp->cd(1);
   hyptar_cent_foil->Draw();
	  TCutG *cutg;     
     TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
      gPad->Update();
      Double_t xlo[7],xhi[7];
    if (!tempg) cout << " no cut" << endl;
    if (tempg)		  {
          cutg=(TCutG*)(tempg->Clone());
      		  cutg->SetName("cuttemp");
     		  //cutg->Print();
		  Int_t npts=cutg->GetN();
		  Double_t ydum;
		  Int_t nc=0;
                  for (int n=0;n<14;n++) {
		  cutg->GetPoint(n++,xlo[nc],ydum);
		  cutg->GetPoint(n,xhi[nc++],ydum);
		  cout << nc-1 << " xlo " << xlo[nc-1] << " xhi " << xhi[nc-1] << endl;
		  }
     }
		  cout << endl;
                 for (int n=0;n<7;n++) {
		  cout << xlo[n] << "," ;
		  }
		  cout << endl;
                 for (int n=0;n<7;n++) {
		  cout << xhi[n] << "," ;
		  }
		  cout << endl;
                   
     

   //
}
