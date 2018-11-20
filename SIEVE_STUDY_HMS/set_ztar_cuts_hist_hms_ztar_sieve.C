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

void set_ztar_cuts_hist_hms_ztar_sieve(TString basename) {
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
   TH2F *hztar_yptar = (TH2F*)fhistroot->Get("hztar_yptar");
   //
   TCanvas *cyp = new TCanvas("cyp","hztar_yptar",700,700);
   cyp->Divide(1,1);
   cyp->cd(1);
   hztar_yptar->Draw();
	  TCutG *cutg;     
     TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
      gPad->Update();
      Double_t xlo[4],xhi[4];
      Double_t ylo[4],yhi[4];
    if (!tempg) cout << " no cut" << endl;
    if (tempg)		  {
          cutg=(TCutG*)(tempg->Clone());
      		  cutg->SetName("cuttemp");
     		  //cutg->Print();
		  Int_t npts=cutg->GetN();
		  Double_t ydum;
		  Int_t nc=0;
		  cutg->GetPoint(0,xlo[0],ylo[0]);
		  cutg->GetPoint(1,xlo[1],ylo[1]);
		  cutg->GetPoint(2,xhi[1],yhi[1]);
		  cutg->GetPoint(3,xhi[0],yhi[0]);
                  cout << " slope hi = " << (yhi[0]-yhi[1])/(xhi[0]-xhi[1]) << " " << yhi[0] -xhi[0]*(yhi[0]-yhi[1])/(xhi[0]-xhi[1]) << endl;

                  cout << " slope lo = " << (ylo[0]-ylo[1])/(xlo[0]-xlo[1]) << " " << ylo[0] -xlo[0]*(ylo[0]-ylo[1])/(xlo[0]-xlo[1]) << endl;

    }
                   
     

   //
}
