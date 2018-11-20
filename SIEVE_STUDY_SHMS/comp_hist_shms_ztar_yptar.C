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

void comp_hist_shms_ztar_yptar(TString basename="") {
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
 outputpdf="plots/"+basename+"_comp_shms_ztar_yptar_hist.pdf";
 //
 static const Int_t nftot=6;
  TString flab[nftot]={"4749","4751","4743","4756","4735","4741"};
 TString label[nftot]={"Run 4749 P=6.3 Q1=1.","Run 4749 P=6.3 Q1=1.003","Run 4743 P=6.842 Q1=1.","Run 4756 P=6.842 Q1=1.003","Run 4743 P=8.035 Q1=1.","Run 4743 P=8.035 Q1=1.004"};
 //
 /*
 static const Int_t nftot=7;
   TString flab[nftot]={"4763","4761","4760","4758","4749","4743","4735"};
 TString label[nftot]={"Run 4763 P=4.0 Q1=1.0","Run 4761 P=4.5 Q1=1.","Run 4760 P=5.0 Q1=1.","Run 4583 P=5.5 Q1=1.","Run 4749 P=6.3 Q1=1.","Run 4743 P=6.842 Q1=1.","Run 4743 P=8.035 Q1=1."};
 */
  TString inputroot[nftot];
   TFile *fhistroot[nftot];
    TH2F *fhist[nftot];
      TH1F *zhist[nftot][9];
       TH1F *yphist[nftot][9];
       Double_t ztar_cent[nftot][9];
       Double_t yptar_cent[nftot][9];
       Double_t ztar_err[nftot][9];
       Double_t yptar_err[nftot][9];
  for (Int_t nf=0;nf<nftot;nf++) {
     inputroot[nf]="hist/shms_replay_matrixopt_"+flab[nf]+"_-1_shms_ztar_sieve_hist.root";
     cout << " infile root = " << inputroot[nf] << endl;
   fhistroot[nf] =  new TFile(inputroot[nf]);
   fhist[nf] = (TH2F*)fhistroot[nf]->Get("hztar_yptar");
       for (Int_t n=0;n<9;n++) {
	 zhist[nf][n] = (TH1F*)fhistroot[nf]->Get(Form("hztar_%d",n));
	 yphist[nf][n] = (TH1F*)fhistroot[nf]->Get(Form("hyptar_%d",n));
	 ztar_cent[nf][n]=zhist[nf][n]->GetMean();
	 yptar_cent[nf][n]=yphist[nf][n]->GetMean();
	 ztar_err[nf][n]=zhist[nf][n]->GetRMS();
	 yptar_err[nf][n]=yphist[nf][n]->GetRMS();
       }
    }
  //  
     TGraphErrors *gr[nftot];
  TF1 *linfit[nftot];
 for (Int_t nf=0;nf<nftot;nf++) {
     gr[nf] = new TGraphErrors(9,yptar_cent[nf],ztar_cent[nf],yptar_err[nf],ztar_err[nf]);
    linfit[nf] = new TF1("linfit","pol1",-1.,1.);
     gr[nf]->Fit("linfit");    
  }
  //
      TCanvas *c1 = new TCanvas("c1","ztar v yp", 700,1200);
      c1->Divide(2,4);
  for (Int_t nf=0;nf<nftot;nf++) {
    c1->cd(nf+1);
gPad->SetGridx();
gPad->SetGridy();
    fhist[nf]->Draw("colz");
    fhist[nf]->SetTitle(label[nf]);
         linfit[nf]->Draw("same");   
  }     
  //
}
