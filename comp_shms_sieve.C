#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TCutG.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TPolyLine.h>
#include <TGraph.h>
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

void comp_shms_sieve(TString basename) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/"+basename+"_comp_shms_sieve.pdf";
const Int_t nftot=4;  
   TFile *fhistroot[nftot];
  TString inputroot;
  //  TString fname[nftot]={"shms_replay_matrixopt_1703_newmatrix","shms_replay_matrixopt_1703_newmatrix_orig_survey"};
  //  TString fname[nftot]={"shms_replay_matrixopt_1778_newmatrix","shms_replay_matrixopt_1796_newmatrix"};
  TString fname[nftot]={"shms_replay_matrixopt_1703_newmatrix","shms_replay_matrixopt_1778_newmatrix","shms_replay_matrixopt_1796_newmatrix","shms_replay_matrixopt_2258_newmatrix"};
  TString nrun[nftot]={"1703 ","1778 ","1796 ","2258 "};
const Int_t nplots=2;  
 TString hname[nplots]={"hys_xscent_foil_1","hxs_ys_foil_1"};
   TH2F *fhist[nftot][nplots];
    for (Int_t nf=0;nf<nftot;nf++) {  
  inputroot="hist/"+fname[nf]+"_hist.root";
  fhistroot[nf] =  new TFile(inputroot);
    for (Int_t nh=0;nh<nplots;nh++) {  

     fhist[nf][nh] = (TH2F*)fhistroot[nf]->Get(hname[nh]);
    }
    }
     //
   //
     TCanvas *ctime2[nplots];
    for (Int_t nh=0;nh<nplots;nh++) {  
     ctime2[nh] = new TCanvas(Form("ctime2_%d",nh),hname[nh], 700,700);
     ctime2[nh]->Divide(2,2);
    for (Int_t nf=0;nf<nftot;nf++) {  
     ctime2[nh]->cd(nf+1);
          gPad->SetLogz();
     fhist[nf][nh]->Draw("colz");
     fhist[nf][nh]->SetMinimum(5);
    }
     if (nh==0) ctime2[nh]->Print(outputpdf+"(");
     if (nh!=0&&nh!=nplots-1) ctime2[nh]->Print(outputpdf);
     if (nh==nplots-1) ctime2[nh]->Print(outputpdf+")");
    }

    //
}
