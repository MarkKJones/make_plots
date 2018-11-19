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
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TLine.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_xs_ys_fit(TString basename=""){
  //
   if (basename=="") {
     cout << " Input the basename of the root file" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="hist/"+basename+".root";
   TString pdffile;
  TFile *fhistroot;
   fhistroot =  new TFile(inputroot);
    //
    TH2F *hxs_ys_old[18];
    TH2F *hxs_ys_new[18];

    for (int nd=0;nd<18;nd++) {
       hxs_ys_old[nd] = (TH2F*)fhistroot->Get(Form("hxs_ys_%d_old",nd));
       hxs_ys_new[nd] = (TH2F*)fhistroot->Get(Form("hxs_ys_%d_new",nd));
    if (!hxs_ys_new[nd] ) cout << " no hist = " << Form("hxs_ys_%d_new",nd) << endl;
     }
   pdffile="plots/"+basename+"_xs_ys.pdf";
      Double_t max_ny;
      TCanvas *cplot3[18];
      TLine *nxline[11];
      TLine *nyline[11];
      for (int nx=0;nx<11;nx++) {
	    nxline[nx] = new TLine(-10,(5-nx)*2.5,10,(5-nx)*2.5);
	    nyline[nx] = new TLine((5-nx)*1.64,-15,(5-nx)*1.64,15);
	  }
    for (int nd=0;nd<18;nd++) {
      cplot3[nd] = new TCanvas(Form("cplot3_%d",nd),Form("Delta run =%d",-12+2*nd), 900,700);
      cplot3[nd]->Divide(2,1);
      cplot3[nd]->cd(1);
      hxs_ys_old[nd]->Draw("colz");
      for (int nx=0;nx<11;nx++) {
       	nxline[nx]->Draw();
       	nyline[nx]->Draw();
      }
      cplot3[nd]->cd(2);
      hxs_ys_new[nd]->Draw("colz");
      for (int nx=0;nx<11;nx++) {
       	nxline[nx]->Draw();
       	nyline[nx]->Draw();
      }
     if (nd==0) cplot3[nd]->Print(pdffile+"(");
    if (nd>0&&nd<17) cplot3[nd]->Print(pdffile);
    if (nd==17) cplot3[nd]->Print(pdffile+")");
    cplot3[nd]->Close();
    }  
    //
    // end bracket
}
