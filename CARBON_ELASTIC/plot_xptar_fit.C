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
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_xptar_fit(TString basename=""){
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
    TH1F *hxptardiffold[18][11][11];
    TH1F *hxptardiffnew[18][11][11];
    for (int nd=0;nd<18;nd++) {
    for (int nnx=0;nnx<11;nnx++) {
    for (int nny=0;nny<11;nny++) {
      hxptardiffnew[nd][nnx][nny] = (TH1F*)fhistroot->Get(Form("hxptardiffnew_%d_%d_%d",nd,nnx,nny));
    if (!hxptardiffnew[nd][nnx][nny] ) cout << " no hist = " << Form("hxptardiffnew_%d_%d_%d",nd,nnx,nny) << endl;
      hxptardiffold[nd][nnx][nny] = (TH1F*)fhistroot->Get(Form("hxptardiffold_%d_%d_%d",nd,nnx,nny));
    }
    }
    }
    //
   Double_t meannew[18][11][11];
   Double_t rmsnew[18][11][11];
   Double_t meanold[18][11][11];
   Double_t rmsold[18][11][11];
   Double_t xmeannew[18][11][11];
   Double_t xrmsnew[18][11][11];
   Double_t xmeanold[18][11][11];
   Double_t xrmsold[18][11][11];
   Double_t y[18][11][11];
   Double_t x[18][11][11];
    for (int nd=0;nd<18;nd++) {
    for (int nnx=0;nnx<11;nnx++) {
    for (int nny=0;nny<11;nny++) {
      if (hxptardiffnew[nd][nnx][nny]->Integral()>50) {
      meannew[nd][nnx][nny]=hxptardiffnew[nd][nnx][nny]->GetMean();
      rmsnew[nd][nnx][nny]=hxptardiffnew[nd][nnx][nny]->GetRMS();
      meanold[nd][nnx][nny]=hxptardiffold[nd][nnx][nny]->GetMean();
      rmsold[nd][nnx][nny]=hxptardiffold[nd][nnx][nny]->GetRMS();
      xmeannew[nd][nny][nnx]=hxptardiffnew[nd][nnx][nny]->GetMean();
      xrmsnew[nd][nny][nnx]=hxptardiffnew[nd][nnx][nny]->GetRMS();
      xmeanold[nd][nny][nnx]=hxptardiffold[nd][nnx][nny]->GetMean();
      xrmsold[nd][nny][nnx]=hxptardiffold[nd][nnx][nny]->GetRMS();
      } else {
      meannew[nd][nnx][nny]=-100.;
      rmsnew[nd][nnx][nny]=0.;
      meanold[nd][nnx][nny]=-100.;
      rmsold[nd][nnx][nny]=0.;
      xmeannew[nd][nny][nnx]=-100.;
      xrmsnew[nd][nny][nnx]=0.;
      xmeanold[nd][nny][nnx]=-100.;
      xrmsold[nd][nny][nnx]=0.;
      }
      y[nd][nnx][nny]=nny;
      x[nd][nny][nnx]=nnx;
    }}}
    //
   pdffile="plots/"+basename+"_xptardiff_nx.pdf";
      TCanvas *cplot[18];
      TGraphErrors *gr[18][11];
      TGraphErrors *grold[18][11];
      TMultiGraph *mgr[18][11];
    for (int nd=0;nd<18;nd++) {
      cplot[nd] = new TCanvas(Form("cplot_%d",nd),Form("Delta run =%d",-12+2*nd), 900,700);
      cplot[nd]->Divide(3,4);
    for (int nnx=0;nnx<11;nnx++) {
      cplot[nd]->cd(nnx+1);
      mgr[nd][nnx] = new TMultiGraph();
      gr[nd][nnx] = new TGraphErrors(11,y[nd][nnx],meannew[nd][nnx],0,rmsnew[nd][nnx]);
      gr[nd][nnx]->SetTitle(Form("Delta run = %d, Sieve Nx = %d; Ny sieve",-12+2*nd,nnx));
     gr[nd][nnx]->SetMarkerStyle(22);
     gr[nd][nnx]->SetMarkerColor(2);
      grold[nd][nnx] = new TGraphErrors(11,y[nd][nnx],meanold[nd][nnx],0,rmsold[nd][nnx]);
      grold[nd][nnx]->SetTitle(Form("Delta run = %d, Sieve Nx = %d; Ny sieve",-12+2*nd,nnx));
     grold[nd][nnx]->SetMarkerStyle(21);
     grold[nd][nnx]->SetMarkerColor(1);
     mgr[nd][nnx]->SetMinimum(-.01);
     mgr[nd][nnx]->SetMaximum(.01);
     mgr[nd][nnx]->Add(gr[nd][nnx]);
     mgr[nd][nnx]->Add(grold[nd][nnx]);
     mgr[nd][nnx]->SetTitle(Form("Delta run= %d, Sieve Nx = %d; Ny sieve; Xptar diff",-12+2*nd,nnx));
     mgr[nd][nnx]->Draw("AP");
    }
    if (nd==0) cplot[nd]->Print(pdffile+"(");
    if (nd>0&&nd<17) cplot[nd]->Print(pdffile);
    if (nd==17) cplot[nd]->Print(pdffile+")");
    cplot[nd]->Close();
    }  
   //
   pdffile="plots/"+basename+"_xptardiff_ny.pdf";
    //
      TCanvas *cplot2[18];
      TGraphErrors *grx[18][11];
      TGraphErrors *grxold[18][11];
      TMultiGraph *mgrx[18][11];
    for (int nd=0;nd<18;nd++) {
      cplot2[nd] = new TCanvas(Form("cplot2_%d",nd),Form("Delta run =%d",-12+2*nd), 900,700);
      cplot2[nd]->Divide(3,4);
    for (int nny=0;nny<11;nny++) {
      cplot2[nd]->cd(nny+1);
      mgrx[nd][nny] = new TMultiGraph();
      grx[nd][nny] = new TGraphErrors(11,x[nd][nny],xmeannew[nd][nny],0,xrmsnew[nd][nny]);
      grx[nd][nny]->SetTitle(Form("Delta run = %d, Sieve Ny = %d; Nx sieve",-12+2*nd,nny));
     grx[nd][nny]->SetMarkerStyle(22);
     grx[nd][nny]->SetMarkerColor(2);
      grxold[nd][nny] = new TGraphErrors(11,x[nd][nny],xmeanold[nd][nny],0,xrmsold[nd][nny]);
      grxold[nd][nny]->SetTitle(Form("Delta run = %d, Sieve Ny = %d; Nx sieve",-12+2*nd,nny));
     grxold[nd][nny]->SetMarkerStyle(21);
     grxold[nd][nny]->SetMarkerColor(1);
     mgrx[nd][nny]->SetMinimum(-.02);
     mgrx[nd][nny]->SetMaximum(.02);
     mgrx[nd][nny]->Add(grx[nd][nny]);
     mgrx[nd][nny]->Add(grxold[nd][nny]);
     mgrx[nd][nny]->SetTitle(Form("Delta run = %d, Sieve Ny = %d; Nx sieve; Xptar diff",-12+2*nd,nny));
     mgrx[nd][nny]->Draw("AP");
    }
    if (nd==0) cplot2[nd]->Print(pdffile+"(");
    if (nd>0&&nd<17) cplot2[nd]->Print(pdffile);
    if (nd==17) cplot2[nd]->Print(pdffile+")");
    cplot2[nd]->Close();
    }  
   //
   pdffile="plots/"+basename+"_xptardiff.pdf";
      Double_t max_ny;
      TCanvas *cplot3[18];
    for (int nd=0;nd<18;nd++) {
      cplot3[nd] = new TCanvas(Form("cplot3_%d",nd),Form("Delta run =%d",-12+2*nd), 900,700);
      cplot3[nd]->Divide(3,4);
    for (int nny=0;nny<11;nny++) {
      cplot3[nd]->cd(nny+1);
      max_ny=0;
      for (int nnx=0;nnx<11;nnx++) {
	if (hxptardiffnew[nd][nnx][nny]->GetMaximum()>max_ny) max_ny = hxptardiffnew[nd][nnx][nny]->GetMaximum(); 
	hxptardiffnew[nd][nnx][nny]->SetTitle(Form("Delta run =%d Ny = %d ",-12+2*nd,nny));
      }
      //gPad->SetLogy();
      // cout << nd << " " << nny << " " << max_ny << " " << endl;
      for (int nnx=0;nnx<11;nnx++) {
	if ( hxptardiffnew[nd][nnx][nny]->Integral()>0) {
      hxptardiffnew[nd][nnx][nny]->SetMaximum(max_ny*1.1); 
      hxptardiffnew[nd][nnx][nny]->SetMinimum(0); 
      hxptardiffnew[nd][nnx][nny]->SetAxisRange(-.5,.5,"X"); 
      hxptardiffnew[nd][nnx][nny]->Draw("same"); 
      hxptardiffnew[nd][nnx][nny]->SetLineColor(nnx); 
	}
      }
    }
    if (nd==0) cplot3[nd]->Print(pdffile+"(");
    if (nd>0&&nd<17) cplot3[nd]->Print(pdffile);
    if (nd==17) cplot3[nd]->Print(pdffile+")");
    cplot3[nd]->Close();
    }  
    //
    // end bracket
}
