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
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_shms_dipole(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
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
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_shms_dipole_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf="plots/"+basename+"_shms_dipole.pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
 Double_t  sumnpe;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
 Double_t  hgsumnpe;
   tsimc->SetBranchAddress("P.hgcer.npeSum",&hgsumnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etracknorm);
 Double_t  etotnorm;
   tsimc->SetBranchAddress("P.cal.etotnorm",&etotnorm);
   Double_t  ntrack;
   tsimc->SetBranchAddress("P.dc.ntrack",&ntrack);
  Double_t  PruneSelect;
   tsimc->SetBranchAddress("P.tr.PruneSelect",&PruneSelect);
  Double_t  InsideDipoleExit;
   tsimc->SetBranchAddress("P.dc.InsideDipoleExit",&InsideDipoleExit);
  Double_t  xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
 Double_t  yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
   // Define histograms
  Double_t cer_cut=2.;
  Double_t ep_cut=.8;
  Double_t zoff=-307.;
  //
   TH1F *hPruneSelect = new TH1F("hPruneSelect", Form(" Run %d ; PruneSelect ; Counts",nrun), 15, 0,15);
   HList.Add(hPruneSelect);
   TH1F *hetotnorm = new TH1F("hetotnorm", Form(" Run %d ; Cal E/p ; Counts",nrun), 200, 0.0,1.5);
   HList.Add(hetotnorm);
  TH2F *hsumnpe_etracknorm = new TH2F("hsumnpe_etracknorm", Form("Run %d ; Track  E/p ; NG Sumnpe",nrun), 150, 0.0,1.5,60,0,30);
   HList.Add(hsumnpe_etracknorm);
  TH2F *hxfp_yfp_proj;
  hxfp_yfp_proj = new TH2F(Form("hxfp_yfp_proj"), Form(" Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff), 200,-40.,40.,200,-40.,40.);
   HList.Add(hxfp_yfp_proj);
  TH2F *hxfp_yfp_insidedip_proj;
  hxfp_yfp_insidedip_proj = new TH2F(Form("hxfp_yfp_insidedip_proj"), Form("Inside dipole exit Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff), 200,-40.,40.,200,-40.,40.);
   HList.Add(hxfp_yfp_insidedip_proj);
  TH2F *hxfp_yfp_outsidedip_proj;
  hxfp_yfp_outsidedip_proj = new TH2F(Form("hxfp_yfp_outsidedip_proj"), Form("Outside dipole exit Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff), 200,-40.,40.,200,-40.,40.);
   HList.Add(hxfp_yfp_outsidedip_proj);
  TH2F *hxfp_yfp_ntr1_proj;
  hxfp_yfp_ntr1_proj = new TH2F(Form("hxfp_yfp_ntr1_proj"), Form("Ntrack==1 Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff), 200,-40.,40.,200,-40.,40.);
   HList.Add(hxfp_yfp_ntr1_proj);
  TH2F *hxfp_yfp_insidedip_ntr1_proj;
  hxfp_yfp_insidedip_ntr1_proj = new TH2F(Form("hxfp_yfp_insidedip_ntr1_proj"), Form("Ntrack==1 Inside dipole exit Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff), 200,-40.,40.,200,-40.,40.);
   HList.Add(hxfp_yfp_insidedip_ntr1_proj);
  TH2F *hxfp_yfp_outsidedip_ntr1_proj;
  hxfp_yfp_outsidedip_ntr1_proj = new TH2F(Form("hxfp_yfp_outsidedip_ntr1_proj"), Form("Ntrack==1 Outside dipole exit Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff), 200,-40.,40.,200,-40.,40.);
   HList.Add(hxfp_yfp_outsidedip_ntr1_proj);
   TH2F *hxfp_yfp_ntr2_proj;
  hxfp_yfp_ntr2_proj = new TH2F(Form("hxfp_yfp_ntr2_proj"), Form("Ntrack>1 Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff), 200,-40.,40.,200,-40.,40.);
   HList.Add(hxfp_yfp_ntr2_proj);
  TH2F *hxfp_yfp_insidedip_ntr2_proj;
  hxfp_yfp_insidedip_ntr2_proj = new TH2F(Form("hxfp_yfp_insidedip_ntr2_proj"), Form("Ntrack>1 Inside dipole exit Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff), 200,-40.,40.,200,-40.,40.);
   HList.Add(hxfp_yfp_insidedip_ntr2_proj);
  TH2F *hxfp_yfp_outsidedip_ntr2_proj;
  hxfp_yfp_outsidedip_ntr2_proj = new TH2F(Form("hxfp_yfp_outsidedip_ntr2_proj"), Form("Ntrack>1 Outside dipole exit Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff), 200,-40.,40.,200,-40.,40.);
   HList.Add(hxfp_yfp_outsidedip_ntr2_proj);
// loop over entries
Long64_t nentries = tsimc->GetEntries();
// nentries=300000;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (i%50000==0) cout << " Entry = " << i << endl;
		if (ntrack>0) {
		  hxfp_yfp_proj->Fill(xfp+xpfp*zoff,yfp+ypfp*zoff); 
		  if (ntrack>1) hPruneSelect->Fill(PruneSelect);
		  if (InsideDipoleExit==1) hxfp_yfp_insidedip_proj->Fill(xfp+xpfp*zoff,yfp+ypfp*zoff); 
		  if (InsideDipoleExit==0) hxfp_yfp_outsidedip_proj->Fill(xfp+xpfp*zoff,yfp+ypfp*zoff); 
		  if (ntrack==1) {
		  hxfp_yfp_ntr1_proj->Fill(xfp+xpfp*zoff,yfp+ypfp*zoff); 
		  if (InsideDipoleExit==1) hxfp_yfp_insidedip_ntr1_proj->Fill(xfp+xpfp*zoff,yfp+ypfp*zoff); 
		  if (InsideDipoleExit==0) hxfp_yfp_outsidedip_ntr1_proj->Fill(xfp+xpfp*zoff,yfp+ypfp*zoff); 
		  }
		  if (ntrack>1) {
		  hxfp_yfp_ntr2_proj->Fill(xfp+xpfp*zoff,yfp+ypfp*zoff); 
		  if (InsideDipoleExit==1) hxfp_yfp_insidedip_ntr2_proj->Fill(xfp+xpfp*zoff,yfp+ypfp*zoff); 
		  if (InsideDipoleExit==0) hxfp_yfp_outsidedip_ntr2_proj->Fill(xfp+xpfp*zoff,yfp+ypfp*zoff); 
		  }
		}
	}
// Define exit window at dipole exit ( shape is rectangle with semi-circle sides)
    const int numcirpt=12;
    Double_t angstep=180./numcirpt;
    Double_t ang;
    const int numpt=2*numcirpt+5;
    int ctpt=0;
    Double_t test_offset = 2.;
    Double_t crad=23.81; // radius of semicircle
    Double_t voffset= crad-26.035+test_offset;
    Double_t hwid=11.549/2.;
    Double_t xwin[numpt]; // vertical
    Double_t ywin[numpt]; // horizontal
    xwin[ctpt]= crad+voffset ;
    ywin[ctpt]= hwid;
    xwin[numpt-1]= crad+voffset ;
    ywin[numpt-1]= hwid;
    Double_t xcent=hwid;
    for (int i=0 ; i < numcirpt ; i++) {
    ctpt++;
    ang = (90.-(i+1)*angstep)*(3.14159/180.);
    xwin[ctpt]= crad*TMath::Sin(ang)+voffset ;
    ywin[ctpt]= hwid+crad*TMath::Cos(ang);
    }
    ctpt++;
    xwin[ctpt]= -crad+voffset ;
    ywin[ctpt]= hwid;
    ctpt++;
    xwin[ctpt]= -crad+voffset ;
    ywin[ctpt]= -hwid;
     for (int i=0 ; i < numcirpt ; i++) {
    ctpt++;
    ang = (-90.+(i+1)*angstep)*(3.14159/180.);
    xwin[ctpt]= crad*TMath::Sin(ang)+voffset ;
    ywin[ctpt]= -hwid-crad*TMath::Cos(ang);
    }
   ctpt++;
    xwin[ctpt]= crad+voffset ;
    ywin[ctpt]= -hwid;
    for (Int_t j=0;j<numpt;j++) {
      //      cout << xwin[j] << " " << ywin[j] <<endl;
    }
    TPolyLine *exitwindow = new TPolyLine(numpt,xwin,ywin);
    exitwindow->SetFillColor(0);
    exitwindow->SetLineColor(kRed);
    exitwindow->SetLineWidth(3);
	// plot data
    TCanvas *cr4;
   cr4 = new TCanvas(Form("cr4"), Form("Z=%5.2f Xfp v Yfp",zoff), 700,700);
   cr4->Divide(2,2);
   cr4->cd(1);
   hxfp_yfp_proj->Draw("colz");
   exitwindow->Draw();
   cr4->cd(2);
   hxfp_yfp_insidedip_proj->Draw("colz");
   exitwindow->Draw();
   cr4->cd(3);
   hxfp_yfp_outsidedip_proj->Draw("colz");
   exitwindow->Draw();
     cr4->Print(outputpdf+"(");
	// plot data
    TCanvas *cr_ntr1;
   cr_ntr1 = new TCanvas(Form("cr_ntr1"), Form("Z=%5.2f NTr=1 Xfp v Yfp",zoff), 700,700);
   cr_ntr1->Divide(2,2);
   cr_ntr1->cd(1);
   hxfp_yfp_ntr1_proj->Draw("colz");
   exitwindow->Draw();
   cr_ntr1->cd(2);
   hxfp_yfp_insidedip_ntr1_proj->Draw("colz");
   exitwindow->Draw();
   cr_ntr1->cd(3);
   hxfp_yfp_outsidedip_ntr1_proj->Draw("colz");
   exitwindow->Draw();
     cr_ntr1->Print(outputpdf);
	// plot data
    TCanvas *cr_ntr2;
   cr_ntr2 = new TCanvas(Form("cr_ntr2"), Form("Z=%5.2f NTr1> Xfp v Yfp",zoff), 700,700);
   cr_ntr2->Divide(2,2);
   cr_ntr2->cd(1);
   hxfp_yfp_ntr2_proj->Draw("colz");
   exitwindow->Draw();
   cr_ntr2->cd(2);
   hxfp_yfp_insidedip_ntr2_proj->Draw("colz");
   exitwindow->Draw();
   cr_ntr2->cd(3);
   hxfp_yfp_outsidedip_ntr2_proj->Draw("colz");
   exitwindow->Draw();
     cr_ntr2->Print(outputpdf);
     //
    TCanvas *cr5;
   cr5 = new TCanvas(Form("cr5"),"PruneSelect", 700,700);
   cr5->Divide(1,1);
   hPruneSelect->Draw();
     cr5->Print(outputpdf+")");
 //
 TFile hsimc(outputhist,"recreate");
        HList.Write();
}
