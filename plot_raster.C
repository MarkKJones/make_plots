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
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_raster(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="../../ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t raster_ay;
   tsimc->SetBranchAddress("H.rb.raster.fr_ya",&raster_ay);
 Double_t raster_ax;
   tsimc->SetBranchAddress("H.rb.raster.fr_xa",&raster_ax);
 Double_t raster_raw_ay;
   tsimc->SetBranchAddress("H.rb.raster.frxaRawAdc",&raster_raw_ay);
 Double_t raster_raw_ax;
   tsimc->SetBranchAddress("H.rb.raster.fryaRawAdc",&raster_raw_ax);
 Double_t raster_raw_by;
   tsimc->SetBranchAddress("H.rb.raster.frxbRawAdc",&raster_raw_by);
 Double_t raster_raw_bx;
   tsimc->SetBranchAddress("H.rb.raster.frybRawAdc",&raster_raw_bx);
 Double_t  sumnpe;
   tsimc->SetBranchAddress("H.cer.npeSum",&sumnpe);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
   // Define histograms
   TH1F *hrastx = new TH1F("hrastx", Form("Run %d ;Raster X (cm) ; Counts",nrun), 2000, -.3,.3);
   TH1F *hrasty = new TH1F("hrasty",Form("Run %d ;Raster Y (cm) ;Counts",nrun), 2000, -.3,.3);
   TH1F *hrastrawxa = new TH1F("hrastrawxa", Form("Run %d ;Raster Xa ADC;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hrastrawya = new TH1F("hrastrawya",Form("Run %d ;Raster Ya ADC ;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hrastrawxb = new TH1F("hrastrawxb", Form("Run %d ;Raster Xb ADC;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hrastrawyb = new TH1F("hrastrawyb",Form("Run %d ;Raster Yb ADC ;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hytar = new TH1F("hytar",Form("Run %d ; Y_tar (cm);Counts",nrun), 200, -5.,5.);
   TH2F *hytar_yrast = new TH2F("hytar_yrast",Form("Run %d ; Y_tar (cm); Raster Y ADC (mV)",nrun), 200, -5.,5.,2000, -.3,.3);
   TH2F *hytar_xrast = new TH2F("hytar_xrast",Form("Run %d ; Y_tar (cm); Raster X ADC (mV)",nrun), 200, -5.,5.,2000, -.3,.3);
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (sumnpe > 2. &&  delta>-8 &&  delta<8 ) {
		  hrasty->Fill(raster_ay);
		  hrastx->Fill(raster_ax);
		  hrastrawya->Fill(raster_raw_ay);
		  hrastrawxa->Fill(raster_raw_ax);
		  hrastrawyb->Fill(raster_raw_by);
		  hrastrawxb->Fill(raster_raw_bx);
		  hytar->Fill(ytar);
		  hytar_yrast->Fill(ytar,raster_ay);
		  hytar_xrast->Fill(ytar,raster_ax);
		}
	}
	// plot data
TCanvas *c = new TCanvas("c", "raster", 800, 1200);
c->Divide(1,2);
c->cd(1);
hrastx->Draw();
c->cd(2);
hrasty->Draw();
    outputpdf=basename+"_raster.pdf";
c->Print(outputpdf);
//
TCanvas *c3 = new TCanvas("c3", "raw raster", 800, 1200);
c3->Divide(2,2);
c3->cd(1);
hrastrawxa->Draw();
c3->cd(2);
hrastrawya->Draw();
c3->cd(3);
hrastrawxb->Draw();
c3->cd(4);
hrastrawyb->Draw();
    outputpdf=basename+"_raster_raw.pdf";
c->Print(outputpdf);
//
TCanvas *c1 = new TCanvas("c1", "Ytar", 800, 1200);
c1->Divide(1,1);
c1->cd(1);
hytar->Draw();
    outputpdf=basename+"_ytar.pdf";
c1->Print(outputpdf);
//
TCanvas *c2 = new TCanvas("c2", "Ytar v rast", 800, 1200);
c2->Divide(1,2);
c2->cd(1);
hytar_yrast->Draw("colz");
c2->cd(2);
hytar_xrast->Draw("colz");
    outputpdf=basename+"_ytar_rast.pdf";
c2->Print(outputpdf);


 //
}
