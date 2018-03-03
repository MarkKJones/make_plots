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

void plot_shms(TString basename="",Int_t nrun=2043){
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
   tsimc->SetBranchAddress("P.rb.raster.fr_ya",&raster_ay);
 Double_t raster_ax;
   tsimc->SetBranchAddress("P.rb.raster.fr_xa",&raster_ax);
 Double_t raster_raw_ay;
   tsimc->SetBranchAddress("P.rb.raster.frxaRawAdc",&raster_raw_ay);
 Double_t raster_raw_ax;
   tsimc->SetBranchAddress("P.rb.raster.fryaRawAdc",&raster_raw_ax);
 Double_t raster_raw_by;
   tsimc->SetBranchAddress("P.rb.raster.frxbRawAdc",&raster_raw_by);
 Double_t raster_raw_bx;
   tsimc->SetBranchAddress("P.rb.raster.frybRawAdc",&raster_raw_bx);
 Double_t  sumnpe;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
 Double_t  delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etracknorm",&etracknorm);
 Double_t  W;
   tsimc->SetBranchAddress("P.kin.W",&W);
 Double_t  ytar;
   tsimc->SetBranchAddress("P.gtr.y",&ytar);
   // Define histograms
   TH1F *hW = new TH1F("hW", Form("Run %d ; W (GeV) ; Counts",nrun), 200, -0.8,1.5);
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (sumnpe > 2. ) {
		  hW->Fill(W);
		}
	}
	// plot data
TCanvas *c = new TCanvas("c", "pos raster", 800, 1200);
c->Divide(1,2);
c->cd(1);
hW->Draw();
    outputpdf=basename+"_shms_kin.pdf";
c->Print(outputpdf);

 //
}
