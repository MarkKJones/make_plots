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

void plot_beta(TString basename="",Int_t nrun=2043){
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
 Double_t raster_ay_adc;
   tsimc->SetBranchAddress("H.rb.raster.frya_adc",&raster_ay_adc);
 Double_t raster_ax_adc;
   tsimc->SetBranchAddress("H.rb.raster.frxa_adc",&raster_ax_adc);
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
 Double_t  betanotrack;
   tsimc->SetBranchAddress("H.hod.betanotrack",&betanotrack);
 Double_t  betatrack;
   tsimc->SetBranchAddress("H.gtr.beta",&betatrack);
 Double_t  edtmtime;
   tsimc->SetBranchAddress("T.hms.hEDTM_tdcTime",&edtmtime);
 Double_t  gtrIndex;
   tsimc->SetBranchAddress("H.gtr.index",&gtrIndex);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
   // Define histograms
   TH1F *hbetanotrack = new TH1F("hbetanotrack", Form("Run %d ;Beta No track ; Counts",nrun), 340, -1.2,1.6);
   TH1F *hbetatrack = new TH1F("hbetatrack", Form("Run %d ;Beta track ; Counts",nrun), 340, -.1,1.6);
   TH2F *hbetatrack_xfp = new TH2F("hbetatrack_xfp", Form("Run %d ;Beta track ; x_fp",nrun), 340, -.1,1.6,300,-50,50);
   TH2F *hbetatrack_yfp = new TH2F("hbetatrack_yfp", Form("Run %d ;Beta track ; y_fp",nrun), 340, -.1,1.6,300,-50,50);

// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (edtmtime==0) {
		if (sumnpe > 2. &&  delta>-8 &&  delta<8 ) {
		hbetanotrack->Fill(betanotrack);
		if (gtrIndex>-1) {
                    hbetatrack->Fill(betatrack);
		    hbetatrack_xfp->Fill(betatrack,xfp);
		    hbetatrack_yfp->Fill(betatrack,yfp);
		}
		}
	}
	}
	// plot data
TCanvas *c = new TCanvas("c", "beta", 800, 1200);
c->Divide(2,2);
c->cd(1);
 hbetanotrack->Draw();
c->cd(2);
 hbetatrack->Draw();
c->cd(3);
 hbetatrack_xfp->Draw("colz");
   outputpdf=basename+"_beta_.pdf";
c->Print(outputpdf);
//


 //
}
