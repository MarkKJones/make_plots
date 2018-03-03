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
#include "TF1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TCutG.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void comp_hms_sieve(TString basename="",Int_t nrun=2043,Int_t nfoils=3){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
   tsimc->SetBranchAddress("H.cer.npeSum",&sumnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t h_extcor_delta_xptar;
   tsimc->SetBranchAddress("H.extcor.delta_xptar",&h_extcor_delta_xptar);
 Double_t  reactx;
   tsimc->SetBranchAddress("H.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("H.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("H.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
  Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("H.gtr.x",&xtar);
Double_t  yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("H.gtr.th",&xptar);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp);
   //
   TH1F *hytar=new TH1F("hytar", Form("Run %d ; Ytar ; Counts",nrun), 240,-6.,6.);
   TH1F *hreactz=new TH1F("hreactz", Form("Run %d ; Z_beam (cm) ; Counts",nrun), 600,-15.,15.);
   TH1F *hreactx=new TH1F("hreactx", Form("Run %d ; X_beam (cm) ; Counts",nrun), 600,-.5,.5);
   TH1F *hreacty=new TH1F("hreacty", Form("Run %d ; Y_beam (cm) ; Counts",nrun), 600,-.5,.5);
   //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (sumnpe>2&&etracknorm > .6 ) {
		  hytar->Fill(ytar);
		  hreactz->Fill(reactz);
		  hreactx->Fill(reactx);
		  hreacty->Fill(reacty);
		}
	}
	//
TCanvas *c3b = new TCanvas("c3b", "react", 900,800);
 c3b->Divide(2,2);
 c3b->cd(1);
 hytar->Draw();
 c3b->cd(2);
 hreactx->Draw();
 c3b->cd(3);
 hreacty->Draw();
 c3b->cd(4);
 hreactz->Draw();
outputpdf="plots/"+basename+"_comp_hms_sieve.pdf";
c3b->Print(outputpdf);

   //
}
