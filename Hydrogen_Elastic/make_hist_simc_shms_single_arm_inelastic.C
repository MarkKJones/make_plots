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
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
using namespace std;

void make_hist_simc_shms_single_arm_inelastic(TString basename="",Int_t nrun=1272){
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
   outputhist= "hist/"+basename+"_simc_shms_ep_inelastic_hist.root";
 TObjArray HList(0);
   //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("h1411");
// Define branches
 Float_t  W;
   tsimc->SetBranchAddress("ysieve",&W);
 Float_t  Weight;
   tsimc->SetBranchAddress("xsnum",&Weight);
 Float_t  e_xfp;
   tsimc->SetBranchAddress("psxfp",&e_xfp);
 Float_t  e_yfp;
   tsimc->SetBranchAddress("psyfp",&e_yfp);
 Float_t  e_xpfp;
   tsimc->SetBranchAddress("psxpfp",&e_xpfp);
 Float_t  e_ypfp;
   tsimc->SetBranchAddress("psypfp",&e_ypfp);
 Float_t  e_xptar;
   tsimc->SetBranchAddress("psxptar",&e_xptar);
 Float_t  e_yptar;
   tsimc->SetBranchAddress("psyptar",&e_yptar);
 Float_t  e_delta;
   tsimc->SetBranchAddress("psdelta",&e_delta);
   // Define histograms
   TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 100, 2.,3.);
    HList.Add(hW);
      TH1F *hxfp = new TH1F("hxfp",Form("Run %d ; SHMS X_fp;Counts",nrun), 100, -40.,40.);
    HList.Add(hxfp);
    TH1F *hyfp = new TH1F("hyfp",Form("Run %d ; SHMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(hyfp);
    TH1F *hxpfp = new TH1F("hxpfp",Form("Run %d ; SHMS Xp_fp;Counts",nrun), 100, -.1,.1);
    HList.Add(hxpfp);
    TH1F *hypfp = new TH1F("hypfp",Form("Run %d ; SHMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(hypfp);
    TH1F *hxptar = new TH1F("hxptar",Form("Run %d ; SHMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hxptar);
    TH1F *hyptar = new TH1F("hyptar",Form("Run %d ; SHMS Yp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hyptar);
     TH1F *hdelta = new TH1F("hdelta",Form("Run %d ; SHMS Delta;Counts",nrun), 100,-15.,25.);
    HList.Add(hdelta);
  //
    Double_t Exp_current=1.0; // assume SIMC charge is 20uA
    Double_t Exp_eff=1.0; // Livetime*Track_eff*Hodo_eff
    Double_t simc_fac=1.;
      Exp_charge= 44611.452;
      Exp_eff =.99;
      simc_fac = Exp_charge/1.*Exp_eff;
      //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (e_delta > -10. && e_delta < 22. && Weight >0) {
		hW->Fill(W,Weight*simc_fac);
 		  hxptar->Fill(e_xptar,Weight*simc_fac);		  
		  hyptar->Fill(e_yptar,Weight*simc_fac);		  
		  hdelta->Fill(e_delta,Weight*simc_fac);		  
		  hxfp->Fill(e_xfp,Weight*simc_fac);		  
		  hyfp->Fill(e_yfp,Weight*simc_fac);		  
		  hxpfp->Fill(e_xpfp,Weight*simc_fac);		  
		  hypfp->Fill(e_ypfp,Weight*simc_fac);
		}		  
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
