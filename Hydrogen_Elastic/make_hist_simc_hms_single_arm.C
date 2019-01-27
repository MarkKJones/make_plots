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

void make_hist_simc_hms_single_arm(TString basename="",Int_t nrun=1272){
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
   outputhist= "hist/"+basename+"_simc_hms_ep_elastic_hist.root";
 TObjArray HList(0);
   //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("h666");
// Define branches
 Float_t  W;
   tsimc->SetBranchAddress("W",&W);
 Float_t  Weight;
   tsimc->SetBranchAddress("Weight",&Weight);
 Float_t  e_xfp;
   tsimc->SetBranchAddress("hsxfp",&e_xfp);
 Float_t  e_yfp;
   tsimc->SetBranchAddress("hsyfp",&e_yfp);
 Float_t  e_xpfp;
   tsimc->SetBranchAddress("hsxpfp",&e_xpfp);
 Float_t  e_ypfp;
   tsimc->SetBranchAddress("hsypfp",&e_ypfp);
 Float_t  e_xptar;
   tsimc->SetBranchAddress("hsxptar",&e_xptar);
 Float_t  e_yptar;
   tsimc->SetBranchAddress("hsyptar",&e_yptar);
 Float_t  e_delta;
   tsimc->SetBranchAddress("hsdelta",&e_delta);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 100, 0.75,1.2);
    HList.Add(hW);
      TH1F *hxfp = new TH1F("hxfp",Form("Run %d ; HMS X_fp;Counts",nrun), 100, -40.,40.);
    HList.Add(hxfp);
    TH1F *hyfp = new TH1F("hyfp",Form("Run %d ; HMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(hyfp);
    TH1F *hxpfp = new TH1F("hxpfp",Form("Run %d ; HMS Xp_fp;Counts",nrun), 100, -.1,.1);
    HList.Add(hxpfp);
    TH1F *hypfp = new TH1F("hypfp",Form("Run %d ; HMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(hypfp);
    TH1F *hxptar = new TH1F("hxptar",Form("Run %d ; HMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hxptar);
    TH1F *hyptar = new TH1F("hyptar",Form("Run %d ; HMS Yp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hyptar);
     TH1F *hdelta = new TH1F("hdelta",Form("Run %d ; HMS Delta;Counts",nrun), 100,-10.,10.);
    HList.Add(hdelta);
  //
    Double_t Normfac;
    Double_t Nent_simc;
    Double_t Exp_charge; // assume SIMC charge is 1 mC
    Double_t Exp_eff; // Livetime*Track_eff*Hodo_eff
    Double_t simc_fac=1.;
    if (nrun ==1329) {
      Nent_simc=100000.;
      Normfac = 0.107289E+08;
      Exp_charge= 1428.397 /1000.;
      Exp_eff = 23.47/100.*0.9245*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==1272) {
      Nent_simc=100000.;
      Normfac = 0.107289E+08;
      Exp_charge= 1356.930/1000.;
      Exp_eff = 5.6348/100.*0.9856*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==1161) {
      Nent_simc=100000.;
      Normfac = 0.107191E+08;
      Exp_charge= 148.989/1000.;
      Exp_eff = 2.6679/100.*0.9653*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==4793) {
      Nent_simc=200000.;
      Normfac = 0.169261E+08;
      Exp_charge= 96046.593/1000.;
      Exp_eff = 100./100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==4788) {
      Nent_simc=200000.;
      Normfac = 0.169261E+08;
      Exp_charge= 32119.659/1000.;
      Exp_eff = 96./100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==4798) {
      Nent_simc=200000.;
      Normfac = 0.169261E+08;
      Exp_charge= 61659.676/1000.;
      Exp_eff = 100./100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==4816) {
      Nent_simc=200000.;
      Normfac = 0.171086E+08;
      Exp_charge= 35168.773/1000.;
      Exp_eff = 77.0/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==4811) {
      Nent_simc=200000.;
      Normfac = 0.171700E+08;
      Exp_charge= 15417.161/1000.;
      Exp_eff = 68.7/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==4821) {
      Nent_simc=200000.;
      Normfac = 0.165910E+08;
      Exp_charge= 26387.430/1000.;
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==6611) {
      Nent_simc=100000.;
      Normfac = 0.165892E+08;
      Exp_charge= 5812.897/1000.;
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6595) {
      Nent_simc=100000.;
      Normfac = 0.169575E+08;
      Exp_charge= 7285.215/1000./17.;
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6601) {
      Nent_simc=100000.;
      Normfac = 0.165659E+08;
      Exp_charge= 6523.945/1000./9.; // do not use ps factor=9 in report
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6602) {
      Nent_simc=100000.;
      Normfac = 0.165694E+08;
      Exp_charge= 6734.147/1000./3.;
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6609) {
      Nent_simc=100000.;
      Normfac = 0.165835E+08;
      Exp_charge= 5477.710/1000./2.;
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6871) {
      Nent_simc=100000.;
      Normfac = 0.173024E+08;
      Exp_charge= 1120.581/1000./33.;
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6875) {
      Nent_simc=100000.;
      Normfac = 0.172127E+08;
      Exp_charge= 4783.038/1000./65.;
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6876) {
      Nent_simc=100000.;
      Normfac = 0.171463E+08;
      Exp_charge= 3722.930/1000./9.;
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6879) {
      Nent_simc=100000.;
      Normfac = 0.164582E+08;
      Exp_charge= 12170.401/1000./1.;
      Exp_eff = 99.8/100.*0.99*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hW->Fill(W,Weight*simc_fac);
                if (W<1.075 && TMath::Abs(e_delta)<8.) {
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
