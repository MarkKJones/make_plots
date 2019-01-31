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

void make_hist_simc_shms_single_arm(TString basename="",Int_t nrun=1272){
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
   outputhist= "hist/"+basename+"_simc_shms_ep_elastic_hist.root";
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
   tsimc->SetBranchAddress("ssxfp",&e_xfp);
 Float_t  e_yfp;
   tsimc->SetBranchAddress("ssyfp",&e_yfp);
 Float_t  e_xpfp;
   tsimc->SetBranchAddress("ssxpfp",&e_xpfp);
 Float_t  e_ypfp;
   tsimc->SetBranchAddress("ssypfp",&e_ypfp);
 Float_t  e_xptar;
   tsimc->SetBranchAddress("ssxptar",&e_xptar);
 Float_t  e_yptar;
   tsimc->SetBranchAddress("ssyptar",&e_yptar);
 Float_t  e_delta;
   tsimc->SetBranchAddress("ssdelta",&e_delta);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 100, 0.75,1.2);
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
    Double_t Normfac=1.0;
    Double_t Nent_simc=1.0;
    Double_t Exp_charge=1.0; // assume SIMC charge is 1 mC
    Double_t Exp_eff=1.0; // Livetime*Track_eff*Hodo_eff
    Double_t simc_fac=1.;
     if (nrun ==6871) {
      Nent_simc=100000.;
      Normfac = 0.118724E+08;
      Exp_charge= 1157.288/1000./513;
      Exp_eff = 0.934;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6875) {
      Nent_simc=100000.;
      Normfac = 0.117555E+08;
      Exp_charge= 4783.038/1000./257;
      Exp_eff = 100./100.*0.94*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6876) {
      Nent_simc=100000.;
      Normfac = 0.116965E+08;
      Exp_charge= 3722.930/1000./65.;
      Exp_eff = 100./100.*0.97*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6879) {
      Nent_simc=100000.;
      Normfac = 0.115178E+08;
      Exp_charge= 12343.869/1000./1.0;
      Exp_eff = 0.982;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6621) {
      Nent_simc=100000.;
      Normfac = 0.116489E+08;
      Exp_charge= 7688.217/1000./17.;
      Exp_eff = 0.9756;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6623) {
      Nent_simc=100000.;
      Normfac = 0.116091E+08 ;
      Exp_charge= 5402.358/1000./5.0;
      Exp_eff =.9803;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6627) {
      Nent_simc=100000.;
      Normfac = 0.115876E+08;
      Exp_charge= 8522.162/1000./5.0;
      Exp_eff =.9817;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6629) {
      Nent_simc=100000.;
      Normfac =  0.115866E+08 ;
      Exp_charge=8420.190/1000./3.0;
      Exp_eff =.9823;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==6633) {
      Nent_simc=100000.;
      Normfac = 0.115915E+08 ;
      Exp_charge= 11485.874/1000./3.0;
      Exp_eff =.9833;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hW->Fill(W,Weight*simc_fac);
                if (W<1.075 &&  e_delta > -10. && e_delta < 20. ) {
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