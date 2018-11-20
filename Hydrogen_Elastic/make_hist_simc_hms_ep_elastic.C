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

void make_hist_simc_hms_ep_elastic(TString basename="",Int_t nrun=1272){
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
 Float_t  e_xfp;
   tsimc->SetBranchAddress("hsxfp",&e_xfp);
 Float_t  e_yfp;
   tsimc->SetBranchAddress("hsyfp",&e_yfp);
 Float_t  p_xfp;
   tsimc->SetBranchAddress("ssxfp",&p_xfp);
 Float_t  p_yfp;
   tsimc->SetBranchAddress("ssyfp",&p_yfp);
 Float_t  e_xpfp;
   tsimc->SetBranchAddress("hsxpfp",&e_xpfp);
 Float_t  e_ypfp;
   tsimc->SetBranchAddress("hsypfp",&e_ypfp);
 Float_t  p_xpfp;
   tsimc->SetBranchAddress("ssxpfp",&p_xpfp);
 Float_t  p_ypfp;
   tsimc->SetBranchAddress("ssypfp",&p_ypfp);
 Float_t  e_xptar;
   tsimc->SetBranchAddress("hsxptar",&e_xptar);
 Float_t  e_yptar;
   tsimc->SetBranchAddress("hsyptar",&e_yptar);
 Float_t  e_delta;
   tsimc->SetBranchAddress("hsdelta",&e_delta);
 Float_t  p_xptar;
   tsimc->SetBranchAddress("ssxptar",&p_xptar);
 Float_t  p_yptar;
   tsimc->SetBranchAddress("ssyptar",&p_yptar);
 Float_t  p_delta;
   tsimc->SetBranchAddress("ssdelta",&p_delta);
 Float_t  emiss;
   tsimc->SetBranchAddress("Em",&emiss);
 Float_t  pmiss;
   tsimc->SetBranchAddress("Pm",&pmiss);
 Float_t  pmissx;
   tsimc->SetBranchAddress("Pmper",&pmissx);
 Float_t  pmissy;
   tsimc->SetBranchAddress("Pmoop",&pmissy);
 Float_t  pmissz;
   tsimc->SetBranchAddress("Pmpar",&pmissz);
 Float_t  Weight;
   tsimc->SetBranchAddress("Weight",&Weight);
   // Define histograms
    TH1F *hxptar = new TH1F("hxptar",Form("Run %d ; HMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hxptar);
    TH1F *hyptar = new TH1F("hyptar",Form("Run %d ; HMS Yp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hyptar);
    TH1F *pxptar = new TH1F("pxptar",Form("Run %d ; SHMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(pxptar);
    TH1F *pyptar = new TH1F("pyptar",Form("Run %d ; SHMS Yp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(pyptar);
    TH1F *pdelta = new TH1F("pdelta",Form("Run %d ; SHMS Delta;Counts",nrun), 100,-10.,10.);
    HList.Add(pdelta);
    TH1F *hdelta = new TH1F("hdelta",Form("Run %d ; HMS Delta;Counts",nrun), 100,-10.,10.);
    HList.Add(hdelta);
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 100, 0.8,1.2);
    HList.Add(hW);
     TH1F *hEmiss = new TH1F("hEmiss",Form("Run %d ; Emiss (GeV);Counts",nrun), 100, -.2,.2);
    HList.Add(hEmiss);
    TH1F *hPmiss = new TH1F("hPmiss",Form("Run %d ; Pmiss (GeV);Counts",nrun), 100, -.2,.2);
    HList.Add(hPmiss);
    TH1F *hPmissx = new TH1F("hPmissx",Form("Run %d ; Pmissx (GeV);Counts",nrun), 100, -.2,.2);
    HList.Add(hPmissx);
    TH1F *hPmissy = new TH1F("hPmissy",Form("Run %d ; Pmissy (GeV);Counts",nrun), 100, -.2,.2);
     HList.Add(hPmissy);
   TH1F *hPmissz = new TH1F("hPmissz",Form("Run %d ; Pmissz (GeV);Counts",nrun), 100, -.2,.2);
    HList.Add(hPmissz);
     TH1F *hxfp = new TH1F("hxfp",Form("Run %d ; HMS X_fp;Counts",nrun), 100, -40.,40.);
    HList.Add(hxfp);
    TH1F *hyfp = new TH1F("hyfp",Form("Run %d ; HMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(hyfp);
    TH1F *hxpfp = new TH1F("hxpfp",Form("Run %d ; HMS Xp_fp;Counts",nrun), 100, -.1,.1);
    HList.Add(hxpfp);
    TH1F *hypfp = new TH1F("hypfp",Form("Run %d ; HMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(hypfp);
    TH1F *pxfp = new TH1F("pxfp",Form("Run %d ; SHMS X_fp;Counts",nrun), 100, -50.,50.);
    HList.Add(pxfp);
    TH1F *pyfp = new TH1F("pyfp",Form("Run %d ; SHMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(pyfp);
    TH1F *pxpfp = new TH1F("pxpfp",Form("Run %d ; SHMS Xp_fp;Counts",nrun), 100, -.12,.12);
    HList.Add(pxpfp);
    TH1F *pypfp = new TH1F("pypfp",Form("Run %d ; SHMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(pypfp);
  //
    Float_t Normfac;
    Float_t Nent_simc;
    Float_t Exp_charge; // assume SIMC charge is 1 mC
    Float_t Exp_eff; // Livetime*Track_eff*Hodo_eff
    Float_t simc_fac=1.;
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
    if (nrun ==4858) {
      Nent_simc=200000.;
      Normfac = 0.594325E+07;
      Exp_charge= 138.926; //mC
      Exp_eff = 100./100.*0.98*0.99;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==6009) {
      Nent_simc=200000.;
      Normfac = 0.840095E+07;
      Exp_charge= 36.765; //mC
      Exp_eff = 100./100.*0.99*0.90;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==6010) {
      Nent_simc=200000.;
      Normfac = 0.719275E+07;
      Exp_charge= 63.217; //mC
      Exp_eff = 100./100.*0.99*0.91;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
                if (TMath::Abs(emiss)<0.1 &&	TMath::Abs(e_delta) < 8  ) {	
        	hW->Fill(W,Weight*simc_fac);
		  hEmiss->Fill(emiss,Weight*simc_fac);		  
		  hPmiss->Fill(pmiss,Weight*simc_fac);		  
		  hPmissx->Fill(pmissx,Weight*simc_fac);		  
		  hPmissy->Fill(pmissy,Weight*simc_fac);		  
		  hPmissz->Fill(pmissz,Weight*simc_fac);		  
		  hxptar->Fill(e_xptar,Weight*simc_fac);		  
		  pxptar->Fill(p_xptar,Weight*simc_fac);		  
		  hyptar->Fill(e_yptar,Weight*simc_fac);		  
		  pyptar->Fill(p_yptar,Weight*simc_fac);		  
		  hdelta->Fill(e_delta,Weight*simc_fac);		  
		  pdelta->Fill(p_delta,Weight*simc_fac);		  
		  hxfp->Fill(e_xfp,Weight*simc_fac);		  
		  hyfp->Fill(e_yfp,Weight*simc_fac);		  
		  hxpfp->Fill(e_xpfp,Weight*simc_fac);		  
		  hypfp->Fill(e_ypfp,Weight*simc_fac);		  
		  pxfp->Fill(p_xfp,Weight*simc_fac);		  
		  pyfp->Fill(p_yfp,Weight*simc_fac);		  
		  pxpfp->Fill(p_xpfp,Weight*simc_fac);		  
		  pypfp->Fill(p_ypfp,Weight*simc_fac);		  
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}