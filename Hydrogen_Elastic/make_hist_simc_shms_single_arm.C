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
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 200, 0.9,1.1);
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
     if (nrun ==4785) {
      Nent_simc=100000.;
      Normfac = 0.112592E+08 ;
      Exp_charge= 36242.993/1000./1.0;
      Exp_eff =0.9821;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==4788) {
      Nent_simc=100000.;
      Normfac = .113035E+08 ;
      Exp_charge= 32103.140/1000./1.0;
      Exp_eff =.9806;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==4791) {
      Nent_simc=100000.;
      Normfac = 0.106707E+08 ;
      Exp_charge= 89083.534/1000./1.0;
      Exp_eff =0.979;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
      if (nrun ==4793) {
      Nent_simc=100000.;
      Normfac = 0.113566E+08 ;
      Exp_charge= 95849.709/1000./1.0;
      Exp_eff =0.982;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==4797) {
      Nent_simc=100000.;
      Normfac =0.936453E+07  ;
      Exp_charge= 62420.635/1000./9.0;
      Exp_eff =0.976;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==4805) {
      Nent_simc=100000.;
      Normfac =  0.113736E+08;
      Exp_charge= 30080.535/1000./3.0;
      Exp_eff =0.9811;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==7166) {
      Nent_simc=100000.;
      Normfac =  0.922484E+07;
      Exp_charge= 41343.595/1000./1.0;
      Exp_eff =88.24/100.*0.9932;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==7167) {
      Nent_simc=100000.;
      Normfac =   0.112456E+08 ;
      Exp_charge= 89105.875/1000./1.0;
      Exp_eff =99/100.*.99;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     if (nrun ==7168) {
      Nent_simc=100000.;
      Normfac =  0.113501E+08  ;
      Exp_charge= 169163.460/1000./1.0;
      Exp_eff =99.80/100.*0.993;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
     }
       if (nrun ==7190) {
      Nent_simc=100000.;
      Normfac =  0.115064E+08 ;
      Exp_charge= 6725.374/1000./1.0;
      Exp_eff =0.23*0.9939;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
  }
       if (nrun ==7191) {
      Nent_simc=100000.;
      Normfac =  0.115064E+08 ;
      Exp_charge= 7962.344/1000./2.0;
      Exp_eff =56.7/100.*0.9939;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
  }
       if (nrun ==8482) {
      Nent_simc=100000.;
      Normfac =  0.119638E+08 ;
      Exp_charge= 250.5/1000./257.0;
      Exp_eff =99./100.;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
  }
       if (nrun ==8486) {
      Nent_simc=100000.;
      Normfac =  0.119638E+08 ;
      Exp_charge= 689.074/1000./257.0;
      Exp_eff =99./100.;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
  }
       if (nrun ==8487) {
      Nent_simc=100000.;
      Normfac =  0.119638E+08 ;
      Exp_charge= 1697.681/1000./257.0;
      Exp_eff =99./100.;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
  }
       if (nrun ==8557) {
      Nent_simc=100000.;
      Normfac =  0.117513E+08;
      Exp_charge= 3134.803/1000./9.0;
      Exp_eff =100./100.;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
  }
       if (nrun ==8564) {
      Nent_simc=100000.;
      Normfac =  0.116765E+08;
      Exp_charge= 4414.566/1000./5.0;
      Exp_eff =100./100.;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
  }
       if (nrun ==8566) {
      Nent_simc=100000.;
      Normfac =  0.118243E+08;
      Exp_charge= 1561.090/1000./5.0;
      Exp_eff =10.4441/100.;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
  }
       if (nrun ==9329) {
      Nent_simc=200000.;
      Normfac =  0.114180E+08;
      Exp_charge= 994.650/1000./1.;
      Exp_eff =100./100.;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
  }
      //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hW->Fill(W,Weight*simc_fac);
                if (W<1.075 &&  e_delta > -10. && e_delta < 30. ) {
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
