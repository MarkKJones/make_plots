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

void make_hist_simc_shms_ep_elastic(TString basename="",Int_t nrun=3288){
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
 Float_t  p_xfp;
   tsimc->SetBranchAddress("hsxfp",&p_xfp);
 Float_t  p_yfp;
   tsimc->SetBranchAddress("hsyfp",&p_yfp);
 Float_t  e_xfp;
   tsimc->SetBranchAddress("ssxfp",&e_xfp);
 Float_t  e_yfp;
   tsimc->SetBranchAddress("ssyfp",&e_yfp);
 Float_t  p_xpfp;
   tsimc->SetBranchAddress("hsxpfp",&p_xpfp);
 Float_t  p_ypfp;
   tsimc->SetBranchAddress("hsypfp",&p_ypfp);
 Float_t  e_xpfp;
   tsimc->SetBranchAddress("ssxpfp",&e_xpfp);
 Float_t  e_ypfp;
   tsimc->SetBranchAddress("ssypfp",&e_ypfp);
 Float_t  p_xptar;
   tsimc->SetBranchAddress("hsxptar",&p_xptar);
 Float_t  p_ytar;
   tsimc->SetBranchAddress("hsytar",&p_ytar);
  Float_t  p_yptar;
   tsimc->SetBranchAddress("hsyptar",&p_yptar);
Float_t  p_delta;
   tsimc->SetBranchAddress("hsdelta",&p_delta);
 Float_t  e_xptar;
   tsimc->SetBranchAddress("ssxptar",&e_xptar);
 Float_t  e_yptar;
   tsimc->SetBranchAddress("ssyptar",&e_yptar);
 Float_t  e_ytar;
   tsimc->SetBranchAddress("ssytar",&e_ytar);
 Float_t  e_delta;
   tsimc->SetBranchAddress("ssdelta",&e_delta);
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
    TH1F *h_eytar = new TH1F("h_eytar",Form("Run %d ; SHMS Ytar (cm);Counts",nrun), 100, -5,5.);
    HList.Add(h_eytar);
    TH1F *h_hytar = new TH1F("h_hytar",Form("Run %d ; HMS Ytar (cm);Counts",nrun), 100, -5,5.);
    HList.Add(h_hytar);
    TH1F *h_pxptar = new TH1F("h_pxptar",Form("Run %d ; HMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(h_pxptar);
    TH1F *h_pyptar = new TH1F("h_pyptar",Form("Run %d ; HMS Yp_tar;Counts",nrun), 100, -.04,.04);
    HList.Add(h_pyptar);
    TH1F *h_exptar = new TH1F("h_exptar",Form("Run %d ; SHMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(h_exptar);
    TH1F *h_eyptar = new TH1F("h_eyptar",Form("Run %d ; SHMS Yp_tar;Counts",nrun), 100, -.04,.04);
    HList.Add(h_eyptar);
    TH1F *h_edelta = new TH1F("h_edelta",Form("Run %d ; SHMS Delta;Counts",nrun), 100,-10.,20.);
    HList.Add(h_edelta);
    TH1F *h_pdelta = new TH1F("h_pdelta",Form("Run %d ; HMS Delta;Counts",nrun), 100,-10.,10.);
    HList.Add(h_pdelta);
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 125, 0.8,1.3);
    HList.Add(hW);
    TH1F *hW_2 = new TH1F("hW_2",Form("Run %d ; W (GeV) (p_delta < 10);Counts",nrun), 250, 0.8,1.3);
    HList.Add(hW_2);
    TH1F *hW_3 = new TH1F("hW_3",Form("Run %d ; W (GeV) (p_delta < 8);Counts",nrun), 250, 0.8,1.3);
    HList.Add(hW_3);
     TH1F *hEmiss = new TH1F("hEmiss",Form("Run %d ; Emiss (GeV);Counts",nrun), 200, -.05,.1);
    HList.Add(hEmiss);
    TH1F *hPmiss = new TH1F("hPmiss",Form("Run %d ; Pmiss (GeV);Counts",nrun), 200, -.1,.1);
    HList.Add(hPmiss);
    TH1F *hPmissx = new TH1F("hPmissx",Form("Run %d ; Pmissx (GeV);Counts",nrun), 200, -.1,.1);
    HList.Add(hPmissx);
    TH1F *hPmissy = new TH1F("hPmissy",Form("Run %d ; Pmissy (GeV);Counts",nrun), 200, -.1,.1);
     HList.Add(hPmissy);
   TH1F *hPmissz = new TH1F("hPmissz",Form("Run %d ; Pmissz (GeV);Counts",nrun), 200, -.1,.1);
    HList.Add(hPmissz);
    TH1F *hprot_mom_calc = new TH1F("hprot_mom_calc",Form("Run %d ; (pcalc-p)/p ;Counts",nrun), 100,-0.1,.1 );
    HList.Add(hprot_mom_calc);
     TH1F *h_pxfp = new TH1F("h_pxfp",Form("Run %d ; HMS X_fp;Counts",nrun), 100, -40.,40.);
    HList.Add(h_pxfp);
    TH1F *h_pyfp = new TH1F("h_pyfp",Form("Run %d ; HMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(h_pyfp);
    TH1F *h_pxpfp = new TH1F("h_pxpfp",Form("Run %d ; HMS Xp_fp;Counts",nrun), 100, -.1,.1);
    HList.Add(h_pxpfp);
    TH1F *h_pypfp = new TH1F("h_pypfp",Form("Run %d ; HMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(h_pypfp);
    TH1F *h_exfp = new TH1F("h_exfp",Form("Run %d ; SHMS X_fp;Counts",nrun), 100, -50.,50.);
    HList.Add(h_exfp);
    TH1F *h_eyfp = new TH1F("h_eyfp",Form("Run %d ; SHMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(h_eyfp);
   TH2F *h_eyfp_xfp = new TH2F("h_eyfp_xfp",Form("Run %d ; SHMS Y_fp; X_fp",nrun), 100, -20.,20., 100, -50.,50.);
    HList.Add(h_eyfp_xfp);
   TH2F *h_expfp_xfp = new TH2F("h_expfp_xfp",Form("Run %d ; SHMS Xp_fp; X_fp",nrun), 100, -.1,.1, 100, -50.,50.);
    HList.Add(h_expfp_xfp);
   TH2F *h_eypfp_xfp = new TH2F("h_eypfp_xfp",Form("Run %d ; SHMS Yp_fp; X_fp",nrun), 100, -.05,.05, 100, -50.,50.);
    HList.Add(h_eypfp_xfp);
    TH1F *h_expfp = new TH1F("h_expfp",Form("Run %d ; SHMS Xp_fp;Counts",nrun), 100, -.12,.12);
    HList.Add(h_expfp);
    TH1F *h_eypfp = new TH1F("h_eypfp",Form("Run %d ; SHMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(h_eypfp);
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
    if (nrun ==3288) {
      Nent_simc=200000.;
      Normfac =0.247006E+07;
      Exp_charge=144.502 ;
      Exp_eff = .98*.985*.986;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==3371) {
      Nent_simc=200000.;
      Normfac =0.375730E+07;
      Exp_charge=49.335 ;
      Exp_eff = .947*.986*.984;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==3374) {
      Nent_simc=200000.;
      Normfac =0.141957E+07;
      Exp_charge=49.943 ;
      Exp_eff = .937*.988*.967;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==3377) {
      Nent_simc=200000.;
      Normfac =0.940077E+06;
      Exp_charge=39.343 ;
      Exp_eff = .806*.988*.981*1.0;
      Exp_charge=3.377 ; //50000 events
      Exp_eff = .838*.985*.935*1.0;
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
    if (nrun ==6634) {
      Nent_simc=100000.;
      Normfac = 0.133479E+08;
      Exp_charge= 9.701; //mC
      Exp_eff = 100./100.*0.9915*0.9815;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==6881) {
      Nent_simc=100000.;
      Normfac = 0.781056E+07;
      Exp_charge= 3.331; //mC
      Exp_eff = 100./100.*0.9919*0.9783;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==2452) {
      Nent_simc=100000.;
      Normfac = 0.157892E+08;
      Exp_charge= 173.521; //mC
      Exp_eff = 100./100.*0.96*0.984;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==7000) {
      Nent_simc=200000.;
      Normfac = 0.786524E+07;
      Exp_charge= 173.521; //mC
      Exp_eff = 100./100.*0.96*0.984;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==7214) {
      Nent_simc=200000.;
      Normfac = 0.670510E+07 ;
      Exp_charge= 80.473; //mC
      Exp_eff = 100./100.*0.96*0.984;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    if (nrun ==7215) {
      Nent_simc=200000.;
      Normfac = 0.672480E+07 ;
      Exp_charge= 28.529; //mC
      Exp_eff = 100./100.*0.96*0.984;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
   }
   if (nrun ==7217) {
      Nent_simc=200000.;
      Normfac = 0.143824E+07  ;
      Exp_charge=1.034; //mC
      Exp_eff = 100./100.*0.94*0.98/2;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
   if (nrun ==8825) {
      Nent_simc=200000.;
      Normfac = 0.143824E+07 ;
      Exp_charge= 1.034; //mC
      Exp_eff = 100./100./3.;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
   if (nrun ==8900) {
      Nent_simc=200000.;
      Normfac = 0.136112E+07 ;
      Exp_charge= 1.106; //mC
      Exp_eff = 84.1/100./3.*.99*.99;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
   if (nrun ==9369) {
      Nent_simc=200000.;
      Normfac = 0.119135E+07;
      Exp_charge= 2.661; //mC
      Exp_eff = 100./100./5.*.99*.99;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
     Double_t th_cent=25.6;
  Double_t Mp = .93827;
   Double_t Ei=10.600;
   Double_t pcent=4.9302;
   if (nrun==2452) {
   th_cent=17.83;
   Ei=10.600;
   }
  Double_t cos_ts=TMath::Cos(th_cent/180*3.14159);
   Double_t sin_ts=TMath::Sin(th_cent/180*3.14159);
   //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		  if (e_delta>-10. && e_delta<22. && p_delta>-8. && p_delta<8.) {
         	hW->Fill(W,Weight*simc_fac);
         	if (p_delta>-10. && p_delta<10.) hW_2->Fill(W,Weight*simc_fac);
         	if (p_delta>-8. && p_delta<8.) hW_3->Fill(W,Weight*simc_fac);
                if (W>0.85 && W < 1.05  ) {	
 		      Double_t theta_hms = TMath::ACos((cos_ts + p_yptar*sin_ts) / TMath::Sqrt( 1. + p_xptar*p_xptar + p_yptar * p_yptar ));
		      Double_t pcalc=2*Mp*Ei*(Ei+Mp)*cos(theta_hms)/(Mp*Mp+2*Mp*Ei+Ei*Ei*sin(theta_hms)*sin(theta_hms));
		      Double_t pmom = pcent*(1+p_delta/100.);
		      hprot_mom_calc->Fill((pcalc-pmom)/pmom,Weight*simc_fac);
 		  hEmiss->Fill(emiss,Weight*simc_fac);		  
		  hPmiss->Fill(pmiss,Weight*simc_fac);		  
		  hPmissx->Fill(pmissx,Weight*simc_fac);		  
		  hPmissy->Fill(pmissy,Weight*simc_fac);		  
		  hPmissz->Fill(pmissz,Weight*simc_fac);		  
		  h_exptar->Fill(e_xptar,Weight*simc_fac);		  
		  h_pxptar->Fill(p_xptar,Weight*simc_fac);		  
		  h_eyptar->Fill(e_yptar,Weight*simc_fac);		  
		  h_eytar->Fill(e_ytar,Weight*simc_fac);		  
		  h_pyptar->Fill(p_yptar,Weight*simc_fac);		  
		  h_hytar->Fill(p_ytar,Weight*simc_fac);		  
		  h_edelta->Fill(e_delta,Weight*simc_fac);		  
		  h_pdelta->Fill(p_delta,Weight*simc_fac);		  
		  h_exfp->Fill(e_xfp,Weight*simc_fac);		  
		  h_eyfp->Fill(e_yfp,Weight*simc_fac);		  
		  h_eyfp_xfp->Fill(e_yfp,e_xfp,Weight*simc_fac);		  
		  h_expfp_xfp->Fill(e_xpfp,e_xfp,Weight*simc_fac);		  
		  h_eypfp_xfp->Fill(e_ypfp,e_xfp,Weight*simc_fac);		  
		  h_expfp->Fill(e_xpfp,Weight*simc_fac);		  
		  h_eypfp->Fill(e_ypfp,Weight*simc_fac);		  
		  h_pxfp->Fill(p_xfp,Weight*simc_fac);		  
		  h_pyfp->Fill(p_yfp,Weight*simc_fac);		  
		  h_pxpfp->Fill(p_xpfp,Weight*simc_fac);		  
		  h_pypfp->Fill(p_ypfp,Weight*simc_fac);		  
		}
		  }
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
