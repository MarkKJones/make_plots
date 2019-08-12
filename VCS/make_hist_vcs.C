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

void make_hist_vcs(TString basename="",Int_t nrun=3288,TString kintype="kin1b" ){
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
      inputroot="Ana_ROOTfiles/"+basename+".root";
   //  inputroot="Ana_files_mboer/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_vcs_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot);
TTree *tsimc = (TTree*) fsimc->Get("HallCTree");
//
 Float_t ctime;
tsimc->SetBranchAddress("time_roc2",&ctime);
/*
 Float_t ein_4vec[4];
tsimc->SetBranchAddress("ALV_el_in_data",ein_4vec);
 Float_t eout_4vec[4];
tsimc->SetBranchAddress("ALV_el_out_data",eout_4vec);
 Float_t g_4vec[4];
tsimc->SetBranchAddress("ALV_gamma_out_data",g_4vec);
 Float_t p_4vec[4];
tsimc->SetBranchAddress("ALV_proton_out_data",p_4vec);
*/
 Float_t W;
tsimc->SetBranchAddress("W_data",&W);
 Float_t Q2;
tsimc->SetBranchAddress("Q2_data",&Q2);
 Float_t MM2;
tsimc->SetBranchAddress("M2miss_data",&MM2);
 Float_t Thcm;
tsimc->SetBranchAddress("ThCM_vcs_data",&Thcm);
 Float_t Cur;
tsimc->SetBranchAddress("SHMS_B2_cur_cut",&Cur);
 Int_t Nrun;
tsimc->SetBranchAddress("runindex",&Nrun);
 Float_t RunTime;
tsimc->SetBranchAddress("SHMS_act_time",&RunTime);
//
    TH1F *hCtime = new TH1F("hCtime",Form("Run %d ; Ctime (ns);Counts",nrun), 800, 20,50);
    HList.Add(hCtime);
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (Gev);Counts",nrun), 300, .9,1.5);
    HList.Add(hW);
    TH1F *hThcm_vcs = new TH1F("hThcm_vcs",Form("Run %d ; Thcm_vcs (deg);Counts",nrun), 90, 90,180);
    HList.Add(hThcm_vcs);
    TH1F *hThcm_vcs_cut = new TH1F("hThcm_vcs_cut",Form("Run %d Cuts ; Thcm_vcs (deg);Counts",nrun), 90, 90,180);
    HList.Add(hThcm_vcs_cut);
   TH1F *hQ2 = new TH1F("hQ2",Form("Run %d ; Q2 (Gev2);Counts",nrun), 300, .0,.5);
    HList.Add(hQ2);
    TH1F *hMM2 = new TH1F("hMM2",Form("Run %d ; MM2 (Gev2);Counts",nrun), 100, -.05,.05);
    HList.Add(hMM2);
   TH1F *hMM2_vcs_cut = new TH1F("hMM2_vcs_cut",Form("Run %d ; MM2_vcs_cut (Gev2);Counts",nrun), 100, -.05,.05);
    HList.Add(hMM2_vcs_cut);
   TH1F *hMM2_vcs_hiQ2_cut = new TH1F("hMM2_vcs_hiQ2_cut",Form("Run %d ; MM2_vcs_hiQ2_cut (Gev2);Counts",nrun), 50, -.05,.05);
    HList.Add(hMM2_vcs_hiQ2_cut);
    TH1F *hMM2_Wcut = new TH1F("hMM2_Wcut",Form("Run %d ; MM2_Wcut (Gev2);Counts",nrun), 100, -.05,.05);
    HList.Add(hMM2_Wcut);
   TH1F *hMM2_ran = new TH1F("hMM2_ran",Form("Run %d ; MM2 random (Gev2);Counts",nrun), 100, -.05,.05);
    HList.Add(hMM2_ran);
   TH1F *hMM2_ran_vcs_hiQ2_cut = new TH1F("hMM2_ran_vcs_hiQ2_cut",Form("Run %d ; MM2 random_vcs_hiQ2_cut (Gev2);Counts",nrun), 50, -.05,.05);
    HList.Add(hMM2_ran_vcs_hiQ2_cut);
    TH1F *hW_pi0 = new TH1F("hW_pi0",Form("Run %d ; W _pi0 (Gev);Counts",nrun), 300, .9,1.5);
    HList.Add(hW_pi0);
    TH1F *hW_vcs = new TH1F("hW_vcs",Form("Run %d ; W _vcs (Gev);Counts",nrun), 300, .9,1.5);
    HList.Add(hW_vcs);
    TH2F *hMM2_W = new TH2F("hMM2_W",Form("Run %d ; MM2 (Gev2);W",nrun), 100, -.05,.05,100,.9,1.5);
    HList.Add(hMM2_W);
    TH2F *hMM2_thcm = new TH2F("hMM2_thcm",Form("Run %d ; MM2 (Gev2); thcm",nrun), 100, -.05,.05,100,90,180);
    HList.Add(hMM2_thcm);
    TH2F *hMM2_thcm_vcs_hiQ2_cut = new TH2F("hMM2_thcm_vcs_hiQ2_cut",Form("Run %d ; MM2 (Gev2)_vcs_hiQ2_cut; thcm",nrun), 100, -.05,.05,100,90,180);
    HList.Add(hMM2_thcm_vcs_hiQ2_cut);
    //
    Double_t Ctime_ran_cent=32.62;
    Double_t Ctime_ran_sig=0.5;
    Double_t Ctime_cent=38.52;
    Double_t Ctime_sig=0.4;
    if (kintype =="kin1a") Ctime=39.01;
    if (kintype =="kin1b") Ctime=38.82;
    if (kintype =="kin2a") Ctime=38.75;
    if (kintype =="kin2b") Ctime=38.53;
    if (kintype =="kin2b_pass1") Ctime=38.53;
     if (kintype =="kin3b") Ctime=37.75;
     Double_t pi_cent=0.02;
    Double_t pi_sig =0.01;
    Int_t CurRun=0;
    Double_t TotCharge=0;
    Double_t vcs_cent=0.002;
    Double_t vcs_sig =0.002;
    Double_t W_cent=1.232;
    Double_t W_sig=0.005;
    Double_t Q2_cent=0.33;
    Double_t Q2_sig=0.05;
    Double_t Q2_cent_2=0.42;
    Double_t Q2_sig_2=0.04;
    Double_t thcm_cent=141;
    Double_t thcm_cent_2=155;
    Double_t thcm_sig=4;
    if (kintype =="kin1a") thcm_cent=155;
    if (kintype =="kin1b") thcm_cent=155;
    if (kintype =="kin2b") thcm_cent=141;
    if (kintype =="kin2a") thcm_cent=141;
    if (kintype =="kin2b_pass1") thcm_cent=141;
    if (kintype =="kin3b") thcm_cent=120.;
    if (kintype =="kin1a") thcm_cent_2=145;
    if (kintype =="kin1b") thcm_cent_2=145;
    if (kintype =="kin2b") thcm_cent_2=134;
    if (kintype =="kin2a") thcm_cent_2=134;
    if (kintype =="kin2b_pass1") thcm_cent_2=134;
    if (kintype =="kin3b") thcm_cent_2=112.;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (Nrun!=CurRun) {
		  CurRun=Nrun;
		  TotCharge+=Cur*RunTime/1000;
		  cout << " Run = " << CurRun << " time = " << RunTime << " " << Cur << " " << Cur*RunTime << " Total charge = " << TotCharge << " mC" << endl;
		}
                if (i%50000==0) cout << " Entry = " << i << endl;
		
		hCtime->Fill(ctime);
		if (abs(ctime-Ctime_cent)<3*Ctime_sig) {
		  hW->Fill(W);
		  hQ2->Fill(Q2);
		  hMM2_W->Fill(MM2,W); 
		  hMM2->Fill(MM2); 
		  if (abs(MM2-pi_cent)<pi_sig) hW_pi0->Fill(W);
		    if (abs(W-W_cent)<W_sig && abs(Q2-Q2_cent)<Q2_sig && abs(Thcm-thcm_cent)<thcm_sig ) {
		           hMM2_vcs_cut->Fill(MM2); 
		    }
		    if (abs(W-W_cent)<W_sig && abs(Q2-Q2_cent_2)<Q2_sig_2 && abs(Thcm-thcm_cent_2)<thcm_sig ) {
		           hMM2_vcs_hiQ2_cut->Fill(MM2); 
		    }
		    if (abs(MM2-vcs_cent)<vcs_sig) {
		    hW_vcs->Fill(W);
                    hThcm_vcs->Fill(Thcm);
		    if (abs(W-W_cent)<W_sig && abs(Q2-Q2_cent)<Q2_sig && abs(Thcm-thcm_cent)<thcm_sig ) {
                           hThcm_vcs_cut->Fill(Thcm);
		    }
		  }
		    if (abs(W-W_cent)<W_sig && abs(Q2-Q2_cent)<Q2_sig ) {
		      hMM2_Wcut->Fill(MM2);
		      hMM2_thcm->Fill(MM2,Thcm);
		    } 
		    if (abs(W-W_cent)<W_sig && abs(Q2-Q2_cent_2)<Q2_sig_2 ) {
		      hMM2_thcm_vcs_hiQ2_cut->Fill(MM2,Thcm);
		    } 
		}
		if (abs(ctime-Ctime_ran_cent)<6*Ctime_ran_sig) {
		  if (abs(W-W_cent)<W_sig && abs(Q2-Q2_cent)<Q2_sig && abs(Thcm-thcm_cent)<thcm_sig ) hMM2_ran->Fill(MM2); 
		  if (abs(W-W_cent)<W_sig && abs(Q2-Q2_cent_2)<Q2_sig_2 && abs(Thcm-thcm_cent_2)<thcm_sig ) hMM2_ran_vcs_hiQ2_cut->Fill(MM2); 
		}
	}
//
	cout << kintype << " INtegral = " << hThcm_vcs_cut->Integral() << " Yield cnt/C  = " <<  hThcm_vcs_cut->Integral()/(TotCharge/1000.) << endl;
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
} 
