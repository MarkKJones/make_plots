
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
using namespace std;

void make_hist_hms_cer(TString basename="",Int_t nrun=2043){
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
   inputroot="hallc_online-ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_cer_hist.root";
 TObjArray HList(0);
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t cer_pamp[2];
 tsimc->SetBranchAddress("H.cer.goodAdcPulseAmp",cer_pamp) ;
 Double_t cer_mult[2];
 tsimc->SetBranchAddress("H.cer.goodAdcMult",cer_mult) ;
 Double_t cer_pint[2];
 tsimc->SetBranchAddress("H.cer.goodAdcPulseInt",cer_pint) ;
 Double_t cer_npe[2];
 tsimc->SetBranchAddress("H.cer.npe",cer_npe) ;
 Double_t cer_npesum;
 tsimc->SetBranchAddress("H.cer.npeSum",&cer_npesum) ;
 Double_t evtyp;
 tsimc->SetBranchAddress("g.evtyp",&evtyp) ;
 Double_t xfp;
 tsimc->SetBranchAddress("H.dc.x_fp",&xfp) ;
 Double_t xpfp;
 tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp) ;
 Double_t el_real;
 tsimc->SetBranchAddress("T.hms.hEL_REAL_tdcTime",&el_real) ;
 Double_t delta;
 tsimc->SetBranchAddress("H.gtr.dp",&delta) ;
 Double_t etracknorm;
 tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm) ;
   // Define histograms
 TH1F *cer1_pamp = new TH1F("cer1_pamp","; PMT1 Pulse Amp; counts ",100,0.,100.);
 HList.Add(cer1_pamp);
 TH1F *cer1_pamp_nocut = new TH1F("cer1_pamp_nocut","; PMT1 Pulse Amp; counts ",100,0.,200.);
 HList.Add(cer1_pamp_nocut);
 TH1F *cer2_pamp = new TH1F("cer2_pamp"," PMT2 Pulse Amp; counts",100,0.,100.);
 HList.Add(cer2_pamp);
 TH1F *cer2_pamp_nocut = new TH1F("cer2_pamp_nocut"," PMT2 Pulse Amp; counts",100,0.,200.);
 HList.Add(cer2_pamp_nocut);
 TH1F *cer1_pint = new TH1F("cer1_pint","; PMT1 Pulse Int; counts ",100,0.,100.);
 HList.Add(cer1_pint);
 TH1F *cer2_pint = new TH1F("cer2_pint"," PMT2 Pulse Int; counts",100,0.,100.);
 HList.Add(cer2_pint);
 TH1F *cer1_npe = new TH1F("cer1_npe","; PMT1 Npe; counts ",80,0.,20.);
 TH1F *cer1_npe_cut = new TH1F("cer1_npe_cut","; PMT1 Npe; counts ",80,0.,20.);
 TH2F *cer1_npe_xfp = new TH2F("cer1_npe_xfp","; PMT1 Npe; X at Z=230 ",80,0.,20.,80,-100,100);
 TH2F *cer1_npe_delta = new TH2F("cer1_npe_delta","; PMT1 Npe; Delta ",80,0.,20.,96,-15,15);
 TH2F *cer2_npe_xfp = new TH2F("cer2_npe_xfp","; PMT2 Npe; X at Z=230 ",80,0.,20.,80,-100,100);
 TH2F *cer2_npe_delta = new TH2F("cer2_npe_delta","; PMT2 Npe; Delta ",80,0.,20.,80,-15,15);
 TH1F *cer2_npe = new TH1F("cer2_npe"," PMT2 Npe; counts",80,0.,20.);
 TH1F *cer2_npe_cut = new TH1F("cer2_npe_cut"," PMT2 Npe; counts",80,0.,20.);
 TH1F *cer1_pamp_elcut = new TH1F("cer1_pamp_elcut","; PMT1 Pulse Amp (EL_REAL); counts ",200,0.,100.);
 TH1F *cer2_pamp_elcut = new TH1F("cer2_pamp_elcut"," PMT2 Pulse Amp (EL_REAL) ; counts",200,0.,100.);
 TH1F *h_etotnorm = new TH1F("h_etotnorm"," Etot track norm; counts",200,0.,1.5);
 TH1F *h_elreal = new TH1F("h_elreal"," EL REAL time; counts",200,-100.,1000.);
			 //
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
                if (etracknorm>.8 && TMath::Abs(delta) <15.) {
		  if (    cer_pamp[0]>0 && TMath::Abs(xfp+xpfp*230-5) < 5) cer1_pamp->Fill(cer_pamp[0]);
		  if (   cer_pamp[1]>0 && TMath::Abs(xfp+xpfp*230-10) <5)cer2_pamp->Fill(cer_pamp[1]);
		  if (  cer_pamp[0]>0 && TMath::Abs(xfp+xpfp*230-5) < 5) cer1_pint->Fill(cer_pint[0]);
		  if (  cer_pamp[1]>0 && TMath::Abs(xfp+xpfp*230-10) <5)cer2_pint->Fill(cer_pint[1]);
                if (  cer_npe[0]>0 ) cer1_npe->Fill(cer_npe[0]);
                if (  cer_npe[0]>0 ) cer1_pamp_nocut->Fill(cer_pamp[0]);
                if (  cer_npe[1]>0 ) cer2_pamp_nocut->Fill(cer_pamp[1]);
                if (  cer_npe[0]>0 && TMath::Abs(xfp+xpfp*230-5) < 5)cer1_npe_cut->Fill(cer_npe[0]);
                if ( cer_npe[0]>0)cer1_npe_xfp->Fill(cer_npe[0],xfp+xpfp*230);
                if (  cer_npe[0]>0)cer1_npe_delta->Fill(cer_npe[0],delta);
                if (  cer_npe[1]>0)cer2_npe_delta->Fill(cer_npe[1],delta);
                if (   cer_npe[1]>0)cer2_npe_xfp->Fill(cer_npe[1],xfp+xpfp*230);
                if (  cer_npe[1]>0 ) cer2_npe->Fill(cer_npe[1]);
                if (  cer_npe[1]>0 && TMath::Abs(xfp+xpfp*230-10) <5)cer2_npe_cut->Fill(cer_npe[1]);
		h_elreal->Fill(el_real);
		if (el_real >0) {
                if (  cer_pamp[0]>0) cer1_pamp_elcut->Fill(cer_pamp[0]);
                if (  cer_pamp[1]>0) cer2_pamp_elcut->Fill(cer_pamp[1]);
		}
		}
		h_etotnorm->Fill(etracknorm);
  	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
