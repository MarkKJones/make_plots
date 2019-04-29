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

void make_hist_shms_cal(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_ngcer_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t preshow_adcint[2][14];
   tsimc->SetBranchAddress("P.cal.pr.goodNegAdcPulseInt",&preshow_adcint[0]) ;
   tsimc->SetBranchAddress("P.cal.pr.goodPosAdcPulseInt",&preshow_adcint[1]) ;
 Double_t preshow_adcamp[2][14];
   tsimc->SetBranchAddress("P.cal.pr.goodNegAdcPulseAmp",&preshow_adcamp[0]) ;
   tsimc->SetBranchAddress("P.cal.pr.goodPosAdcPulseAmp",&preshow_adcamp[1]) ;
 Double_t preshow_energy[2][14];
   tsimc->SetBranchAddress("P.cal.pr.eneg",&preshow_energy[0]) ;
   tsimc->SetBranchAddress("P.cal.pr.epos",&preshow_energy[1]) ;
Double_t npeSum;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&npeSum) ;
 Double_t etottracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etottracknorm) ;
Double_t xAtCer;
   tsimc->SetBranchAddress("P.ngcer.xAtCer",&xAtCer) ;
 Double_t yAtCer;
   tsimc->SetBranchAddress("P.ngcer.yAtCer",&yAtCer) ;
Double_t xAtCal;
   tsimc->SetBranchAddress("P.cal.xtrack",&xAtCal) ;
 Double_t yAtCal;
   tsimc->SetBranchAddress("P.cal.ytrack",&yAtCal) ;
 Double_t ntrack;
   tsimc->SetBranchAddress("P.gtr.index",&ntrack) ;
   //
   TString temp;
   //
   const char* side[2]={"Neg","Pos"};
   TH1F *h_preshow_adcint[2][14];
   TH1F *h_preshow_adcint_pion[2][14];
   TH1F *h_preshow_adcamp[2][14];
   TH1F *h_preshow_adcamp_pion[2][14];
   TH1F *h_preshow_energy[2][14];
   TH1F *h_preshow_energy_pion[2][14];
   for (Int_t ns=0;ns<2;ns++) {
   for (Int_t nr=0;nr<14;nr++) {
     temp=Form("Run %d ; %s Row %d Pulse Int ; Counts",nrun,side[ns],nr+1);
     h_preshow_adcint[ns][nr] = new TH1F(Form("h_preshow_adcint_%d_%d",ns,nr),temp,100,0,100.0);
     temp=Form("Run %d Pion ; %s Row %d Pulse Int ; Counts",nrun,side[ns],nr+1);
     h_preshow_adcint_pion[ns][nr] = new TH1F(Form("h_preshow_adcint_pion_%d_%d",ns,nr),temp,100,0,100.0);
     temp=Form("Run %d ; %s Row %d Pulse Amp ; Counts",nrun,side[ns],nr+1);
     h_preshow_adcamp[ns][nr] = new TH1F(Form("h_preshow_adcamp_%d_%d",ns,nr),temp,100,0,100.0);
     temp=Form("Run %d Pion ; %s Row %d Pulse Amp ; Counts",nrun,side[ns],nr+1);
     h_preshow_adcamp_pion[ns][nr] = new TH1F(Form("h_preshow_adcamp_pion_%d_%d",ns,nr),temp,100,0,100.0);
     temp=Form("Run %d ; %s Row %d Energy ; Counts",nrun,side[ns],nr+1);
     h_preshow_energy[ns][nr] = new TH1F(Form("h_preshow_energy_%d_%d",ns,nr),temp,100,0,.5);
     temp=Form("Run %d Pion ; %s Row %d Energy ; Counts",nrun,side[ns],nr+1);
     h_preshow_energy_pion[ns][nr] = new TH1F(Form("h_preshow_energy_pion_%d_%d",ns,nr),temp,100,0,.5);
     HList.Add(h_preshow_energy_pion[ns][nr]);
   }}
   temp=Form("Run %d ; Etottracknorm  ; Counts",nrun);
   TH1F *hetotnorm = new TH1F("hetotnorm",temp,100,0,2.0);
   HList.Add(hetotnorm);
   temp=Form("Run %d ; Etottracknorm (npe cut)  ; Counts",nrun);
   TH1F *hetotnorm_npecut = new TH1F("hetotnorm_npecut",temp,100,0,2.0);
   HList.Add(hetotnorm_npecut);
   temp=Form("Run %d ; NpeSUm  ; Counts",nrun);
   TH1F *hcernpeSum = new TH1F("hcernpeSum",temp,160,0,40.);
   HList.Add(hcernpeSum);
   temp=Form("Run %d ; Cal  Y ; Cal X ",nrun);
   TH2F *h_calXY = new TH2F("h_calXY",temp,80,-80,80.,80,-80,80.0);
   TH2F *h_calXY_now= new TH2F("h_calXY_now",temp,80,-80,80.,80,-80,80.0);
   HList.Add(h_calXY);
   HList.Add(h_calXY_now);
   temp=Form("Run %d ; Cer  Y ; Cer X ",nrun);
   TH2F *h_cerXY = new TH2F("h_cerXY",temp,50,-50,50.,50,-50,50.0);
   TH2F *h_cerXY_now= new TH2F("h_cerXY_now",temp,50,-50,50.,50,-50,50.0);
   HList.Add(h_cerXY);
   HList.Add(h_cerXY_now);
  //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		if (1==1) {
		hetotnorm->Fill(etottracknorm);
		hcernpeSum->Fill(npeSum);
                for (Int_t ns=0;ns<2;ns++) {
                for (Int_t nr=0;nr<14;nr++) {
		  h_preshow_adcint[ns][nr]->Fill(preshow_adcint[ns][nr]);
		  h_preshow_adcamp[ns][nr]->Fill(preshow_adcamp[ns][nr]);
		  h_preshow_energy[ns][nr]->Fill(preshow_energy[ns][nr]);
		  if  (TMath::Abs(etottracknorm-.25) < .1 && npeSum==0)  {
                     h_preshow_adcint_pion[ns][nr]->Fill(preshow_adcint[ns][nr]);
                     h_preshow_adcamp_pion[ns][nr]->Fill(preshow_adcamp[ns][nr]);
		     h_preshow_energy_pion[ns][nr]->Fill(preshow_energy[ns][nr]);
		  }
		}}
		if (npeSum>2.) {
                hetotnorm_npecut->Fill(etottracknorm);
		h_calXY->Fill(yAtCal,xAtCal,etottracknorm);
		h_calXY_now->Fill(yAtCal,xAtCal);
		}
		if (etottracknorm>.8) {
		h_cerXY->Fill(yAtCer,xAtCer,npeSum);
		h_cerXY_now->Fill(yAtCer,xAtCer);
		}
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
        HList.Write();
	//
}
