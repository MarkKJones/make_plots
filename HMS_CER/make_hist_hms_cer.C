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

void make_hist_hms_cer(TString basename="",Int_t nrun=1272){
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
   outputhist= "hist/"+basename+"_hms_cer_hist.root";
 TObjArray HList(0);
   //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 static const Int_t plnum=4;
 static const Int_t ns=2;
 const char* sname[ns]={"Neg","Pos"};
 const char* sname2[ns]={"neg","pos"};
 const char* plname[plnum]={"1pr","2ta","3ta","4ta"};
 Int_t  e_nhits[plnum][ns];
 Double_t  energy[plnum][ns][13];
  for (Int_t ip=0;ip<plnum;ip++) {
 for (Int_t is=0;is<ns;is++) {
   tsimc->SetBranchAddress(Form("Ndata.H.cal.%s.e%s",plname[ip],sname2[is]),&e_nhits[ip][is]) ;   
   tsimc->SetBranchAddress(Form("H.cal.%s.e%s",plname[ip],sname2[is]),&energy[ip][is]) ;   
 }}
 Double_t npeSum;
   tsimc->SetBranchAddress("H.cer.npeSum",&npeSum) ;
 Double_t etottracknorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&etottracknorm) ;
   //
   TString temp=Form("Run %d ; NpeSUm  ; Counts",nrun);
   TH1F *hcernpeSum = new TH1F("hcernpeSum",temp,160,0,40.);
   HList.Add(hcernpeSum);
   temp=Form("Run %d ; Etottracknorm  ; Counts",nrun);
   TH1F *hEtottracknorm = new TH1F("Etottracknorm",temp,150,0,1.5);
   HList.Add(hEtottracknorm);
   TH1F *HCalEnergy[plnum][ns][13];
   TH1F *HCalEnergy_elec[plnum][ns][13];
  for (Int_t ip=0;ip<plnum;ip++) {
 for (Int_t is=0;is<ns;is++) {
  for (Int_t ih=0;ih<13;ih++) {
    HCalEnergy[ip][is][ih] = new TH1F(Form("cal_pion_%s_%s_%d",plname[ip],sname2[is],ih)," ; Pion Energy",200,0.,1.);
    HCalEnergy_elec[ip][is][ih] = new TH1F(Form("cal_elec_%s_%s_%d",plname[ip],sname2[is],ih)," ; Elec Energy",200,0.,1.);
  }}}   //
   //
Long64_t nentries = tsimc->GetEntries();
// nentries=10;
	for (int i = 0; i < nentries; i++) {
	  //	for (int i = 0; i < 10; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hcernpeSum->Fill(npeSum);
		if (npeSum>2.) {
		  hEtottracknorm->Fill(etottracknorm);
		}
		 
 for (Int_t is=0;is<ns;is++) {
  for (Int_t ip=0;ip<plnum;ip++) {
  for (Int_t ih=0;ih<e_nhits[ip][is];ih++) {
    if (npeSum==0. && energy[ip][is][ih] > 0) HCalEnergy[ip][is][ih]->Fill(energy[ip][is][ih]);
    if (npeSum>2. && etottracknorm>1. && energy[ip][is][ih] > 0) HCalEnergy_elec[ip][is][ih]->Fill(energy[ip][is][ih]);
  }}   //
		  
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
