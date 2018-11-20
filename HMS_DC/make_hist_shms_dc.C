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

void make_hist_shms_tw(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 static const Int_t plnum=4;
 static const Int_t iside=2;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 const char* sidename[isde]={"Neg","Pos"};
 static const Int_t npad[plnum]={16,10,16,10};
 //
 Double_t tw_corr[plnum][iside][16];
 Double_t pulseamp[plnum][iside[16];
 for (Int_t ipl=0;ipl<plnum;ipl++) {
 for (Int_t is=0;is<iside;is++) {
 for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
   tsimc->SetBranchAddress(Form("P.hod.%s.Good%sAdcTimeWalkCorr",plname[ipl],sidename[is],&tw_corr[ipl][is][ipad]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Good%sAdcPulseAmp",plname[ipl],sidename[is],&pulseamp[ipl][is][ipad]) ;
 }}  
 }
 //
 TH2F *hTW_adc[plnum][is][ipad];
 for (Int_t ipl=0;ipl<plnum;ipl++) {
 for (Int_t is=0;is<iside;is++) {
 for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
   hTW_adc[ipl][is][ipad]= new TH2F(Form("tw_adc_%s_%s_pad_%d",plname[ipl],sidename[is],ipad),Form("%s %spad_%d; Adc Amp; TDc TW corr",plname[ipl],sidename[is],ipad),100,0,200.,100,-50,50);
   HList.Add(Form("tw_adc_%s_%s_pad_%d",plname[ipl],sidename[is],ipad));
 }}  
 }
 //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
 for (Int_t ipl=0;ipl<plnum;ipl++) {
 for (Int_t is=0;is<iside;is++) {
 for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
   hTW_adc[ipl][is][ipad]->Fill(pulseamp[ipl][is][ipad],tw_corr[ipl][is][ipad]);
 }}}		
	}
//
}
