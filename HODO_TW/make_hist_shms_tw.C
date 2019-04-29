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
 const char* sidename[iside]={"Neg","Pos"};
 const char* sidename2[iside]={"neg","pos"};
 static const Int_t npad[plnum]={13,13,14,21};
 //
 Double_t plhits[plnum];
 Int_t pladchits[plnum][iside];
 Int_t pltdchits[plnum][iside];
 Double_t pltdcpad[plnum][iside][100];
 Double_t pladcpad[plnum][iside][100];
 Double_t tw_corr[plnum][iside][21];
 Double_t tw_uncorr[plnum][iside][21];
 Double_t pulseamp[plnum][iside][21];
  Double_t pcal_etrkNorm;
  Double_t pngcer_npeSum;
  Double_t pdc_ntrack;
  Double_t beta;
  TString npcal_etrkNorm = "P.cal.etracknorm";
  TString npngcer_npeSum = "P.hgcer.npeSum";
  TString npdc_ntrack = "P.dc.ntrack";
  TString nbeta = "P.hod.beta";
  Double_t etrknrm_low_cut = 0.7;
  Double_t npngcer_npeSum_low_cut = 0.7;
  Double_t betanotrack_low_cut = 0.5;
  Double_t betanotrack_hi_cut = 1.5;
  TString nTrackXPos;
  TString nTrackYPos;
  Double_t TrackXPos[plnum];
  Double_t TrackYPos[plnum];
  Double_t xfp;
  Double_t xpfp;
  Double_t yfp;
  Double_t ypfp;
  Double_t z0[plnum] = {52.1, 61.7, 271.4, 282.4};     //zentral z dist from hod planes to focal plane
  tsimc->SetBranchAddress("P.dc.x_fp", &xfp);
  tsimc->SetBranchAddress("P.dc.y_fp", &yfp);
  tsimc->SetBranchAddress("P.dc.xp_fp", &xpfp);
  tsimc->SetBranchAddress("P.dc.yp_fp", &ypfp);
  tsimc->SetBranchAddress(npcal_etrkNorm, &pcal_etrkNorm);
  tsimc->SetBranchAddress(npngcer_npeSum, &pngcer_npeSum);
  tsimc->SetBranchAddress(npdc_ntrack, &pdc_ntrack);
  tsimc->SetBranchAddress(nbeta, &beta);
 for (Int_t ipl=0;ipl<plnum;ipl++) {
   tsimc->SetBranchAddress(Form("P.hod.%s.nhits",plname[ipl]),&plhits[ipl]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.TrackXPos",plname[ipl]),&TrackXPos[ipl]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.TrackYPos",plname[ipl]),&TrackYPos[ipl]) ;
 for (Int_t is=0;is<iside;is++) {
   tsimc->SetBranchAddress(Form("P.hod.%s.Good%sTdcTimeWalkUnCorr",plname[ipl],sidename[is]),&tw_uncorr[ipl][is]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Good%sTdcTimeWalkCorr",plname[ipl],sidename[is]),&tw_corr[ipl][is]) ;
   tsimc->SetBranchAddress(Form("Ndata.P.hod.%s.%sAdcCounter",plname[ipl],sidename2[is]),&pladchits[ipl][is]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.%sAdcCounter",plname[ipl],sidename2[is]),&pladcpad[ipl][is]) ;
   tsimc->SetBranchAddress(Form("Ndata.P.hod.%s.%sTdcCounter",plname[ipl],sidename2[is]),&pltdchits[ipl][is]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.%sTdcCounter",plname[ipl],sidename2[is]),&pltdcpad[ipl][is]) ;
     tsimc->SetBranchAddress(Form("P.hod.%s.Good%sAdcPulseAmp",plname[ipl],sidename[is]),&pulseamp[ipl][is]) ;
 }}  
 //
 TH1F *hPlane_hits[plnum];
 TH1F *hPlane_padhits[plnum];
 TH1F *hPlane_Xpos[plnum];
 TH1F *hPlane_Ypos[plnum];
 TH1F *hPlane_padhits_cut1[plnum];
 TH1F *hPlane_padhits_cut2[plnum];
 TH2F *hTW_adc[plnum][iside][21];
 for (Int_t ipl=0;ipl<plnum;ipl++) {
   hPlane_hits[ipl]= new TH1F(Form("Nhits_%s",plname[ipl]),Form("%s ; Nhits; Counts",plname[ipl]),10,0,10);
   HList.Add(hPlane_hits[ipl]);
   hPlane_Xpos[ipl]= new TH1F(Form("Xpos_%s",plname[ipl]),Form("%s ; Xpos; Counts",plname[ipl]),100,-300,300);
   HList.Add(hPlane_Xpos[ipl]);
   hPlane_Ypos[ipl]= new TH1F(Form("Ypos_%s",plname[ipl]),Form("%s ; Ypos; Counts",plname[ipl]),100,-100,100);
   HList.Add(hPlane_Ypos[ipl]);
   hPlane_padhits[ipl]= new TH1F(Form("padhits_%s",plname[ipl]),Form("%s ; Paddle; Counts",plname[ipl]),npad[ipl],0,npad[ipl]);
   HList.Add(hPlane_padhits[ipl]);
   hPlane_padhits_cut1[ipl]= new TH1F(Form("padhits_cut1_%s",plname[ipl]),Form("Cut1 %s ; Paddle; Counts",plname[ipl]),npad[ipl],0,npad[ipl]);
   HList.Add(hPlane_padhits_cut1[ipl]);
   hPlane_padhits_cut2[ipl]= new TH1F(Form("padhits_cut2_%s",plname[ipl]),Form("Cut1 %s ; Paddle; Counts",plname[ipl]),npad[ipl],0,npad[ipl]);
   HList.Add(hPlane_padhits_cut2[ipl]);
 for (Int_t is=0;is<iside;is++) {
 for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
   //hTW_adc[ipl][is][ipad]= new TH2F(Form("tw_adc_%s_%s_pad_%d",plname[ipl],sidename[is],ipad),Form("%s %spad_%d; Adc Amp; TDc TW corr",plname[ipl],sidename[is],ipad),100,0,200.,100,-50,50);
   //HList.Add(hTW_adc[ipl][is][ipad]);
 }}  
 }
 //
Long64_t nentries = tsimc->GetEntries();
// nentries=1000;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		Bool_t hodTrk = TMath::Abs(xfp+xpfp*z0[0])<200&&TMath::Abs(yfp+ypfp*z0[0])<200&&TMath::Abs(xfp+xpfp*z0[1])<200&&TMath::Abs(yfp+ypfp*z0[1])<200&&TMath::Abs(xfp+xpfp*z0[2])<200&&TMath::Abs(yfp+ypfp*z0[2])<200&&TMath::Abs(xfp+xpfp*z0[3])<200&&TMath::Abs(yfp+ypfp*z0[3])<200;
      Bool_t pcal = pcal_etrkNorm>etrknrm_low_cut;
      Bool_t pngcer = pngcer_npeSum>npngcer_npeSum_low_cut;
      Bool_t pdctrk = pdc_ntrack>0.0;
      Bool_t betaCut = beta>betanotrack_low_cut&& beta<betanotrack_hi_cut;
      Bool_t  pid_pelec = pcal&&pngcer&&pdctrk;
		if (plhits[0]==1&&plhits[1]==1&&plhits[2]==1&&plhits[3]==1) {
 for (Int_t ipl=0;ipl<plnum;ipl++) {
   hPlane_hits[ipl]->Fill(plhits[ipl]);
   hPlane_Xpos[ipl]->Fill(xfp+xpfp*z0[ipl]);
   hPlane_Ypos[ipl]->Fill(yfp+ypfp*z0[ipl]);
 for (Int_t is=0;is<iside;is++) {
   for (Int_t nh=0;nh<pltdchits[ipl][is];nh++) {
     hPlane_padhits[ipl]->Fill(pltdcpad[ipl][is][nh]);
     if(pid_pelec&&betaCut) hPlane_padhits_cut1[ipl]->Fill(pltdcpad[ipl][is][nh]);
     if(pid_pelec) hPlane_padhits_cut2[ipl]->Fill(pltdcpad[ipl][is][nh]);
   }
   /*
 for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
   //hTW_adc[ipl][is][ipad]->Fill(pulseamp[ipl][is][ipad],tw_corr[ipl][is][ipad]);
 }
   */
}}		
		}
	}
//
 }
