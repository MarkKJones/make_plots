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

void make_hist_coin_edtm(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_edtm_hist.root";
 TObjArray HList(0);
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Int_t  nedtm;
   tsimc->SetBranchAddress("Ndata.D.pedtm",&nedtm);
 Double_t  edtm[5];
   tsimc->SetBranchAddress("D.pedtm",edtm);
 Int_t  nptrigref1;
   tsimc->SetBranchAddress("Ndata.D.ptrigref1",&nptrigref1);
 Double_t  ptrigref1[5];
   tsimc->SetBranchAddress("D.ptrigref1",ptrigref1);
 Int_t  nptrigref2;
   tsimc->SetBranchAddress("Ndata.D.ptrigref2",&nptrigref2);
 Double_t  ptrigref2[5];
   tsimc->SetBranchAddress("D.ptrigref2",ptrigref2);
 Int_t  npdcref1;
   tsimc->SetBranchAddress("Ndata.D.pdcref1",&npdcref1);
 Double_t  pdcref1[5];
   tsimc->SetBranchAddress("D.pdcref1",pdcref1);
 Int_t  npdcref10;
   tsimc->SetBranchAddress("Ndata.D.pdcref10",&npdcref10);
 Double_t  pdcref10[5];
   tsimc->SetBranchAddress("D.pdcref10",pdcref10);
   //
 Int_t  nhdcref1;
   tsimc->SetBranchAddress("Ndata.D.hdcref1",&nhdcref1);
 Int_t  nhdcref4;
   tsimc->SetBranchAddress("Ndata.D.hdcref4",&nhdcref4);
 Double_t  hdcref1[5];
   tsimc->SetBranchAddress("D.hdcref1",hdcref1);
 Double_t  hdcref4[5];
   tsimc->SetBranchAddress("D.hdcref4",hdcref4);
  Int_t  nhtrigref1;
   tsimc->SetBranchAddress("Ndata.D.htrigref1",&nhtrigref1);
 Double_t  htrigref1[5];
   tsimc->SetBranchAddress("D.htrigref1",htrigref1);
 Int_t  nhtrigref2;
   tsimc->SetBranchAddress("Ndata.D.htrigref2",&nhtrigref2);
 Double_t  htrigref2[5];
   tsimc->SetBranchAddress("D.htrigref2",htrigref2);

   // Define histograms
   TH1F *hdiffpref1[3];
   TH1F *hdiffpref2[3];
   TH1F *hdiffhref1[3];
   TH1F *hdiffhref2[3];
   TH1F *hdiffhdcref1[3];
   TH1F *hdiffhdcref4[3];
   TH1F *hdiffpdcref1[3];
   TH1F *hdiffpdcref10[3];
   for (Int_t i=0;i<2;i++) {
     hdiffpref1[i]= new TH1F(Form("hdiffpref1_%d",i),Form(" SHMS Trig ref 1; Time diff pulse %d - 1  (ns); Counts ",i+2),3000,0,300);
   HList.Add( hdiffpref1[i]);
     hdiffpdcref1[i]= new TH1F(Form("hdiffpdcref1_%d",i),Form(" SHMS DC ref 1; Time diff pulse %d - 1  (ns); Counts ",i+2),3000,0,300);
   HList.Add( hdiffpdcref1[i]);
     hdiffpdcref10[i]= new TH1F(Form("hdiffpdcref10_%d",i),Form(" SHMS DC ref 10; Time diff pulse %d - 1  (ns); Counts ",i+2),3000,0,300);
   HList.Add( hdiffpdcref10[i]);
     hdiffpref2[i]= new TH1F(Form("hdiffpref2_%d",i),Form(" SHMS Trig Ref 2; Time diff Pulse %d -1 (ns); Counts ",i+2),3000,0,300);
   HList.Add( hdiffpref2[i]);
     hdiffhref1[i]= new TH1F(Form("hdiffhref1_%d",i),Form("HMS Trig Ref 1; Time diff Pulse %d-1 (ns); Counts ",i+2),3000,0,300);
   HList.Add( hdiffhref1[i]);
     hdiffhdcref1[i]= new TH1F(Form("hdiffhdcref1_%d",i),Form("HMS DC Ref 1; Time diff Pulse %d-1 (ns); Counts ",i+2),3000,0,300);
   HList.Add( hdiffhdcref1[i]);
     hdiffhdcref4[i]= new TH1F(Form("hdiffhdcref4_%d",i),Form("HMS DC Ref 4; Time diff Pulse %d-1 (ns); Counts ",i+2),3000,0,300);
   HList.Add( hdiffhdcref4[i]);
     hdiffhref2[i]= new TH1F(Form("hdiffhref2_%d",i),Form("HMS Trig Ref 2; Time diff  Pulse %d-1 (ns); Counts ",i+2),4000,0,400);
   HList.Add( hdiffhref2[i]);
   }
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		if (nptrigref1>=2) hdiffpref1[0]->Fill(.1*(ptrigref1[1]-ptrigref1[0]));
		if (nptrigref1>=3)hdiffpref1[1]->Fill(.1*(ptrigref1[2]-ptrigref1[0]));
		if (npdcref1>=2) hdiffpdcref1[0]->Fill(.1*(pdcref1[1]-pdcref1[0]));
		if (npdcref1>=3)hdiffpdcref1[1]->Fill(.1*(pdcref1[2]-pdcref1[0]));
		if (npdcref10>=2) hdiffpdcref10[0]->Fill(.1*(pdcref10[1]-pdcref10[0]));
		if (npdcref10>=3)hdiffpdcref10[1]->Fill(.1*(pdcref10[2]-pdcref10[0]));
		if (nptrigref2>=2)hdiffpref2[0]->Fill(.1*(ptrigref2[1]-ptrigref2[0]));
		if (nptrigref2>=3)hdiffpref2[1]->Fill(.1*(ptrigref2[2]-ptrigref2[0]));
		if (nhtrigref1>=2)hdiffhref1[0]->Fill(.1*(htrigref1[1]-htrigref1[0]));
		if (nhtrigref1>=3)hdiffhref1[1]->Fill(.1*(htrigref1[2]-htrigref1[0]));
		if (nhtrigref2>=2)hdiffhref2[0]->Fill(.1*(htrigref2[1]-htrigref2[0]));
		if (nhtrigref2>=3)hdiffhref2[1]->Fill(.1*(htrigref2[2]-htrigref2[0]));
		if (nhdcref1>=2) hdiffhdcref1[0]->Fill(.1*(hdcref1[1]-hdcref1[0]));
		if (nhdcref1>=3)hdiffhdcref1[1]->Fill(.1*(hdcref1[2]-hdcref1[0]));
		if (nhdcref4>=2) hdiffhdcref4[0]->Fill(.1*(hdcref4[1]-hdcref4[0]));
		if (nhdcref4>=3)hdiffhdcref4[1]->Fill(.1*(hdcref4[2]-hdcref4[0]));
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
