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
#include<vector>
using namespace std;

void calib_bcm(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_calib_bcm.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("TSP");
//
 Double_t Unserrate ;
   tsimc->SetBranchAddress("P.Unser.scalerRate",&Unserrate);
 Double_t BCM1rate ;
   tsimc->SetBranchAddress("P.BCM1.scalerRate",&BCM1rate);
 Double_t BCM2rate ;
   tsimc->SetBranchAddress("P.BCM2.scalerRate",&BCM2rate);
 Double_t BCM4Arate ;
   tsimc->SetBranchAddress("P.BCM4A.scalerRate",&BCM4Arate);
 Double_t BCM4Brate ;
   tsimc->SetBranchAddress("P.BCM4B.scalerRate",&BCM4Brate);
 Double_t BCM4Crate ;
   tsimc->SetBranchAddress("P.BCM4C.scalerRate",&BCM4Crate);

//
 vector<double> vUnserrate;
 vector<double> vBCM1rate;
 vector<double> vBCM2rate;
 vector<double> vBCM4Arate;
 vector<double> vBCM4Brate;
 vector<double> vBCM4Crate;
 //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		vUnserrate.push_back(Unserrate);
		vBCM1rate.push_back(BCM1rate);
		vBCM2rate.push_back(BCM2rate);
		vBCM4Arate.push_back(BCM4Arate);
		vBCM4Brate.push_back(BCM4Brate);
		vBCM4Crate.push_back(BCM4Crate);
                if (i%50000==0) cout << " Entry = " << i << endl;
	}
	//
	for (int i = 0; i <vUnserrate.size() ; i++) {
	  cout << i << "," << vUnserrate[i] << "," << vBCM1rate[i] << "," << vBCM2rate[i] << "," << vBCM4Arate[i] << "," << vBCM4Brate[i] << "," << vBCM4Crate[i] << endl;
	}

}
