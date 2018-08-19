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

void make_hist_simc_hms_single_arm(TString basename="",Int_t nrun=1272){
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
   outputhist= "hist/"+basename+"_simc_hms_ep_elastic_hist.root";
 TObjArray HList(0);
   //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("h666");
// Define branches
 Float_t  W;
   tsimc->SetBranchAddress("W",&W);
 Float_t  Weight;
   tsimc->SetBranchAddress("Weight",&Weight);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 300, 0.9,1.2);
    HList.Add(hW);
    //
    Double_t Normfac;
    Double_t Nent_simc;
    Double_t Exp_charge; // assume SIMC charge is 1 mC
    Double_t Exp_eff; // Livetime*Track_eff*Hodo_eff
    Double_t simc_fac=1.;
    if (nrun ==1272) {
      Nent_simc=100000.;
      Normfac = 0.107289E+08;
      Exp_charge= 1356.930/1000.;
      Exp_eff = 5.6348/100.*0.9856*1.0;
      simc_fac = Normfac*Exp_charge*Exp_eff/Nent_simc;
    }
    //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hW->Fill(W,Weight*simc_fac);
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
