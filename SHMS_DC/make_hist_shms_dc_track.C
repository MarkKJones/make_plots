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

void make_hist_shms_dc_track(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_dc_track_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 static const Int_t maxtracks=10;
 Int_t ntracks;
 Double_t flag[maxtracks];
 Double_t trackn[maxtracks];
 Double_t trackchi[maxtracks];
 tsimc->SetBranchAddress("P.tr.flag",&flag) ;
   tsimc->SetBranchAddress("Ndata.P.tr.flag",&ntracks) ;
   tsimc->SetBranchAddress("P.tr.n",&trackn) ;
   tsimc->SetBranchAddress("P.tr.chi2",&trackchi) ;
 Double_t gindex;
   tsimc->SetBranchAddress("P.gtr.index",&gindex) ;
 Double_t etotnorm;
   tsimc->SetBranchAddress("P.cal.etotnorm",&etotnorm) ;
 Double_t starttime;
   tsimc->SetBranchAddress("P.hod.starttime",&starttime) ;
 Double_t npesum;
    tsimc->SetBranchAddress("P.hgcer.npeSum",&npesum);
 TH1F* hetotnorm;
 TH2F* hetotnorm_npesum;
 TString temp;
 TH1F* hstarttime;
   temp=Form("Run %d ; Hod Starttime ; Counts",nrun);
   hstarttime = new TH1F("hstarttime",temp,140,0,70);
   temp=Form("Run %d ; Etotnorm  ; Counts",nrun);
   hetotnorm = new TH1F("hetotnorm",temp,100,0,2.0);
   HList.Add(hetotnorm);
   if (nrun > 1721) {
   temp=Form("Run %d ; Etotnorm  ; HG npeSum",nrun);
   hetotnorm_npesum = new TH2F("hetotnorm_npesum",temp,150,0,1.5,250,0,50);
   } else {
   temp=Form("Run %d ; Etotnorm  ; HG npeSum",nrun);
   hetotnorm_npesum = new TH2F("hetotnorm_npesum",temp,150,0,1.5,60,0,15);
   }
   HList.Add(hetotnorm_npesum);
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		if (ntracks >1) {
		  cout << " ntracks = " << ntracks << " " << gindex << endl;
		for (Int_t ii=0; ii<ntracks ;ii++) {
		  cout << "track = " << ii << " " << flag[ii] << " " << trackchi[ii] << endl;
		}
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
