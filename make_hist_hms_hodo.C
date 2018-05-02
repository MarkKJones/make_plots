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

void make_hist_hms_hodo(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_hodo_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 static const Int_t plnum=4;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 static const Int_t npad[plnum]={16,10,16,10};
 Double_t negtimediff[plnum][16];
 Double_t postimediff[plnum][16];
 Double_t negtime[plnum][16];
 Double_t postime[plnum][16];
 Double_t negadcmult[plnum][16];
 Double_t posadcmult[plnum][16];
 for (Int_t ip=0;ip<plnum;ip++) {
   tsimc->SetBranchAddress(Form("H.hod.%s.GoodNegAdcTdcDiffTime",plname[ip]),&negtimediff[ip]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.GoodPosAdcTdcDiffTime",plname[ip]),&postimediff[ip]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.GoodNegTdcTimeUnCorr",plname[ip]),&negtime[ip]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.GoodPosTdcTimeUnCorr",plname[ip]),&postime[ip]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.GoodNegAdcMult",plname[ip]),&negadcmult[ip]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.GoodPosAdcMult",plname[ip]),&posadcmult[ip]) ;
 }
   // Define histograms
 /*
 TH1F* neghist[plnum][21];
 TH1F* poshist[plnum][21];
 for (Int_t ip=0;ip<plnum;ip++) {
 for (Int_t ipd=0;ipd<npad[ip];ipd++) {
   neghist[ip][ipd] = new TH1F(Form("hist_%s_neg_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Neg ADC Time-Tdc Time ns; COunts",plname[ip],ipd+1),200,-30,20);
   poshist[ip][ipd] = new TH1F(Form("hist_%s_pos_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Pos ADC Time-Tdc Time ns; COunts",plname[ip],ipd+1),200,-30,20);
 } }
 */
 TH2D *neg2dhist[plnum];
 TH2D *pos2dhist[plnum];
 TH2D *neg2dtimehist[plnum];
 TH2D *pos2dtimehist[plnum];
 for (Int_t ip=0;ip<plnum;ip++) {
   neg2dhist[ip] = new  TH2D(Form("hist_%s_neg_pad",plname[ip]),Form("Plane %s ; Neg Pad ; ADC Time-Tdc Time ns",plname[ip]),npad[ip],0,npad[ip]+1,320,-80,-40);
   HList.Add(neg2dhist[ip]);
   pos2dhist[ip] = new  TH2D(Form("hist_%s_pos_pad",plname[ip]),Form("Plane %s ; Pos Pad ; ADC Time-Tdc Time ns",plname[ip]),npad[ip],0,npad[ip]+1,320,-80,-40);
   HList.Add(pos2dhist[ip]);
   neg2dtimehist[ip] = new  TH2D(Form("hist_%s_neg_time_pad",plname[ip]),Form("Plane %s ; Neg Pad ;Tdc Time ns",plname[ip]),npad[ip],0,npad[ip]+1,300,-100,100);
   HList.Add(neg2dtimehist[ip]);
   pos2dtimehist[ip] = new  TH2D(Form("hist_%s_pos_time_pad",plname[ip]),Form("Plane %s ; Pos Pad ; Tdc Time ns",plname[ip]),npad[ip],0,npad[ip]+1,300,-100,100);
   HList.Add(pos2dtimehist[ip]);
 }
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
 for (Int_t ip=0;ip<plnum;ip++) {
 for (Int_t ipd=0;ipd<npad[ip];ipd++) {
   neg2dhist[ip]->Fill(float(ipd+1),negtimediff[ip][ipd]);
   pos2dhist[ip]->Fill(float(ipd+1),postimediff[ip][ipd]);
   neg2dtimehist[ip]->Fill(float(ipd+1),negtime[ip][ipd]);
   pos2dtimehist[ip]->Fill(float(ipd+1),postime[ip][ipd]);
 } 
}		
  	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
