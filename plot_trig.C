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
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_trig(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="../../ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  edtm_time;
   tsimc->SetBranchAddress("T.shms.pEDTM_tdcTime",&edtm_time);
 Double_t  p1x_time;
   tsimc->SetBranchAddress("T.shms.p1X_tdcTime",&p1x_time);
 Double_t  p1y_time;
   tsimc->SetBranchAddress("T.shms.p1Y_tdcTime",&p1y_time);
 Double_t  p2x_time;
   tsimc->SetBranchAddress("T.shms.p2X_tdcTime",&p2x_time);
 Double_t  p2y_time;
   tsimc->SetBranchAddress("T.shms.p2Y_tdcTime",&p2y_time);
   //
   TH1F *hp1xtime = new TH1F("hp1xtime",Form("Run %d ; SHMS 1X time (ns); Counts",nrun),600,-100,200);
  TH1F *hp1ytime = new TH1F("hp1ytime",Form("Run %d ; SHMS 1Y time (ns); Counts",nrun),600,-100,200);
   TH1F *hp2xtime = new TH1F("hp2xtime",Form("Run %d ; SHMS 1X time (ns); Counts",nrun),600,-100,200);
  TH1F *hp2ytime = new TH1F("hp2ytime",Form("Run %d ; SHMS 1Y time (ns); Counts",nrun),600,-100,200);
  TH1F *hp1xtime_allgood = new TH1F("hp1xtime_allgood",Form("All Scin Good Run %d ; SHMS 1X time (ns); Counts",nrun),600,-100,200);
  TH1F *hp1ytime_allgood = new TH1F("hp1ytime_allgood",Form("All Scin Good Run %d ; SHMS 1Y time (ns); Counts",nrun),600,-100,200);
   TH1F *hp2xtime_allgood = new TH1F("hp2xtime_allgood",Form("All Scin Good Run %d ; SHMS 2X time (ns); Counts",nrun),600,-100,200);
  TH1F *hp2ytime_allgood = new TH1F("hp2ytime_allgood",Form("All Scin Good Run %d ; SHMS 2Y time (ns); Counts",nrun),600,-100,200);
  TH1F *hp1ytime_notgoodpx1 = new TH1F("hp1ytime_notgoodpx1",Form("no X1, Other 3 scin good Run %d ; SHMS 1Y time (ns); Counts",nrun),600,-100,200);
   TH1F *hp2xtime_notgoodpx1 = new TH1F("hp2xtime_notgoodpx1",Form("no X1, Other 3 scin good Run %d ; SHMS 2X time (ns); Counts",nrun),600,-100,200);
  TH1F *hp2ytime_notgoodpx1 = new TH1F("hp2ytime_notgoodpx1",Form("no X1, Other 3 scin good Run %d ; SHMS 2Y time (ns); Counts",nrun),600,-100,200);
  TH1F *hp1xtime_notgoodpy1 = new TH1F("hp1xtime_notgoodpy1",Form("no Y1, Other 3 scin good Run %d ; SHMS 1X time (ns); Counts",nrun),600,-100,200);
   TH1F *hp2xtime_notgoodpy1 = new TH1F("hp2xtime_notgoodpy1",Form("no Y1, Other 3 scin good Run %d ; SHMS 2X time (ns); Counts",nrun),600,-100,200);
  TH1F *hp2ytime_notgoodpy1 = new TH1F("hp2ytime_notgoodpy1",Form("no Y1, Other 3 scin good Run %d ; SHMS 2Y time (ns); Counts",nrun),600,-100,200);
  TH1F *hp1xtime_notgoodpx2 = new TH1F("hp1xtime_notgoodpx2",Form("no X2, Other 3 scin good Run %d ; SHMS 1X time (ns); Counts",nrun),600,-100,200);
   TH1F *hp1ytime_notgoodpx2 = new TH1F("hp1ytime_notgoodpx2",Form("no X2, Other 3 scin good Run %d ; SHMS 1Y time (ns); Counts",nrun),600,-100,200);
  TH1F *hp2ytime_notgoodpx2 = new TH1F("hp2ytime_notgoodpx2",Form("no X2, Other 3 scin good Run %d ; SHMS 2Y time (ns); Counts",nrun),600,-100,200);
  TH1F *hp1xtime_notgoodpy2 = new TH1F("hp1xtime_notgoodpy2",Form("no Y2, Other 3 scin good Run %d ; SHMS 1X time (ns); Counts",nrun),600,-100,200);
   TH1F *hp1ytime_notgoodpy2 = new TH1F("hp1ytime_notgoodpy2",Form("no Y2, Other 3 scin good Run %d ; SHMS 1Y time (ns); Counts",nrun),600,-100,200);
  TH1F *hp2xtime_notgoodpy2 = new TH1F("hp2xtime_notgoodpy2",Form("no Y2, Other 3 scin good Run %d ; SHMS 2Y time (ns); Counts",nrun),600,-100,200);
  //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (edtm_time==0) {
                  Bool_t good_p1x = p1x_time!=0;
                  Bool_t good_p1y = p1y_time!=0;
                  Bool_t good_p2x = p2x_time!=0;
                  Bool_t good_p2y = p2y_time!=0;
                  hp1xtime->Fill(p1x_time);
                  hp1ytime->Fill(p1y_time);
                  hp2xtime->Fill(p2x_time);
                  hp2ytime->Fill(p2y_time);
		  if (good_p1x&&good_p1y&&good_p2x&&good_p2y) {
                  hp1xtime_allgood->Fill(p1x_time);
                  hp1ytime_allgood->Fill(p1y_time);
                  hp2xtime_allgood->Fill(p2x_time);
                  hp2ytime_allgood->Fill(p2y_time);
		  }
		  if (!good_p1x&&good_p1y&&good_p2x&&good_p2y)  {
                   hp1ytime_notgoodpx1->Fill(p1y_time);
                   hp2xtime_notgoodpx1->Fill(p2x_time);
                   hp2ytime_notgoodpx1->Fill(p2y_time);
		  }
		  if (good_p1x&&!good_p1y&&good_p2x&&good_p2y)  {
                  hp1xtime_notgoodpy1->Fill(p1x_time);
                   hp2xtime_notgoodpy1->Fill(p2x_time);
                  hp2ytime_notgoodpy1->Fill(p2y_time);
		  }
		  if (good_p1x&&good_p1y&&!good_p2x&&good_p2y)  {
                  hp1xtime_notgoodpx2->Fill(p1x_time);
                   hp1ytime_notgoodpx2->Fill(p1y_time);
                  hp2ytime_notgoodpx2->Fill(p2y_time);
		  }
		  if (good_p1x&&good_p1y&&good_p2x&&!good_p2y)  {
                  hp1xtime_notgoodpy2->Fill(p1x_time);
                   hp1ytime_notgoodpy2->Fill(p1y_time);
                  hp2xtime_notgoodpy2->Fill(p2x_time);
		  }
		}
	}
	//
//
TCanvas *c1 = new TCanvas("c1", " S1X ", 900,800);
c1->Divide(2,2);
 c1->cd(1);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp1xtime_allgood->Draw();
 c1->cd(2);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp1xtime_notgoodpy1->Draw();
 //hp1xtime_notgoodpy1->SetLineColor(2);
 c1->cd(3);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp1xtime_notgoodpx2->Draw();
 //hp1xtime_notgoodpx2->SetLineColor(3);
 c1->cd(4);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp1xtime_notgoodpy2->Draw();
 //hp1xtime_notgoodpy2->SetLineColor(4);
    outputpdf="plots/"+basename+"_trig.pdf(";
c1->Print(outputpdf);

//
TCanvas *c2 = new TCanvas("c2", " S1Y ", 900,800);
c2->Divide(2,2);
 c2->cd(1);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp1ytime_allgood->Draw();
 c2->cd(2);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp1ytime_notgoodpx1->Draw();
 //hp1ytime_notgoodpy1->SetLineColor(2);
 c2->cd(3);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp1ytime_notgoodpx2->Draw();
 //hp1ytime_notgoodpx2->SetLineColor(3);
 c2->cd(4);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp1ytime_notgoodpy2->Draw();
 //hp1ytime_notgoodpy2->SetLineColor(4);
    outputpdf="plots/"+basename+"_trig.pdf";
c2->Print(outputpdf);
//
TCanvas *c3 = new TCanvas("c3", " S2X ", 900,800);
c3->Divide(2,2);
 c3->cd(1);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp2xtime_allgood->Draw();
 c3->cd(2);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp2xtime_notgoodpx1->Draw();
 //hp2xtime_notgoodpy1->SetLineColor(2);
 c3->cd(3);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp2xtime_notgoodpy1->Draw();
 //hp2xtime_notgoodpx2->SetLineColor(3);
 c3->cd(4);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp2xtime_notgoodpy2->Draw();
 //hp2xtime_notgoodpy2->SetLineColor(4);
    outputpdf="plots/"+basename+"_trig.pdf";
c3->Print(outputpdf);
//
TCanvas *c4 = new TCanvas("c4", " S2Y ", 900,800);
c4->Divide(2,2);
 c4->cd(1);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp2ytime_allgood->Draw();
 c4->cd(2);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp2ytime_notgoodpx1->Draw();
 //hp2ytime_notgoodpy1->SetLineColor(2);
 c4->cd(3);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp2ytime_notgoodpy1->Draw();
 //hp2ytime_notgoodpx2->SetLineColor(3);
 c4->cd(4);
 gPad->SetGridy();
 gPad->SetGridx();
 gPad->SetLogy();
 hp2ytime_notgoodpx2->Draw();
 //hp2ytime_notgoodpy2->SetLineColor(4);
    outputpdf="plots/"+basename+"_trig.pdf)";
c4->Print(outputpdf);

//
}
