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
#include "TF1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TCutG.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_hms_sieve(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
   tsimc->SetBranchAddress("H.cer.npeSum",&sumnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("H.gtr.x",&xtar);
 Double_t h_extcor_delta_xptar;
   tsimc->SetBranchAddress("H.extcor.delta_xptar",&h_extcor_delta_xptar);
 Double_t  reactx;
   tsimc->SetBranchAddress("H.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("H.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("H.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
 Double_t  yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("H.gtr.th",&xptar);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp);
   // Define histograms
   TH2F *hxfp_yfp_foil[9];
   TH1F *hys_foil[9];
   //   Double_t xptar_cut_lo[9]={-0.07,-0.046,-0.031,-0.014,-0.0,0.015,0.031,0.046,0.061};
   //Double_t xptar_cut_hi[9]={-0.046,-0.031,-0.014,-0.0,0.015,0.031,0.046,0.061,0.080};
   // old matrix cuts
   /*if (nrun==1136 ) {
   Double_t xptar_cut_lo[9]={-12.,-8.8,-6.5,-4.25,-1.25,1.25,4.,6.5,9.};
   Double_t xptar_cut_hi[9]={-8.8,-6.5,-4.25,-1.25,1.25,4.,6.5,9.,12};
   Double_t ys5_cut_lo[7]={-4.5,-3.1,-1.6,0.0,1.5,3.0,4.6};
   Double_t ys5_cut_hi[7]={-4.5,-3.1,-1.6,0.0,1.5,3.0,4.6};
   */
   //if (nrun==1337) {
   Double_t xptar_cut_lo[9]={-13.,-9.8,-7.5,-5.0,-2.5,0,2.5,5.,7.5};
   Double_t xptar_cut_hi[9]={-9.8,-7.5,-5.0,-2.5,0,2.5,5.,7.5,10.5};
   Double_t ys5_cut_lo[7]= {-5.0,-3.6,-2.1,-0.8,0.8,2.1,3.8};
   Double_t ys5_cut_hi[7]={-3.6,-2.1,-0.8,0.8,2.1,3.8,4.5};
   //
   /*if (nrun==1528) {
   Double_t xptar_cut_lo[9]={-13.,-9.8,-7.5,-5.0,-2.5,0,2.5,5.,7.5};
   Double_t xptar_cut_hi[9]={-9.8,-7.5,-5.0,-2.5,0,2.5,5.,7.5,10.5};
   Double_t ys5_cut_lo[7]= {-5.0,-3.6,-2.0,-0.8,0.8,2.,3.5};
   Double_t ys5_cut_hi[7]={-3.6,-2.0,-0.8,0.8,2.,3.5,4.5};
   */
   TH2F *hxfp_yfp_foil5_ys[7];
   for (Int_t i=0;i<7;i++) {
    hxfp_yfp_foil5_ys[i] = new TH2F(Form("hxfp_yfp_foil5_ys_%d",i+1), Form("Run %d Ys cut %5.3f to %5.3f ; X_fp ; Y_fp",nrun,ys5_cut_lo[i]-.25,ys5_cut_hi[i]+.25), 40,-10,10,40,-10,10);
   }
   for (Int_t i=0;i<9;i++) {
    hxfp_yfp_foil[i] = new TH2F(Form("hxfp_yfp_%d",i+1), Form("Run %d Xptar cut %5.3f to %5.3f ; X_fp ; Y_fp",nrun,xptar_cut_lo[i],xptar_cut_hi[i]), 200,-40,40,200,-25,25);
    hys_foil[i] = new TH1F(Form("hys_%d",i+1), Form("Run %d Cut %d; Y_sieve ; Counts",nrun,i+1), 200, -10,10);
   }
   TH2F *hxfp_yfp = new TH2F("hxfp_yfp", Form("Run %d ; X_fp ; Y_fp",nrun), 80,-40,40,60,-30,30);
   TH2F *hytar_yptar = new TH2F("hytar_yptar", Form("Run %d ; Y_tar ; Yp_tar",nrun), 80,-2.,2.,60,-.03,.03);
   TH2F *hxs_ys = new TH2F("hxs_ys", Form("Run %d ; Y_s ; X_s",nrun), 200,-6,6,200,-15,15);
   TH1F *hxptar = new TH1F("hxptar", Form("Run %d ; Xp_tar ; Counts",nrun), 200, -.1,.1);
   TH1F *hreactx = new TH1F("hreactx", Form("Run %d ; React X (cm) ; Counts",nrun), 100, -1,1);
   TH1F *hreacty = new TH1F("hreacty", Form("Run %d ;  React Y (cm) ; Counts",nrun), 100, -1,1);
   TH1F *hyptar = new TH1F("hyptar", Form("Run %d ; Yp_tar ; Counts",nrun), 200, -.04,.04);
   TH1F *hxtar = new TH1F("hxtar", Form("Run %d ; X_tar ; Counts",nrun), 200, -.5,.5);
    TH1F *hdeltaxptar = new TH1F("hdeltaxptar", Form("Run %d ; Delta(Xptar) ; Counts",nrun), 200, -.010,.010);
   Double_t xs;
   Double_t ys;
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (sumnpe>0.5&&etracknorm > .6 && ytar > -1 && ytar < 1) {
		  hxfp_yfp->Fill(xfp,yfp);
		  hreactx->Fill(reactx);
		  hreacty->Fill(reacty);
		  hxptar->Fill(xptar);
		  hyptar->Fill(yptar);
		  hytar_yptar->Fill(ytar,yptar);
		  hxtar->Fill(xtar);
		  hdeltaxptar->Fill(h_extcor_delta_xptar);
		  xs=xptar*168.;
		  ys=yptar*168.;
                  hxs_ys->Fill(ys,xs);
                  for (Int_t i=0;i<9;i++) {
                    if (xs>xptar_cut_lo[i]&&xs<=xptar_cut_hi[i]) {
  		      hys_foil[i]->Fill(ys);
		      hxfp_yfp_foil[i]->Fill(xfp,yfp);
                        for (Int_t j=0;j<7;j++) {
			  if (ys>(ys5_cut_lo[j]-.25)&&ys<=(ys5_cut_hi[j]+.25)) hxfp_yfp_foil5_ys[j]->Fill(xfp,yfp);
                        }
                    }
		  }
		} // electrons
	} // entries
	//
TCanvas *c1 = new TCanvas("c1", "Divide FP", 900,800);
 TF1 *linfit[7];
c1->Divide(2,4);
for (Int_t j=0;j<7;j++) {
c1->cd(j+1);
 gPad->SetGridy();
 gPad->SetGridx();
 linfit[j] = new TF1(Form("linfit_%d",j),"pol1",-10,10);
 hxfp_yfp_foil5_ys[j]->Draw("colz");
 hxfp_yfp_foil5_ys[j]->Fit(Form("linfit_%d",j),"Q");
  cout << linfit[j]->GetParameter(0) << " " << linfit[j]->GetParError(0) << " " <<linfit[j]->GetParameter(1) << " " <<linfit[j]->GetParError(1) <<endl;
 }
outputpdf="plots/"+basename+"_hms_sieve.pdf(";
c1->Print(outputpdf);
//
TCanvas *c2 = new TCanvas("c2", "FP", 900,800);
 c2->Divide(1,1);
 TH2F *htemp = (TH2F*) hxfp_yfp_foil5_ys[0]->Clone();
 htemp->Draw("axis");
 htemp->SetTitle(Form("Run = %d, Plot all fits",nrun));
 gPad->SetGridy(1);
 gPad->SetGridx(1);
for (Int_t j=0;j<7;j++) {
  linfit[j]->Draw("same");
  linfit[j]->SetLineWidth(3);
 }
outputpdf="plots/"+basename+"_hms_sieve.pdf";
c2->Print(outputpdf);
//
TCanvas *c3 = new TCanvas("c3", "Ysieve", 900,800);
 c3->Divide(3,3);
for (Int_t j=0;j<9;j++) {
c3->cd(j+1);
 hys_foil[j]->Draw();
  }
outputpdf="plots/"+basename+"_hms_sieve.pdf";
c3->Print(outputpdf);
//
//
TCanvas *c3a = new TCanvas("c3a", "Xfp v Yfp", 900,800);
 c3a->Divide(1,1);
 gPad->SetGridy(1);
 gPad->SetGridx(1);
 hxfp_yfp->Draw("colz");
outputpdf="plots/"+basename+"_hms_sieve.pdf";
c3a->Print(outputpdf);
//
TCanvas *c4 = new TCanvas("c4", "Ysieve", 900,800);
 c4->Divide(1,1);
 c4->cd(1);
 gPad->SetGridy(1);
 gPad->SetGridx(1);
 hxs_ys->Draw("colz");
outputpdf="plots/"+basename+"_hms_sieve.pdf)";
c4->Print(outputpdf);
//
TCanvas *c3b = new TCanvas("c3b", "react", 900,800);
 c3b->Divide(2,3);
 c3b->cd(1);
 hxptar->Draw();
 c3b->cd(2);
 hreactx->Draw();
 c3b->cd(3);
 hreacty->Draw();
 c3b->cd(4);
 hxtar->Draw();
 c3b->cd(5);
 hdeltaxptar->Draw();
outputpdf="plots/"+basename+"_hms_react.pdf";
c3b->Print(outputpdf);
//
	//
}
