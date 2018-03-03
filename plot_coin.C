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

void plot_coin(TString basename="",Int_t nrun=2043){
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
 Double_t h_reactx;
   tsimc->SetBranchAddress("H.react.x",&h_reactx);
 Double_t h_reacty;
   tsimc->SetBranchAddress("H.react.y",&h_reacty);
 Double_t h_reactz;
   tsimc->SetBranchAddress("H.react.z",&h_reactz);
 Double_t p_reactx;
   tsimc->SetBranchAddress("P.react.x",&p_reactx);
 Double_t p_reacty;
   tsimc->SetBranchAddress("P.react.y",&p_reacty);
 Double_t p_reactz;
   tsimc->SetBranchAddress("P.react.z",&p_reactz);
 Double_t h_extcor_delta_dp;
   tsimc->SetBranchAddress("H.extcor.delta_dp",&h_extcor_delta_dp);
 Double_t h_extcor_delta_xptar;
   tsimc->SetBranchAddress("H.extcor.delta_xptar",&h_extcor_delta_xptar);
 Double_t h_extcor_x;
   tsimc->SetBranchAddress("H.extcor.x",&h_extcor_x);
 Double_t h_extcor_y;
   tsimc->SetBranchAddress("H.extcor.y",&h_extcor_y);
 Double_t p_extcor_x;
   tsimc->SetBranchAddress("P.extcor.x",&p_extcor_x);
 Double_t p_extcor_y;
   tsimc->SetBranchAddress("P.extcor.y",&p_extcor_y);
 Double_t raster_ay;
   tsimc->SetBranchAddress("H.rb.raster.fr_ya",&raster_ay);
 Double_t raster_ax;
   tsimc->SetBranchAddress("H.rb.raster.fr_xa",&raster_ax);
 Double_t raster_ay_adc;
   tsimc->SetBranchAddress("H.rb.raster.frya_adc",&raster_ay_adc);
 Double_t raster_ax_adc;
   tsimc->SetBranchAddress("H.rb.raster.frxa_adc",&raster_ax_adc);
 Double_t raster_raw_ay;
   tsimc->SetBranchAddress("H.rb.raster.fryaRawAdc",&raster_raw_ay);
 Double_t raster_raw_ax;
   tsimc->SetBranchAddress("H.rb.raster.frxaRawAdc",&raster_raw_ax);
 Double_t raster_raw_by;
   tsimc->SetBranchAddress("H.rb.raster.frybRawAdc",&raster_raw_by);
 Double_t raster_raw_bx;
   tsimc->SetBranchAddress("H.rb.raster.frxbRawAdc",&raster_raw_bx);
 Double_t  sumnpe;
   tsimc->SetBranchAddress("H.cer.npeSum",&sumnpe);
 Double_t  hdelta;
   tsimc->SetBranchAddress("H.gtr.dp",&hdelta);
 Double_t  pdelta;
   tsimc->SetBranchAddress("P.gtr.dp",&pdelta);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t  hxfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&hxfp);
 Double_t  hyfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&hyfp);
 Double_t  hxpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&hxpfp);
 Double_t  hypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&hypfp);
 Double_t  pxfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&pxfp);
 Double_t  pyfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&pyfp);
 Double_t  pxpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&pxpfp);
 Double_t  pypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&pypfp);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  hmom;
   tsimc->SetBranchAddress("H.gtr.p",&hmom);
 Double_t  yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("H.gtr.th",&xptar);
 Double_t  W;
    tsimc->SetBranchAddress("H.kin.primary.W",&W);
 Double_t  Q2;
    tsimc->SetBranchAddress("H.kin.primary.Q2",&Q2);
 Double_t  scat_ang;
    tsimc->SetBranchAddress("H.kin.primary.scat_ang_deg",&scat_ang);
 Double_t  emiss;
    tsimc->SetBranchAddress("P.kin.secondary.emiss",&emiss);
 Double_t  erecoil;
    tsimc->SetBranchAddress("P.kin.secondary.Erecoil",&erecoil);
 Double_t  pmiss;
    tsimc->SetBranchAddress("P.kin.secondary.pmiss",&pmiss);
 Double_t  mrecoil;
    tsimc->SetBranchAddress("P.kin.secondary.Mrecoil",&mrecoil);
 // tsimc->SetBranchAddress("H.kin.W",&W);
   // Define histograms
   TH1F *hemiss = new TH1F("hemiss", Form("Run %d ; Emiss (GeV) ; Counts",nrun), 400, -.5,.5);
   TH2F *hemiss_xfp = new TH2F("hemissxfp", Form("Run %d ; Emiss (GeV) ; X_fp",nrun), 400, -.5,.5,100,-40,40);
   TH1F *hpmiss = new TH1F("hpmiss", Form("Run %d ; Pmiss (GeV) ; Counts",nrun), 400, -.5,.5);
   TH1F *hmrecoil = new TH1F("hmrecoil", Form("Run %d ; Mrecoil (GeV) ; Counts",nrun), 400, -.5,.5);
   TH1F *hreactx = new TH1F("hreactx", Form("Run %d ;HMS React X (cm) ; Counts",nrun), 400, -.5,.5);
   TH1F *hreacty = new TH1F("hreacty", Form("Run %d ;HMS React Y (cm) ; Counts",nrun), 400, -.5,.5);
   TH1F *hreactz = new TH1F("hreactz", Form("Run %d ;HMS React Z (cm) ; Counts",nrun), 400, -6.,6.);
  TH1F *preactx = new TH1F("preactx", Form("Run %d ;SHMS React X (cm) ; Counts",nrun), 400, -.5,.5);
   TH1F *preacty = new TH1F("preacty", Form("Run %d ;SHMS React Y (cm) ; Counts",nrun), 400, -.5,.5);
   TH1F *preactz = new TH1F("preactz", Form("Run %d ;SHMS React Z (cm) ; Counts",nrun), 400, -6.,6.);
   TH2F *preactz_hreactz = new TH2F("preactz_hreactz", Form("Run %d ;SHMS React Z (cm) ;HMS React Z (cm) ",nrun), 100, -6.,6., 100, -6.,6.);
   TH1F *hextcorx = new TH1F("hextcorx", Form("Run %d ;HMS Extcor X (cm) ; Counts",nrun), 400, -2.,2.);
   TH2F *hextcorx_hreactz = new TH2F("hextcorx_reactz", Form("Run %d ;HMS Extcor X (cm) ; HMS React Z (cm)",nrun), 100, -1.,1.,100, -6.,6.);
   TH2F *hextcorx_deltadp = new TH2F("hextcorx_deltadp", Form("Run %d ;HMS Extcor X (cm) ; Delta(dp)",nrun), 100, -2.,2.,100, -1.,1.);
   TH2F *hextcorx_deltaxptar = new TH2F("hextcorx_deltaxptar", Form("Run %d ;HMS Extcor X (cm) ; Delta(xptar)",nrun), 100, -2.,2.,100, -.01,.01);
   TH1F *hextcory = new TH1F("hextcory", Form("Run %d ;HMS Extcor Y (cm) ; Counts",nrun), 400, -5.,5.);
   TH1F *pextcorx = new TH1F("pextcorx", Form("Run %d ;SHMS Extcor X (cm) ; Counts",nrun), 400, -2.,2.);
   TH1F *pextcory = new TH1F("pextcory", Form("Run %d ;SHMS Extcor Y (cm) ; Counts",nrun), 400, -5.,5.);
   TH1F *hrastx = new TH1F("hrastx", Form("Run %d ;Raster X (cm) ; Counts",nrun), 600, -.3,.3);
   TH1F *hrasty = new TH1F("hrasty",Form("Run %d ;Raster Y (cm) ;Counts",nrun), 600, -.3,.3);
   TH1F *hrastxcalc = new TH1F("hrastxcalc", Form("Run %d ;calc Raster X (cm) ; Counts",nrun), 600, -.3,.3);
   TH1F *hrastycalc = new TH1F("hrastycalc",Form("Run %d ;calc Raster Y (cm) ;Counts",nrun), 600, -.3,.3);
   TH1F *hrastrawxa = new TH1F("hrastrawxa", Form("Run %d ;Raster Xa Raw ADC;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hrastrawya = new TH1F("hrastrawya",Form("Run %d ;Raster Ya Raw ADC ;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hrastxa_adc = new TH1F("hrastxa_adc", Form("Run %d ;Raster Xa ADC;Counts",nrun), 2000, -50000.,50000.);
   TH1F *hrastya_adc = new TH1F("hrastya_adc",Form("Run %d ;Raster Ya ADC ;Counts",nrun), 2000, -50000.,50000.);
   TH1F *hrastrawxb = new TH1F("hrastrawxb", Form("Run %d ;Raster Xb ADC;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hrastrawyb = new TH1F("hrastrawyb",Form("Run %d ;Raster Yb ADC ;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hytar = new TH1F("hytar",Form("Run %d ; Y_tar (cm);Counts",nrun), 200, -5.,5.);
   TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 100, 0.85,1.050);
   TH1F *hWcalc = new TH1F("hWcalc",Form("Run %d ; W Calc (GeV);Counts",nrun), 100, 0.85,1.050);
   TH2F *hWdiff = new TH2F("hWdiff",Form("Run %d ; W dfff (GeV);Counts",nrun), 100, -.1,.1, 100, 0.85,1.050);
   TH2F *hthdiff = new TH2F("hthdiff",Form("Run %d ; TH dfff (GeV);yptar",nrun), 100, -1,1.,100, 20.,25.);
   TH2F *hytar_yrast = new TH2F("hytar_yrast",Form("Run %d ; Y_tar (cm); Raster Y ADC (mV)",nrun), 200, -5.,5.,2000, -.3,.3);
   TH2F *hytar_xrast = new TH2F("hytar_xrast",Form("Run %d ; Y_tar (cm); Raster X ADC (mV)",nrun), 200, -5.,5.,2000, -.3,.3);
   TH2F *hW_yrast = new TH2F("hW_yrast",Form("Run %d ; W (GeV) ; Raster Y ADC (mV)",nrun),100, 0.85,1.050 ,200, -.3,.3);
   TH2F *hW_xrast = new TH2F("hW_xrast",Form("Run %d ; W (GeV)  ; Raster X ADC (mV)",nrun),100, 0.85,1.050,200, -.3,.3);
   TH2F *hW_yptar = new TH2F("hW_yptar",Form("Run %d ; W (GeV)  ; Yptar",nrun),50, 0.85,1.050,30, -.03,.03);
   TH2F *hW_xptar = new TH2F("hW_xptar",Form("Run %d ; W (GeV)  ; Xptar",nrun),50, 0.85,1.050,30, -.06,.06);
   TH2F *hW_xpfp = new TH2F("hW_xpfp",Form("Run %d ; W (GeV)  ; Xpfp",nrun),50, 0.85,1.050,40, -.03,.03);
   TH2F *hW_ypfp = new TH2F("hW_ypfp",Form("Run %d ; W (GeV)  ; Ypfp",nrun),50, 0.85,1.050,40, -.03,.03);
    TH2F *hW_xfp = new TH2F("hW_xfp",Form("Run %d ; W (GeV)  ; Xfp",nrun),50, 0.85,1.050,40, -30.,30.);
   TH2F *hW_yfp = new TH2F("hW_yfp",Form("Run %d ; W (GeV)  ; Yfp",nrun),50, 0.85,1.050,40, -30.,30.);
   TH2F *hW_xpcor = new TH2F("hW_xpcor",Form("Run %d ; W (GeV)  ; Xp_cor",nrun),50, 0.85,1.050,40, -.01,.01);
// loop over entries
   Double_t theta,Q2calc,nu,Wcalc,W2calc,Ef;
   Double_t cos_ts=TMath::Cos(22.01/180*3.14159);
   Double_t sin_ts=TMath::Sin(22.01/180*3.14159);
   Double_t Mp = .93827;
   Double_t Ei=6.405;
   Double_t p_spec=4.284;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (sumnpe > 2. &&  hdelta>-8 &&  hdelta<8 ) {
		  Ef = p_spec * (1.)* (1.0 + 0.01*hdelta);
		  theta = TMath::ACos((cos_ts + yptar * sin_ts) / TMath::Sqrt( 1. + (xptar+.005) * (xptar+.005) + yptar * yptar )); // polar 			scattering angle relative to the beam line //rad
		Q2calc = 4.0 * Ei * Ef * (TMath::Sin(theta / 2.0) * TMath::Sin(theta / 2.0)); //GeV^2
		nu = Ei - Ef; //GeV
		W2calc = -Q2calc + Mp * Mp + 2.0 * Mp * nu; // GeV^2
                Wcalc=0.;
		if (W2calc > 0) Wcalc = TMath::Sqrt(-Q2calc + Mp * Mp + 2.0 * Mp * nu); //Ge		  
		hemiss->Fill(emiss);
		hemiss_xfp->Fill(emiss,pxfp);
		hpmiss->Fill(pmiss);
		hWcalc->Fill(Wcalc);
 		hWdiff->Fill(Wcalc-W,W);
 		hthdiff->Fill(theta*180/3.14159-scat_ang,scat_ang);
                 hextcorx_deltadp->Fill(h_extcor_x,h_extcor_delta_dp);
                  hextcorx_deltaxptar->Fill(h_extcor_x,h_extcor_delta_xptar);
		  hextcorx->Fill(h_extcor_x);
		  hextcorx_hreactz->Fill(h_extcor_x,h_reactz);
		  hextcory->Fill(h_extcor_y);
		  pextcorx->Fill(p_extcor_x);
		  pextcory->Fill(p_extcor_y);
		  hreactx->Fill(h_reactx);
		  hreacty->Fill(h_reacty);
		  hreactz->Fill(h_reactz);
		  preactx->Fill(p_reactx);
		  preacty->Fill(p_reacty);
		  preactz->Fill(p_reactz);
		  preactz_hreactz->Fill(p_reactz,h_reactz);
		  hrasty->Fill(raster_ay);
		  hrastx->Fill(raster_ax);
		  hrastycalc->Fill((raster_raw_ay-68016)/101130.*(6.4/Ei));
		  hrastxcalc->Fill((raster_raw_ax-64251)/119215.*(6.4/Ei));
		  hrastrawya->Fill(raster_raw_ay);
		  hrastrawxa->Fill(raster_raw_ax);
		  hrastya_adc->Fill(raster_ay_adc);
		  hrastxa_adc->Fill(raster_ax_adc);
		  hrastrawyb->Fill(raster_raw_by);
		  hrastrawxb->Fill(raster_raw_bx);
		  hytar->Fill(ytar);
		  hytar_yrast->Fill(ytar,raster_ay);
		  hytar_xrast->Fill(ytar,raster_ax);
		  hW->Fill(W);
		  hW_yrast->Fill(W,raster_ay);
		  hW_xrast->Fill(W,raster_ax);
		  hW_yptar->Fill(W,yptar);
		  hW_xptar->Fill(W,xptar);
		  hW_xpfp->Fill(W,hxpfp);
		  hW_ypfp->Fill(W,hypfp);
		  hW_xfp->Fill(W,hxfp);
		  hW_yfp->Fill(W,hyfp);
		  hW_xpcor->Fill(W,h_extcor_delta_xptar);
		}
	}
	// plot data
	//
TCanvas *c = new TCanvas("c", "pos raster", 900,700);
c->Divide(1,2);
c->cd(1);
hrastx->Draw();
hrastxcalc->Draw("same");
hrastxcalc->SetLineColor(2);
c->cd(2);
hrasty->Draw();
hrastycalc->Draw("same");
hrastycalc->SetLineColor(2);
    outputpdf="plots/"+basename+"_raster_pos.pdf";
c->Print(outputpdf);
//
TCanvas *c111 = new TCanvas("c111", "Emiss", 900,700);
c111->Divide(1,2);
c111->cd(1);
hemiss->Draw();
c111->cd(2);
hemiss_xfp->Draw("colz");
outputpdf="plots/"+basename+"_emiss.pdf";
c111->Print(outputpdf);
//
 if (1==-1) {
TCanvas *c11 = new TCanvas("c11", "A adc raster", 900,700);
c11->Divide(1,2);
c11->cd(1);
hrastxa_adc->Draw();
c11->cd(2);
hrastya_adc->Draw();
    outputpdf="plots/"+basename+"_rastera_adc.pdf";
c11->Print(outputpdf);
//
TCanvas *c3 = new TCanvas("c3", "raw raster", 900,700);
c3->Divide(2,2);
c3->cd(1);
hrastrawxa->Draw();
c3->cd(2);
hrastrawya->Draw();
c3->cd(3);
hrastrawxb->Draw();
c3->cd(4);
hrastrawyb->Draw();
    outputpdf="plots/"+basename+"_raster_raw.pdf";
c3->Print(outputpdf);
//
TCanvas *c1 = new TCanvas("c1", "Ytar", 900,700);
c1->Divide(1,1);
c1->cd(1);
hytar->Draw();
    outputpdf="plots/"+basename+"_ytar.pdf";
c1->Print(outputpdf);
	//
 }
TCanvas *cW = new TCanvas("cW", "W", 900,700);
cW->Divide(1,1);
cW->cd(1);
hW->Draw();
//hWcalc->Draw("same");
// hWcalc->SetLineColor(2);
    outputpdf="plots/"+basename+"_W.pdf";
cW->Print(outputpdf);
//
 if (1==-1) {
TCanvas *c2 = new TCanvas("c2", "Ytar v rast", 900,700);
c2->Divide(1,2);
c2->cd(1);
hytar_yrast->Draw("colz");
c2->cd(2);
hytar_xrast->Draw("colz");
    outputpdf="plots/"+basename+"_ytar_rast.pdf";
c2->Print(outputpdf);
 }
//
TCanvas *c4 = new TCanvas("c4", "W v rast", 900,700);
c4->Divide(1,2);
c4->cd(1);
hW_yrast->Draw("colz");
c4->cd(2);
hW_xrast->Draw("colz");
    outputpdf="plots/"+basename+"_W_rast.pdf";
c4->Print(outputpdf);

//
TCanvas *c5 = new TCanvas("c5", "W v angles", 900,700);
c5->Divide(1,2);
c5->cd(1);
hW_yptar->Draw("colz");
c5->cd(2);
hW_xptar->Draw("colz");
    outputpdf="plots/"+basename+"_W_angles.pdf";
c5->Print(outputpdf);


 //
//
 if (1==-1) {
TCanvas *cextcor = new TCanvas("cextcor", " Ext Cor ", 900,700);
cextcor->Divide(2,2);
cextcor->cd(1);
hextcorx->Draw();
cextcor->cd(2);
hextcory->Draw();
cextcor->cd(3);
pextcorx->Draw();
cextcor->cd(4);
pextcory->Draw();
    outputpdf="plots/"+basename+"_extcor.pdf";
cextcor->Print(outputpdf);
 }
//
// end scripts
}
