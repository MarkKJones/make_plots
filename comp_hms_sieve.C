#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
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

Int_t npeaks=30;
Double_t fpeaks(Double_t *x,Double_t *par) {
  Double_t result = par[0] +par[1]*x[0];
  for (Int_t p=0;p<npeaks;p++) {
    Double_t norm = par[3*p+2];
    Double_t mean = par[3*p+3];
    Double_t sigma = par[3*p+4];
    result +=norm*TMath::Gaus(x[0],mean,sigma);
  }
  return result;
}

void comp_hms_sieve(TString basename="",Int_t nrun=2043,Int_t nfoils=3){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1111);
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
 Double_t h_extcor_delta_xptar;
   tsimc->SetBranchAddress("H.extcor.delta_xptar",&h_extcor_delta_xptar);
 Double_t raster_raw_ay;
   tsimc->SetBranchAddress("H.rb.raster.fryaRawAdc",&raster_raw_ay);
 Double_t raster_raw_ax;
   tsimc->SetBranchAddress("H.rb.raster.frxaRawAdc",&raster_raw_ax);
 Double_t raster_raw_by;
   tsimc->SetBranchAddress("H.rb.raster.frybRawAdc",&raster_raw_by);
 Double_t raster_raw_bx;
   tsimc->SetBranchAddress("H.rb.raster.frxbRawAdc",&raster_raw_bx);
 Double_t raster_ay;
   tsimc->SetBranchAddress("H.rb.raster.fr_ya",&raster_ay);
 Double_t raster_ax;
   tsimc->SetBranchAddress("H.rb.raster.fr_xa",&raster_ax);
 Double_t  reactx;
   tsimc->SetBranchAddress("H.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("H.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("H.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
  Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("H.gtr.x",&xtar);
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
   //
	Double_t fcent[3]={-10.,0.,10.};
	Double_t fwid[3]={2.,2.,2.};
   TH1F *hytar=new TH1F("hytar", Form("Run %d ; Ytar ; Counts",nrun), 240,-6.,6.);
   TH1F *hreactz=new TH1F("hreactz", Form("Run %d ; Z_beam (cm) ; Counts",nrun), 600,-15.,15.);
   TH1F *hreactx=new TH1F("hreactx", Form("Run %d ; X_beam (cm) ; Counts",nrun), 600,-.5,.5);
   TH1F *hreacty=new TH1F("hreacty", Form("Run %d ; Y_beam (cm) ; Counts",nrun), 600,-.5,.5);
   TH1F *hrastx = new TH1F("hrastx", Form("Run %d ;Raster X (cm) ; Counts",nrun), 2000, -.3,.3);
   TH2F *hxs_ys = new TH2F("hxs_ys", Form("Run %d ;Y_sieve (cm) ; X_sieve",nrun), 200, -7.,7.,200, -12.,12.);
   TH2F *hxs_ys_foil[3];
   TH1F *hxs_foil[3];
   TH1F *hys_foil[3];
 for (Int_t nf=0;nf<nfoils;nf++) {
   hxs_ys_foil[nf] = new TH2F(Form("hxs_ys_foil_%d",nf), Form("Foil %5.2f Run %d ;Y_sieve (cm) ; X_sieve",fcent[nf],nrun), 200, -7.,7.,200, -12.,12.);
   hxs_foil[nf] = new TH1F(Form("hxs_foil_%d",nf), Form("Foil %5.2f Run %d ;X_sieve (cm) ; Counts",fcent[nf],nrun), 200, -12.,12.);
   hys_foil[nf] = new TH1F(Form("hys_foil_%d",nf), Form("Foil %5.2f Run %d ;Y_sieve (cm) ; Counts",fcent[nf],nrun), 200, -7.,7.);
 }
   TH2F *hrastx_ytar = new TH2F("hrastx_ytar", Form("Run %d ;Raster X (cm) ; Ytar ",nrun), 600, -.3,.3,240, -6.,6);
   TH2F *hrastx_zb = new TH2F("hrastx_zb", Form("Run %d ;Raster X (cm) ; Z_beam ",nrun), 600, -.3,.3,600,-15.,15.);
   TH1F *hrasty = new TH1F("hrasty",Form("Run %d ;Raster Y (cm) ;Counts",nrun), 2000, -.3,.3);
   TH1F *hrastrawxa = new TH1F("hrastrawxa", Form("Run %d ;Raster Xa ADC;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hrastrawya = new TH1F("hrastrawya",Form("Run %d ;Raster Ya ADC ;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hrastrawxb = new TH1F("hrastrawxb", Form("Run %d ;Raster Xb ADC;Counts",nrun), 2000, 50000.,90000.);
   TH1F *hrastrawyb = new TH1F("hrastrawyb",Form("Run %d ;Raster Yb ADC ;Counts",nrun), 2000, 50000.,90000.);
    //
   Double_t zsieve=168.;
   Double_t xs,ys;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (sumnpe>2&&etracknorm > .6 && TMath::Abs(delta)<10) {
		  hytar->Fill(ytar);
		  hreactz->Fill(reactz);
		  hreactx->Fill(reactx);
		  hreacty->Fill(reacty);
		  xs=xtar+xptar*zsieve;
		  ys=ytar+yptar*zsieve;
		  hxs_ys->Fill(ys,xs);
                  for (Int_t nf=0;nf<nfoils;nf++) {
		    if (TMath::Abs(reactz-fcent[nf])<fwid[nf]) hxs_ys_foil[nf]->Fill(ys,xs);
		    if (TMath::Abs(reactz-fcent[nf])<fwid[nf]) hxs_foil[nf]->Fill(xs);
		    if (TMath::Abs(reactz-fcent[nf])<fwid[nf]) hys_foil[nf]->Fill(ys);
                  }
		  hrastx_ytar->Fill(raster_ax,ytar);
		  hrastx_zb->Fill(raster_ax,reactz);
		  hrastrawya->Fill(raster_raw_ay);
		  hrastrawxa->Fill(raster_raw_ax);
		  hrastrawyb->Fill(raster_raw_by);
		  hrastrawxb->Fill(raster_raw_bx);
		  hrasty->Fill(raster_ay);
		  hrastx->Fill(raster_ax);
		}
	}
	//
	TF1 *foilfit[3];
TCanvas *c3b = new TCanvas("c3b", "react", 700,700);
 c3b->Divide(2,2);
 c3b->cd(1);
 hytar->SetStats(0);
 hytar->Draw();
 c3b->cd(2);
 hreactx->Draw();
 c3b->cd(3);
 hreacty->Draw();
 c3b->cd(4);
 hreactz->Draw();
 for (Int_t nf=0;nf<nfoils;nf++) {
   foilfit[nf] = new TF1(Form("foilfit_%d",nf),"gaus",fcent[nf]-fwid[nf],fcent[nf]+fwid[nf]);
   foilfit[nf]->SetParameter(1,fcent[nf]);
   foilfit[nf]->SetParameter(2,fwid[nf]);
   foilfit[nf]->SetParameter(0,1000.);
   hreactz->Fit(Form("foilfit_%d",nf),"Q","SAME",fcent[nf]-fwid[nf],fcent[nf]+fwid[nf]);
   cout << foilfit[nf]->GetParameter(1) << " " << foilfit[nf]->GetParError(1) << " " <<foilfit[nf]->GetParameter(2) << " " <<foilfit[nf]->GetParError(2) <<endl;
 }
outputpdf="plots/"+basename+"_comp_hms_sieve.pdf";
c3b->Print(outputpdf);
//
TCanvas *c3 = new TCanvas("c3", "raw raster", 700, 700);
c3->Divide(2,2);
c3->cd(1);
hrastrawxa->Draw();
c3->cd(2);
hrastrawya->Draw();
c3->cd(3);
hrastrawxb->Draw();
c3->cd(4);
hrastrawyb->Draw();
    outputpdf="plots/"+basename+"_raster.pdf(";
c3->Print(outputpdf);
//
TCanvas *craster2d = new TCanvas("craster2d", "raster", 700, 700);
craster2d->Divide(1,2);
craster2d->cd(1);
hrastx_ytar->Draw("colz");
craster2d->cd(2);
hrastx_zb->Draw("colz");
    outputpdf="plots/"+basename+"_raster.pdf";
craster2d->Print(outputpdf);
//
//
TCanvas *craster = new TCanvas("craster", "raster", 700, 700);
craster->Divide(1,2);
craster->cd(1);
hrastx->Draw();
craster->cd(2);
hrasty->Draw();
    outputpdf="plots/"+basename+"_raster.pdf)";
craster->Print(outputpdf);
//
//
       Double_t xdist=2.54;
       Double_t ydist=1.524;
       Double_t zdist=168.;
       Int_t nyholes=9.;
       Int_t nxholes=9.;
       TLine* yline;
TCanvas *csieve = new TCanvas("csieve", "sieve", 700, 700);
csieve->Divide(2,2);
csieve->cd(1);
 hxs_ys->SetStats(0);
hxs_ys->Draw("colz");
 for (Int_t nf=0;nf<nfoils;nf++) {
csieve->cd(nf+2);
 hxs_ys_foil[nf]->SetStats(0);
 hxs_ys_foil[nf]->Draw("colz");
       Double_t minx= hxs_ys_foil[nf]->GetXaxis()->GetXmin();
       Double_t maxx= hxs_ys_foil[nf]->GetXaxis()->GetXmax();
       Double_t miny= hxs_ys_foil[nf]->GetYaxis()->GetXmin();
       Double_t maxy= hxs_ys_foil[nf]->GetYaxis()->GetXmax();
      for (Int_t nh=0;nh<nyholes;nh++) {
         Double_t ypos=(nh-(nyholes-1)/2.)*ydist;
	 yline= new TLine(ypos,miny,ypos,maxy);
	 yline->Draw();
	 yline->SetLineColor(2);
       }
       for (Int_t nh=0;nh<nxholes;nh++) {
         Double_t xpos=(nh-(nxholes-1)/2.)*xdist;
	 yline= new TLine(minx,xpos,maxx,xpos);
	 yline->Draw();
	 yline->SetLineColor(2);
       }
 }
    outputpdf="plots/"+basename+"hms_sieve.pdf";
csieve->Print(outputpdf);
//
 TSpectrum *sfoil[3];
 Double_t xs_cent_measured[3][9];
 Double_t xs_cent_true[3][9];
TCanvas *cxsieve = new TCanvas("cxsieve", "xsieve", 700, 700);
 cxsieve->Divide(2,nfoils);
 for (Int_t nf=0;nf<nfoils;nf++) {
cxsieve->cd(2*nf+1);
 hxs_foil[nf]->Draw();
 TH1F *h2 = (TH1F*)hxs_foil[nf]->Clone("h2");
 sfoil[nf] = new TSpectrum(2*nxholes);
 Int_t nfound= sfoil[nf]->Search(hxs_foil[nf],2,"",.1);
 TH1 *hb = sfoil[nf]->Background(hxs_foil[nf],20,"same");
 Double_t par[100];
 if (hb) cxsieve->Update();
 if (nfound>0) {
   TF1 *fline = new TF1("fline","pol1",-12,12.);
   hxs_foil[nf]->Fit("fline","qn");
   par[0]= fline->GetParameter(0);
   par[1]= fline->GetParameter(1);
   Double_t *xpeaks;
   xpeaks = sfoil[nf]->GetPositionX();
   npeaks=0;
   for (Int_t p=0;p<nfound;p++) {
     Double_t xp=xpeaks[p];
     Int_t bin= hxs_foil[nf]->GetXaxis()->FindBin(xp);
     Double_t yp=hxs_foil[nf]->GetBinContent(bin);
     par[3*npeaks+2] = yp;
     par[3*npeaks+3] = xp;
     par[3*npeaks+4] = .1;
     npeaks++;
   }
   TF1 *fit = new TF1("fit",fpeaks,-12,12,2+3*npeaks);
   TVirtualFitter::Fitter(h2,10+3*npeaks);
   fit->SetParameters(par);
   cxsieve->cd(2*nf+2);
   h2->Fit("fit","qn");
   fit->SetParameters(par);
   h2->Fit("fit","q");
 }
 // nfoils
 }
//
 TSpectrum *ysfoil[3];
 Double_t ys_cent_measured[3][9];
 Double_t ys_cent_true[3][9];
TCanvas *cysieve = new TCanvas("cysieve", "ysieve", 700, 700);
 cysieve->Divide(2,nfoils);
 for (Int_t nf=0;nf<nfoils;nf++) {
cysieve->cd(2*nf+1);
 hys_foil[nf]->Draw();
 TH1F *h2 = (TH1F*)hys_foil[nf]->Clone("h2");
 ysfoil[nf] = new TSpectrum(2*nxholes);
 Int_t nfound= ysfoil[nf]->Search(hys_foil[nf],2,"",.05);
 TH1 *hb = ysfoil[nf]->Background(hys_foil[nf],20,"same");
 Double_t par[100];
 if (hb) cysieve->Update();
 if (nfound>0) {
   TF1 *fline = new TF1("fline","pol1",-12,12.);
   hys_foil[nf]->Fit("fline","qn");
   par[0]= fline->GetParameter(0);
   par[1]= fline->GetParameter(1);
   Double_t *xpeaks;
   xpeaks = ysfoil[nf]->GetPositionX();
   npeaks=0;
   for (Int_t p=0;p<nfound;p++) {
     Double_t xp=xpeaks[p];
     Int_t bin= hys_foil[nf]->GetXaxis()->FindBin(xp);
     Double_t yp=hys_foil[nf]->GetBinContent(bin);
     par[3*npeaks+2] = yp;
     par[3*npeaks+3] = xp;
     par[3*npeaks+4] = .1;
     npeaks++;
   }
   TF1 *fit = new TF1("fit",fpeaks,-12,12,2+3*npeaks);
   TVirtualFitter::Fitter(h2,10+3*npeaks);
   fit->SetParameters(par);
   cysieve->cd(2*nf+2);
   h2->Fit("fit","qn");
   fit->SetParameters(par);
   h2->Fit("fit","q");
 }
 // nfoils
 }

   //
}
