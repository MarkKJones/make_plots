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
#include <TLine.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_hist_hms_sieve(TString basename="",Double_t fTheta_lab=15.0, Int_t nf = 3){
  //
   if (basename=="") {
     cout << " Input the basename of the root file" << endl;
     cin >> basename;
   }
   Double_t zoff=0.1;
     Double_t zf[3]={10.+zoff,0.+zoff,-10.+zoff};
   if (nf==2) {
     zf[0] = 5.+zoff;
     zf[1] = -5.+zoff;
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
   inputroot="hist/"+basename+"_hms_sieve_hist.root";
   TString pdffile;
  TFile *fhistroot;
   fhistroot =  new TFile(inputroot);
   TH1F *hytar_cent[3];
   TH1F *hztar_cent[3];
   TH1F *hys[3];
   TH1F *hxs[3];
   TH1F *hyptar[3];
   TH1F *hxptar[3];
   TH2F *hxptar_yptar[3];
   TH2F *hxs_ys[3];
   Double_t fMispointing_y = 0.1*(0.52-0.012*TMath::Abs(fTheta_lab)+0.002*TMath::Abs(fTheta_lab)*TMath::Abs(fTheta_lab));
   Double_t fMispointing_x = 0.1*(2.37-0.086*TMath::Abs(fTheta_lab)+0.0012*TMath::Abs(fTheta_lab)*TMath::Abs(fTheta_lab));
   //   if (fTheta_lab<22.) fMispointing_y = 0.15;
   Double_t ztar_exp=0.1;
   for (int nd=0;nd<3;nd++) {
        hxptar_yptar[nd] = (TH2F*)fhistroot->Get(Form("hxptar_yptar_%d",nd));
        hxs_ys[nd] = (TH2F*)fhistroot->Get(Form("hxs_ys_%d",nd));
        hys[nd] = (TH1F*)fhistroot->Get(Form("hys_%d",nd));
        hxs[nd] = (TH1F*)fhistroot->Get(Form("hxs_%d",nd));
        hyptar[nd] = (TH1F*)fhistroot->Get(Form("hyptar_%d",nd));
        hxptar[nd] = (TH1F*)fhistroot->Get(Form("hxptar_%d",nd));
        hytar_cent[nd] = (TH1F*)fhistroot->Get(Form("hytar_cent_%d",nd));
        hztar_cent[nd] = (TH1F*)fhistroot->Get(Form("hztar_cent_%d",nd));
  }
TH1F* hbeamy =(TH1F*)fhistroot->Get("hbeamy");
TH1F* hbeamx =(TH1F*)fhistroot->Get("hbeamx");
   //
//
     TF1 *gfit_ztar_cent = new TF1("gaus_ztar_cent","gaus",-.05,.05);
     TF1 *gfit_ytar_cent = new TF1("gaus_ytar_cent","gaus",-.05,.05);
     TF1 *gfit_yptar_cent = new TF1("gaus_yptar_cent","gaus",-.05,.05);
     TF1 *gfit_ys_cent = new TF1("gaus_ys_cent","gaus",-.05,.05);
     TF1 *gfit_xptar_cent = new TF1("gaus_xptar_cent","gaus",-.05,.05);
     TF1 *gfit_xs_cent = new TF1("gaus_xs_cent","gaus",-.05,.05);
     TCanvas *cplot = new TCanvas("cplot","plots", 900,700);
      cplot->Divide(2,2);
      cplot->cd(1);
      hbeamx->Draw();
      Double_t beamx=hbeamx->GetMean();
      cplot->cd(2);
      hbeamy->Draw();
       Double_t beamy=hbeamy->GetMean();
  Double_t xtarver[3];
   Double_t ytarver[3];
   Double_t ztarver[3];
   Double_t sinTheta=TMath::Sin(fTheta_lab/57.3);
   Double_t cosTheta=TMath::Cos(fTheta_lab/57.3);
   for ( Int_t n=0 ; n<nf ;n++) {
     xtarver[n] = -beamy;
     ytarver[n] = zf[n]*sinTheta + beamx*cosTheta -fMispointing_y ;
     ztarver[n] = zf[n]*cosTheta - beamx*sinTheta  ;
   }
      Double_t ytar_exp = ztar_exp*TMath::Sin(fTheta_lab/57.3) + beamx*TMath::Cos(fTheta_lab/57.3)-fMispointing_y ;
      Double_t yptar_exp = ytar_exp/(166-ztarver[1]);
       cplot->cd(3);
      hztar_cent[1]->Draw();
      hztar_cent[1]->Fit("gaus_ztar_cent","Q","",hztar_cent[1]->GetMean()-3*hztar_cent[1]->GetRMS(),hztar_cent[1]->GetMean()+3*hztar_cent[1]->GetRMS());
       cplot->cd(4);
      hytar_cent[1]->Draw();
      hytar_cent[1]->Fit("gaus_ytar_cent","Q","",hytar_cent[1]->GetMean()-3*hytar_cent[1]->GetRMS(),hytar_cent[1]->GetMean()+3*hytar_cent[1]->GetRMS());
//
      TCanvas *cplot1 = new TCanvas("cplot1","ang plots", 900,700);
      cplot1->Divide(2,2);
      cplot1->cd(1);
      hxptar[1]->Draw();
     hxptar[1]->Fit("gaus_xptar_cent","Q","",-.005,.005);
      cplot1->cd(2);
      hyptar[1]->Draw();
     hyptar[1]->Fit("gaus_yptar_cent","Q","",-.005,.005);
       cplot1->cd(3);
      hxs[1]->Draw();
      hxs[1]->Fit("gaus_xs_cent","Q","",-.8,.8);
       cplot1->cd(4);
      hys[1]->Draw();
      hys[1]->Fit("gaus_ys_cent","Q","",-.8,.8);
      cout << " Theta = " << fTheta_lab <<endl;
      cout << " Mis_y = " << fMispointing_y << " beam_x = " << beamx << endl;
      cout << " Ztar_expect  = " << ztar_exp                 << " Ytar expect = " << ytar_exp<< endl;
      cout << " Ztar meas    = " << gfit_ztar_cent->GetParameter(1) << " Ytar meas   = " << gfit_ytar_cent->GetParameter(1)<<endl;
      cout << "yptar_exp = " << yptar_exp << " yptar_meas = " << gfit_yptar_cent->GetParameter(1) << endl;
      cout << " Mis_x = " << fMispointing_x << " beam_y = " << beamy << endl;
      Double_t xptar_exp = (fMispointing_x+beamy)/168.;
      cout << "xptar_exp = " << xptar_exp << " xptar_meas = " << gfit_xptar_cent->GetParameter(1) << endl;
     cout << " xs_meas = " << gfit_xs_cent->GetParameter(1) << " ys_meas = " << gfit_ys_cent->GetParameter(1) << endl;



     cout << " " << fTheta_lab << " " << ztar_exp << " " <<fMispointing_y << " " <<fMispointing_x << " " <<beamx << " " << beamy << " " << gfit_ztar_cent->GetParameter(1) << " " <<gfit_ytar_cent->GetParameter(1) << " " <<  yptar_exp << " " <<gfit_yptar_cent->GetParameter(1) << " " <<  xptar_exp <<" " <<gfit_xptar_cent->GetParameter(1) << " " << gfit_xs_cent->GetParameter(1) << " " <<gfit_ys_cent->GetParameter(1) << endl;
     // plot xs v ys for each foil
      TCanvas *cplot3[3];
      TCanvas *cplot2[3];
      TLine *nxline[9];
      TLine *nyline[9];
      for (int nx=0;nx<9;nx++) {
	    nxline[nx] = new TLine(-15,(4-nx)*1.524,15,(4-nx)*1.524);
	    nyline[nx] = new TLine((4-nx)*2.54,-8,(4-nx)*2.54,8);
	  }
    for (int nd=0;nd<nf;nd++) {
      cplot3[nd] = new TCanvas(Form("cplot3_%d",nd),Form("Ztar cut =%d",nd), 900,700);
      cplot3[nd]->Divide(1,1);
      cplot3[nd]->cd(1);
      hxs_ys[nd]->Draw("colz");
      for (int nx=0;nx<9;nx++) {
       	nxline[nx]->Draw();
       	nyline[nx]->Draw();
      }
    }
    //
      TLine *nxline2[3][9];
      TLine *nyline2[3][9];
      for (int n=0;n<nf;n++) {
	Double_t zdis= 168.-ztarver[n];
      for (int nx=0;nx<9;nx++) {
	Double_t yptemp=((4-nx)*1.524-ytarver[n])/zdis;
	Double_t xptemp=((4-nx)*2.54-xtarver[n])/zdis;
	nxline2[n][nx] = new TLine(-.1,yptemp,.1,yptemp);
	nyline2[n][nx] = new TLine(xptemp,-.05,xptemp,.05);
	  }
      }
    for (int nd=0;nd<nf;nd++) {
      cplot2[nd] = new TCanvas(Form("cplot2_%d",nd),Form("xp v yp Ztar cut =%d",nd), 900,700);
      cplot2[nd]->Divide(1,1);
      cplot2[nd]->cd(1);
      hxptar_yptar[nd]->Draw("colz");
      for (int nx=0;nx<9;nx++) {
       	nxline2[nd][nx]->Draw();
       	nyline2[nd][nx]->Draw();
      }
    }
       //
}
