#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void fit_xfp_yfp_hist_hms_ztar_sieve(TString basename,TString label) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
 //
      TString foilname[3]={"Z=-10cm","Z=0","Z=10cm"};
 outputpdf="plots/"+basename+"_fit_xfp_yfp_hms_ztar_sieve_hist.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_hms_ztar_sieve_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   TH2F* hyfp_yxfp_foil[3];
   TH2F* hys_xs_foil[3];
    TH2F* hyfp_yxfp_foil_ypcut[3][7];
       TH1F *yphist[3][7];
  	for (Int_t nf = 0; nf < 3; nf++) {
            hyfp_yxfp_foil[nf]  = (TH2F*)fhistroot->Get(Form("hyfp_yxfp_foil_%d",nf));
            hys_xs_foil[nf]  = (TH2F*)fhistroot->Get(Form("hys_xs_foil_%d",nf));
	for (Int_t iz = 0; iz < 7; iz++) {
	  hyfp_yxfp_foil_ypcut[nf][iz] = (TH2F*)fhistroot->Get(Form("hyfp_yxfp_foil_%d_ypcut_%d",nf,iz));
  	  yphist[nf][iz] = (TH1F*)fhistroot->Get(Form("hyptar_%d_%d",nf,iz));
	}
	}
	//
       Double_t yptar_cent[7];
	  for (Int_t n=0;n<7;n++) {
	    yptar_cent[n] = yphist[1][n]->GetMean();
	  }

	    //
 TF1 *linfit[7];
 Double_t xfp_cent[7];
 Double_t xfp_cent_err[7];
 Double_t yfp_cent[7];
 Double_t yfp_cent_err[7];
	    TCanvas *cfit = new TCanvas("cfit"," Ind fits",900,700);
	  cfit->Divide(4,2);
	  for (Int_t n=0;n<7;n++) {
            linfit[n] = new TF1(Form("linfit_%d",n),"pol1",-20,20);
            linfit[n]->SetLineColor(2);
            linfit[n]->SetLineWidth(2);
	    cfit->cd(n+1);
	    hyfp_yxfp_foil_ypcut[1][n]->Draw("colz");
            hyfp_yxfp_foil_ypcut[1][n]->Fit(Form("linfit_%d",n),"Q");
            Double_t err0=linfit[n]->GetParError(0)/linfit[n]->GetParameter(0);
            Double_t err1=linfit[n]->GetParError(1)/linfit[n]->GetParameter(1);
	    Double_t rat=-linfit[n]->GetParameter(0)/linfit[n]->GetParameter(1);
	    Double_t raterr = rat*TMath::Sqrt(err0*err0+err1*err1);
	    yfp_cent[n]=yptar_cent[n];
	    yfp_cent_err[n]=0.01;
	    xfp_cent[n]=rat;
	    xfp_cent_err[n]=raterr;
	    cout << yptar_cent[n] << " " << linfit[n]->GetParameter(0) << " " << linfit[n]->GetParError(0) << " " <<linfit[n]->GetParameter(1) << " " <<linfit[n]->GetParError(1) << "  Xfp (yfp=0) " << rat << " +/- " << raterr<< endl;
	  }
	    cfit->Print(outputpdf+"(");
	  //
	    TCanvas *call = new TCanvas("call","yfp v xfp",700,700);
	    call->Divide(1,1);
	    call->cd(1);
	    hyfp_yxfp_foil[1]->Draw("colz");
            TString nlabel = label+" "+foilname[1];
	    hyfp_yxfp_foil[1]->SetTitle(nlabel);
	  for (Int_t n=0;n<7;n++) {
            linfit[n]->Draw("same");
	  }
	    call->Print(outputpdf);
	//
	    TGraphErrors *gr_xfp = new TGraphErrors(7,yfp_cent,xfp_cent,yfp_cent_err,xfp_cent_err);
	    TCanvas *cgr = new TCanvas("cgr"," xfp cent",700,700);
	    cgr->Divide(1,1);
	    gr_xfp->SetMinimum(-5.);
	    gr_xfp->SetMaximum(5.);
	    gr_xfp->Draw();
	    gr_xfp->Fit("pol0");
	  //
      TLine *nxline[9];
      TLine *nyline[9];
      for (int nx=0;nx<9;nx++) {
	nxline[nx] = new TLine((4-nx)*1.524-.25,-15,(4-nx)*1.524-.25,15);
	nyline[nx] = new TLine(-8,(4-nx)*2.54,8,(4-nx)*2.54);
	  }
	    TCanvas *cs = new TCanvas("cs","ys v xs",1300,700);
	    cs->Divide(1,1);
	    for (int nf=1;nf<2;nf++) {
	      //	    cs->cd(nf+1);
	    hys_xs_foil[nf]->SetMinimum(10);
	    hys_xs_foil[nf]->Draw("colz");
            TString nlabel = label+" "+foilname[nf];
	    hys_xs_foil[nf]->SetTitle(nlabel);
      for (int nx=0;nx<9;nx++) {
       	nxline[nx]->Draw();
       	nyline[nx]->Draw();
      }
	    }
	    cs->Print(outputpdf+")");
}
