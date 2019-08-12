#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
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
void comp_hist_vcs( TString fn1 , TString fn2 ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/comp_hist_vcs_"+fn1+".pdf";
    //
    Double_t dummy_ch=1;
    Double_t data_ch=1;
    Double_t tot_ch=0;
    if (fn1=="all_kin3B") {
    dummy_ch=61.27;
    data_ch=363.969;
    tot_ch=4.4;
    }
    if (fn1=="all_kin2A") {
    dummy_ch=83.4759;
    data_ch=227.264;
    tot_ch=1.0;
    }
    if (fn1=="all_kin2B") {
    dummy_ch=61.49;
    data_ch=397.627;
    tot_ch=3.71;
    }
    if (fn1=="all_kin2B_pass1") {
    dummy_ch=135.244;
    data_ch=175.71;
    tot_ch=1;
    }
    if (fn1=="all_kin1B") {
    dummy_ch=59.63;
    data_ch=348.19;
    tot_ch=4.25;
    }
    if (fn1=="all_kin1A") {
    dummy_ch=64.82;
    data_ch=293.73;
    tot_ch=2.92;
    }
    Double_t dummy_wall=0.1816+0.1815;
    Double_t data_wall=0.104+0.133;
    Double_t Ctime_ran_sig=0.5*6;
    Double_t Ctime_sig=0.4*3;
    Double_t dummy_scal_fac =1;
    dummy_scal_fac = (data_ch/dummy_ch)/4.;
    cout << data_ch << " " << dummy_ch << endl;
const UInt_t nftot=2;
 Int_t colind[nftot]={1,2};
 TString lname[nftot]={"data","dummy"};
 lname[0]=fn1+" Q2 = 0.33 ";
 lname[1]=lname[1]+" Q2 = 0.42 ";
 TString lname_hiQ2[nftot]={"data","dummy"};
 lname_hiQ2[0]=fn1+" Q2 = 0.42 ";
 lname_hiQ2[1]=lname_hiQ2[1]+" Q2 = 0.42 ";
   TFile *fhistroot[nftot];
const UInt_t nplots=2;
 TString h1name[nplots]={"hMM2_vcs_cut","hMM2_ran"};
 TString h2name[nplots]={"hMM2_vcs_cut_scale","hMM2_ran_scale"};
 TString h1name_hiQ2[nplots]={"hMM2_vcs_hiQ2_cut","hMM2_ran_vcs_hiQ2_cut"};
 TString h2name_hiQ2[nplots]={"hMM2_vcs_hiQ2_cut_scale","hMM2_ran_vcs_hiQ2_cut_scale"};
  TH1F *fhist[nftot][nplots];
  TH1F *fhist_scale[nftot][nplots];
  TH1F *fhist_hiQ2[nftot][nplots];
  TH1F *fhist_hiQ2_scale[nftot][nplots];
  TH1F *fhist_acc_sub[nftot];
  TH1F *fhist_hiQ2_acc_sub[nftot];
  TH2F *fhist_MM2_thcm[nftot];
  TH1F *fhist_dum_sub;
  TH1F *fhist_hiQ2_dum_sub;
 TString inputroot[2];
 inputroot[0] = "hist/"+fn1+"_vcs_hist.root";
 inputroot[1] = "hist/"+fn2+"_vcs_hist.root";
 Double_t scale_fac[nftot][nplots]={1.,Ctime_sig/Ctime_ran_sig,1.,Ctime_sig/Ctime_ran_sig};
   for (UInt_t nf=0;nf<nftot;nf++) {
     cout << " infile root = " << inputroot[nf] << endl;
   fhistroot[nf] =  new TFile(inputroot[nf]);
     fhist_MM2_thcm[nf] = (TH2F*)fhistroot[nf]->Get("hMM2_thcm");
     fhist_MM2_thcm[nf]->SetTitle(lname[nf]);
   for (UInt_t nh=0;nh<nplots;nh++) {
     fhist[nf][nh] = (TH1F*)fhistroot[nf]->Get(h1name[nh]);
     fhist_scale[nf][nh] = (TH1F*)(fhist[nf][nh]->Clone(h2name[nh]));
     cout << scale_fac[nf][nh] << endl;
     fhist_scale[nf][nh]->Scale(scale_fac[nf][nh]);
     fhist[nf][nh]->SetTitle(lname[nf]);
      fhist_hiQ2[nf][nh] = (TH1F*)fhistroot[nf]->Get(h1name_hiQ2[nh]);
     fhist_hiQ2_scale[nf][nh] = (TH1F*)(fhist_hiQ2[nf][nh]->Clone(h2name_hiQ2[nh]));
    fhist_hiQ2_scale[nf][nh]->Scale(scale_fac[nf][nh]);
     fhist_hiQ2[nf][nh]->SetTitle(lname_hiQ2[nf]);
   }
   }
    for (UInt_t nf=0;nf<nftot;nf++) {
       fhist_acc_sub[nf] = (TH1F*)(fhist_scale[nf][0]->Clone("hMM2_Wcut_acc_sub"));
       fhist_acc_sub[nf]->Add(fhist_scale[nf][1],-1.);
    }
    cout << " Dummy scale fac = " << dummy_scal_fac<< endl;
       fhist_dummy_sub = (TH1F*)(fhist_acc_sub[0]->Clone("hMM2_Wcut_dum_sub"));
       fhist_dummy_sub->Add(fhist_acc_sub[1],-dummy_scal_fac);
       //       
    for (UInt_t nf=0;nf<nftot;nf++) {
       fhist_hiQ2_acc_sub[nf] = (TH1F*)(fhist_hiQ2_scale[nf][0]->Clone("hMM2_Wcut_acc_sub"));
       fhist_hiQ2_acc_sub[nf]->Add(fhist_hiQ2_scale[nf][1],-1.);
    }
       fhist_hiQ2_dummy_sub = (TH1F*)(fhist_hiQ2_acc_sub[0]->Clone("hMM2_Wcut_dum_sub"));
       fhist_hiQ2_dummy_sub->Add(fhist_hiQ2_acc_sub[1],-dummy_scal_fac);
       
     //
       TCanvas *cW;
    TLegend *lW[2];
     cW = new TCanvas("cW","W", 1000,700);
     cW->Divide(1,2);
     for (UInt_t nf=0;nf<nftot;nf++) {
       cW->cd(nf+1);
       lW[nf] = new TLegend(.79,.65,.99,.95,"");
       fhist[nf][0]->Draw();
       fhist_scale[nf][1]->Draw("same");
       fhist_scale[nf][1]->SetLineColor(2);
       fhist_acc_sub[nf]->Draw("same");
       fhist_acc_sub[nf]->SetLineColor(3);
       lW[nf]->AddEntry(fhist[nf][0],"Data");
       lW[nf]->AddEntry(fhist_scale[nf][1],"Acc coincidence");
       lW[nf]->AddEntry(fhist_acc_sub[nf],"Data -Acc");
       lW[nf]->Draw();
    }
     cW->Print(outputpdf+"(");
     //
      TCanvas *c2d;
      c2d = new TCanvas("c2d","2d", 1000,700);
      c2d->Divide(1,1);
      c2d->cd(1);
      gPad->SetLogz(1);
      fhist_MM2_thcm[0]->SetMaximum(200);
      fhist_MM2_thcm[0]->Draw("colz");
      c2d->Print(outputpdf);
      //     
      TCanvas *cSub;
      TF1* fMM2;
     cSub = new TCanvas("cSub","Sub", 1000,700);
     cSub->Divide(1,1);
     fhist_dummy_sub->Draw();
     fhist_dummy_sub->SetTitle(" Data - dummy ");
     Double_t MM2_Max = fhist_dummy_sub->GetBinCenter(fhist_dummy_sub->GetMaximumBin());
     fMM2 = new TF1("fMM2","gaus",MM2_Max-.006,MM2_Max+.006);
     fhist_dummy_sub->Fit("fMM2","QR");
     Double_t MM2_peak=fMM2->GetParameter(1);
     Double_t MM2_sig=fMM2->GetParameter(2);
     Double_t cnts=fhist_dummy_sub->Integral(fhist_dummy_sub->FindBin(MM2_peak-3*MM2_sig),fhist_dummy_sub->FindBin(MM2_peak+3*MM2_sig));
     Double_t dummy_cnts=fhist_acc_sub[1]->Integral(fhist_acc_sub[1]->FindBin(MM2_peak-3*MM2_sig),fhist_acc_sub[1]->FindBin(MM2_peak+3*MM2_sig));
     Double_t data_cnts=fhist_acc_sub[0]->Integral(fhist_acc_sub[0]->FindBin(MM2_peak-3*MM2_sig),fhist_acc_sub[0]->FindBin(MM2_peak+3*MM2_sig));
     cout << fn1 << " Integral = " << cnts << " Yield/C = " << cnts/(data_ch/1000.) << " tot charge = " << tot_ch << " total vcs counts = " << tot_ch*cnts/(data_ch/1000.)<< endl;
     cout << fn1 << "dummy cnts = " << dummy_cnts << " Ratio dummy/data = " << dummy_cnts/data_cnts << endl;
     cout << MM2_peak << " " << MM2_peak-3*MM2_sig << " " << MM2_peak+3*MM2_sig <<endl;
     cSub->Print(outputpdf);
    //
     cout << " HI Q2 " << endl;
     //
       TCanvas *cW_hiQ2;
    TLegend *lW_hiQ2[2];
     cW_hiQ2 = new TCanvas("cW_hiQ2","W_hiQ2", 1000,700);
     cW_hiQ2->Divide(1,2);
     for (UInt_t nf=0;nf<nftot;nf++) {
       cW_hiQ2->cd(nf+1);
       lW_hiQ2[nf] = new TLegend(.79,.65,.99,.95,"");
       fhist_hiQ2[nf][0]->Draw();
       fhist_hiQ2_scale[nf][1]->Draw("same");
       fhist_hiQ2_scale[nf][1]->SetLineColor(2);
       fhist_hiQ2_acc_sub[nf]->Draw("same");
       fhist_hiQ2_acc_sub[nf]->SetLineColor(3);
       lW_hiQ2[nf]->AddEntry(fhist_hiQ2[nf][0],"Data");
       lW_hiQ2[nf]->AddEntry(fhist_hiQ2_scale[nf][1],"Acc coincidence");
       lW_hiQ2[nf]->AddEntry(fhist_hiQ2_acc_sub[nf],"Data -Acc");
       lW_hiQ2[nf]->Draw();
    }
     cW_hiQ2->Print(outputpdf);
      TCanvas *cSub_hiQ2;
      TF1* fMM2_hiQ2;
     cSub_hiQ2 = new TCanvas("cSub_hiQ2","Sub_hiQ2", 1000,700);
     cSub_hiQ2->Divide(1,1);
     fhist_hiQ2_dummy_sub->Draw();
     fhist_hiQ2_dummy_sub->SetTitle(" Data - dummy Q2 = 0.43");
     MM2_Max = fhist_hiQ2_dummy_sub->GetBinCenter(fhist_hiQ2_dummy_sub->GetMaximumBin());
     MM2_Max = .002;
     fMM2_hiQ2 = new TF1("fMM2_hiQ2","gaus",MM2_Max-.01,MM2_Max+.01);
     fhist_hiQ2_dummy_sub->Fit("fMM2_hiQ2","QR");
     MM2_peak=fMM2_hiQ2->GetParameter(1);
     MM2_sig=fMM2_hiQ2->GetParameter(2);
     cnts=fhist_hiQ2_dummy_sub->Integral(fhist_hiQ2_dummy_sub->FindBin(MM2_peak-3*MM2_sig),fhist_hiQ2_dummy_sub->FindBin(MM2_peak+3*MM2_sig));
     dummy_cnts=fhist_hiQ2_acc_sub[1]->Integral(fhist_hiQ2_acc_sub[1]->FindBin(MM2_peak-3*MM2_sig),fhist_hiQ2_acc_sub[1]->FindBin(MM2_peak+3*MM2_sig));
     data_cnts=fhist_hiQ2_acc_sub[0]->Integral(fhist_hiQ2_acc_sub[0]->FindBin(MM2_peak-3*MM2_sig),fhist_hiQ2_acc_sub[0]->FindBin(MM2_peak+3*MM2_sig));
     cout << fn1 << " Integral = " << cnts << " Yield/C = " << cnts/(data_ch/1000.) << " tot charge = " << tot_ch << " total vcs counts = " << tot_ch*cnts/(data_ch/1000.)<< endl;
     cout << fn1 << "dummy cnts = " << dummy_cnts << " Ratio dummy/data = " << dummy_cnts/data_cnts << endl;
     cout << MM2_peak << " " << MM2_peak-3*MM2_sig << " " << MM2_peak+3*MM2_sig <<endl;
     cSub_hiQ2->Print(outputpdf+")");
    //
}
