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

void plot_shms_trig_test(TString basename,TString label="") {
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
 outputpdf="plots/"+basename+".pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_trig_test_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
 static const Int_t trignum=15;
 const char* trigname[trignum]={
 "D.pdcref1",
 "D.pdcref2",
 "D.pdcref3",
 "D.pdcref4",
 "D.pdcref5",
 "D.pdcref6",
 "D.pdcref7",
 "D.pdcref8",
 "D.pdcref9",
 "D.pdcref10",
 "D.phodoref1",
"D.phodoref2",
 "D.ptrigref1",
 "D.ptrigref2",
"D.pedtm"
};
 TH1F *hist[trignum][5][5];
 TH1F *nhist[trignum];
 //
 for (Int_t ip=0;ip<trignum;ip++) {
       nhist[ip] = (TH1F*)fhistroot->Get(Form("hist_ndata_%s",trigname[ip]));
       if (!nhist[ip]) cout << " no hist = " << Form("hist_ndata_%s",trigname[ip]) << endl;
   for (Int_t ic=0;ic<5;ic++) {
     for (Int_t ih=0;ih<ic+1;ih++) {
       hist[ip][ic][ih]= (TH1F*)fhistroot->Get(Form("hist_%s_m%d_h%d",trigname[ip],ic+1,ih+1));
       if (!hist[ip][ic][ih]) cout << " no hist = " << Form("hist_%s_m%d_h%d",trigname[ip],ic+1,ih+1) << endl;
   }
   }
 }
 static const Int_t dnum=12;
   TH1F *dhist[dnum];
 const char* dname[dnum]={
"hist_D.ptrigref1_mult_2_diff_1_2",
"hist_D.ptrigref1_mult_3_diff_1_2",
"hist_D.ptrigref1_mult_3_diff_2_3",
"hist_D.ptrigref1_mult_4_diff_1_2",
"hist_D.ptrigref1_mult_4_diff_2_3",
"hist_D.ptrigref1_mult_4_diff_3_4",
"hist_cut_D.ptrigref1_mult_2_diff_1_2",
"hist_cut_D.ptrigref1_mult_3_diff_1_2",
"hist_cut_D.ptrigref1_mult_3_diff_2_3",
"hist_cut_D.ptrigref1_mult_4_diff_1_2",
"hist_cut_D.ptrigref1_mult_4_diff_2_3",
"hist_cut_D.ptrigref1_mult_4_diff_3_4",
 };
 for (Int_t ip=0;ip<dnum;ip++) {
       dhist[ip] = (TH1F*)fhistroot->Get(dname[ip]);
       if (!dhist[ip]) cout << " no hist = " << dname[ip] << endl;
 }
 // canvas
 TCanvas *cdiff ;
     cdiff= new TCanvas(Form("diff_%d",0),"Diff Time", 700,700);
 cdiff->Divide(2,3);
     for (Int_t ih=0;ih<6;ih++) {
     cdiff->cd(ih+1);
     gPad->SetLogy();
     if (dhist[ih]) dhist[ih]->Draw();
     }
     //
     TCanvas *cdiff2 ;
     cdiff2= new TCanvas(Form("diff2_%d",0),"Diff Time NDC==0", 700,700);
    cdiff2->Divide(2,3);
     for (Int_t ih=6;ih<dnum;ih++) {
     cdiff2->cd(ih-5);
     gPad->SetLogy();
     if (dhist[ih])     dhist[ih]->Draw();
     }
//
 TCanvas *cmult;
     cmult= new TCanvas(Form("mult_%d",0),"Multiplicity", 700,700);
     cmult->Divide(2,2);
     cmult->cd(1);
     gPad->SetLogy();
     nhist[0]->Draw();
     cmult->cd(2);
     gPad->SetLogy();
     nhist[1]->Draw();
     cmult->cd(3);
     gPad->SetLogy();
     nhist[12]->Draw();
     cmult->cd(4);
     gPad->SetLogy();
     nhist[13]->Draw();
 TCanvas *ctrig[5];
   for (Int_t ic=0;ic<5;ic++) {
     ctrig[ic]= new TCanvas(Form("ctrig_%d",ic),Form("Trig 1  Mult=%d",ic), 700,700);
     if (ic==0) ctrig[ic]->Divide(1,1);
     if (ic==1) ctrig[ic]->Divide(1,2);
     if (ic==2) ctrig[ic]->Divide(1,3);
     if (ic==3) ctrig[ic]->Divide(2,2);
      if (ic==4) ctrig[ic]->Divide(2,3);
     for (Int_t ih=0;ih<ic+1;ih++) {
     ctrig[ic]->cd(ih+1);
     gPad->SetLogy();
     hist[12][ic][ih]->Draw();
     }
   }
  TCanvas *ctime[5];
   for (Int_t ic=0;ic<5;ic++) {
     ctime[ic]= new TCanvas(Form("ctime_%d",ic),Form("DCref1  Mult=%d",ic), 700,700);
     if (ic==0) ctime[ic]->Divide(1,1);
     if (ic==1) ctime[ic]->Divide(1,2);
     if (ic==2) ctime[ic]->Divide(1,3);
     if (ic==3) ctime[ic]->Divide(2,2);
      if (ic==4) ctime[ic]->Divide(2,3);
     for (Int_t ih=0;ih<ic+1;ih++) {
     ctime[ic]->cd(ih+1);
     gPad->SetLogy();
     hist[0][ic][ih]->Draw();
     }
   }


   //
}
