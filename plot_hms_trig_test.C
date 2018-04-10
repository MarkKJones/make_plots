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

void plot_hms_trig_test(TString basename,TString label="") {
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
 outputpdf="plots/"+basename+"_trig_test_hist.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_trig_test_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
  static const Int_t trignum=11;
 const char* trigname[trignum]={
 "D.hdcref1",
 "D.hdcref2",
 "D.hdcref3",
 "D.hdcref4",
 "D.hhodoref1",
 "D.htrigref1",
 "D.htrigref2",
 "D.htrigref3",
 "D.htrig1",
 "D.htrig2",
 "D.htrig3",
};
 TH1F *hist[trignum][5][5];
 TH1F *nhist[trignum];
 //
 for (Int_t ip=0;ip<trignum;ip++) {
       nhist[ip] = (TH1F*)fhistroot->Get(Form("hist_ndata_%s",trigname[ip]));
   for (Int_t ic=0;ic<5;ic++) {
     for (Int_t ih=0;ih<ic+1;ih++) {
       hist[ip][ic][ih]= (TH1F*)fhistroot->Get(Form("hist_%s_m%d_h%d",trigname[ip],ic+1,ih+1));
   }
   }
 }
 TH2F *mhist[3];
 mhist[0] = (TH2F*)fhistroot->Get(Form("hist_mult_%s_%s",trigname[5],trigname[6]));
 mhist[1] = (TH2F*)fhistroot->Get(Form("hist_mult_%s_%s",trigname[5],trigname[0]));
 mhist[2] = (TH2F*)fhistroot->Get(Form("hist_mult_%s_%s",trigname[6],trigname[0]));
 TH2F *tref2hist[3];
 tref2hist[0] = (TH2F*)fhistroot->Get(Form("hist_time_%s_hit12",trigname[6]));
tref2hist[1] = (TH2F*)fhistroot->Get(Form("hist_time_%s_hit23",trigname[6]));
tref2hist[2] = (TH2F*)fhistroot->Get(Form("hist_time_%s_hit13",trigname[6]));
 TH2F *tref1hist[3];
tref1hist[0] = (TH2F*)fhistroot->Get(Form("hist_time_%s_hit12",trigname[5]));
tref1hist[1] = (TH2F*)fhistroot->Get(Form("hist_time_%s_hit23",trigname[5]));
tref1hist[2] = (TH2F*)fhistroot->Get(Form("hist_time_%s_hit13",trigname[5]));
 // canvas
//
 TCanvas *cmult2;
     cmult2= new TCanvas(Form("mult2_%d",0),"2d Multiplicity", 700,700);
     cmult2->Divide(2,2);
      cmult2->cd(1);
    mhist[0]->Draw("colz");
      cmult2->cd(2);
    mhist[1]->Draw("colz");
       cmult2->cd(3);
    mhist[2]->Draw("colz");
    cmult2->Print(outputpdf+"(");
    //
 TCanvas *cmult;
     cmult= new TCanvas(Form("mult_%d",0),"Multiplicity", 700,700);
     cmult->Divide(2,3);
     cmult->cd(1);
     gPad->SetLogy();
     nhist[0]->Draw();
     cmult->cd(2);
     gPad->SetLogy();
     nhist[1]->Draw();
     cmult->cd(3);
     gPad->SetLogy();
     nhist[5]->Draw();
     cmult->cd(4);
     gPad->SetLogy();
     nhist[6]->Draw();
    cmult->cd(5);
     gPad->SetLogy();
     nhist[7]->Draw();
     cmult->cd(6);
     gPad->SetLogy();
     nhist[9]->Draw();
    cmult2->Print(outputpdf);
 TCanvas *ctrig[5];
   for (Int_t ic=0;ic<5;ic++) {
     ctrig[ic]= new TCanvas(Form("ctrig_%d",ic),Form("Trig Ref 1  Mult=%d",ic), 700,700);
     if (ic==0) ctrig[ic]->Divide(1,1);
     if (ic==1) ctrig[ic]->Divide(1,2);
     if (ic==2) ctrig[ic]->Divide(1,3);
     if (ic==3) ctrig[ic]->Divide(2,2);
      if (ic==4) ctrig[ic]->Divide(2,3);
     for (Int_t ih=0;ih<ic+1;ih++) {
     ctrig[ic]->cd(ih+1);
     gPad->SetLogy();
     hist[5][ic][ih]->Draw();
     }
    ctrig[ic]->Print(outputpdf);
   }
 TCanvas *ctrig2[5];
   for (Int_t ic=0;ic<5;ic++) {
     ctrig2[ic]= new TCanvas(Form("ctrig2_%d",ic),Form("Trig Ref 3  Mult=%d",ic), 700,700);
     if (ic==0) ctrig2[ic]->Divide(1,1);
     if (ic==1) ctrig2[ic]->Divide(1,2);
     if (ic==2) ctrig2[ic]->Divide(1,3);
     if (ic==3) ctrig2[ic]->Divide(2,2);
      if (ic==4) ctrig2[ic]->Divide(2,3);
     for (Int_t ih=0;ih<ic+1;ih++) {
     ctrig2[ic]->cd(ih+1);
     gPad->SetLogy();
     hist[7][ic][ih]->Draw();
     }
    ctrig2[ic]->Print(outputpdf);
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
     if (ic!=4) ctime[ic]->Print(outputpdf);
     if (ic==4) ctime[ic]->Print(outputpdf+")");
   }

TCanvas *ctref2;
     ctref2= new TCanvas("tref2"," 2d tref2", 700,700);
        ctref2->Divide(2,2);
    for (Int_t ih=0;ih<3;ih++) {
       ctref2->cd(ih+1);
       gPad->SetLogz();
       tref2hist[ih]->Draw("colz");
     }
TCanvas *ctref1;
     ctref1= new TCanvas("tref1"," 2d tref1", 700,700);
        ctref1->Divide(2,2);
    for (Int_t ih=0;ih<3;ih++) {
       ctref1->cd(ih+1);
       gPad->SetLogz();
       tref1hist[ih]->Draw("colz");
     }
   //
}
