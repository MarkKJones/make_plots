#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TCutG.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TPolyLine.h>
#include <TGraph.h>
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

void fit_shms_sieve(TString basename) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
    outputpdf="plots/"+basename+"_fit_shms_sieve.pdf";
const Int_t nftot=4;  
   TFile *fhistroot[nftot];
  TString inputroot;
  //  TString fname[nftot]={"shms_replay_matrixopt_1703_newmatrix","shms_replay_matrixopt_1703_newmatrix_orig_survey"};
  //  TString fname[nftot]={"shms_replay_matrixopt_1778_newmatrix","shms_replay_matrixopt_1796_newmatrix"};
  TString fname[nftot]={"shms_replay_matrixopt_1703_newmatrix","shms_replay_matrixopt_1778_newmatrix","shms_replay_matrixopt_1796_newmatrix","shms_replay_matrixopt_2258_newmatrix"};
  TString nrun[nftot]={"1703 ","1778 ","1796 ","2258 "};
const Int_t nplots=9;  
 TString hname[nplots];
   TH2F *fhist[nftot][nplots];
    for (Int_t nf=0;nf<nftot;nf++) {  
  inputroot="hist/"+fname[nf]+"_hist.root";
  fhistroot[nf] =  new TFile(inputroot);
    for (Int_t nh=0;nh<nplots;nh++) {  
      hname[nh]=Form("hxfp_yfp_foil1_xscent_yscut_%d",nh);
     fhist[nf][nh] = (TH2F*)fhistroot[nf]->Get(hname[nh]);
    }
    }
     //
     TCanvas *ctime;
      TF1 *linfit[nftot][nplots];
      TF1 *linfith[nftot][nplots];
      TCutG *cutg[nftot][nplots];
      TGraph *gr[nftot][nplots];
     ctime = new TCanvas(Form("ctime_%d",0),"fit", 700,700);
     ctime->Divide(1,1);
     gStyle->SetOptFit(111);
     ctime->cd();
    for (Int_t nh=0;nh<nplots;nh++) {  
    for (Int_t nf=0;nf<nftot;nf++) {  
     //      gPad->SetLogz();
     gPad->SetGridy();
     gPad->SetGridx();
     linfit[nf][nh] = new TF1(Form("linfit_%d_%d",nf,nh),"pol2",-20,20);
     linfith[nf][nh] = new TF1(Form("linfith_%d_%d",nf,nh),"pol2",-20,20);
     //fhist[nf][nh]->SetMinimum(3);
     fhist[nf][nh]->Draw("colz");
     TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
      gPad->Update();
     TString cutname=Form("cut_%d_%d",nf,nh);
     TString grname=Form("gr_%d_%d",nf,nh);
    if (!tempg) cout << " no cut" << endl;
    if (tempg)		  {
          cutg[nf][nh]=(TCutG*)(tempg->Clone());
      		  cutg[nf][nh]->SetName(cutname);
     		  //cutg[nf][nh]->Print();
		  Int_t npts=cutg[nf][nh]->GetN();
		  cutg[nf][nh]->RemovePoint(npts-1);
     		  //cutg[nf][nh]->Print();
     		  cutg[nf][nh]->Fit(Form("linfit_%d_%d",nf,nh),"Q");
     		  cutg[nf][nh]->Draw();
    }
    fhist[nf][nh]->Fit(Form("linfith_%d_%d",nf,nh),"Q");
     gPad->Update();
     //     gPad->WaitPrimitive();
    }
    }
    //
     TCanvas *ctime2[nplots];
    for (Int_t nh=0;nh<nplots;nh++) {  
     ctime2[nh] = new TCanvas(Form("ctime2_%d",nh),hname[nh], 700,700);
     ctime2[nh]->Divide(2,2);
     gStyle->SetOptFit(111);
    for (Int_t nf=0;nf<nftot;nf++) {  
     ctime2[nh]->cd(nf+1);
     //      gPad->SetLogz();
     gPad->SetGridy();
     gPad->SetGridx();
     //fhist[nf][nh]->SetMinimum(3);
     fhist[nf][nh]->Draw("colz");
     linfit[nf][nh]->Draw("same");
     cout << nrun[nf] << linfit[nf][nh]->GetParameter(0) << " " << linfit[nf][nh]->GetParError(0) << " " <<linfit[nf][nh]->GetParameter(1) << " " <<linfit[nf][nh]->GetParError(1) << " " <<linfit[nf][nh]->GetParameter(2) << " " <<linfit[nf][nh]->GetParError(2) <<endl;
    }
     if (nh==0) ctime2[nh]->Print(outputpdf+"(");
     if (nh!=0&&nh!=nplots-1) ctime2[nh]->Print(outputpdf);
     if (nh==nplots-1) ctime2[nh]->Print(outputpdf+")");
    }

    //
}
