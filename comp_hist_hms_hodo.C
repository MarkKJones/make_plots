#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
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

void comp_hist_hms_hodo(TString basename, TString basename1) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
 outputpdf="plots/"+basename+"_hms_trig.pdf";
  TString inputroot;
 static const Int_t nftot=2;
   TFile *fhistroot[nftot];
     inputroot="hist/"+basename+"_hodo_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[0] =  new TFile(inputroot);
     inputroot="hist/"+basename1+"_hodo_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[1] =  new TFile(inputroot);
 static const Int_t plnum=4;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 TH2F *neg2dhist[2][plnum];
 TH2F *pos2dhist[2][plnum];
 const char* tname[nftot]={"After","Before"};
  for (Int_t ifn=0;ifn<2;ifn++) {
  for (Int_t ip=0;ip<plnum;ip++) {
       TString hname= Form("hist_%s_neg_pad",plname[ip]);
       neg2dhist[ifn][ip] = (TH2F*)fhistroot[ifn]->Get(hname);
       neg2dhist[ifn][ip]->SetTitle(hname+" "+tname[ifn]);
       if (!neg2dhist[ifn][ip]) cout << " no hist = " << hname << endl;
       hname= Form("hist_%s_pos_pad",plname[ip]);
       pos2dhist[ifn][ip] = (TH2F*)fhistroot[ifn]->Get(hname);
       pos2dhist[ifn][ip]->SetTitle(hname+" "+tname[ifn]);
       if (!pos2dhist[ifn][ip]) cout << " no hist = " << hname << endl;
  }} 
  //
  /* canvas
    TCanvas *cplot[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplot[nh] = new TCanvas(Form("cplot_%d",nh),Form("neg_%s",plname[nh]), 700,700);
       cplot[nh]->Divide(1,2);
        for (Int_t nft=0;nft<2;nft++) {
	  cplot[nh]->cd(nft+1);
	  gPad->SetLogz();
	  neg2dhist[nft][nh]->Draw("colz"); 
	}
     }
    TCanvas *cplotpos[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplotpos[nh] = new TCanvas(Form("cplotpos_%d",nh),Form("pos_%s",plname[nh]), 700,700);
       cplotpos[nh]->Divide(1,2);
        for (Int_t nft=0;nft<2;nft++) {
	  cplotpos[nh]->cd(nft+1);
	  gPad->SetLogz();
	  pos2dhist[nft][nh]->Draw("colz"); 
	}
     }
  */
Int_t npad_lo[plnum]={0,0,0,0};
Int_t npad_hi[plnum]={16,10,16,10};
 TH1D *hn1[4][16];
 TH1D *hn2[4][16];
    TCanvas *cplotnegpad[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplotnegpad[nh] = new TCanvas(Form("cplotnegpad_%d",nh),Form("padneg_%s",plname[nh]), 900,700);
       if (nh==0||nh==2) cplotnegpad[nh]->Divide(4,4);
       if (nh==1||nh==3) cplotnegpad[nh]->Divide(5,2);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh];np++) {
	  cplotnegpad[nh]->cd(cnt);
	  gPad->SetLogy();
	  hn1[nh][cnt] = neg2dhist[0][nh]->ProjectionY(Form("hist1proj_%s_neg_pad%d",plname[nh],np),np+1,np+1);
	  hn1[nh][cnt]->Draw();
	  hn1[nh][cnt]->SetTitle(Form("Neg %s Pad%d",plname[nh],np+1));
	   hn2[nh][cnt]= neg2dhist[1][nh]->ProjectionY(Form("hist2proj_%s_neg_pad%d",plname[nh],np),np+1,np+1);
	  hn2[nh][cnt]->Draw("same");
	  hn2[nh][cnt]->SetLineColor(2);
          cnt++;
        }      
     }
  //
 TH1D *hp1[4][16];
 TH1D *hp2[4][16];
    TCanvas *cplotpospad[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplotpospad[nh] = new TCanvas(Form("cplotpospad_%d",nh),Form("padpos_%s",plname[nh]), 900,700);
       if (nh==0||nh==2) cplotpospad[nh]->Divide(4,4);
       if (nh==1||nh==3) cplotpospad[nh]->Divide(5,2);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh];np++) {
	  cplotpospad[nh]->cd(cnt);
	  gPad->SetLogy();
	  hp1[nh][cnt] = pos2dhist[0][nh]->ProjectionY(Form("hist1proj_%s_pos_pad%d",plname[nh],np),np+1,np+1);
	  hp1[nh][cnt]->Draw();
	  hp1[nh][cnt]->SetTitle(Form("Pos %s Pad%d",plname[nh],np+1));
	   hp2[nh][cnt]= pos2dhist[1][nh]->ProjectionY(Form("hist2proj_%s_pos_pad%d",plname[nh],np),np+1,np+1);
	  hp2[nh][cnt]->Draw("same");
	  hp2[nh][cnt]->SetLineColor(2);
          cnt++;
        }      
     }
  //

 //
}
 
