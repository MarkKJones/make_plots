#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
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

void comp_hist_shms_hodo_cosmic(TString basename, TString basename1,Int_t plot_pm_diff=-1,Int_t plot_plane_diff=-1,Int_t plot_2hit_diff=1,Int_t plot_adc_diff=-1) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
 outputpdf="plots/"+basename+"_shms_hodo_adc.pdf";
  TString inputroot;
 static const Int_t nftot=2;
   TFile *fhistroot[nftot];
     inputroot="hist/"+basename+"_hodo_cosmic_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[0] =  new TFile(inputroot);
     inputroot="hist/"+basename1+"_hodo_cosmic_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[1] =  new TFile(inputroot);
 static const Int_t plnum=4;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 static const Int_t npad[plnum]={13,13,14,21};
 const char* tname[nftot]={"4373","4375"};
 TH1F* diffhistcent[2][plnum][21];
 TH1F* difftdccorrcent[2][plnum][21];
 TH1F* diffadchistcent[2][plnum][21];
 TH1F* negadchistcent[2][plnum][21];
 TH1F* posadchistcent[2][plnum][21];
 TH1F* diffposhistcent[2][plnum][21];
 TH1F* diffneghistcent[2][plnum][21];
 TH1F* diffpos2hitcent[2][plnum][21];
 TH1F* diffneg2hitcent[2][plnum][21];
  for (Int_t ifn=0;ifn<2;ifn++) {
  for (Int_t ip=0;ip<plnum;ip++) {
    for (Int_t ipd=0;ipd<npad[ip];ipd++) {
      TString hname = Form("histcent_%s_diff_pad_%d",plname[ip],ipd+1);
   diffhistcent[ifn][ip][ipd] = (TH1F*)fhistroot[ifn]->Get(hname);
      hname = Form("histcent_%s_tdccorrdiff_pad_%d",plname[ip],ipd+1);
   difftdccorrcent[ifn][ip][ipd] = (TH1F*)fhistroot[ifn]->Get(hname);
      hname = Form("hist_%s_adcneg_pad_%d",plname[ip],ipd+1);
   negadchistcent[ifn][ip][ipd] = (TH1F*)fhistroot[ifn]->Get(hname);
      hname = Form("hist_%s_adcpos_pad_%d",plname[ip],ipd+1);
   posadchistcent[ifn][ip][ipd] = (TH1F*)fhistroot[ifn]->Get(hname);
      hname = Form("histcent_%s_adcdiff_pad_%d",plname[ip],ipd+1);
   diffadchistcent[ifn][ip][ipd] = (TH1F*)fhistroot[ifn]->Get(hname);
      hname = Form("histcent_%s_diffpos_pad_%d",plname[ip],ipd+1);
   diffposhistcent[ifn][ip][ipd] =(TH1F*)fhistroot[ifn]->Get(hname);
      hname = Form("hdiff_%s_diffpos_2hit_pad_%d",plname[ip],ipd+1);
   diffpos2hitcent[ifn][ip][ipd] =(TH1F*)fhistroot[ifn]->Get(hname);
       hname = Form("hdiff_%s_diffneg_2hit_pad_%d",plname[ip],ipd+1);
   diffneg2hitcent[ifn][ip][ipd] =(TH1F*)fhistroot[ifn]->Get(hname);
      hname = Form("histcent_%s_diffneg_pad_%d",plname[ip],ipd+1);
   diffneghistcent[ifn][ip][ipd] = (TH1F*)fhistroot[ifn]->Get(hname);
    }
  }} 
  //
Int_t npad_lo[plnum]={0,0,0,2};
Int_t npad_hi[plnum]={13,13,14,18};
//
 if (plot_adc_diff==1) {
     TCanvas *cplotdiffpad[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplotdiffpad[nh] = new TCanvas(Form("cplotdiffpad_%d",nh),Form("adcdiff_%s",plname[nh]), 900,700);
       if (nh!=plnum-1) cplotdiffpad[nh]->Divide(5,3);
       if (nh==plnum-1) cplotdiffpad[nh]->Divide(5,4);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh];np++) {
	  cplotdiffpad[nh]->cd(cnt);
          diffadchistcent[0][nh][np]->DrawNormalized();
          diffadchistcent[1][nh][np]->SetLineColor(2);
          diffadchistcent[1][nh][np]->DrawNormalized("same");
	  //	  gPad->SetLogy();
          cnt++;
        }      
	if (nh==0)   cplotdiffpad[nh]->Print(outputpdf+"(");  
	if (nh>0&&  nh<3)   cplotdiffpad[nh]->Print(outputpdf);  
	if (nh==3)   cplotdiffpad[nh]->Print(outputpdf+")");  
     }
  //
  outputpdf="plots/"+basename+"_shms_hodo_adc_pos_neg.pdf";
    TCanvas *cplotpospad[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplotpospad[nh] = new TCanvas(Form("cplotpospad_%d",nh),Form("adc_pos_neg_%s",plname[nh]), 900,700);
       if (nh!=plnum-1) cplotpospad[nh]->Divide(5,3);
       if (nh==plnum-1) cplotpospad[nh]->Divide(5,4);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh];np++) {
	  cplotpospad[nh]->cd(cnt);
          posadchistcent[0][nh][np]->DrawNormalized();
          negadchistcent[0][nh][np]->SetLineColor(2);
          negadchistcent[0][nh][np]->DrawNormalized("same");
	  //	  gPad->SetLogy();
          cnt++;
        }      
	if (nh==0)   cplotpospad[nh]->Print(outputpdf+"(");  
	if (nh>0&&  nh<3)   cplotpospad[nh]->Print(outputpdf);  
	if (nh==3)   cplotpospad[nh]->Print(outputpdf+")");  
     }
  //
 }
 if (plot_2hit_diff==1) {
  //
      Double_t pos_2hit_mean[2][4][21];
     Double_t pos_2hit_padnum[2][4][21];
  outputpdf="plots/"+basename+"_shms_hodo_diff_2hit_pos.pdf";
    TCanvas *cplot2hit[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplot2hit[nh] = new TCanvas(Form("cplot2hit_%d",nh),Form("twohit_pos_%s",plname[nh]), 900,700);
       if (nh!=plnum-1) cplot2hit[nh]->Divide(5,3);
       if (nh==plnum-1) cplot2hit[nh]->Divide(5,4);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh]-1;np++) {
	  cplot2hit[nh]->cd(cnt);
          diffpos2hitcent[0][nh][np]->DrawNormalized();
          diffpos2hitcent[1][nh][np]->SetLineColor(2);
          diffpos2hitcent[1][nh][np]->DrawNormalized("same");
	  pos_2hit_mean[0][nh][cnt-1]=diffpos2hitcent[0][nh][np]->GetMean();
	  pos_2hit_mean[1][nh][cnt-1]=diffpos2hitcent[1][nh][np]->GetMean();
	  pos_2hit_padnum[0][nh][cnt-1]=np+1;
	  pos_2hit_padnum[1][nh][cnt-1]=np+1;
	  //	  gPad->SetLogy();
          cnt++;
        }      
	if (nh==0)   cplot2hit[nh]->Print(outputpdf+"(");  
	if (nh>0&&  nh<3)   cplot2hit[nh]->Print(outputpdf);  
	if (nh==3)   cplot2hit[nh]->Print(outputpdf+")");  
     }
  //
  //
  outputpdf="plots/"+basename+"_shms_hodo_diff_2hit_pos_mean.pdf";
     TGraph *gr_pos_2hit_mean[2][plnum];
     TMultiGraph *mgr_pos_2hit[plnum];
     TCanvas *cplot_pos_2hit_mean;
       cplot_pos_2hit_mean = new TCanvas("cplot_pos_2hit_mean","pos_2hit_mean", 900,700);
       cplot_pos_2hit_mean->Divide(2,2);
     for (Int_t nh=0;nh<plnum;nh++) {
             mgr_pos_2hit[nh] = new TMultiGraph();
       cplot_pos_2hit_mean->cd(nh+1);
       Int_t totp = npad_hi[nh]-npad_lo[nh]-1;
       gr_pos_2hit_mean[0][nh] = new TGraph(totp,pos_2hit_padnum[0][nh],pos_2hit_mean[0][nh]);
       gr_pos_2hit_mean[1][nh] = new TGraph(totp,pos_2hit_padnum[1][nh],pos_2hit_mean[1][nh]);
       gr_pos_2hit_mean[0][nh]->SetMarkerColor(1);
        gr_pos_2hit_mean[1][nh]->SetMarkerColor(2);
       gr_pos_2hit_mean[0][nh]->SetMarkerStyle(22);
     gr_pos_2hit_mean[1][nh]->SetMarkerStyle(21);
      mgr_pos_2hit[nh]->Add(gr_pos_2hit_mean[0][nh]);
       mgr_pos_2hit[nh]->Add(gr_pos_2hit_mean[1][nh]);
       mgr_pos_2hit[nh]->SetTitle(Form("Plane %s  ; Paddle Number ; Mean of Time Diff Pos PMT Adj. pads",plname[nh]));
       mgr_pos_2hit[nh]->Draw("AP");
       cplot_pos_2hit_mean->Print(outputpdf);
     }
 //
  //
      Double_t neg_2hit_mean[2][4][21];
     Double_t neg_2hit_padnum[2][4][21];
  outputpdf="plots/"+basename+"_shms_hodo_diff_2hit_neg.pdf";
    TCanvas *cplot2hitneg[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplot2hitneg[nh] = new TCanvas(Form("cplot2hitneg_%d",nh),Form("twohit_neg_%s",plname[nh]), 900,700);
       if (nh!=plnum-1) cplot2hitneg[nh]->Divide(5,3);
       if (nh==plnum-1) cplot2hitneg[nh]->Divide(5,4);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh];np++) {
	  cplot2hitneg[nh]->cd(cnt);
          diffneg2hitcent[0][nh][np]->DrawNormalized();
          diffneg2hitcent[1][nh][np]->SetLineColor(2);
          diffneg2hitcent[1][nh][np]->DrawNormalized("same");
	  neg_2hit_mean[0][nh][cnt-1]=diffneg2hitcent[0][nh][np]->GetMean();
	  neg_2hit_mean[1][nh][cnt-1]=diffneg2hitcent[1][nh][np]->GetMean();
	  neg_2hit_padnum[0][nh][cnt-1]=np+1;
	  neg_2hit_padnum[1][nh][cnt-1]=np+1;
	  //	  gPad->SetLogy();
          cnt++;
        }      
	if (nh==0)   cplot2hitneg[nh]->Print(outputpdf+"(");  
	if (nh>0&&  nh<3)   cplot2hitneg[nh]->Print(outputpdf);  
	if (nh==3)   cplot2hitneg[nh]->Print(outputpdf+")");  
     }
  //
  outputpdf="plots/"+basename+"_shms_hodo_diff_2hit_neg_mean.pdf";
     TGraph *gr_neg_2hit_mean[2][plnum];
     TMultiGraph *mgr_neg_2hit[plnum];
     TCanvas *cplot_neg_2hit_mean;
       cplot_neg_2hit_mean = new TCanvas("cplot_neg_2hit_mean","neg_2hit_mean", 900,700);
       cplot_neg_2hit_mean->Divide(2,2);
     for (Int_t nh=0;nh<plnum;nh++) {
             mgr_neg_2hit[nh] = new TMultiGraph();
       cplot_neg_2hit_mean->cd(nh+1);
       Int_t totp = npad_hi[nh]-npad_lo[nh]-1;
       gr_neg_2hit_mean[0][nh] = new TGraph(totp,neg_2hit_padnum[0][nh],neg_2hit_mean[0][nh]);
       gr_neg_2hit_mean[1][nh] = new TGraph(totp,neg_2hit_padnum[1][nh],neg_2hit_mean[1][nh]);
       gr_neg_2hit_mean[0][nh]->SetMarkerColor(1);
        gr_neg_2hit_mean[1][nh]->SetMarkerColor(2);
       gr_neg_2hit_mean[0][nh]->SetMarkerStyle(22);
     gr_neg_2hit_mean[1][nh]->SetMarkerStyle(21);
      mgr_neg_2hit[nh]->Add(gr_neg_2hit_mean[0][nh]);
       mgr_neg_2hit[nh]->Add(gr_neg_2hit_mean[1][nh]);
       mgr_neg_2hit[nh]->SetTitle(Form("Plane %s  ; Paddle Number ; Mean of Time Diff Neg PMT Adj. pads",plname[nh]));
       mgr_neg_2hit[nh]->Draw("AP");
       cplot_neg_2hit_mean->Print(outputpdf);
     }
 }
  //
     if ( plot_plane_diff==1) {
       //
  //
  outputpdf="plots/"+basename+"_shms_hodo_diff_pos.pdf";
     Double_t pl_pos_mean[2][4][21];
     Double_t pl_pos_padnum[2][4][21];
     TCanvas *cplotposdiff[plnum];
    for (Int_t nh=0;nh<plnum;nh++) {
       cplotposdiff[nh] = new TCanvas(Form("cplotposdiff_%d",nh),Form("diff_pos_%s",plname[nh]), 900,700);
       if (nh!=plnum-1) cplotposdiff[nh]->Divide(5,3);
       if (nh==plnum-1) cplotposdiff[nh]->Divide(5,4);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh];np++) {
	  cplotposdiff[nh]->cd(cnt);
          diffposhistcent[0][nh][np]->DrawNormalized();
          diffposhistcent[1][nh][np]->SetLineColor(2);
          diffposhistcent[1][nh][np]->DrawNormalized("same");
	  pl_pos_padnum[0][nh][cnt-1]=np+1;
	  pl_pos_padnum[1][nh][cnt-1]=np+1;
	  pl_pos_mean[0][nh][cnt-1]=diffposhistcent[0][nh][np]->GetMean();
	  pl_pos_mean[1][nh][cnt-1]=diffposhistcent[1][nh][np]->GetMean();
	  //	  gPad->SetLogy();
          cnt++;
        }      
	if (nh==0)   cplotposdiff[nh]->Print(outputpdf+"(");  
	if (nh>0&&  nh<3)   cplotposdiff[nh]->Print(outputpdf);  
	if (nh==3)   cplotposdiff[nh]->Print(outputpdf+")");  
     }
  //
  //
  outputpdf="plots/"+basename+"_shms_hodo_diff_neg.pdf";
     Double_t pl_neg_mean[2][4][21];
     Double_t pl_neg_padnum[2][4][21];
     TCanvas *cplotnegdiff[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplotnegdiff[nh] = new TCanvas(Form("cplotnegdiff_%d",nh),Form("diff_neg_%s",plname[nh]), 900,700);
       if (nh!=plnum-1) cplotnegdiff[nh]->Divide(5,3);
       if (nh==plnum-1) cplotnegdiff[nh]->Divide(5,4);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh];np++) {
	  cplotnegdiff[nh]->cd(cnt);
          diffneghistcent[0][nh][np]->DrawNormalized();
          diffneghistcent[1][nh][np]->SetLineColor(2);
          diffneghistcent[1][nh][np]->DrawNormalized("same");
	  pl_neg_padnum[0][nh][cnt-1]=np+1;
	  pl_neg_padnum[1][nh][cnt-1]=np+1;
	  pl_neg_mean[0][nh][cnt-1]=diffneghistcent[0][nh][np]->GetMean();
	  pl_neg_mean[1][nh][cnt-1]=diffneghistcent[1][nh][np]->GetMean();
	  //	  gPad->SetLogy();
          cnt++;
        }      
	if (nh==0)   cplotnegdiff[nh]->Print(outputpdf+"(");  
	if (nh>0&&  nh<3)   cplotnegdiff[nh]->Print(outputpdf);  
	if (nh==3)   cplotnegdiff[nh]->Print(outputpdf+")");  
     }
     // 
     const char*  cross_pl_pad[plnum]= {"S1YP7","S1XP7","S2YP11","S2XP7-8"};
     outputpdf="plots/"+basename+"_shms_hodo_diff_pl_pos_mean.pdf";
     TGraph *gr_pl_pos_mean[2][plnum];
     TMultiGraph *mgr[plnum];
     TCanvas *cplotpl_posmean;
       cplotpl_posmean = new TCanvas("cplotpl_posmean","diff_pl_pos_mean", 900,700);
       cplotpl_posmean->Divide(2,2);
     for (Int_t nh=0;nh<plnum;nh++) {
             mgr[nh] = new TMultiGraph();
       cplotpl_posmean->cd(nh+1);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       gr_pl_pos_mean[0][nh] = new TGraph(totp,pl_pos_padnum[0][nh],pl_pos_mean[0][nh]);
       gr_pl_pos_mean[1][nh] = new TGraph(totp,pl_pos_padnum[1][nh],pl_pos_mean[1][nh]);
       gr_pl_pos_mean[0][nh]->SetMarkerColor(1);
        gr_pl_pos_mean[1][nh]->SetMarkerColor(2);
       gr_pl_pos_mean[0][nh]->SetMarkerStyle(22);
     gr_pl_pos_mean[1][nh]->SetMarkerStyle(21);
      mgr[nh]->Add(gr_pl_pos_mean[0][nh]);
       mgr[nh]->Add(gr_pl_pos_mean[1][nh]);
       mgr[nh]->SetTitle(Form("Plane %s  ; Paddle Number ; Mean of Pos - Pos(%s) Time",plname[nh],cross_pl_pad[nh]));
       mgr[nh]->Draw("AP");
     }
       cplotpl_posmean->Print(outputpdf);
     //
    outputpdf="plots/"+basename+"_shms_hodo_diff_pl_neg_mean.pdf";
     TGraph *gr_pl_neg_mean[2][plnum];
     TMultiGraph *mgr_pl_neg[plnum];
     TCanvas *cplotpl_negmean;
       cplotpl_negmean = new TCanvas("cplotpl_negmean","diff_pl_neg_mean", 900,700);
       cplotpl_negmean->Divide(2,2);
     for (Int_t nh=0;nh<plnum;nh++) {
             mgr_pl_neg[nh] = new TMultiGraph();
       cplotpl_negmean->cd(nh+1);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       gr_pl_neg_mean[0][nh] = new TGraph(totp,pl_neg_padnum[0][nh],pl_neg_mean[0][nh]);
       gr_pl_neg_mean[1][nh] = new TGraph(totp,pl_neg_padnum[1][nh],pl_neg_mean[1][nh]);
       gr_pl_neg_mean[0][nh]->SetMarkerColor(1);
        gr_pl_neg_mean[1][nh]->SetMarkerColor(2);
       gr_pl_neg_mean[0][nh]->SetMarkerStyle(22);
     gr_pl_neg_mean[1][nh]->SetMarkerStyle(21);
      mgr_pl_neg[nh]->Add(gr_pl_neg_mean[0][nh]);
       mgr_pl_neg[nh]->Add(gr_pl_neg_mean[1][nh]);
       mgr_pl_neg[nh]->SetTitle(Form("Plane %s  ; Paddle Number ; Mean of Neg - Pos(%s) Time",plname[nh],cross_pl_pad[nh]));
       mgr_pl_neg[nh]->Draw("AP");
     }
       cplotpl_negmean->Print(outputpdf);
     //
   //
     }
     // Comparing pos - neg on one paddle
     if (plot_pm_diff==1) {
  //
       if ( 1==1) {
     outputpdf="plots/"+basename+"_shms_hodo_compdiff_pm.pdf";
     TCanvas *cplotcomppmdiff[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplotcomppmdiff[nh] = new TCanvas(Form("cplotcomppmdiff_%d",nh),Form("diff_comppm_%s",plname[nh]), 900,700);
       if (nh!=plnum-1) cplotcomppmdiff[nh]->Divide(5,3);
       if (nh==plnum-1) cplotcomppmdiff[nh]->Divide(5,4);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh];np++) {
	  cplotcomppmdiff[nh]->cd(cnt);
          diffhistcent[0][nh][np]->DrawNormalized();
          difftdccorrcent[0][nh][np]->SetLineColor(2);
          difftdccorrcent[0][nh][np]->DrawNormalized("same");
	  //	  gPad->SetLogy();
          cnt++;
        }      
	if (nh==0)   cplotcomppmdiff[nh]->Print(outputpdf+"(");  
	if (nh>0&&  nh<3)   cplotcomppmdiff[nh]->Print(outputpdf);  
	if (nh==3)   cplotcomppmdiff[nh]->Print(outputpdf+")");  
     }
       }
  //
     Double_t pm_mean[2][4][21];
     Double_t pm_mean_out[4][21];
     Double_t pm_padnum[2][4][21];
     for (Int_t nh=0;nh<plnum;nh++) {
      for (Int_t np=0;np<21;np++) {
	pm_mean_out[nh][np]=0.0;
      }}
 outputpdf="plots/"+basename+"_shms_hodo_diff_pm.pdf";
     TCanvas *cplotpmdiff[plnum];
     for (Int_t nh=0;nh<plnum;nh++) {
       cplotpmdiff[nh] = new TCanvas(Form("cplotpmdiff_%d",nh),Form("diff_pm_%s",plname[nh]), 900,700);
       if (nh!=plnum-1) cplotpmdiff[nh]->Divide(5,3);
       if (nh==plnum-1) cplotpmdiff[nh]->Divide(5,4);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       Int_t cnt=1;
        for (Int_t np=npad_lo[nh];np<totp+npad_lo[nh];np++) {
	  cplotpmdiff[nh]->cd(cnt);
          diffhistcent[0][nh][np]->DrawNormalized();
	  pm_padnum[0][nh][cnt-1]=np+1;
	  pm_padnum[1][nh][cnt-1]=np+1;
          diffhistcent[1][nh][np]->SetLineColor(2);
          diffhistcent[1][nh][np]->DrawNormalized("same");
	  pm_mean[0][nh][cnt-1]=diffhistcent[0][nh][np]->GetMean();
	  pm_mean_out[nh][np]=diffhistcent[0][nh][np]->GetMean();
	  pm_mean[1][nh][cnt-1]=diffhistcent[1][nh][np]->GetMean();
	  //	  gPad->SetLogy();
          cnt++;
        }      
	if (nh==0)   cplotpmdiff[nh]->Print(outputpdf+"(");  
	if (nh>0&&  nh<3)   cplotpmdiff[nh]->Print(outputpdf);  
	if (nh==3)   cplotpmdiff[nh]->Print(outputpdf+")");  
     }
      for (Int_t np=0;np<21;np++) {
         for (Int_t nh=0;nh<plnum;nh++) {
	   if (nh<plnum-1) {
	     cout <<  -pm_mean_out[nh][np]/2. << ",";
	   } else {
	     cout <<  -pm_mean_out[nh][np]/2. ;
	   }
         }
	 cout << endl;
      }
  //
     TGraph *gr_pm_mean[2][plnum];
     TMultiGraph *mgr[plnum];
     TCanvas *cplotpmmean;
       cplotpmmean = new TCanvas("cplotpmmean","diff_pm_mean", 900,700);
       cplotpmmean->Divide(2,2);
     for (Int_t nh=0;nh<plnum;nh++) {
             mgr[nh] = new TMultiGraph();
       cplotpmmean->cd(nh+1);
       Int_t totp = npad_hi[nh]-npad_lo[nh];
       gr_pm_mean[0][nh] = new TGraph(totp-1,pm_padnum[0][nh],pm_mean[0][nh]);
       gr_pm_mean[1][nh] = new TGraph(totp-1,pm_padnum[1][nh],pm_mean[1][nh]);
       gr_pm_mean[0][nh]->SetMarkerColor(1);
        gr_pm_mean[1][nh]->SetMarkerColor(2);
       gr_pm_mean[0][nh]->SetMarkerStyle(22);
     gr_pm_mean[1][nh]->SetMarkerStyle(21);
      mgr[nh]->Add(gr_pm_mean[0][nh]);
       mgr[nh]->Add(gr_pm_mean[1][nh]);
       mgr[nh]->SetTitle(Form("Plane %s  ; Paddle Number ; Mean of Pos-Neg Time",plname[nh]));
       mgr[nh]->Draw("AP");
     }
     //
     }
 //
}
 
