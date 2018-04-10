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
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_shms_tracks(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
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
Double_t  pr_neg_tdiff[14];
    tsimc->SetBranchAddress("P.cal.pr.goodNegAdcTdcDiffTime",&pr_neg_tdiff[0]);
 Double_t  ngcer_tdiff[4];
 tsimc->SetBranchAddress("P.ngcer.goodAdcTdcDiffTime",&ngcer_tdiff[0]) ;
 TString hlab[4]={"1x","1y","2x","2y"};
 Int_t hnumpad[4]={13,13,14,21};
 Double_t  hodo_neg_tdiff[4][21];
 Double_t  hodo_pos_tdiff[4][21];
 Double_t  hodo_neg_amp[4][21];
 Double_t  hodo_pos_amp[4][21];
 Double_t  hodo_neg_tdc[4][21];
 Double_t  hodo_pos_tdc[4][21];
 Double_t  hodo_neg_adctime[4][21];
 Double_t  hodo_pos_adctime[4][21];
 Double_t  hodo_neg_mult[4][21];
 Double_t  hodo_pos_mult[4][21];
 TString temp;
 for (Int_t i=0;i<4;i++) {
   temp="P.hod."+hlab[i]+".GoodNegAdcTdcDiffTime";
 tsimc->SetBranchAddress(temp,&hodo_neg_tdiff[i]) ;
 temp="P.hod."+hlab[i]+".GoodPosAdcTdcDiffTime";
 tsimc->SetBranchAddress(temp,&hodo_pos_tdiff[i]) ;
   temp="P.hod."+hlab[i]+".GoodNegAdcPulseAmp";
 tsimc->SetBranchAddress(temp,&hodo_neg_amp[i]) ;
 temp="P.hod."+hlab[i]+".GoodPosAdcPulseAmp";
 tsimc->SetBranchAddress(temp,&hodo_pos_amp[i]) ;
   temp="P.hod."+hlab[i]+".GoodNegAdcPulseTime";
 tsimc->SetBranchAddress(temp,&hodo_neg_adctime[i]) ;
 temp="P.hod."+hlab[i]+".GoodPosAdcPulseTime";
 tsimc->SetBranchAddress(temp,&hodo_pos_adctime[i]) ;
   temp="P.hod."+hlab[i]+".GoodNegTdcTimeUnCorr";
 tsimc->SetBranchAddress(temp,&hodo_neg_tdc[i]) ;
 temp="P.hod."+hlab[i]+".GoodPosTdcTimeUnCorr";
 tsimc->SetBranchAddress(temp,&hodo_pos_tdc[i]) ;
   temp="P.hod."+hlab[i]+".GoodNegAdcMult";
 tsimc->SetBranchAddress(temp,&hodo_neg_mult[i]) ;
 temp="P.hod."+hlab[i]+".GoodPosAdcMult";
 tsimc->SetBranchAddress(temp,&hodo_pos_mult[i]) ;
 }
 Double_t  ngcer_mult[4];
 tsimc->SetBranchAddress("P.ngcer.goodAdcMult",&ngcer_mult[0]) ;
 Double_t  ngcer_hitused[4];
 tsimc->SetBranchAddress("P.ngcer.goodAdcHitUsed",&ngcer_hitused[0]) ;
Double_t  edtm_time;
   tsimc->SetBranchAddress("T.shms.pEDTM_tdcTime",&edtm_time);
 Double_t  trig1_time;
   tsimc->SetBranchAddress("T.shms.pTRIG1_tdcTime",&trig1_time);
 Double_t  trig2_time;
   tsimc->SetBranchAddress("T.shms.pTRIG2_tdcTime",&trig2_time);
 Double_t  trig3_time;
   tsimc->SetBranchAddress("T.shms.pTRIG3_tdcTime",&trig3_time);
 Double_t  el_hi_time;
   tsimc->SetBranchAddress("T.shms.pEL_HI_tdcTime",&el_hi_time);
 Double_t  el_lo_time;
   tsimc->SetBranchAddress("T.shms.pEL_LO_tdcTime",&el_lo_time);
 Double_t  el_lo_lo_time;
   tsimc->SetBranchAddress("T.shms.pEL_LO_LO_tdcTime",&el_lo_lo_time);
 Double_t  ngcer_time;
   tsimc->SetBranchAddress("T.shms.pNGCER_tdcTime",&ngcer_time);
 Double_t  prshow_time;
   tsimc->SetBranchAddress("T.shms.pPSHWR_adcPulseTime",&prshow_time);
 Double_t  stof_time;
   tsimc->SetBranchAddress("T.shms.pSTOF_tdcTime",&stof_time);
 Double_t  prhi_time;
   tsimc->SetBranchAddress("T.shms.pPRHI_tdcTime",&prhi_time);
 Double_t  prlo_time;
   tsimc->SetBranchAddress("T.shms.pPRLO_tdcTime",&prlo_time);
 Double_t  sumnpe;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
 Double_t  delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etracknorm",&etracknorm);
 Double_t  etotnorm;
   tsimc->SetBranchAddress("P.cal.etotnorm",&etotnorm);
 Double_t  W;
   tsimc->SetBranchAddress("P.kin.W",&W);
 Double_t  ytar;
   tsimc->SetBranchAddress("P.gtr.y",&ytar);
 Double_t  track_index;
   tsimc->SetBranchAddress("P.gtr.index",&track_index);
 Double_t  xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
 Double_t  yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
   // Define histograms
  Double_t cer_cut=.5;
  Double_t ep_cut=.8;
  Double_t zoff[4]={-80.,-290.,-320.,-307.};
  Double_t pscin_1x_zpos =  56.9-4.8;
  Double_t pscin_1y_zpos =  56.9+4.8;
  Double_t pscin_2x_zpos =  267.7+3.7;
  Double_t pscin_2y_zpos =  267.7+14.7;
  const Int_t num1x=13;
  Double_t pscin_1x_center[num1x]= {-45.0,-37.5,30.0,22.5,15.0,7.5,0.0,7.5,15.0,22.5,30.0,37.5,45.0};
  //
   TH1F *hetracknorm = new TH1F("hetracknorm", Form(" Run %d ; Track Cal E/p ; Counts",nrun), 200, 0.0,1.5);
   TH1F *hngcer_sum = new TH1F("hngcer_sum", Form(" Run %d ; NG Cer Sum Npe ; Counts",nrun), 200, 0.0,50);
   TH1F *hngcer_sum_elec = new TH1F("hngcer_sum_elec", Form("NG > %3.2f , Cal E/p> %3.2f  Run %d ; NG Cer Sum Npe ; Counts",cer_cut,ep_cut,nrun), 200, 0.0,50);
   TH1F *hngcer_sum_zero = new TH1F("hngcer_sum_zero", Form("E/p==0 Run %d ; NG Cer Sum Npe ; Counts",nrun), 200, 0.0,50);
   TH1F *hetracknorm_win = new TH1F("hetracknorm_win", Form(" X(z=%5.2f) > 14.  Run %d ; Track Cal E/p ; Counts",zoff[3],nrun), 200, 0.0,1.5);
   TH1F *hetotnorm = new TH1F("hetotnorm", Form(" Run %d ; Cal E/p ; Counts",nrun), 200, 0.0,1.5);
   TH2F *hsumnpe_etracknorm = new TH2F("hsumnpe_etracknorm", Form("Run %d ; Track  E/p ; NG Sumnpe",nrun), 150, 0.0,1.5,60,0,30);
   TH2F *hxfp_xpfp_all = new TH2F("hxfp_xpfp_all", Form(" Run %d ; X_fp (cm) ; Xp_fp",nrun), 200,-40.,40.,200,-.1,.1);
   TH2F *hxfp_xpfp_zero = new TH2F("hxfp_xpfp_zero", Form("E/p==0 Run %d ; X_fp (cm) ; Xp_fp",nrun), 200,-40.,40.,200,-.1,.1);
   TH2F *hxfp_xpfp = new TH2F("hxfp_xpfp", Form("NG > %3.2f , Cal E/p < %3.2f Run %d ; X_fp (cm) ; Xp_fp",cer_cut,ep_cut,nrun), 200,-40.,40.,200,-.1,.1);
   TH2F *hxfp_xpfp_win = new TH2F("hxfp_xpfp_win", Form("NG > %3.2f , Cal E/p < %3.2f X(z=%5.2f) > 14. Run %d ; X_fp (cm) ; Xp_fp",cer_cut,ep_cut,zoff[3],nrun), 200,-40.,40.,200,-.1,.1);
  TH2F *hxfp_xpfp_elec = new TH2F("hxfp_xpfp_elec", Form("NG > %3.2f , Cal E/p> %3.2f   Run %d ; X_fp (cm) ; Xp_fp",cer_cut,ep_cut,nrun), 200,-40.,40.,200,-.1,.1);
  TH1F *htime[10];
  htime[0] = new TH1F("htrig2_time", Form("Run %d ; ELREALTime (ns) ; Counts",nrun), 200,150.,250.);
  TH1F *htime_0_elclean = new TH1F("htrig2_time_elclean", Form("Run %d ; ELREAL Time (ns) ; Counts",nrun), 200,150.,250.);
  TH1F *htime_0_prlo_zero = new TH1F("htrig2_time_prlo_zero", Form("Run %d ; ELREAL Time (ns) ; Counts",nrun), 200,150.,250.);
  TH1F *htime_0_prlo = new TH1F("htrig2_time_prlo", Form("Run %d ; ELREAL Time (ns) ; Counts",nrun), 200,150.,250.);
  TH1F *htime_0_ello_zero = new TH1F("htrig2_time_ello_zero", Form("Run %d ; ELREAL Time (ns) ; Counts",nrun), 200,150.,250.);
  TH1F *htime_0_ello = new TH1F("htrig2_time_ello", Form("Run %d ; ELREAL Time (ns) ; Counts",nrun), 200,150.,250.);
  TH1F *htime_0_stof = new TH1F("htrig2_time_stof", Form("Run %d ; ELREAL Time (ns) ; Counts",nrun), 200,150.,250.);
  TH1F *htime_0_elec = new TH1F("htrig2_time_elec", Form("Run %d ; ELREAL Time (ns) ; Counts",nrun), 200,150.,250.);
  TH1F *htime_0_zero = new TH1F("htrig2_time_zero", Form("NG==0 , Cal E/p==0   Run %d ; Trig2 Time (ns) ; Counts",nrun), 200,150,250.);
  htime[1] = new TH1F("htrig3_time", Form("Run %d ; ELCLEAN Time (ns) ; Counts",nrun), 500,0.,250.);
  htime[2] = new TH1F("helhi_time", Form("Run %d ; Elhi Time (ns) ; Counts",nrun), 500,0.,250.);
  htime[3] = new TH1F("hello_time", Form("Run %d ; Ello Time (ns) ; Counts",nrun), 500,0.,250.);
  htime[4] = new TH1F("hellolo_time", Form("Run %d ; Ellolo Time (ns) ; Counts",nrun), 500,0.,250.);
  htime[5] = new TH1F("hngcer_time", Form("Run %d ; Ngcer Time (ns) ; Counts",nrun), 500,0.,250.);
  htime[6] = new TH1F("hstof_time", Form("Run %d ; Stof Time (ns) ; Counts",nrun), 200.,75.,175.);
  TH1F *htime_6_elclean = new TH1F("hstof_time_elclean", Form("Run %d ; STOF Time (ns) ; Counts",nrun), 200.,75.,175.);
  TH1F *htime_6_prlo_zero = new TH1F("hstof_time_prlo_zero", Form("Run %d ; STOF Time (ns) ; Counts",nrun), 200.,75.,175.);
  TH1F *htime_6_prlo = new TH1F("hstof_time_prlo", Form("Run %d ; STOF Time (ns) ; Counts",nrun), 200.,75.,175.);
  TH1F *htime_6_elec = new TH1F("hstof_time_elec", Form("Run %d ; STOF Time (ns) ; Counts",nrun), 200.,75.,175.);
  htime[7] = new TH1F("hprlo_time", Form("Run %d ; Prlo Time (ns) ; Counts",nrun), 200.,0.,200.);
  htime[8] = new TH1F("hprhi_time", Form("Run %d ; Prhi Time (ns) ; Counts",nrun), 200.,0.,200.);
  htime[9] = new TH1F("hprshow_time", Form("Run %d ; Preshower Time (ns) ; Counts",nrun), 200.,-200.,200.);
   TH2F *hxfp_yfp_zero = new TH2F("hxfp_yfp_zero", Form("NG=0,cal E/p=0 Run %d ; X_fp (cm) ; Y_fp",nrun), 200,-40.,40.,200,-40.,40.);
   TH2F *hxfp_yfp = new TH2F("hxfp_yfp", Form("NG > %3.2f , Cal E/p< %3.2f Run %d ; X_fp (cm) ; Y_fp",cer_cut,ep_cut,nrun), 200,-40.,40.,200,-40.,40.);
  TH2F *hxfp_yfp_win = new TH2F("hxfp_yfp_win", Form("NG > %3.2f , Cal E/p< %3.2f X(z=%5.2f) > 14. Run %d ; X_fp (cm) ; Y_fp",cer_cut,ep_cut,zoff[3],nrun), 200,-40.,40.,200,-40.,40.);
   TH2F *hx_y_scin1x = new TH2F("hx_y_scin1x", Form("NG > %3.2f , Cal E/p> %3.2f Run %d Scin 1X ; X (cm) ; Y (cm)",cer_cut,ep_cut,nrun), 80,-40.,40.,80,-40.,40.);
   TH2F *hx_y_scin1x_rat = new TH2F("hx_y_scin1x_rat", Form("NG > %3.2f , Cal E/p> %3.2f Run %d Scin 1X ; X (cm) ; Y (cm)",cer_cut,ep_cut,nrun), 80,-40.,40.,80,-40.,40.);
   TH2F *hx_y_scin1x_cut3 = new TH2F("hx_y_scin1x_cut3", Form("NGcer3 time cut NG > %3.2f , Cal E/p> %3.2f Run %d Scin 1X ; X (cm) ; Y (cm)",cer_cut,ep_cut,nrun), 80,-40.,40.,80,-40.,40.);;
   TH2F *hx_y_scin1x_cut4 = new TH2F("hx_y_scin1x_cut4", Form("NGcer4 time cut NG > %3.2f , Cal E/p> %3.2f Run %d Scin 1X ; X (cm) ; Y (cm)",cer_cut,ep_cut,nrun), 80,-40.,40.,80,-40.,40.);;
   TH2F *hx_y_scin1x_w = new TH2F("hx_y_scin1x_w", Form("NG > %3.2f , Cal E/p> %3.2f Run %d Scin 1X ; X (cm) ; Y (cm)",cer_cut,ep_cut,nrun), 80,-40.,40.,80,-40.,40.);
  TH2F *hxfp_yfp_elec = new TH2F("hxfp_yfp_elec", Form("NG > %3.2f , Cal E/p> %3.2f Run %d ; X_fp (cm) ; Y_fp",cer_cut,ep_cut,nrun), 200,-40.,40.,200,-40.,40.);
  TH2F *hxfp_yfp_proj[4];
  TH2F *hxfp_yfp_proj_elec[4];
  TH2F *hxfp_yfp_proj_zero[4];
	for (int i = 0; i <4; i++) {  
	  hxfp_yfp_proj[i] = new TH2F(Form("hxfp_yfp_proj_%d",i), Form("NG > %3.2f , Cal E/p< %3.2f Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",cer_cut,ep_cut,nrun,zoff[i]), 200,-40.,40.,200,-40.,40.);
	  hxfp_yfp_proj_zero[i] = new TH2F(Form("hxfp_yfp_proj_zero_%d",i), Form("Cal E/p==0 Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",nrun,zoff[i]), 200,-40.,40.,200,-40.,40.);
	  hxfp_yfp_proj_elec[i] = new TH2F(Form("hxfp_yfp_proj_elec_%d",i), Form("NG > %3.2f , Cal E/p> %3.2f Run %d Zoff = %5.2f; X_fp (cm) ; Y_fp",cer_cut,ep_cut,nrun,zoff[i]), 200,-40.,40.,200,-40.,40.);
	}
	TH1F *hhodo_neg_tdiff[4][21];
	TH1F *hhodo_pos_tdiff[4][21];
	TH1F *hhodo_tdc_tdiff[4][21];
	TH1F *hhodo_adc_tdiff[4][21];
	TH2F *hhodo_neg_tdiff_amp[4][21];
	TH2F *hhodo_pos_tdiff_amp[4][21];
	TH1F *hhodo_neg_tdiff_mult1[4][21];
	TH1F *hhodo_pos_tdiff_mult1[4][21];
	TH1F *hhodo_neg_tdiff_mult2[4][21];
	TH1F *hhodo_pos_tdiff_mult2[4][21];
	for (int ip = 0; ip <4; ip++) {
	for (int ipad = 0; ipad <hnumpad[ip]; ipad++) {
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Neg Time diff (ns) ; Amp",ipad+1);
	  hhodo_neg_tdiff_amp[ip][ipad]= new TH2F(Form("hhod_neg_tdiff_amp_%d_%d",ip,ipad),temp, 500,-50,30,500,0,100);
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Pos Time diff (ns) ; Amp (mv)",ipad+1);
	  hhodo_pos_tdiff_amp[ip][ipad]= new TH2F(Form("hhod_pos_tdiff_amp_%d_%d",ip,ipad), temp, 500,-50,30,500,0,100);
	  //
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Neg-Pos TDC Time (ns) ; Counts",ipad+1);
	  hhodo_tdc_tdiff[ip][ipad]= new TH1F(Form("hhod_tdc_tdiff_%d_%d",ip,ipad),temp, 500,-50,50);
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Neg-Pos ADC Time (ns) ; Counts",ipad+1);
	  hhodo_adc_tdiff[ip][ipad]= new TH1F(Form("hhod_adc_tdiff_%d_%d",ip,ipad),temp, 500,-50,50);
	  //
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Neg Time diff (ns) ; Counts",ipad+1);
	  hhodo_neg_tdiff[ip][ipad]= new TH1F(Form("hhod_neg_tdiff_%d_%d",ip,ipad),temp, 500,-50,30);
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Pos Time diff (ns) ; Counts",ipad+1);
	  hhodo_pos_tdiff[ip][ipad]= new TH1F(Form("hhod_pos_tdiff_%d_%d",ip,ipad), temp, 500,-50,30);
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Neg Time diff Mutl=1 (ns) ; Counts",ipad+1);
	  hhodo_neg_tdiff_mult1[ip][ipad]= new TH1F(Form("hhod_neg_tdiff_mult1_%d_%d",ip,ipad),temp, 500,-50,30);
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Pos Time diff (ns) Mult =1 ; Counts",ipad+1);
	  hhodo_pos_tdiff_mult1[ip][ipad]= new TH1F(Form("hhod_pos_tdiff_mult1_%d_%d",ip,ipad), temp, 500,-50,30);
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Neg Time diff Mutl>1 (ns) ; Counts",ipad+1);
	  hhodo_neg_tdiff_mult2[ip][ipad]= new TH1F(Form("hhod_neg_tdiff_mult2_%d_%d",ip,ipad),temp, 500,-50,30);
	  temp=Form("Run %d ; Plane ",nrun)+hlab[ip]+Form(" %d Pos Time diff (ns) Mult >1 ; Counts",ipad+1);
	  hhodo_pos_tdiff_mult2[ip][ipad]= new TH1F(Form("hhod_pos_tdif2_mult1_%d_%d",ip,ipad), temp, 500,-50,30);
	}
	}
	TH1F *hpr_neg_tdiff[14];
	for (int i = 0; i <14; i++) {
	  hpr_neg_tdiff[i]= new TH1F(Form("hpr_neg_tdiff_%d",i), Form("Run %d ;Time diff (ns) ; Coutns",nrun), 500,-30,30);
	}
	TH1F *hngcer_tdiff[4];
	TH1F *hngcer_tdiff_mult1[4];
	TH1F *hngcer_tdiff_mult2[4];
	TH1F *hngcer_tdiff_elec[4];
	TH1F *hngcer_tdiff_zero[4];
	for (int i = 0; i <4; i++) {
	  hngcer_tdiff[i]= new TH1F(Form("hngcer_tdiff_%d",i), Form("Run %d ; NGCEr %d Time Diff (ns) ; Counts",nrun,i+1), 800,-200,200);
	  hngcer_tdiff_elec[i]= new TH1F(Form("hngcer_tdiff_elec_%d",i), Form("NG > %3.2f , Cal E/p> %3.2f  Run %d ; NGCEr %d Time Diff (ns) ; Counts",cer_cut,ep_cut,nrun,i+1), 800,-200,200);
	  hngcer_tdiff_zero[i]= new TH1F(Form("hngcer_tdiff_zero_%d",i), Form("Run %d ; NGCEr %d Time Diff (ns) ; Counts",nrun,i+1), 800,-200,200);
	  hngcer_tdiff_mult1[i]= new TH1F(Form("hngcer_tdiff_mult1_%d",i), Form("NG > %3.2f , Cal E/p> %3.2f Run %d ; NGCEr %d Time Diff (ns) ; Counts",cer_cut,ep_cut,nrun,i+1), 800,-200,200);
	  hngcer_tdiff_mult2[i]= new TH1F(Form("hngcer_tdiff_mult2_%d",i), Form("NG > %3.2f , Cal E/p> %3.2f Run %d ; NGCEr %d Time Diff (ns) ; Counts",cer_cut,ep_cut,nrun,i+1), 800,-200,200);
	}
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (edtm_time==0) {
	for (int ii = 0; ii <14; ii++) {
	  hpr_neg_tdiff[ii]->Fill(pr_neg_tdiff[ii]);
	}
	for (int ii = 0; ii <4; ii++) {
	  hngcer_tdiff[ii]->Fill(ngcer_tdiff[ii]);
	}
	for (int ip = 0; ip <4; ip++) {
	for (int ipad = 0; ipad <hnumpad[ip]; ipad++) {
	  //	    hhodo_neg_tdiff[ip][ipad]->Fill(hodo_neg_tdiff[ip][ipad]);
	  //   hhodo_pos_tdiff[ip][ipad]->Fill(hodo_pos_tdiff[ip][ipad]);
            if (hodo_neg_mult[ip][ipad]==1&&hodo_pos_mult[ip][ipad]==1&&hodo_neg_tdc[ip][ipad]<10000&&hodo_pos_tdc[ip][ipad]<10000) {
	 	    hhodo_neg_tdiff[ip][ipad]->Fill(hodo_neg_tdiff[ip][ipad]);
	        hhodo_pos_tdiff[ip][ipad]->Fill(hodo_pos_tdiff[ip][ipad]);
	      hhodo_tdc_tdiff[ip][ipad]->Fill(hodo_neg_tdc[ip][ipad]-hodo_pos_tdc[ip][ipad]);
	    if (hodo_neg_mult[ip][ipad]==1)hhodo_neg_tdiff_amp[ip][ipad]->Fill(hodo_neg_tdiff[ip][ipad],hodo_neg_amp[ip][ipad]);
	    if (hodo_pos_mult[ip][ipad]==1) hhodo_pos_tdiff_amp[ip][ipad]->Fill(hodo_pos_tdiff[ip][ipad],hodo_pos_amp[ip][ipad]);
	    }
            if (hodo_neg_mult[ip][ipad]==1&&hodo_pos_mult[ip][ipad]==1&&hodo_neg_adctime[ip][ipad]<10000&&hodo_pos_adctime[ip][ipad]<10000) {
	      hhodo_adc_tdiff[ip][ipad]->Fill(hodo_neg_adctime[ip][ipad]-hodo_pos_adctime[ip][ipad]);
	    }
	    if (hodo_neg_mult[ip][ipad]==1)hhodo_neg_tdiff_mult1[ip][ipad]->Fill(hodo_neg_tdiff[ip][ipad]);
	    if (hodo_neg_mult[ip][ipad]>1)hhodo_neg_tdiff_mult2[ip][ipad]->Fill(hodo_neg_tdiff[ip][ipad]);
	    if (hodo_pos_mult[ip][ipad]==1)hhodo_pos_tdiff_mult1[ip][ipad]->Fill(hodo_pos_tdiff[ip][ipad]);
	    if (hodo_pos_mult[ip][ipad]>1)hhodo_pos_tdiff_mult2[ip][ipad]->Fill(hodo_pos_tdiff[ip][ipad]);
	}
	}
	hngcer_sum->Fill(sumnpe);
                     hxfp_xpfp_all->Fill(xfp,xpfp);
		    htime[0]->Fill(trig2_time);		  
		    if(trig3_time!=0) htime_0_elclean->Fill(trig2_time);		  
		    if(prlo_time==0) htime_0_prlo_zero->Fill(trig2_time);		  
		    if(prlo_time>0) htime_0_prlo->Fill(trig2_time);		  
		    if(el_lo_time==0) htime_0_ello_zero->Fill(trig2_time);		  
		    if(el_lo_time>0) htime_0_ello->Fill(trig2_time);		  
		    if(stof_time==0) htime_0_stof->Fill(trig2_time);		  
		    if(trig3_time!=0) htime_6_elclean->Fill(stof_time);		  
		    if(prlo_time==0) htime_6_prlo_zero->Fill(stof_time);		  
		    if(prlo_time>0) htime_6_prlo->Fill(stof_time);		  
		    htime[1]->Fill(trig3_time);		  
		    htime[2]->Fill(el_hi_time);		  
		    htime[3]->Fill(el_lo_time);		  
		    htime[4]->Fill(el_lo_lo_time);		  
		    htime[5]->Fill(ngcer_time);		  
		    htime[6]->Fill(stof_time);		  
		    htime[7]->Fill(prlo_time);		  
		    htime[8]->Fill(prhi_time);		  
		    htime[9]->Fill(prshow_time);		  
		  if  (etotnorm==0) {
                     hxfp_yfp_zero->Fill(xfp,yfp);
		    htime_0_zero->Fill(trig2_time);
                     hxfp_xpfp_zero->Fill(xfp,xpfp);
	            hngcer_sum_zero->Fill(sumnpe);
		  }
		  if ( track_index>-1) {
	        	hsumnpe_etracknorm->Fill(etracknorm,sumnpe);
		  hetracknorm->Fill(etracknorm);
		  if ((xfp+xpfp*zoff[3]) > 14) hetracknorm_win->Fill(etracknorm);
		  hetotnorm->Fill(etotnorm);
		  }
                  if(etracknorm==0) {
                    for (int j = 0; j <4; j++) {  
		       hxfp_yfp_proj_zero[j]->Fill(xfp+xpfp*zoff[j],yfp+ypfp*zoff[j]);
		     }
               	     for (int ii = 0; ii <4; ii++) {
	                 hngcer_tdiff_zero[ii]->Fill(ngcer_tdiff[ii]);
	                }
 		  }
		if (sumnpe > cer_cut && etracknorm!=0) {
                   if  (etracknorm > ep_cut ) {
                	for (int ii = 0; ii <4; ii++) {
	                 hngcer_tdiff_elec[ii]->Fill(ngcer_tdiff[ii]);
	                 if (ngcer_mult[ii]==1) hngcer_tdiff_mult1[ii]->Fill(ngcer_tdiff[ii]);
	                 if (ngcer_mult[ii]>1) hngcer_tdiff_mult2[ii]->Fill(ngcer_tdiff[ii]);
	                }
	            hngcer_sum_elec->Fill(sumnpe);
                    hx_y_scin1x->Fill(xfp+xpfp*pscin_1x_zpos,yfp+ypfp*pscin_1x_zpos);
		                        if (TMath::Abs(ngcer_tdiff[2]+14.5)<2.5) hx_y_scin1x_cut3->Fill(xfp+xpfp*pscin_1x_zpos,yfp+ypfp*pscin_1x_zpos);
		    //                     if (TMath::Abs(ngcer_tdiff[2]+10.5)<2.5) hx_y_scin1x_cut3->Fill(xfp+xpfp*pscin_1x_zpos,yfp+ypfp*pscin_1x_zpos);
                    if (TMath::Abs(ngcer_tdiff[3]+14.5)<2.5) hx_y_scin1x_cut4->Fill(xfp+xpfp*pscin_1x_zpos,yfp+ypfp*pscin_1x_zpos);
                    hx_y_scin1x_w->Fill(xfp+xpfp*pscin_1x_zpos,yfp+ypfp*pscin_1x_zpos,trig2_time);
		    htime_6_elec->Fill(stof_time);
		    htime_0_elec->Fill(trig2_time);
                     hxfp_xpfp_elec->Fill(xfp,xpfp);
                     hxfp_yfp_elec->Fill(xfp,yfp);
                     for (int j = 0; j <4; j++) {  
		       hxfp_yfp_proj_elec[j]->Fill(xfp+xpfp*zoff[j],yfp+ypfp*zoff[j]);
		     }
		  } else {
                    hxfp_xpfp->Fill(xfp,xpfp);
                    if ((xfp+xpfp*zoff[3]) > 14) hxfp_xpfp_win->Fill(xfp,xpfp);
                    if ((xfp+xpfp*zoff[3]) > 14) hxfp_yfp_win->Fill(xfp,yfp);
                     hxfp_yfp->Fill(xfp,yfp);
                     for (int j = 0; j <4; j++) {  
		       hxfp_yfp_proj[j]->Fill(xfp+xpfp*zoff[j],yfp+ypfp*zoff[j]);
		     }
		  }
		}
		}
	}
// Define exit window at dipole exit ( shape is rectangle with semi-circle sides)
    const int numcirpt=12;
    Double_t angstep=180./numcirpt;
    Double_t ang;
    const int numpt=2*numcirpt+5;
    int ctpt=0;
    Double_t test_offset = 2.;
    Double_t crad=23.81; // radius of semicircle
    Double_t voffset= crad-26.035+test_offset;
    Double_t hwid=11.549/2.;
    Double_t xwin[numpt]; // vertical
    Double_t ywin[numpt]; // horizontal
    xwin[ctpt]= crad+voffset ;
    ywin[ctpt]= hwid;
    xwin[numpt-1]= crad+voffset ;
    ywin[numpt-1]= hwid;
    Double_t xcent=hwid;
    for (int i=0 ; i < numcirpt ; i++) {
    ctpt++;
    ang = (90.-(i+1)*angstep)*(3.14159/180.);
    xwin[ctpt]= crad*TMath::Sin(ang)+voffset ;
    ywin[ctpt]= hwid+crad*TMath::Cos(ang);
    }
    ctpt++;
    xwin[ctpt]= -crad+voffset ;
    ywin[ctpt]= hwid;
    ctpt++;
    xwin[ctpt]= -crad+voffset ;
    ywin[ctpt]= -hwid;
     for (int i=0 ; i < numcirpt ; i++) {
    ctpt++;
    ang = (-90.+(i+1)*angstep)*(3.14159/180.);
    xwin[ctpt]= crad*TMath::Sin(ang)+voffset ;
    ywin[ctpt]= -hwid-crad*TMath::Cos(ang);
    }
   ctpt++;
    xwin[ctpt]= crad+voffset ;
    ywin[ctpt]= -hwid;
    for (Int_t j=0;j<numpt;j++) {
      //      cout << xwin[j] << " " << ywin[j] <<endl;
    }
    TPolyLine *exitwindow = new TPolyLine(numpt,xwin,ywin);
    exitwindow->SetFillColor(0);
    exitwindow->SetLineColor(kRed);
    exitwindow->SetLineWidth(3);
	// plot data
    if (1==-1) {
TCanvas *c1 = new TCanvas("c1", "Etrack norm", 700,700);
c1->Divide(1,2);
c1->cd(1);
 gPad->SetLogz();
 hsumnpe_etracknorm->Draw("colz");
c1->cd(2);
 gPad->SetLogy();
 hetracknorm->Draw();
 hetracknorm_win->Draw("same");
 hetracknorm_win->SetLineColor(2);
    outputpdf="plots/"+basename+"_shms_track.pdf(";
    c1->Print(outputpdf);

TCanvas *c2 = new TCanvas("c2", "Xfp v Xpfp", 700,700);
c2->Divide(2,3);
c2->cd(1);
 gPad->SetLogz();
 hxfp_xpfp->Draw("colz");
c2->cd(2);
 gPad->SetLogz();
 hxfp_xpfp_elec->Draw("colz");
c2->cd(3);
 gPad->SetLogz();
 hxfp_xpfp_win->Draw("colz");
c2->cd(4);
 gPad->SetLogz();
 hxfp_xpfp_all->Draw("colz");
c2->cd(5);
 gPad->SetLogz();
 hxfp_xpfp_zero->Draw("colz");
    outputpdf="plots/"+basename+"_shms_track.pdf";
    c2->Print(outputpdf);
    //
TCanvas *c1x = new TCanvas("c1x", "ELREAL X v Y", 700,700);
c1x->Divide(1,2);
c1x->cd(1);
 hx_y_scin1x->Draw("colz");
c1x->cd(2);
 hx_y_scin1x_rat->Divide(hx_y_scin1x_w,hx_y_scin1x,1,1,"B");
 hx_y_scin1x_rat->Draw("colz");
 hx_y_scin1x_rat->SetMinimum(150);
 hx_y_scin1x_rat->SetMaximum(250);
    outputpdf="plots/"+basename+"_shms_track.pdf";
    c1x->Print(outputpdf);
   //
    TCanvas *ctime[10];
    TLegend *ltime[10];
    TString clab[10]={"ELREAL","ELCLEAN","ELHI","ELLO","EL-LO-LO","NGCER","STOF","PRLO","PRHI","PRESHOW"};
 for (int j = 0; j <10; j++) {  
   ctime[j] = new TCanvas(Form("ctime_%d",j),clab[j], 700,700);
   ctime[j]->Divide(1,1);
   ctime[j]->cd(1);
 gPad->SetLogy();
   htime[j]->Draw();
   ltime[j] = new TLegend(.15,.75,.35,.95,"");
   if (j==0) {
     htime_0_elclean->Draw("same");
     htime_0_elclean->SetLineColor(2);
     htime_0_elclean->SetLineWidth(2);
     htime_0_prlo_zero->Draw("same");
     htime_0_prlo_zero->SetLineColor(3);
     htime_0_prlo->Draw("same");
     htime_0_prlo->SetLineColor(6);
     htime_0_elec->Draw("same");
     htime_0_elec->SetLineColor(7);
     //     htime_0_zero->Draw("same");
     //htime_0_zero->SetLineColor(8);
     ltime[j]->AddEntry(htime_0_elclean," ELCLEAN time > 0");
     ltime[j]->AddEntry(htime_0_prlo_zero," PR_LO time = 0");
     ltime[j]->AddEntry(htime_0_prlo," PR_LO time > 0");
     ltime[j]->AddEntry(htime_0_elec," Electron");
     //ltime[j]->AddEntry(htime_0_zero," NG=0 E/p=0");
     ltime[j]->Draw();
   }
   if (j==6) {
     htime_6_elclean->Draw("same");
     htime_6_elclean->SetLineColor(2);
     htime_6_elclean->SetLineWidth(2);
     htime_6_prlo_zero->Draw("same");
     htime_6_prlo_zero->SetLineColor(3);
     htime_6_prlo->Draw("same");
     htime_6_prlo->SetLineColor(6);
     htime_6_elec->Draw("same");
     htime_6_elec->SetLineColor(7);
     ltime[j]->AddEntry(htime_6_elclean," ELCLEAN time > 0");
     ltime[j]->AddEntry(htime_6_prlo_zero," PRLO time = 0");
     ltime[j]->AddEntry(htime_6_prlo," PRLO time > 0");
     ltime[j]->AddEntry(htime_6_elec," Electron");
     ltime[j]->Draw();
   }
   outputpdf="plots/"+basename+"_shms_track.pdf";
     ctime[j]->Print(outputpdf);
 }
 //
    TCanvas *ctime2;
    TLegend *ltime2;
   ctime2 = new TCanvas("ctime2","ELREAL 2", 700,700);
   ctime2->Divide(1,1);
   ctime2->cd(1);
 gPad->SetLogy();
   htime[0]->Draw();
   ltime2 = new TLegend(.15,.75,.35,.95,"");
     htime_0_zero->Draw("same");
     htime_0_zero->SetLineColor(2);
     ltime2->AddEntry(htime_0_zero," Track E/p=0");
     ltime2->Draw();
   outputpdf="plots/"+basename+"_shms_track.pdf";
     ctime2->Print(outputpdf);
 //
    TCanvas *ctime3;
    TLegend *ltime3;
   ctime3 = new TCanvas("ctime3","ELREAL 2", 700,700);
   ctime3->Divide(1,1);
   ctime3->cd(1);
 gPad->SetLogy();
   htime[0]->Draw();
   ltime3 = new TLegend(.15,.75,.35,.95,"");
     htime_0_ello_zero->Draw("same");
     htime_0_ello_zero->SetLineColor(3);
     htime_0_ello->Draw("same");
     htime_0_ello->SetLineColor(6);
     ltime3->AddEntry(htime_0_ello_zero," EL_LO time = 0");
     ltime3->AddEntry(htime_0_ello," EL_LO time > 0");
     ltime3->Draw();
   outputpdf="plots/"+basename+"_shms_track.pdf";
     ctime3->Print(outputpdf);

 //
    TCanvas *cr4[4];
 for (int j = 0; j <4; j++) {  
   cr4[j] = new TCanvas(Form("cr4_%d",j), Form("Z=%5.2f Xfp v Yfp",zoff[j]), 700,700);
   cr4[j]->Divide(2,2);
   cr4[j]->cd(1);
   hxfp_yfp_proj[j]->Draw("colz");
   if (j==3) exitwindow->Draw();
   cr4[j]->cd(2);
   hxfp_yfp_proj_elec[j]->Draw("colz");
   if (j==3) exitwindow->Draw();
   cr4[j]->cd(3);
   hxfp_yfp_proj_zero[j]->Draw("colz");
   if (j==3) exitwindow->Draw();
   outputpdf="plots/"+basename+"_shms_track.pdf";
     cr4[j]->Print(outputpdf);
 }
 //
TCanvas *c3 = new TCanvas("c3", "Xfp v Yfp", 700,700);
c3->Divide(2,2);
c3->cd(1);
 hxfp_yfp->Draw("colz");
c3->cd(2);
 hxfp_yfp_win->Draw("colz");
c3->cd(3);
 hxfp_yfp_elec->Draw("colz");
    outputpdf="plots/"+basename+"_shms_track.pdf)";
    c3->Print(outputpdf);
    //
    } // if 1==-1
 //
    TCanvas *ctime4;
    TLegend *ltime4;
   ctime4 = new TCanvas("ctime4","NGCER", 700,700);
   ctime4->Divide(2,2);
   for (Int_t i=0;i<4;i++) {
   ctime4->cd(i+1);
   gPad->SetLogy();
   hngcer_tdiff[i]->Draw();
   ltime4 = new TLegend(.15,.75,.35,.95,"");
   hngcer_tdiff_elec[i]->Draw("same");
     hngcer_tdiff_elec[i]->SetLineColor(2);
     hngcer_tdiff_zero[i]->Draw("same");
     hngcer_tdiff_zero[i]->SetLineColor(3);
     ltime4->AddEntry(hngcer_tdiff_elec[i]," Electron");
     ltime4->AddEntry(hngcer_tdiff_zero[i]," E/p==0");
     ltime4->Draw();
   }
   outputpdf="plots/"+basename+"_shms_track2.pdf(";
     ctime4->Print(outputpdf);
//
    TCanvas *ctime6;
    TLegend *ltime6;
   ctime6 = new TCanvas("ctime6","NGCER3", 700,700);
   ctime6->Divide(2,2);
   for (Int_t i=0;i<4;i++) {
   ctime6->cd(i+1);
   gPad->SetLogy();
   hngcer_tdiff_elec[i]->DrawCopy();
   ltime6 = new TLegend(.15,.75,.55,.95,"");
     hngcer_tdiff_mult2[i]->Draw("same");
     hngcer_tdiff_mult2[i]->SetLineColor(2);
     hngcer_tdiff_mult1[i]->Draw("same");
     hngcer_tdiff_mult1[i]->SetLineColor(3);
     ltime6->AddEntry(hngcer_tdiff_mult2[i]," Electron ADC Mult >1");
     ltime6->AddEntry(hngcer_tdiff_mult1[i],"Electron Adc Mult=1");
     ltime6->Draw();
   }
   outputpdf="plots/"+basename+"_shms_track2.pdf";
     ctime6->Print(outputpdf);
//
//
    TCanvas *ctime5;
    TLegend *ltime5;
    TH1F *h1[4];
    TH1F *h2[4];
    TH1F *h3[4];
   ctime5 = new TCanvas("ctime5","NGCER2", 700,700);
   ctime5->Divide(2,2);
   for (Int_t i=0;i<4;i++) {
   ctime5->cd(i+1);
   gPad->SetLogy();
   h1[i]= (TH1F*)hngcer_tdiff[i]->Clone(Form("h1_%d",i));
   h2[i]= (TH1F*)hngcer_tdiff_elec[i]->Clone(Form("h1_elec_%d",i));
   h3[i]= (TH1F*)hngcer_tdiff_zero[i]->Clone(Form("h1_zero_%d",i));
   h1[i]->Draw();
   h1[i]->GetXaxis()->SetRangeUser(-30.,0.);
   ltime5 = new TLegend(.15,.75,.35,.95,"");
     h2[i]->Draw("same");
     h2[i]->SetLineColor(2);
     h2[i]->GetXaxis()->SetRangeUser(-30.,0.);
     h3[i]->Draw("same");
     h3[i]->SetLineColor(3);
     h3[i]->GetXaxis()->SetRangeUser(-30.,0.);
     ltime5->AddEntry(h2[i]," Electron");
     ltime5->AddEntry(h3[i]," E/p==0");
     ltime5->Draw();
   }
   outputpdf="plots/"+basename+"_shms_track2.pdf)";
     ctime5->Print(outputpdf);
//
    TCanvas *ctime7[4];
    TLegend *ltime7[4];
    TCanvas *ctime8[4];
    TLegend *ltime8[4];
    for (Int_t ip=0;ip<4;ip++) {
      temp="HODO NEG "+hlab[ip];
      ctime7[ip] = new TCanvas(Form("ctime7_%d",ip),temp, 700,700);
   ctime7[ip]->Divide(2,7);
   Int_t thi=hnumpad[ip];
   Int_t tlo=0;
   if (ip==3) thi=16;
   if (ip==3) tlo=2;
   for (Int_t i=tlo;i<thi;i++) {
   ctime7[ip]->cd(i-tlo+1);
   gPad->SetLogy();
   hhodo_neg_tdiff[ip][i]->DrawCopy();
   ltime7[ip] = new TLegend(.35,.75,.55,.95,"");
     hhodo_neg_tdiff_mult2[ip][i]->Draw("same");
     hhodo_neg_tdiff_mult2[ip][i]->SetLineColor(2);
     hhodo_neg_tdiff_mult1[ip][i]->Draw("same");
     hhodo_neg_tdiff_mult1[ip][i]->SetLineColor(3);
     ltime7[ip]->AddEntry(hhodo_neg_tdiff_mult2[ip][i],"Mult >1");
     ltime7[ip]->AddEntry(hhodo_neg_tdiff_mult1[ip][i],"Mult=1");
     ltime7[ip]->Draw();
   }
   if (ip==0)outputpdf="plots/"+basename+"_shms_track3.pdf(";
   if (ip!=0)outputpdf="plots/"+basename+"_shms_track3.pdf";
     ctime7[ip]->Print(outputpdf);
//
//
      temp="HODO POS "+hlab[ip];
     ctime8[ip] = new TCanvas(Form("ctime8_%d",ip),temp, 700,700);
   ctime8[ip]->Divide(2,7);
   for (Int_t i=tlo;i<thi;i++) {
   ctime8[ip]->cd(i-tlo+1);
   gPad->SetLogy();
   hhodo_pos_tdiff[ip][i]->DrawCopy();
   ltime8[ip] = new TLegend(.35,.75,.55,.95,"");
     hhodo_pos_tdiff_mult2[ip][i]->Draw("same");
     hhodo_pos_tdiff_mult2[ip][i]->SetLineColor(2);
     hhodo_pos_tdiff_mult1[ip][i]->Draw("same");
     hhodo_pos_tdiff_mult1[ip][i]->SetLineColor(3);
     ltime8[ip]->AddEntry(hhodo_pos_tdiff_mult2[ip][i],"Mult >1");
     ltime8[ip]->AddEntry(hhodo_pos_tdiff_mult1[ip][i],"Mult=1");
     ltime8[ip]->Draw();
   }
   if (ip==3) outputpdf="plots/"+basename+"_shms_track3.pdf)";
   if (ip!=3) outputpdf="plots/"+basename+"_shms_track3.pdf";
     ctime8[ip]->Print(outputpdf);
    }
//
    TCanvas *ctime9[4];
    TLegend *ltime9[4];
    TCanvas *ctime10[4];
    TLegend *ltime10[4];
    for (Int_t ip=0;ip<4;ip++) {
      temp="2d HODO NEG "+hlab[ip];
      ctime9[ip] = new TCanvas(Form("ctime9_%d",ip),temp, 700,700);
   ctime9[ip]->Divide(2,7);
   Int_t thi=hnumpad[ip];
   Int_t tlo=0;
   if (ip==3) thi=16;
   if (ip==3) tlo=2;
   for (Int_t i=tlo;i<thi;i++) {
   ctime9[ip]->cd(i-tlo+1);
   gPad->SetLogz();
   hhodo_neg_tdiff_amp[ip][i]->DrawCopy("colz");
   }
   if (ip==0)outputpdf="plots/"+basename+"_shms_track4.pdf(";
   if (ip!=0)outputpdf="plots/"+basename+"_shms_track4.pdf";
     ctime9[ip]->Print(outputpdf);
//
//
      temp="2d HODO POS "+hlab[ip];
     ctime10[ip] = new TCanvas(Form("ctime10_%d",ip),temp, 700,700);
   ctime10[ip]->Divide(2,7);
   for (Int_t i=tlo;i<thi;i++) {
   ctime10[ip]->cd(i-tlo+1);
   gPad->SetLogz();
   hhodo_pos_tdiff_amp[ip][i]->DrawCopy("colz");
   }
   if (ip==3) outputpdf="plots/"+basename+"_shms_track4.pdf)";
   if (ip!=3) outputpdf="plots/"+basename+"_shms_track4.pdf";
     ctime10[ip]->Print(outputpdf);
    }
 //
}
