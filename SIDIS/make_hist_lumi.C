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

void make_hist_lumi(TString basename="",Int_t nrun=2043,Double_t mean_current=2.0){
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
   outputhist= "hist/"+basename+"_lumi_shms_jpsi_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
TTree *tscal = (TTree*) fsimc->Get("TSP");
//
 Double_t  Scal_evNumber;
   tscal->SetBranchAddress("evNumber",&Scal_evNumber);
 Double_t  Scal_BCM4B_charge;
   tscal->SetBranchAddress("P.BCM4B.scalerCharge",&Scal_BCM4B_charge);
 Double_t  Scal_BCM4B_current;
   tscal->SetBranchAddress("P.BCM4B.scalerCurrent",&Scal_BCM4B_current);
 Double_t  Scal_time;
   tscal->SetBranchAddress("P.1MHz.scalerTime",&Scal_time);
 Double_t  Scal_EDTM;
   tscal->SetBranchAddress("P.EDTM.scaler",&Scal_EDTM);
Double_t  Scal_TRIG2;
   tscal->SetBranchAddress("P.pTRIG2.scaler",&Scal_TRIG2);
 Double_t  Scal_TRIG3;
   tscal->SetBranchAddress("P.pTRIG3.scaler",&Scal_TRIG3);
 Double_t  Scal_TRIG1;
   tscal->SetBranchAddress("P.pTRIG1.scaler",&Scal_TRIG1);
 Double_t  Scal_TRIG4;
   tscal->SetBranchAddress("P.pTRIG4.scaler",&Scal_TRIG4);
 Double_t  Scal_TRIG5;
   tscal->SetBranchAddress("P.pTRIG5.scaler",&Scal_TRIG5);
 Double_t  Scal_Splane[4];
   tscal->SetBranchAddress("P.S1X.scaler",&Scal_Splane[0]);
   tscal->SetBranchAddress("P.S1Y.scaler",&Scal_Splane[1]);
   tscal->SetBranchAddress("P.S2X.scaler",&Scal_Splane[2]);
   tscal->SetBranchAddress("P.S2Y.scaler",&Scal_Splane[3]);
//loop through scalers
     Int_t nscal_reads=0;
     Int_t nscal_reads_cut=0;
     Double_t prev_read=-1;
     Double_t ave_current=0;
     Double_t ave_current_cut=0;
     Double_t charge_sum=0;
     Double_t charge_sum_cut=0;
     Double_t prev_charge=0;
     Double_t event_flag[10000];
     Double_t scal_event_number[10000];
     Double_t tot_scal_EDTM=0;
     Double_t tot_scal_cut_EDTM=0;
     Double_t prev_EDTM=0;
     Double_t tot_scal_TRIG2=0;
     Double_t tot_scal_TRIG3=0;
     Double_t prev_TRIG2=0;
     Double_t prev_TRIG3=0;
     Double_t tot_scal_cut_TRIG2=0;
     Double_t tot_scal_cut_TRIG3=0;
     Double_t tot_scal_TRIG1=0;
     Double_t tot_scal_TRIG4=0;
     Double_t prev_TRIG1=0;
     Double_t prev_TRIG4=0;
     Double_t tot_scal_cut_TRIG1=0;
     Double_t tot_scal_cut_TRIG4=0;
     Double_t tot_scal_cut_time=0;
     Double_t tot_scal_TRIG5=0;
     Double_t prev_TRIG5=0;
    Double_t tot_scal_cut_TRIG5=0;
    Double_t threshold_cut=3.;
			   //
    Double_t tot_scal_Splane[4]={0,0,0,0};
     Double_t prev_Splane[4]={0,0,0,0};
    Double_t tot_scal_cut_Splane[4]={0,0,0,0};
			   //
      Double_t tot_scal_time=0;
     Double_t prev_time=0;
     //
Long64_t scal_entries = tscal->GetEntries();
 cout << " scal ent = " << scal_entries << endl;
 Double_t nlast=float(scal_entries);
 TH1F *h_cur_entry = new TH1F("h_cur_entry","; ENtry;current",nlast,0,nlast);
Long64_t data_entries = tsimc->GetEntries();
	for (int i = 0; i < scal_entries; i++) {
      		tscal->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
          event_flag[nscal_reads] = 0;
             scal_event_number[nscal_reads] = Scal_evNumber;
          ave_current+=Scal_BCM4B_current;
	  h_cur_entry->Fill(float(i),Scal_BCM4B_current);
	  Bool_t   special_skip=kFALSE;
	  if (i ==178 && nrun==7572) special_skip=kTRUE;
	  if (TMath::Abs(Scal_BCM4B_current-mean_current) < threshold_cut && !special_skip) {
             event_flag[nscal_reads] = 1;
             ave_current_cut+=Scal_BCM4B_current;
 	     tot_scal_cut_time+=(Scal_time-prev_time);
 	     tot_scal_cut_EDTM+=(Scal_EDTM-prev_EDTM);
 	     tot_scal_cut_TRIG2+=(Scal_TRIG2-prev_TRIG2);
 	     tot_scal_cut_TRIG3+=(Scal_TRIG3-prev_TRIG3);
 	     tot_scal_cut_TRIG1+=(Scal_TRIG1-prev_TRIG1);
 	     tot_scal_cut_TRIG4+=(Scal_TRIG4-prev_TRIG4);
 	     tot_scal_cut_TRIG5+=(Scal_TRIG5-prev_TRIG5);
	     for (Int_t s=0;s<4;s++) tot_scal_cut_Splane[s]+=(Scal_Splane[s]-prev_Splane[s]);
	     //	     cout << i << " " << tot_scal_cut_Splane[0] << " " << Scal_Splane[0] << " " << prev_Splane[0] << endl;

             charge_sum_cut+=(Scal_BCM4B_charge-prev_charge);
             nscal_reads_cut++;
	  }
	  prev_charge = Scal_BCM4B_charge;
	  prev_time = Scal_time;
	  prev_EDTM = Scal_EDTM;
	  prev_TRIG2 = Scal_TRIG2;
	  prev_TRIG3 = Scal_TRIG3;
	  prev_TRIG1 = Scal_TRIG1;
	  prev_TRIG4 = Scal_TRIG4;
	  prev_TRIG5 = Scal_TRIG5;
	     for (Int_t s=0;s<4;s++) prev_Splane[s]=Scal_Splane[s];
	  // cout <<  nscal_reads <<  " " << Scal_BCM4B_current << " " << event_flag[nscal_reads] << " " << Scal_TRIG1 << endl;
          nscal_reads++;
          charge_sum=Scal_BCM4B_charge;
	  tot_scal_EDTM=Scal_EDTM;
	  tot_scal_TRIG2=Scal_TRIG2;
	  tot_scal_TRIG3=Scal_TRIG3;
	  tot_scal_TRIG1=Scal_TRIG1;
	  tot_scal_TRIG4=Scal_TRIG4;
	  tot_scal_TRIG5=Scal_TRIG5;
          tot_scal_time=Scal_time;
	  for (Int_t s=0;s<4;s++) tot_scal_Splane[s]=Scal_Splane[s];
	}
   //T.coin.pEDTM_tdcTimeRaw
	Double_t EDTM_timeraw;
	tsimc->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw",&EDTM_timeraw);
	Double_t ctime;
	tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC2",&ctime);
	Double_t gevtyp;
	tsimc->SetBranchAddress("g.evtyp",&gevtyp);
	Double_t gevnum;
	tsimc->SetBranchAddress("g.evnum",&gevnum);
	Double_t ntrack;
	tsimc->SetBranchAddress("P.dc.ntrack",&ntrack);
	Double_t goodscinhit;
	tsimc->SetBranchAddress("P.hod.goodscinhit",&goodscinhit);
	Double_t starttime;
	tsimc->SetBranchAddress("P.hod.starttime",&starttime);
 	Double_t betanotrack;
	tsimc->SetBranchAddress("P.hod.betanotrack",&betanotrack);
  Double_t ngcer_npeSum;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&ngcer_npeSum);
  Double_t aero_npeSum;
   tsimc->SetBranchAddress("P.aero.npeSum",&aero_npeSum);
  Double_t Hhgcer_npeSum;
   tsimc->SetBranchAddress("H.cer.npeSum",&Hhgcer_npeSum);
   Double_t delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
   Double_t etotnorm;
   tsimc->SetBranchAddress("P.cal.etotnorm",&etotnorm);
   Double_t Hetotnorm;
   tsimc->SetBranchAddress("H.cal.etotnorm",&Hetotnorm);
   Double_t etot;
   tsimc->SetBranchAddress("P.cal.etot",&etot);
 static const Int_t plnum=4;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 Double_t ScinTransPos[plnum];
 Double_t ScinLongPos[plnum];
 for (Int_t ip=0;ip<plnum;ip++) {
    tsimc->SetBranchAddress(Form("P.hod.%s.ScintTranversePos",plname[ip]),&ScinTransPos[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.ScintLongPos",plname[ip]),&ScinLongPos[ip]) ;
 }
	//
	TH1F* h_EDTM_CUT = new TH1F("h_EDTM"," ; EDTM (beam cut)",1000,0.,10000.);
	TH1F* h_ctime = new TH1F("h_ctime"," ; epi cointtime",220,15.,85.);
	TH1F* h_aero = new TH1F("h_aero"," ; Aero npesum",120,0.,30.);
	TH1F* h_ev1 = new TH1F("h_ev1"," ; hEv1",10,0.,10.);
	TH1F* h_ev4 = new TH1F("h_ev4"," ; hEv4",10,0.,10.);
	TH1F* h_starttime = new TH1F("h_starttime"," ; Startime (elec good scinhit)",60,0.,120.);
	TH1F* h_starttime_track = new TH1F("h_starttime_track"," ; Startime (elec good scinhit track)",60,0.,120.);
	TH1F* h_starttime_notrack = new TH1F("h_starttime_notrack"," ; Startime (elec good scinhit notrack)",60,0.,120.);
	TH1F* h_ch1_nhit_goodscin_track = new TH1F("h_ch1_nhit_goodscin_track"," ; CH1 nhits",50,0.,50.);
	TH1F* h_ch2_nhit_goodscin_track = new TH1F("h_ch2_nhit_goodscin_track"," ; CH2 nhits",50,0.,50.);
	TH1F* h_ch1_nhit_goodscin = new TH1F("h_ch1_nhit_goodscin "," ; CH1 nhits",50,0.,50.);
	TH1F* h_ch2_nhit_goodscin  = new TH1F("h_ch2_nhit_goodscin "," ; CH2 nhits",50,0.,50.);
	TH1F* h_ch1_nhit_goodscin_notrack = new TH1F("h_ch1_nhit_goodscin_notrack "," ; CH1 nhits",50,0.,50.);
	TH1F* h_ch2_nhit_goodscin_notrack  = new TH1F("h_ch2_nhit_goodscin_notrack "," ; CH2 nhits",50,0.,50.);
	TH1F* h_etotnorm_all = new TH1F("h_etotnorm_all"," All events; etot  norm",100,0.,3.);
	TH1F* h_etottracknorm_trig1 = new TH1F("h_etottracknorm_trig1"," Trig 1 ; etot track norm",100,0.,3.);
	TH1F* h_etottracknorm_trig4 = new TH1F("h_etottracknorm_trig4"," Trig 4 ; etot track norm",100,0.,3.);
	TH1F* h_etottracknorm_all = new TH1F("h_etottracknorm_all"," All events; etot track norm",100,0.,3.);
	TH1F* h_hgcernpeSum_all = new TH1F("h_ngcernpeSum_all"," All events; HG cer npe sum",120,0.,30.);
	TH1F* h_hgcernpeSum_all_etotcut = new TH1F("h_hgcernpeSum_all_etotcut"," All events; HG cer npe sum",120,0.,30.);
	TH1F* h_etottracknorm_curcut = new TH1F("h_etottracknorm_curcut"," Current cut; etot track norm",100,0.,3.);
	//
   Int_t nscal_reads_2=0;
   prev_read=-1;
   Int_t ps2=1;
   Double_t trackeff=1.0;
   if (nrun==7570) ps2=9;
   if (nrun==7570) trackeff=0.9915;
	for (int i = 0; i < data_entries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
	        if (event_flag[nscal_reads_2]==1&&gevtyp==1) h_ev1->Fill(gevtyp);
	        if (event_flag[nscal_reads_2]==1&&gevtyp==4) h_ev4->Fill(gevtyp);
                if (event_flag[nscal_reads_2]==1&&gevtyp==4 && EDTM_timeraw>0) h_EDTM_CUT->Fill(EDTM_timeraw);
		h_aero->Fill(aero_npeSum);
                if (event_flag[nscal_reads_2]==1&&gevtyp==4 && Hetotnorm > 0.6&&Hhgcer_npeSum!=0&&aero_npeSum>0.25&&betanotrack>0) {
                      h_ctime->Fill(ctime);
		}
		h_hgcernpeSum_all->Fill(Hhgcer_npeSum);
		if (Hetotnorm > 0.6) h_hgcernpeSum_all_etotcut->Fill(Hhgcer_npeSum);
 		if (Hetotnorm > 0.6&&Hhgcer_npeSum!=0) h_etotnorm_all->Fill(Hetotnorm);
		if (gevnum>scal_event_number[nscal_reads_2]) {
                   nscal_reads_2++;
		   //cout << " new scal read" << endl;
		}
		if (goodscinhit==1 && abs(ScinLongPos[0])<30 && abs(ScinLongPos[1])<30 ) h_starttime->Fill(starttime);
	        if (goodscinhit==1  && abs(ScinLongPos[0])<30 && abs(ScinLongPos[1])<30 && ntrack>0)h_starttime_track->Fill(starttime);
	}
	//
	Double_t ct_peak=43.3;
	Double_t ct_sig=1.2;
	Double_t ct_bg_low_peak=2.72650e+01;
	Double_t ct_bg_hi_peak=5.94827e+01;
	Double_t ct_peak_integral = h_ctime->Integral( h_ctime->FindBin(ct_peak-3*ct_sig), h_ctime->FindBin(ct_peak+3*ct_sig));
	Double_t ct_bg_low_peak_integral = h_ctime->Integral( h_ctime->FindBin(ct_bg_low_peak-3*ct_sig), h_ctime->FindBin(ct_bg_low_peak+3*ct_sig));
	Double_t ct_bg_hi_peak_integral = h_ctime->Integral( h_ctime->FindBin(ct_bg_hi_peak-3*ct_sig), h_ctime->FindBin(ct_bg_hi_peak+3*ct_sig));
	Double_t ct_bg_low = h_ctime->Integral( h_ctime->FindBin(17.5), h_ctime->FindBin(37.5));
	Double_t ct_bg_hi = h_ctime->Integral( h_ctime->FindBin(49.5), h_ctime->FindBin(69.5));
	//
	Double_t calc_treff= float(h_starttime_track->Integral())/float(h_starttime->Integral()) ;
	Int_t nev4 = h_ev4->Integral();
	Int_t nEDTM = h_EDTM_CUT->Integral();
	Int_t good_ev4 = h_ev4->Integral();
	Double_t good_ev4_err = TMath::Sqrt(good_ev4);
	Double_t good_ev = (good_ev4)/calc_treff;
	Double_t good_ev_err = good_ev4_err/calc_treff;
	//
        Double_t err_trig2 = 1./TMath::Sqrt(tot_scal_cut_TRIG2);
        Double_t err_trig3 = 1./TMath::Sqrt(tot_scal_cut_TRIG3);
        Double_t err_trig1 = 1./TMath::Sqrt(tot_scal_cut_TRIG1);
        Double_t err_trig4 = 1./TMath::Sqrt(tot_scal_cut_TRIG4);
        Double_t err_trig5 = 1./TMath::Sqrt(tot_scal_cut_TRIG5);

	Double_t err_ct_peak_integral = TMath::Sqrt(ct_peak_integral);
	Double_t err_ct_bg_low_peak_integral = TMath::Sqrt(ct_bg_low_peak_integral);
	Double_t err_ct_bg_hi_peak_integral = TMath::Sqrt(ct_bg_hi_peak_integral);
	Double_t err_ct_bg_low = TMath::Sqrt(ct_bg_low);
	Double_t err_ct_bg_hi = TMath::Sqrt(ct_bg_hi);
	cout << " Coin ev = " <<ct_peak_integral/charge_sum_cut << " "  << err_ct_peak_integral/charge_sum_cut << " "  << ct_bg_low_peak_integral/charge_sum_cut/charge_sum_cut  << " "  << err_ct_bg_low_peak_integral/charge_sum_cut/charge_sum_cut << " " << ct_bg_hi_peak_integral/charge_sum_cut/charge_sum_cut << " "  << err_ct_bg_hi_peak_integral/charge_sum_cut/charge_sum_cut << endl;
	Double_t lt=  nEDTM/tot_scal_cut_EDTM;
	Double_t coin_acc_sub = ct_peak_integral-(ct_bg_low_peak_integral);
	Double_t err_coin_acc_sub = TMath::Sqrt(ct_peak_integral+ct_bg_low_peak_integral);
	cout <<  nrun<< " "  << ct_bg_low/ct_bg_hi << " " << ct_bg_low/charge_sum_cut/charge_sum_cut/lt/calc_treff  << " "  << err_ct_bg_low/charge_sum_cut/charge_sum_cut/lt/calc_treff << " " << ct_bg_hi/charge_sum_cut/charge_sum_cut/lt/calc_treff << " "  << err_ct_bg_hi/charge_sum_cut/charge_sum_cut/lt/calc_treff << endl;
	cout << nrun << " " << charge_sum_cut/tot_scal_cut_time << " " << coin_acc_sub/charge_sum_cut/lt/calc_treff << " " << err_coin_acc_sub/charge_sum_cut/lt/calc_treff << " " << lt << " " << calc_treff<< endl;
        cout << "EV4  Livetime = " << (nev4-nEDTM)/(tot_scal_cut_TRIG5-tot_scal_cut_EDTM) << " " <<  nev4  << " " <<  nEDTM << " " <<tot_scal_cut_TRIG5  << " " <<tot_scal_cut_EDTM << endl;
        cout << "EDTM  Livetime = " << nEDTM/tot_scal_cut_EDTM << " " <<  nEDTM << " " << tot_scal_cut_EDTM << endl;
			   cout << nrun << " " << mean_current << " "<< charge_sum_cut/tot_scal_cut_time << " " << tot_scal_cut_TRIG1/tot_scal_cut_time << " " << tot_scal_cut_TRIG4/tot_scal_cut_time << " " << good_ev/tot_scal_cut_time  << " " << tot_scal_cut_TRIG1/charge_sum_cut << " " << err_trig1*tot_scal_cut_TRIG1/charge_sum_cut << " " << tot_scal_cut_TRIG4/charge_sum_cut << " " << err_trig4*tot_scal_cut_TRIG4/charge_sum_cut<< " " << good_ev/charge_sum_cut<< " " << good_ev_err/charge_sum_cut << " " << calc_treff<< endl;
			   

// loop through data
Long64_t nentries = tsimc->GetEntries();
//
}
