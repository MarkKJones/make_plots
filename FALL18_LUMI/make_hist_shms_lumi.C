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

void make_hist_shms_lumi(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_lumi_hist.root";
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
 Double_t  Scal_BCM1_charge;
   tscal->SetBranchAddress("P.BCM1.scalerCharge",&Scal_BCM1_charge);
 Double_t  Scal_BCM1_current;
   tscal->SetBranchAddress("P.BCM1.scalerCurrent",&Scal_BCM1_current);
 Double_t  Scal_time;
   tscal->SetBranchAddress("P.1MHz.scalerTime",&Scal_time);
 Double_t  Scal_EDTM;
   tscal->SetBranchAddress("P.EDTM.scaler",&Scal_EDTM);
Double_t  Scal_pEL_CLEAN;
   tscal->SetBranchAddress("P.pEL_CLEAN.scaler",&Scal_pEL_CLEAN);
Double_t  Scal_pEL_REAL;
   tscal->SetBranchAddress("P.pEL_REAL.scaler",&Scal_pEL_REAL);
 Double_t  Scal_pHGCER;
   tscal->SetBranchAddress("P.HCER.scaler",&Scal_pHGCER);
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
 Double_t  Scal_TRIG6;
   tscal->SetBranchAddress("P.pTRIG6.scaler",&Scal_TRIG6);
 Double_t  Scal_Splane[4];
   tscal->SetBranchAddress("P.S1X.scaler",&Scal_Splane[0]);
   tscal->SetBranchAddress("P.S1Y.scaler",&Scal_Splane[1]);
   tscal->SetBranchAddress("P.S2X.scaler",&Scal_Splane[2]);
   tscal->SetBranchAddress("P.S2Y.scaler",&Scal_Splane[3]);
//loop through scalers
     Int_t nscal_reads=0;
     Int_t nscal_reads_cut=0;
     Double_t prev_read=-1;
     Double_t charge_sum=0;
     Double_t charge_sum_cut=0;
     Double_t prev_charge=0;
     Double_t charge_sum_corr=0;
     Double_t charge_sum_cut_corr=0;
     Double_t prev_charge_corr=0;
     Double_t event_flag[10000];
     Double_t scal_event_number[10000];
     Double_t tot_scal_EDTM=0;
     Double_t tot_scal_cut_EDTM=0;
     Double_t prev_EDTM=0;
     Double_t tot_scal_pEL_CLEAN=0;
     Double_t tot_scal_cut_pEL_CLEAN=0;
     Double_t prev_pEL_CLEAN=0;
     Double_t tot_scal_pEL_REAL=0;
     Double_t tot_scal_cut_pEL_REAL=0;
     Double_t prev_pEL_REAL=0;
     Double_t tot_scal_pHGCER=0;
     Double_t tot_scal_cut_pHGCER=0;
     Double_t prev_pHGCER=0;
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
     Double_t tot_scal_TRIG6=0;
     Double_t prev_TRIG6=0;
    Double_t tot_scal_cut_TRIG6=0;
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
 TH1F *h_cur = new TH1F("h_cur","; Current ;",200,0,100);
Long64_t data_entries = tsimc->GetEntries();
	for (int i = 0; i < scal_entries; i++) {
      		tscal->GetEntry(i);
	  h_cur_entry->Fill(float(i),Scal_BCM1_current);
	  if (Scal_BCM1_current > 3) h_cur->Fill(Scal_BCM1_current);
	}
 Double_t peak_current = h_cur->GetBinCenter(h_cur->GetMaximumBin());
 cout << " Peak current = " << peak_current  <<" " <<  h_cur->GetMaximumBin() << endl;
 Double_t Scal_BCM1_charge_corr=0;
	for (int i = 0; i < scal_entries; i++) {
      		tscal->GetEntry(i);
          event_flag[nscal_reads] = 0;
             scal_event_number[nscal_reads] = Scal_evNumber;
	  Double_t BCM1_correction=1.;
	  if (Scal_BCM1_current >2.) {
	  if (Scal_BCM1_current <= 60) {
	    BCM1_correction =1.0 + 0.045* ( log(60.)-log(Scal_BCM1_current))/( log(60.)-log(2.) );
	  } else {
	    BCM1_correction =1.0 + 0.010*(Scal_BCM1_current-60)/25.;
	  } 
	}
	  Scal_BCM1_charge_corr+=Scal_BCM1_current*(Scal_time-prev_time)*BCM1_correction;
	  //cout << Scal_BCM1_charge << " "  << Scal_BCM1_charge_corr << " " << BCM1_correction << endl;
	  if (TMath::Abs(Scal_BCM1_current-peak_current) < threshold_cut) {
             event_flag[nscal_reads] = 1;
 	     tot_scal_cut_time+=(Scal_time-prev_time);
 	     tot_scal_cut_EDTM+=(Scal_EDTM-prev_EDTM);
 	     tot_scal_cut_pEL_CLEAN+=(Scal_pEL_CLEAN-prev_pEL_CLEAN);
	     tot_scal_cut_pEL_REAL+=(Scal_pEL_REAL-prev_pEL_REAL);
	     tot_scal_cut_pHGCER+=(Scal_pHGCER-prev_pHGCER);
 	     tot_scal_cut_TRIG2+=(Scal_TRIG2-prev_TRIG2);
 	     tot_scal_cut_TRIG3+=(Scal_TRIG3-prev_TRIG3);
 	     tot_scal_cut_TRIG1+=(Scal_TRIG1-prev_TRIG1);
 	     tot_scal_cut_TRIG4+=(Scal_TRIG4-prev_TRIG4);
 	     tot_scal_cut_TRIG5+=(Scal_TRIG5-prev_TRIG5);
 	     tot_scal_cut_TRIG6+=(Scal_TRIG6-prev_TRIG6);
	     for (Int_t s=0;s<4;s++) tot_scal_cut_Splane[s]+=(Scal_Splane[s]-prev_Splane[s]);
	     //	     cout << i << " " << tot_scal_cut_Splane[0] << " " << Scal_Splane[0] << " " << prev_Splane[0] << endl;

             charge_sum_cut+=(Scal_BCM1_charge-prev_charge);
             charge_sum_cut_corr+=(Scal_BCM1_charge_corr-prev_charge_corr);
             nscal_reads_cut++;
	  }
	  prev_charge = Scal_BCM1_charge;
	  prev_charge_corr = Scal_BCM1_charge_corr;
	  prev_time = Scal_time;
	  prev_EDTM = Scal_EDTM;
	  prev_pEL_CLEAN = Scal_pEL_CLEAN;
	  prev_pEL_REAL = Scal_pEL_REAL;
	  prev_pHGCER = Scal_pHGCER;
	  prev_TRIG2 = Scal_TRIG2;
	  prev_TRIG3 = Scal_TRIG3;
	  prev_TRIG1 = Scal_TRIG1;
	  prev_TRIG4 = Scal_TRIG4;
	  prev_TRIG5 = Scal_TRIG5;
	  prev_TRIG6 = Scal_TRIG6;
	     for (Int_t s=0;s<4;s++) prev_Splane[s]=Scal_Splane[s];
	  // cout <<  nscal_reads <<  " " << Scal_BCM4B_current << " " << event_flag[nscal_reads] << " " << Scal_TRIG1 << endl;
          nscal_reads++;
          charge_sum=Scal_BCM1_charge;
          charge_sum_corr=Scal_BCM1_charge_corr;
	  tot_scal_EDTM=Scal_EDTM;
	  tot_scal_pEL_CLEAN=Scal_pEL_CLEAN;
	  tot_scal_pEL_REAL=Scal_pEL_REAL;
	  tot_scal_pHGCER=Scal_pHGCER;
	  tot_scal_TRIG2=Scal_TRIG2;
	  tot_scal_TRIG3=Scal_TRIG3;
	  tot_scal_TRIG1=Scal_TRIG1;
	  tot_scal_TRIG4=Scal_TRIG4;
	  tot_scal_TRIG5=Scal_TRIG5;
	  tot_scal_TRIG6=Scal_TRIG6;
          tot_scal_time=Scal_time;
	  for (Int_t s=0;s<4;s++) tot_scal_Splane[s]=Scal_Splane[s];
	}
	cout << "nscal_reads_cut " << nscal_reads_cut<< endl;
	//
	Double_t EDTM_timeraw;
	tsimc->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw",&EDTM_timeraw);
	Double_t gevtyp;
	tsimc->SetBranchAddress("g.evtyp",&gevtyp);
	Double_t gevnum;
	tsimc->SetBranchAddress("g.evnum",&gevnum);
	Double_t hntrack;
	tsimc->SetBranchAddress("P.dc.ntrack",&hntrack);
	Double_t hgoodscinhit;
	tsimc->SetBranchAddress("P.hod.goodscinhit",&hgoodscinhit);
	Double_t hgoodstarttime;
	tsimc->SetBranchAddress("P.hod.goodstarttime",&hgoodstarttime);
	Double_t hfptime;
	tsimc->SetBranchAddress("P.hod.fpHitsTime",&hfptime);
	Double_t PhgcerRefTime;
	tsimc->SetBranchAddress("P.hgcer.RefTime",&PhgcerRefTime);
 static const Int_t pDCplnum=12;
 const char* Pplname[pDCplnum]={"u1","u2","x1","x2","v1","v2","v2","v1","x2","x1","u2","u1"};
 const char* Pchname[pDCplnum]={"1","1","1","1","1","1","2","2","2","2","2","2"};
            Double_t pDCtdcRefTime[pDCplnum];
           for (Int_t ip=0;ip<pDCplnum;ip++) {
   tsimc->SetBranchAddress(Form("P.dc.%s%s.RefTime",Pchname[ip],Pplname[ip]),&pDCtdcRefTime[ip]) ;
 	    }
	//
	//
          static const Int_t plnum=4;
          static const Int_t ns=2;
          const char* plname[plnum]={"1x","1y","2x","2y"};
           static const Int_t npad[plnum]={13,13,14,21};
            const char* sname[ns]={"Neg","Pos"};
             const char* rname[ns]={"neg","pos"};
            Double_t tdcRefTime[plnum][ns];
            Double_t adcRefTime[plnum][ns];
            Double_t tdcRefDiffTime[plnum][ns];
            Double_t adcRefDiffTime[plnum][ns];
            for (Int_t ip=0;ip<plnum;ip++) {
            for (Int_t is=0;is<ns;is++) {
       tsimc->SetBranchAddress(Form("P.hod.%s.%sTdcRefTime",plname[ip],sname[is]),&tdcRefTime[ip][is]) ;   
       tsimc->SetBranchAddress(Form("P.hod.%s.%sAdcRefTime",plname[ip],sname[is]),&adcRefTime[ip][is]) ;   
       tsimc->SetBranchAddress(Form("P.hod.%s.%sTdcRefDiffTime",plname[ip],sname[is]),&tdcRefDiffTime[ip][is]) ;   
       tsimc->SetBranchAddress(Form("P.hod.%s.%sAdcRefDiffTime",plname[ip],sname[is]),&adcRefDiffTime[ip][is]) ;   
                }}
	Double_t hstarttime;
	tsimc->SetBranchAddress("P.hod.starttime",&hstarttime);
	Double_t hTimeHist_StartTime_Peak;
	tsimc->SetBranchAddress("P.hod.TimeHist_StartTime_Peak",&hTimeHist_StartTime_Peak);
 	Double_t hTimeHist_StartTime_NumPeaks;
	tsimc->SetBranchAddress("P.hod.TimeHist_StartTime_NumPeaks",&hTimeHist_StartTime_NumPeaks);
 	Double_t hTimeHist_StartTime_Hits;
	tsimc->SetBranchAddress("P.hod.TimeHist_StartTime_Hits",&hTimeHist_StartTime_Hits);
 	Double_t hbetatrack;
		tsimc->SetBranchAddress("P.gtr.beta",&hbetatrack);
   Double_t hdelta;
   tsimc->SetBranchAddress("P.gtr.dp",&hdelta);
   Double_t Hetotnorm;
   tsimc->SetBranchAddress("P.cal.etotnorm",&Hetotnorm);
  Double_t Hhgcer_npeSum;
   tsimc->SetBranchAddress("P.hgcer.npeSum",&Hhgcer_npeSum);
   //
 TH1F* h_PhgcerRefTime = new TH1F("h_PhgcerRefTime",";SHMS HG Cer DC ADC Ref Time (chan)",2000,0.,20000.);
 TH1F* h_PDCtdcRefTime;
	      h_PDCtdcRefTime = new TH1F(Form("h_PDCtdcRefTime_%d",0),Form("; Pl %d  SHMS DC TDC Ref Time (chan)",0),2000,0.,20000.);
//
 TH1F* h_tdcRefTime;
 TH1F* h_adcRefTime;
 TH1F* h_tdcRefDiffTime;
 TH1F* h_adcRefDiffTime;
	      h_tdcRefTime = new TH1F(Form("h_tdcRefTime_%d_%s",0,sname[0]),Form("; Pl %d %s  SHMS Hodo TDC Ref Time (chan)",0,sname[0]),2000,0.,20000.);
	      h_adcRefTime= new TH1F(Form("h_adcRefTime_%d_%s",0,sname[0]),Form("; Pl %d %s  SHMS Hodo ADC ref time (chan)",0,sname[0]),2000,0.,20000.);
	      h_tdcRefDiffTime = new TH1F(Form("h_tdcRefDiffTime_%d_%s",0,sname[0]),Form("; Pl %d %s  SHMS Hodo TDC Ref Diff Time (ns)",0,sname[0]),500,0.,500.);
	      h_adcRefDiffTime = new TH1F(Form("h_adcRefDiffTime_%d_%s",0,sname[0]),Form("; Pl %d %s  SHMS Hodo ADC ref Diff time (ns)",0,sname[0]),500,0.,500.);
   //
	TH1F* h_EDTM_CUT = new TH1F("h_EDTM"," ; EDTM (beam cut)",1000,0.,10000.);
	TH1F* h_EDTM_CUT2 = new TH1F("h_EDTM2"," ; EDTM (beam cut2)",1000,0.,10000.);
	TH1F* h_hStartTime = new TH1F("h_hStartTime",";HMS Starttime",200,0,200.);
	TH1F* h_hStartTime_track = new TH1F("h_hStartTime_track",";HMS Starttime",200,0,200.);
	TH1F* h_hFpTime = new TH1F("h_hFpTime",";HMS Fptime",200,0,200.);
	TH1F* h_etotnorm = new TH1F("h_etotnorm"," ; Etot norm",100,0.,2.);
	TH1F* h_npeSum = new TH1F("h_npeSum"," ; Cer Npe",100,0.,20.);
	TH1F* h_ev2 = new TH1F("h_ev2"," ; hEv2",10,0.,10.);
	TH1F* h_goodev2 = new TH1F("h_goodev2"," ; hgoodEv2",10,0.,10.);
	TH1F* h_goodev = new TH1F("h_goodev"," ; Good ev2 with PID cuts",10,0.,10.);
	TH1F* h_goodev_goodreftime = new TH1F("h_goodev_goodreftime"," ; Good ev2 with PID cuts_goodreftime",10,0.,10.);
	TH1F* h_goodevTrack = new TH1F("h_goodevTrack"," ; Good ev2 with PID/track cuts",10,0.,10.);
	
  Int_t nscal_reads_2=0;
  Int_t cnts_goodev=0;
  Int_t cnts_goodev_goodreftime=0;
	for (int i = 0; i < data_entries; i++) {
      		tsimc->GetEntry(i);
		//
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (event_flag[nscal_reads_2]==1&& (gevtyp==1 ||gevtyp==3)  && nscal_reads_2<scal_entries-1) {
		  Bool_t GoodEDTM = EDTM_timeraw>0;
		  if (GoodEDTM) h_EDTM_CUT->Fill(EDTM_timeraw); 
		  if (EDTM_timeraw>0) h_EDTM_CUT2->Fill(EDTM_timeraw); 
		    h_ev2->Fill(gevtyp);
		    if (!GoodEDTM) h_goodev2->Fill(gevtyp);
		    h_etotnorm->Fill(Hetotnorm);
		    if (Hetotnorm > .6) h_npeSum->Fill(Hhgcer_npeSum);
                 h_PhgcerRefTime->Fill(PhgcerRefTime);
	      h_tdcRefTime->Fill(tdcRefTime[0][0]);
	      h_adcRefTime->Fill(adcRefTime[0][0]);
	      if (tdcRefDiffTime[0][0] !=0) h_tdcRefDiffTime->Fill(tdcRefDiffTime[0][0]*.1);
	      if (adcRefDiffTime[0][0] !=0) h_adcRefDiffTime->Fill(adcRefDiffTime[0][0]*.0625);
              h_PDCtdcRefTime->Fill(pDCtdcRefTime[0]);
	      Bool_t GoodAdcRefTime = adcRefTime[0][0]>5500 &&  adcRefTime[0][0]<5800;
	      Bool_t GoodTdcRefTime = tdcRefTime[0][0]>4600 &&  tdcRefTime[0][0]<5300;
	      Bool_t GoodTdcDiffRefTime = tdcRefDiffTime[0][0]*.1>210. || tdcRefDiffTime[0][0]*.1<170.   ;
	      Bool_t GoodAdcDiffRefTime = adcRefDiffTime[0][0]*.0625>200. || adcRefDiffTime[0][0]*.0625<160.   ;
	      Bool_t GoodRefTime = GoodAdcRefTime && GoodTdcRefTime && GoodTdcDiffRefTime && GoodAdcDiffRefTime ;
		   cnts_goodev++;
		   if (GoodRefTime) cnts_goodev_goodreftime++;
	 if (Hetotnorm > 0.6&&Hhgcer_npeSum>1&& hgoodstarttime==1&& hstarttime!=-1000 &&!GoodEDTM ) {
		if (hgoodscinhit==1) h_hStartTime->Fill(hstarttime);
	        if (hgoodscinhit==1 && hntrack>0) h_hStartTime_track->Fill(hstarttime);
		   h_goodev->Fill(gevtyp);
		   if (GoodRefTime) h_goodev_goodreftime->Fill(gevtyp);
		   if (hntrack>0 && hdelta> -10 && hdelta<25.) h_goodevTrack->Fill(gevtyp);
		 
	 }
		    //
	 }
	 if (gevnum>scal_event_number[nscal_reads_2])  nscal_reads_2++;
	}
	//
        Double_t ps = 1.;
	if (nrun == 6482) ps=513.;
	if (nrun == 6483) ps=513.;
	if (nrun == 6484) ps=257.;
	if (nrun == 6485) ps=257.;
	Int_t nev2 = h_ev2->Integral();
	Int_t good_ev2 = h_goodev2->Integral();
	Double_t good_ev2_err = TMath::Sqrt(good_ev2);
	Int_t nEDTM = h_EDTM_CUT->Integral();
	Int_t nEDTM2 = h_EDTM_CUT2->Integral();
	Double_t clt = (good_ev2)/(tot_scal_cut_TRIG2);
	Double_t clt_err = sqrt( (1-clt)*clt/tot_scal_cut_TRIG2);
	if ( nrun >= 6482 &&  nrun <= 6485) {
	clt = (good_ev2)/(tot_scal_cut_TRIG1);
	clt_err = sqrt( (1-clt)*clt/tot_scal_cut_TRIG1);
	}
	Double_t calc_treff= float(h_hStartTime_track->Integral())/float(h_hStartTime->Integral()) ;
	Double_t lt=  nEDTM/tot_scal_cut_EDTM;
	Double_t lt_err = sqrt( (1-lt)*lt/tot_scal_cut_EDTM);
	Double_t GoodRat = float(cnts_goodev)/float(cnts_goodev_goodreftime);
	Double_t good_ev = float(h_goodev->Integral())/clt;
	Double_t good_ev_err = TMath::Sqrt(h_goodev->Integral())/clt;
	Double_t good_ev_goodreftime = float(h_goodev_goodreftime->Integral())/clt;
	Double_t good_ev_goodreftime_err = TMath::Sqrt(h_goodev_goodreftime->Integral())/clt;
	Double_t good_evTrack = float(h_goodevTrack->Integral())/clt/calc_treff;
	Double_t good_evTrack_err = TMath::Sqrt(h_goodevTrack->Integral())/clt/calc_treff;
 	//
	cout << nEDTM << " " << nEDTM2 << " " << tot_scal_cut_EDTM<< " EDTM rate = " << tot_scal_cut_EDTM/tot_scal_cut_time << endl;
	//
       Double_t err_pEL_CLEAN = 1./TMath::Sqrt(tot_scal_cut_pEL_CLEAN);
        Double_t err_pEL_REAL = 1./TMath::Sqrt(tot_scal_cut_pEL_REAL);
        Double_t err_pHGCER = 1./TMath::Sqrt(tot_scal_cut_pHGCER);
        Double_t err_trig2 = 1./TMath::Sqrt(tot_scal_cut_TRIG2);
        Double_t err_trig3 = 1./TMath::Sqrt(tot_scal_cut_TRIG3);
        Double_t err_trig1 = 1./TMath::Sqrt(tot_scal_cut_TRIG1);
        Double_t err_trig4 = 1./TMath::Sqrt(tot_scal_cut_TRIG4);
        Double_t err_trig5 = 1./TMath::Sqrt(tot_scal_cut_TRIG5);
        Double_t err_trig6 = 1./TMath::Sqrt(tot_scal_cut_TRIG6);
        cout << " data " << endl;
	cout << nrun << " "<< charge_sum_cut/tot_scal_cut_time  << " "<< charge_sum_cut_corr/tot_scal_cut_time << " " << clt*ps << " " << clt_err*ps<< " " << lt<< " " << lt_err  << " " << calc_treff << " " << good_ev2/tot_scal_cut_time << " " << good_ev/charge_sum_cut << " " <<  good_ev_err/charge_sum_cut << " " << good_evTrack/charge_sum_cut << " " <<  good_evTrack_err/charge_sum_cut  << " " << tot_scal_cut_pEL_CLEAN/charge_sum_cut << " " << err_pEL_CLEAN*tot_scal_cut_pEL_CLEAN/charge_sum_cut << " " << tot_scal_cut_pEL_REAL/charge_sum_cut << " " << err_pEL_REAL*tot_scal_cut_pEL_REAL/charge_sum_cut << " " <<tot_scal_cut_pEL_REAL/tot_scal_cut_time<< endl;
        cout << " data " << endl;
	cout << nrun << " "<< charge_sum_cut/tot_scal_cut_time  << " "<< charge_sum_cut_corr/tot_scal_cut_time << " " << clt << " " << lt << " " << calc_treff << " " << good_ev2/tot_scal_cut_time << " " << good_ev/charge_sum_cut << " " <<  good_ev_err/charge_sum_cut  << " " << good_evTrack/charge_sum_cut << " " <<  good_evTrack_err/charge_sum_cut << " " << tot_scal_cut_pEL_CLEAN/charge_sum_cut << " " << err_pEL_CLEAN*tot_scal_cut_pEL_CLEAN/charge_sum_cut << " " << tot_scal_cut_pEL_REAL/charge_sum_cut << " " << err_pEL_REAL*tot_scal_cut_pEL_REAL/charge_sum_cut << " " <<tot_scal_cut_TRIG1/tot_scal_cut_time<< " " <<tot_scal_cut_pHGCER/tot_scal_cut_time<< " " << good_ev_goodreftime/charge_sum_cut << " " <<  good_ev_goodreftime_err/charge_sum_cut << " " << GoodRat << endl;
        cout << " Scalers " << endl;
	cout << nrun << " "<< charge_sum_cut/tot_scal_cut_time<< " "<< charge_sum_cut_corr/tot_scal_cut_time << " " << tot_scal_cut_pEL_CLEAN/charge_sum_cut << " " << err_pEL_CLEAN*tot_scal_cut_pEL_CLEAN/charge_sum_cut << " " << tot_scal_cut_pEL_REAL/charge_sum_cut << " " << err_pEL_REAL*tot_scal_cut_pEL_REAL/charge_sum_cut << " " <<tot_scal_cut_pEL_CLEAN/tot_scal_cut_time << " " <<tot_scal_cut_pEL_REAL/tot_scal_cut_time  << endl;
	cout <<  nrun << " "<<tot_scal_cut_TRIG1/tot_scal_cut_time<< " "<<tot_scal_cut_Splane[0]/tot_scal_cut_time<< " "<<tot_scal_cut_Splane[1]/tot_scal_cut_time<< " "<<tot_scal_cut_Splane[2]/tot_scal_cut_time<< " "<<tot_scal_cut_Splane[3]/tot_scal_cut_time << endl;
			   Double_t prob1=tot_scal_cut_Splane[0]/tot_scal_cut_time*50e-9;
			   Double_t prob2=tot_scal_cut_Splane[1]/tot_scal_cut_time*50e-9;
			   Double_t prob3=tot_scal_cut_Splane[2]/tot_scal_cut_time*50e-9;
			   Double_t prob4=tot_scal_cut_Splane[3]/tot_scal_cut_time*50e-9;
			   Double_t prob=prob1*prob2+prob1*prob3+prob1*prob4+prob2*prob3+prob2*prob4+prob3*prob4;
			   cout << nrun << " " << prob1 << " " << prob2 << " " << prob3 << " " << prob4 << " " << prob << " prob trig = " << tot_scal_cut_TRIG1/tot_scal_cut_time*60e-9 << "TRGI1 rate = " << tot_scal_cut_TRIG1/tot_scal_cut_time << endl; 

	//
}
