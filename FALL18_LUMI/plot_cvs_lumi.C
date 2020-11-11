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
#include <TText.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
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

void plot_cvs_lumi() {
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(0);
 gStyle->SetTitleOffset(1.2,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.16);
   TString outputpdf;
   outputpdf= "plots/cvs_lumi_-2.pdf";
   Bool_t DumSub= kTRUE;
   //
   //const Int_t ntar=3;
   TString DataFile;
   DataFile = "FALL18_LUMI/shms_cvs_lumi_-2.dat";
   //
   vector<Double_t>    Cur;
   vector<Double_t>    CurErr;
   vector<Double_t>    CurCorr;
   vector<Double_t>    CurCorrErr;
   vector<Double_t>    CompLT;
   vector<Double_t>    EDTMLT;
   vector<Double_t>    YPid;
   vector<Double_t>    YPidErr;
   vector<Double_t>    YTrack;
   vector<Double_t>    YTrackErr;
   vector<Double_t>    YElReal;
   vector<Double_t>    YElRealErr;
   vector<Double_t>    YElClean;
   vector<Double_t>    YElCleanErr;
   vector<Double_t>    TrigRate;
   vector<Double_t>    TrigRateErr;
   vector<Double_t>    YPidGoodRtime;
   vector<Double_t>    YPidGoodRtimeErr;
   vector<Double_t>    HGC_Rate;
   vector<Double_t>    GoodRat;
   //
  // nrun curr currcorr clt lt treff ev2rate yield_pid err yield_track err Yield_EL_CLEAN err Yield_EL_REAL err TrigRate
   Double_t RatCur=1.0;
   Int_t nrun=0;
     ifstream file_Data(DataFile.Data());
     RatCur=1.0;
     nrun=0;
      if (file_Data.is_open()) {
	Int_t nr;
	file_Data >> nr;
	
	for (Int_t ny=0;ny<nr;ny++)  {
		  Double_t temp;
		  for (Int_t i=0;i<20;i++) {
		    file_Data >> temp;
		    if (i==1) Cur.push_back(temp);
		    if (i==1) CurErr.push_back(0.0001);
		    if (i==2) CurCorr.push_back(temp);
		    if (i==2) CurCorrErr.push_back(0.0001);
		    //if (i==3) RatCur=Cur[nrun]/CurCorr[nrun];
		    if (i==3) CompLT.push_back(temp);
		    if (i==4) EDTMLT.push_back(temp);
		    if (i==7) YPid.push_back(RatCur*temp);
		    if (i==8) YPidErr.push_back(RatCur*temp);
		    if (i==9) YTrack.push_back(RatCur*temp);
		    if (i==10) YTrackErr.push_back(RatCur*temp);
		    if (i==11) YElClean.push_back(RatCur*temp);
		    if (i==12) YElCleanErr.push_back(RatCur*temp);
		    if (i==13) YElReal.push_back(RatCur*temp);
		    if (i==14) YElRealErr.push_back(RatCur*temp);
		    if (i==15) TrigRate.push_back(RatCur*temp/1000.);
		    if (i==15) TrigRateErr.push_back(RatCur*0.0001);
		    if (i==16) HGC_Rate.push_back(temp);
		    if (i==17) YPidGoodRtime.push_back(temp);
		    if (i==18) YPidGoodRtimeErr.push_back(temp);
		    if (i==19) GoodRat.push_back(temp);
		    if (i==19) {
		      Double_t drate = HGC_Rate[nrun]-TrigRate[nrun]*1000;
		      // drate = HGC_Rate[nrun];
		      YPid[nrun]*=1.0/(1-TrigRate[nrun]*50e-9*1000)/(1-(HGC_Rate[nrun])*100e-9);;
		      YPidErr[nrun]*=1./(1-TrigRate[nrun]*50e-9*1000)/(1-(HGC_Rate[nrun])*100e-9);;
		      YPidGoodRtime[nrun]*=GoodRat[nrun]/(1-TrigRate[nrun]*50e-9*1000)/(1-(drate)*100e-9);;
		      YPidGoodRtimeErr[nrun]*=GoodRat[nrun]/(1-TrigRate[nrun]*50e-9*1000)/(1-(drate)*100e-9);;
		      		      YElReal[nrun]=YElReal[nrun]/(1-TrigRate[nrun]*50e-9*1000);
		        YElRealErr[nrun]=YElRealErr[nrun]/(1-TrigRate[nrun]*50e-9*1000);
			YElClean[nrun]=YElClean[nrun]/(1-TrigRate[nrun]*50e-9*1000);
		       YElCleanErr[nrun]=YElCleanErr[nrun]/(1-TrigRate[nrun]*50e-9*1000);
		         YElClean[nrun]=YElClean[nrun]/(1-(drate)*50e-9);
		       YElCleanErr[nrun]=YElCleanErr[nrun]/(1-(drate)*50e-9);
		    }
		  }
		  nrun++;
	}
	file_Data.close();
      } else {
	cout << DataFile.Data() << endl;
	return;
      }
   // Use the first run of subtract dummy from the LD2 and LH2
   // 
      //
      TGraphErrors *gYPid; 
      TGraphErrors *gYPidGoodRF; 
      TGraphErrors *gYTrack; 
      TGraphErrors *gYElReal; 
       TGraphErrors *gYElClean; 
       TGraph *gCompLT; 
       TGraph *gEDTMLT; 
       TGraph *gRatCompEDTM; 
       TString TarLabel = " CVS SHMS lumi";
	cout << TarLabel << " Nrun  = " << Cur.size() << endl;
	gYPid= new TGraphErrors(Cur.size(), &Cur[0], &YPid[0], &CurErr[0], &YPidErr[0]);
	gYPid->SetMarkerStyle(22);
	gYPid->SetTitle(TarLabel);
	gYPid->GetXaxis()->SetTitle("Current (uA)");
	gYPid->GetYaxis()->SetTitle("Yield (cnts/uC) PID only");
	gYPidGoodRF= new TGraphErrors(Cur.size(), &Cur[0], &YPidGoodRtime[0], &CurErr[0], &YPidGoodRtimeErr[0]);
	gYPidGoodRF->SetMarkerStyle(22);
	gYPidGoodRF->SetTitle(TarLabel);
	gYPidGoodRF->GetXaxis()->SetTitle("Current (uA)");
	gYPidGoodRF->GetYaxis()->SetTitle("Yield (cnts/uC) Good RF+PID");
	gYTrack= new TGraphErrors(Cur.size(), &Cur[0], &YTrack[0], &CurErr[0], &YTrackErr[0]);
	gYTrack->SetMarkerStyle(23);
	gYTrack->SetTitle(TarLabel);
	gYTrack->GetXaxis()->SetTitle("Current (uA)");
	gYTrack->GetYaxis()->SetTitle("Yield (cnts/uC) PID+Track");
	gYElClean= new TGraphErrors(Cur.size(), &Cur[0], &YElClean[0], &CurErr[0], &YElCleanErr[0]);
	gYElClean->SetMarkerStyle(23);
	gYElClean->SetTitle(TarLabel);
	gYElClean->GetXaxis()->SetTitle("Current (uA)");
	gYElClean->GetYaxis()->SetTitle("Scaler Yield (cnts/uC) ElClean");
	gYElReal= new TGraphErrors(Cur.size(), &Cur[0], &YElReal[0], &CurErr[0], &YElRealErr[0]);
	gYElReal->SetMarkerStyle(23);
	gYElReal->SetTitle(TarLabel);
	gYElReal->GetXaxis()->SetTitle("Current (uA)");
	gYElReal->GetYaxis()->SetTitle("Scaler Yield (cnts/uC) ElReal");
	//
      //
      TCanvas *cTar;
      Double_t Ypidslope;
      Double_t YpidslopeErr;
      Double_t Ypidint;
      Double_t YpidGoodRFslope;
      Double_t YpidGoodRFslopeErr;
      Double_t YpidGoodRFint;
     Double_t YTrackslope;
     Double_t YTrackslopeErr;
      Double_t YTrackint;
     Double_t YElCleanslope;
     Double_t YElCleanslopeErr;
      Double_t YElCleanint;
     Double_t YElRealslope;
     Double_t YElRealslopeErr;
      Double_t YElRealint;
	 cTar = new TCanvas(Form("Ctar_%d",0),Form("Ctarl_%d",0), 700,700);
	 cTar->Divide(2,3);
	 cTar->cd(1);
	 gYPid->Draw("AP");
	 gYPid->Fit("pol1","Q","",0,60);
	 TF1 *fit = gYPid->GetFunction("pol1");
	 Ypidslope = fit->GetParameter(1);
	 YpidslopeErr = fit->GetParError(1);
	 Ypidint = fit->GetParameter(0);
	 Double_t normSlope=Ypidslope/Ypidint*100*100;
	 Double_t normSlopeErr= YpidslopeErr/Ypidint*100*100;
	 TString ttemp=Form(" Slope ($\%$/100uA) = %3.2f +/- %3.2f",normSlope,normSlopeErr);
	 gYPid->SetTitle(TarLabel+ttemp);
	 cTar->cd(5);
	 gYPidGoodRF->Draw("AP");
	 gYPidGoodRF->Fit("pol1","Q","",0,60);
	 fit = gYPidGoodRF->GetFunction("pol1");
	 YpidGoodRFslope = fit->GetParameter(1);
	 YpidGoodRFslopeErr = fit->GetParError(1);
	 YpidGoodRFint = fit->GetParameter(0);
	 normSlope=YpidGoodRFslope/YpidGoodRFint*100*100;
	 normSlopeErr= YpidGoodRFslopeErr/YpidGoodRFint*100*100;
	 ttemp=Form(" Slope (\%/100uA) = %3.2f +/- %3.2f",normSlope,normSlopeErr);
	 gYPidGoodRF->SetTitle(TarLabel+ttemp);
	 cTar->cd(2);
	 gYTrack->Draw("AP");
	 gYTrack->Fit("pol1","Q","",0,60);
	 fit = gYTrack->GetFunction("pol1");
	 YTrackslope = fit->GetParameter(1);
	 YTrackslopeErr = fit->GetParError(1);
	 YTrackint = fit->GetParameter(0);
	 normSlope=YTrackslope/YTrackint*100*100;
	 normSlopeErr= YTrackslopeErr/YTrackint*100*100;
	 ttemp=Form(" Slope ($\%$/100uA) = %3.2f +/- %3.2f",normSlope,normSlopeErr);
	 gYTrack->SetTitle(TarLabel+ttemp);
	 cTar->cd(3);
	 gYElClean->Draw("AP");
	 gYElClean->Fit("pol1","Q","",0,60);
	 fit = gYElClean->GetFunction("pol1");
	 YElCleanslope = fit->GetParameter(1);
	 YElCleanslopeErr = fit->GetParError(1);
	 YElCleanint = fit->GetParameter(0);
	 normSlope=YElCleanslope/YElCleanint*100*100;
	 normSlopeErr= YElCleanslopeErr/YElCleanint*100*100;
	 ttemp=Form(" Slope ($\%$/100uA) = %3.2f +/- %3.2f",normSlope,normSlopeErr);
	 gYElClean->SetTitle(TarLabel+ttemp);
	 cTar->cd(4);
	 gYElReal->Draw("AP");
	 gYElReal->Fit("pol1","Q","",0,60);
	 fit = gYElReal->GetFunction("pol1");
	 YElRealslope = fit->GetParameter(1);
	 YElRealslopeErr = fit->GetParError(1);
	 YElRealint = fit->GetParameter(0);
	 normSlope=YElRealslope/YElRealint*100*100;
	 normSlopeErr= YElRealslopeErr/YElRealint*100*100;
	 ttemp=Form(" Slope ($\%$/100uA) = %3.2f +/- %3.2f",normSlope,normSlopeErr);
	 gYElReal->SetTitle(TarLabel+ttemp);
	 cout << TarLabel << " Pid slope (%/100uA) = " << Ypidslope/Ypidint*100*100 << " +/- " << YpidslopeErr/Ypidint*100*100 << " Track slope (%/100uA) = " << YTrackslope/YTrackint*100*100 << " +/- "<< YTrackslopeErr/YTrackint*100*100 <<endl;
	 cTar->Print(outputpdf);
       //
  //
}
