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

void plot_kaonlt_yield(Int_t ntar=3, TString spec="HMS") {
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(0);
 gStyle->SetTitleOffset(1.2,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.16);
   TString outputpdf;
   outputpdf= "plots/fall18_"+spec+"_kaonlt_lumi.pdf";
   Bool_t DumSub= kTRUE;
   //
   //const Int_t ntar=3;
   TString DataFile[ntar];
   TString TarLabel[ntar];
   TString TarNames[3]={"LH2","Carbon","Dummy"};
   for (Int_t i=0;i<ntar;i++) {
     TarLabel[i]=spec+" "+TarNames[i];
      DataFile[i] = "FALL18_LUMI/"+spec+"_"+TarNames[i]+"_kaonlt_data.dat";
   }
   //
   vector<vector<Double_t> >   Cur;
   vector<vector<Double_t> >   CurErr;
   vector<vector<Double_t> >   CurCorr;
   vector<vector<Double_t> >   CurCorrErr;
   vector<vector<Double_t> >   CompLT;
   vector<vector<Double_t> >   EDTMLT;
   vector<vector<Double_t> >   YPid;
   vector<vector<Double_t> >   YPidErr;
   vector<vector<Double_t> >   YTrack;
   vector<vector<Double_t> >   YTrackErr;
   vector<vector<Double_t> >   YElReal;
   vector<vector<Double_t> >   YElRealErr;
   vector<vector<Double_t> >   YElClean;
   vector<vector<Double_t> >   YElCleanErr;
   vector<vector<Double_t> >   TrigRate;
   vector<vector<Double_t> >   TrigRateErr;
   //
   Cur.resize(ntar);
   CurCorr.resize(ntar);
   CompLT.resize(ntar);
   EDTMLT.resize(ntar);
    CurErr.resize(ntar);
    CurCorrErr.resize(ntar);
   YPid.resize(ntar);
   YPidErr.resize(ntar);
   YTrack.resize(ntar);
   YTrackErr.resize(ntar);
   YElClean.resize(ntar);
   YElCleanErr.resize(ntar);
   YElReal.resize(ntar);
   YElRealErr.resize(ntar);
   TrigRate.resize(ntar);
   TrigRateErr.resize(ntar);
  // nrun curr currcorr clt lt treff ev2rate yield_pid err yield_track err Yield_EL_CLEAN err Yield_EL_REAL err TrigRate
   Double_t RatCur=1.0;
   Int_t nrun=0;
   for (Int_t nt=0;nt<ntar;nt++) {
     ifstream file_Data(DataFile[nt].Data());
     RatCur=1.0;
     nrun=0;
      if (file_Data.is_open()) {
	Int_t nr;
	file_Data >> nr;
	
	for (Int_t ny=0;ny<nr;ny++)  {
		  Double_t temp;
		  for (Int_t i=0;i<16;i++) {
		    file_Data >> temp;
		    if (i==1) Cur[nt].push_back(temp);
		    if (i==1) CurErr[nt].push_back(0.0001);
		    if (i==2) CurCorr[nt].push_back(temp);
		    if (i==2) CurCorrErr[nt].push_back(0.0001);
		    //if (i==3) RatCur=Cur[nt][nrun]/CurCorr[nt][nrun];
		    if (i==3) CompLT[nt].push_back(temp);
		    if (i==4) EDTMLT[nt].push_back(temp);
		    if (i==7) YPid[nt].push_back(RatCur*temp);
		    if (i==8) YPidErr[nt].push_back(RatCur*temp);
		    if (i==9) YTrack[nt].push_back(RatCur*temp);
		    if (i==10) YTrackErr[nt].push_back(RatCur*temp);
		    if (i==11) YElClean[nt].push_back(RatCur*temp);
		    if (i==12) YElCleanErr[nt].push_back(RatCur*temp);
		    if (i==13) YElReal[nt].push_back(RatCur*temp);
		    if (i==14) YElRealErr[nt].push_back(RatCur*temp);
		    if (i==15) TrigRate[nt].push_back(RatCur*temp/1000.);
		    if (i==15) TrigRateErr[nt].push_back(RatCur*0.0001);
		  }
		  cout << nt << nrun << " cur ratio = " << RatCur<< endl;
		  nrun++;
	}
	file_Data.close();
      } else {
	cout << DataFile[nt].Data() << endl;
	return;
      }
   }
   // Use the first run of subtract dummy from the LD2 and LH2
   // 
   if (DumSub) {
   Int_t DumNt=2;
   Double_t dumscal =  (.015+.019)*2.7/.363;
   for (Int_t nn=0;nn<Cur[0].size();nn++) {
     cout << " nn = " << Cur[0][nn] << " ratio = " << YPid[DumNt][0]*dumscal/YPid[0][nn] << endl;
     YPid[0][nn] = YPid[0][nn] - YPid[DumNt][0]*dumscal;
     YPidErr[0][nn] = sqrt(YPidErr[0][nn]*YPidErr[0][nn]+YPidErr[DumNt][0]*dumscal*YPidErr[DumNt][0]*dumscal);
     YTrack[0][nn] = YTrack[0][nn] - YTrack[DumNt][0]*dumscal;
     YTrackErr[0][nn] = sqrt(YTrackErr[0][nn]*YTrackErr[0][nn]+YTrackErr[DumNt][0]*dumscal*YTrackErr[DumNt][0]*dumscal);
     YElClean[0][nn] = YElClean[0][nn] - YElClean[DumNt][0]*dumscal;
     YElCleanErr[0][nn] = sqrt(YElCleanErr[0][nn]*YElCleanErr[0][nn]+YElCleanErr[DumNt][0]*dumscal*YElClean[DumNt][0]*dumscal);
     YElReal[0][nn] = YElReal[0][nn] - YElReal[DumNt][0]*dumscal;
     YElRealErr[0][nn] = sqrt(YElRealErr[0][nn]*YElRealErr[0][nn]+YElRealErr[DumNt][0]*dumscal*YElRealErr[DumNt][0]*dumscal);
   }
   if (ntar==4) {
   for (Int_t nn=0;nn<Cur[3].size();nn++) {
     cout << " ratio = " << YPid[DumNt][0]*dumscal/YPid[3][nn] << endl;
     YPid[3][nn] = YPid[3][nn] - YPid[DumNt][0]*dumscal;
     YPidErr[3][nn] = sqrt(YPidErr[3][nn]*YPidErr[3][nn]+YPidErr[DumNt][0]*dumscal*YPidErr[DumNt][0]*dumscal);
     YTrack[3][nn] = YTrack[3][nn] - YTrack[DumNt][0]*dumscal;
     YTrackErr[3][nn] = sqrt(YTrackErr[3][nn]*YTrackErr[3][nn]+YTrackErr[DumNt][0]*dumscal*YTrackErr[DumNt][0]*dumscal);
     YElClean[3][nn] = YElClean[3][nn] - YElClean[DumNt][0]*dumscal;
     YElCleanErr[3][nn] = sqrt(YElCleanErr[3][nn]*YElCleanErr[3][nn]+YElCleanErr[DumNt][0]*dumscal*YElCleanErr[DumNt][0]*dumscal);
     YElReal[3][nn] = YElReal[3][nn] - YElReal[DumNt][0]*dumscal;
     YElRealErr[3][nn] = sqrt(YElRealErr[3][nn]*YElRealErr[3][nn]+YElRealErr[DumNt][0]*dumscal*YElRealErr[DumNt][0]*dumscal);
   }
   }
   }
      //
      TGraphErrors *gYPid[ntar]; 
      TGraphErrors *gYTrack[ntar]; 
      TGraphErrors *gYElReal[ntar]; 
       TGraphErrors *gYElClean[ntar]; 
       TGraph *gCompLT[ntar]; 
       TGraph *gEDTMLT[ntar]; 
       TGraph *gRatCompEDTM[ntar]; 
     for (Int_t nt=0;nt<ntar;nt++) {
	cout << TarLabel[nt] << " Nrun  = " << Cur[nt].size() << endl;
	gYPid[nt]= new TGraphErrors(Cur[nt].size(), &Cur[nt][0], &YPid[nt][0], &CurErr[nt][0], &YPidErr[nt][0]);
	gYPid[nt]->SetMarkerStyle(22);
	gYPid[nt]->SetTitle(TarLabel[nt]);
	gYPid[nt]->GetXaxis()->SetTitle("Current (uA)");
	gYPid[nt]->GetYaxis()->SetTitle("Yield (cnts/uC) PID only");
	gYTrack[nt]= new TGraphErrors(Cur[nt].size(), &Cur[nt][0], &YTrack[nt][0], &CurErr[nt][0], &YTrackErr[nt][0]);
	gYTrack[nt]->SetMarkerStyle(23);
	gYTrack[nt]->SetTitle(TarLabel[nt]);
	gYTrack[nt]->GetXaxis()->SetTitle("Current (uA)");
	gYTrack[nt]->GetYaxis()->SetTitle("Yield (cnts/uC) PID+Track");
	gYElClean[nt]= new TGraphErrors(Cur[nt].size(), &Cur[nt][0], &YElClean[nt][0], &CurErr[nt][0], &YElCleanErr[nt][0]);
	gYElClean[nt]->SetMarkerStyle(23);
	gYElClean[nt]->SetTitle(TarLabel[nt]);
	gYElClean[nt]->GetXaxis()->SetTitle("Current (uA)");
	gYElClean[nt]->GetYaxis()->SetTitle("Scaler Yield (cnts/uC) ElClean");
	gYElReal[nt]= new TGraphErrors(Cur[nt].size(), &Cur[nt][0], &YElReal[nt][0], &CurErr[nt][0], &YElRealErr[nt][0]);
	gYElReal[nt]->SetMarkerStyle(23);
	gYElReal[nt]->SetTitle(TarLabel[nt]);
	gYElReal[nt]->GetXaxis()->SetTitle("Current (uA)");
	gYElReal[nt]->GetYaxis()->SetTitle("Scaler Yield (cnts/uC) ElReal");
	//
	gCompLT[nt]= new TGraph(TrigRate[nt].size(), &TrigRate[nt][0], &CompLT[nt][0]);
	gCompLT[nt]->SetMarkerStyle(23);
	gCompLT[nt]->SetTitle(TarLabel[nt]);
	gCompLT[nt]->GetXaxis()->SetTitle("Trig Rate (kHz)");
	gCompLT[nt]->GetYaxis()->SetTitle("Computer LT");
	//
	gEDTMLT[nt]= new TGraph(TrigRate[nt].size(), &TrigRate[nt][0], &EDTMLT[nt][0]);
	gEDTMLT[nt]->SetMarkerStyle(23);
	gEDTMLT[nt]->SetTitle(TarLabel[nt]);
	gEDTMLT[nt]->GetXaxis()->SetTitle("Trig Rate (kHz)");
	gEDTMLT[nt]->GetYaxis()->SetTitle("EDTM LT");
	//
	vector<Double_t > RatCompEDTM;
	for (UInt_t i=0;i<EDTMLT[nt].size();i++) {
	  RatCompEDTM.push_back(CompLT[nt][i]/EDTMLT[nt][i]);
	  }
	gRatCompEDTM[nt]= new TGraph(TrigRate[nt].size(), &TrigRate[nt][0], &RatCompEDTM[0]);
	gRatCompEDTM[nt]->SetMarkerStyle(23);
	gRatCompEDTM[nt]->SetTitle(TarLabel[nt]);
	gRatCompEDTM[nt]->GetXaxis()->SetTitle("Trig Rate (kHz)");
	gRatCompEDTM[nt]->GetYaxis()->SetTitle("Comp LT/EDTM LT");
	RatCompEDTM.clear();
	//
	      }
      //
      TCanvas *cTar[ntar];
      Double_t Ypidslope[ntar];
      Double_t YpidslopeErr[ntar];
      Double_t Ypidint[ntar];
     Double_t YTrackslope[ntar];
     Double_t YTrackslopeErr[ntar];
      Double_t YTrackint[ntar];
     Double_t YElCleanslope[ntar];
     Double_t YElCleanslopeErr[ntar];
      Double_t YElCleanint[ntar];
     Double_t YElRealslope[ntar];
     Double_t YElRealslopeErr[ntar];
      Double_t YElRealint[ntar];
       for (Int_t nt=0;nt<ntar;nt++) {
	 cTar[nt] = new TCanvas(Form("Ctar_%d",nt),Form("Ctarl_%d",nt), 700,700);
	 cTar[nt]->Divide(2,2);
	 cTar[nt]->cd(1);
	 gYPid[nt]->Draw("AP");
	 gYPid[nt]->Fit("pol1","Q","",0,60);
	 TF1 *fit = gYPid[nt]->GetFunction("pol1");
	 Ypidslope[nt] = fit->GetParameter(1);
	 YpidslopeErr[nt] = fit->GetParError(1);
	 Ypidint[nt] = fit->GetParameter(0);
	 Double_t normSlope=Ypidslope[nt]/Ypidint[nt]*100*100;
	 Double_t normSlopeErr= YpidslopeErr[nt]/Ypidint[nt]*100*100;
	 TString ttemp=Form(" Slope (\%/100uA) = %3.2f +/- %3.2f",normSlope,normSlopeErr);
	 gYPid[nt]->SetTitle(TarLabel[nt]+ttemp);
	 cTar[nt]->cd(2);
	 gYTrack[nt]->Draw("AP");
	 gYTrack[nt]->Fit("pol1","Q","",0,60);
	 fit = gYTrack[nt]->GetFunction("pol1");
	 YTrackslope[nt] = fit->GetParameter(1);
	 YTrackslopeErr[nt] = fit->GetParError(1);
	 YTrackint[nt] = fit->GetParameter(0);
	 normSlope=YTrackslope[nt]/YTrackint[nt]*100*100;
	 normSlopeErr= YTrackslopeErr[nt]/YTrackint[nt]*100*100;
	 ttemp=Form(" Slope (\%/100uA) = %3.2f +/- %3.2f",normSlope,normSlopeErr);
	 gYTrack[nt]->SetTitle(TarLabel[nt]+ttemp);
	 cTar[nt]->cd(3);
	 gYElClean[nt]->Draw("AP");
	 gYElClean[nt]->Fit("pol1","Q","",0,60);
	 fit = gYElClean[nt]->GetFunction("pol1");
	 YElCleanslope[nt] = fit->GetParameter(1);
	 YElCleanslopeErr[nt] = fit->GetParError(1);
	 YElCleanint[nt] = fit->GetParameter(0);
	 normSlope=YElCleanslope[nt]/YElCleanint[nt]*100*100;
	 normSlopeErr= YElCleanslopeErr[nt]/YElCleanint[nt]*100*100;
	 ttemp=Form(" Slope (\%/100uA) = %3.2f +/- %3.2f",normSlope,normSlopeErr);
	 gYElClean[nt]->SetTitle(TarLabel[nt]+ttemp);
	 cTar[nt]->cd(4);
	 gYElReal[nt]->Draw("AP");
	 gYElReal[nt]->Fit("pol1","Q","",0,60);
	 fit = gYElReal[nt]->GetFunction("pol1");
	 YElRealslope[nt] = fit->GetParameter(1);
	 YElRealslopeErr[nt] = fit->GetParError(1);
	 YElRealint[nt] = fit->GetParameter(0);
	 normSlope=YElRealslope[nt]/YElRealint[nt]*100*100;
	 normSlopeErr= YElRealslopeErr[nt]/YElRealint[nt]*100*100;
	 ttemp=Form(" Slope (\%/100uA) = %3.2f +/- %3.2f",normSlope,normSlopeErr);
	 gYElReal[nt]->SetTitle(TarLabel[nt]+ttemp);
	 cout << TarLabel[nt] << " Pid slope (%/uA) = " << Ypidslope[nt]/Ypidint[nt]*100*100 << " +/- " << YpidslopeErr[nt]/Ypidint[nt]*100*100 << " Track slope (%/uA) = " << YTrackslope[nt]/YTrackint[nt]*100*100 << " +/- "<< YTrackslopeErr[nt]/YTrackint[nt]*100*100 <<endl;
	 if (nt==0) cTar[nt]->Print(outputpdf+"(");
	 if (nt!=0) cTar[nt]->Print(outputpdf);
       }
       //
       //
       vector<Double_t > CLT_AllRuns;
       vector<Double_t > EDTMLT_AllRuns;
       vector<Double_t > TrigRate_AllRuns;
      TCanvas *cTarRate[ntar];
       for (Int_t nt=0;nt<ntar;nt++) {
	 for (Int_t nr=0;nr<TrigRate[nt].size();nr++) {
	   CLT_AllRuns.push_back(CompLT[nt][nr]);
	   EDTMLT_AllRuns.push_back(EDTMLT[nt][nr]);
	   TrigRate_AllRuns.push_back(TrigRate[nt][nr]);
	 }
	 cTarRate[nt] = new TCanvas(Form("CTarRate_%d",nt),Form("CTarRatel_%d",nt), 700,700);
	 cTarRate[nt]->Divide(2,2);
	 cTarRate[nt]->cd(1);
	 gCompLT[nt]->Draw("AP");
	 cTarRate[nt]->cd(2);
	 gEDTMLT[nt]->Draw("AP");
	 cTarRate[nt]->cd(3);
	 gRatCompEDTM[nt]->Draw("AP");
	 if (nt!=ntar-1) cTarRate[nt]->Print(outputpdf);
	 if (nt==ntar-1) cTarRate[nt]->Print(outputpdf+")");
       }
       //
       TF1* fLT = new TF1("fLT","1./(1+[0]*(x-[1]))");
       TF1* fEDTMLT = new TF1("fEDTMLT","1./(1+[0]*(x-[1]))");
       fEDTMLT->SetLineColor(3);
	TGraph *gCompLT_AllRuns= new TGraph(TrigRate_AllRuns.size(), &TrigRate_AllRuns[0], &CLT_AllRuns[0]);
	gCompLT_AllRuns->SetMarkerStyle(23);
	gCompLT_AllRuns->SetTitle(Form("All %s runs",spec.Data()));
	gCompLT_AllRuns->GetXaxis()->SetTitle("Trig Rate (kHz)");
	gCompLT_AllRuns->GetYaxis()->SetTitle("Computer LT");
	TGraph *gEDTMLT_AllRuns= new TGraph(TrigRate_AllRuns.size(), &TrigRate_AllRuns[0], &EDTMLT_AllRuns[0]);
	gEDTMLT_AllRuns->SetMarkerStyle(22);
	gEDTMLT_AllRuns->SetMarkerColor(3);
	gEDTMLT_AllRuns->SetTitle(Form("All %s runs",spec.Data()));
	gEDTMLT_AllRuns->GetXaxis()->SetTitle("Trig Rate (kHz)");
	gEDTMLT_AllRuns->GetYaxis()->SetTitle("EDTM LT");
      TCanvas *cTarRate_AllRuns;
      cTarRate_AllRuns = new TCanvas("cTarRate_AllRuns","cTarRate_AllRuns",700,700);
      cTarRate_AllRuns->Divide(1,1);
      cTarRate_AllRuns->cd(1);
      gCompLT_AllRuns->Draw("AP");
      gEDTMLT_AllRuns->Draw("P same");
      if (spec.CompareTo("HMS")==0) gCompLT_AllRuns->Fit("fLT","","",5,20);
      if (spec.CompareTo("SHMS")==0) gCompLT_AllRuns->Fit("fLT","","",3,20);
      if (spec.CompareTo("HMS")==0) gEDTMLT_AllRuns->Fit("fEDTMLT","","",5,20);
      if (spec.CompareTo("SHMS")==0) gEDTMLT_AllRuns->Fit("fEDTMLT","","",3,20);
  //
}
