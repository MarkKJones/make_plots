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

void make_hist_shms_hodo_raw(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_hodo_raw_hist.root";
 TObjArray HList(0);
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 static const Int_t plnum=4;
 static const Int_t ns=2;
 static const Int_t clmax=200;
 static const Int_t adcsigtypes=2;
 static const Int_t tdcsigtypes=2;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 static const Int_t npad[plnum]={13,13,14,21};
 const char* sname[ns]={"neg","pos"};
 const char* adcname[adcsigtypes]={"AdcCounter","AdcPulseTime"};
 const char* tdcname[tdcsigtypes]={"TdcCounter","TdcTime"};
 Int_t adc_nhits[plnum][ns];
 Int_t tdc_nhits[plnum][ns];
 Double_t adc[plnum][ns][adcsigtypes][clmax];
 Double_t tdc[plnum][ns][tdcsigtypes][clmax];
 for (Int_t ip=0;ip<plnum;ip++) {
 for (Int_t is=0;is<ns;is++) {
   tsimc->SetBranchAddress(Form("Ndata.P.hod.%s.%s%s",plname[ip],sname[is],adcname[0]),&adc_nhits[ip][is]) ;   
   tsimc->SetBranchAddress(Form("Ndata.P.hod.%s.%s%s",plname[ip],sname[is],tdcname[0]),&tdc_nhits[ip][is]) ;   
 for (Int_t it=0;it<adcsigtypes;it++) {
   tsimc->SetBranchAddress(Form("P.hod.%s.%s%s",plname[ip],sname[is],adcname[it]),&adc[ip][is][it]) ;   
 }
 for (Int_t it=0;it<tdcsigtypes;it++) {
   tsimc->SetBranchAddress(Form("P.hod.%s.%s%s",plname[ip],sname[is],tdcname[it]),&tdc[ip][is][it]) ;   
 }
 }}
 Double_t starttime;
 tsimc->SetBranchAddress("P.hod.starttime",&starttime);
 //
 TH1F* h_nhits[plnum][ns];
TH1F* h_adcpulsetimediff[plnum][ns][16];
TH1F* h_adcpulsetime[plnum][ns][16];
TH1F* h_tdctime[plnum][ns][16];
  for (Int_t ip=0;ip<plnum;ip++) {
 for (Int_t is=0;is<ns;is++) {
   h_nhits[ip][is] = new TH1F(Form("h_nhits_%d_%d",ip,is),Form("; %s %s ",plname[ip],sname[is]),20,0,20);
 for (Int_t ipad=0;ipad<npad[ip];ipad++) {
   h_adcpulsetimediff[ip][is][ipad] = new TH1F(Form("h_adcpulsetimediff_%d_%d_%d",ip,is,ipad),Form("; %s %s Pad %d TDc time - ADC time ",plname[ip],sname[is],ipad+1),200.,-50.,0.);   
   HList.Add(h_adcpulsetimediff[ip][is][ipad]);
   h_adcpulsetime[ip][is][ipad] = new TH1F(Form("h_adcpulsetime_%d_%d_%d",ip,is,ipad),Form("; %s %s Pad %d ADC time",plname[ip],sname[is],ipad+1),300.,50.,150.);   
   HList.Add(h_adcpulsetime[ip][is][ipad]);
   h_tdctime[ip][is][ipad] = new TH1F(Form("h_tdctime_%d_%d_%d",ip,is,ipad),Form("; %s %s Pad %d TDC time",plname[ip],sname[is],ipad+1),200.,10.,110.);   
   HList.Add(h_tdctime[ip][is][ipad]);
 }}}

 //
  Double_t adcpad[plnum][ns][16],tdcpad[plnum][ns][16];
 Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		// look at only 1x pos side
		Int_t ip=0;
		Int_t is=1;
                for (Int_t ip=0;ip<plnum;ip++) {
                for (Int_t is=0;is<ns;is++) {
		for (Int_t np=0;np<16;np++) {
		  adcpad[ip][is][np]=-1000;
		  tdcpad[ip][is][np]=-1000;
		}}}
		//
                for (Int_t ip=0;ip<plnum;ip++) {
                for (Int_t is=0;is<ns;is++) {
		for (Int_t nh=0;nh<adc_nhits[ip][is];nh++) {
		  Int_t np= adc[ip][is][0][nh]-1;
		  adcpad[ip][is][np]=adc[ip][is][1][nh];
		}
		for (Int_t nh=0;nh<tdc_nhits[ip][is];nh++) {
		  Int_t np= tdc[ip][is][0][nh]-1;
		  tdcpad[ip][is][np]=(tdc[ip][is][1][nh]+2000)*0.09776;
		}
		for (Int_t np=0;np<16;np++) {
		  if (adcpad[ip][is][np]!=-1000 && tdcpad[ip][is][np]!=-1000) {
                   h_adcpulsetimediff[ip][is][np]->Fill(tdcpad[ip][is][np]-adcpad[ip][is][np]);
                   h_adcpulsetime[ip][is][np]->Fill(adcpad[ip][is][np]);
                   h_tdctime[ip][is][np]->Fill(tdcpad[ip][is][np]);
		  }
		}
		//
		}}
	}
//
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();

//
}
