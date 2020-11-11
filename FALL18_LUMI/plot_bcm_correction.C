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

void plot_bcm_correction(){
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(0);
 gStyle->SetTitleOffset(1.2,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.16);
  TF1 *Pb_corr = new TF1("pb_corr","1.0+0.045*( log(60)-log(x))/(log(60)-log(2))",1,60.);
  TF1 *Pb_corr2 = new TF1("pb_corr2","1.0+0.01*(x-60)/25.",60.,80.);
  TF1 *dm_corr = new TF1("dm_corr","1+0.00175*(x-50)*(x-50)/x",50.,80.);
  TF1 *offp_corr = new TF1("oofp_corr","1+253./5712./x",1.,60.);
  TF1 *offp_corr2 = new TF1("oofp_corr2","1+2*253./5712./x",1.,60.);
  TGraph *pb_gr = new TGraph(Pb_corr);
  TGraph *pb_gr2 = new TGraph(Pb_corr2);
  TGraph *dm_gr = new TGraph(dm_corr);
  TGraph *offp_gr = new TGraph(offp_corr);
  TGraph *offp_gr2 = new TGraph(offp_corr2);
  TMultiGraph *mgr=new TMultiGraph();
  mgr->SetTitle("Compare BMC corrections;I uA;Correction ");
  TCanvas *cCan = new TCanvas("cCan","BCM comp",700,700);
  TLegend *cLeg = new TLegend(.45,.55,.99,.85,"");
  cCan->Divide(1,1);
  cCan->cd(1);
  mgr->Add(pb_gr);
  mgr->Add(pb_gr2);
  mgr->Add(dm_gr);
  mgr->Add(offp_gr);
  mgr->Add(offp_gr2);
  pb_gr2->SetLineColor(3);
  cLeg->AddEntry(pb_gr,"Peter Bosted Log");
  cLeg->AddEntry(pb_gr2,"Peter Bosted Saturation");
  cLeg->AddEntry(dm_gr,"Dave Mack Saturation");
  cLeg->AddEntry(offp_gr,"Offset off by 0.1%");
  cLeg->AddEntry(offp_gr2,"Offset off by 0.2%");
  dm_gr->SetLineColor(4);
  offp_gr->SetLineColor(6);
  offp_gr2->SetLineColor(7);
  mgr->Draw("AL");
  mgr->GetHistogram()->GetXaxis()->SetRangeUser(0.,80.);
  cLeg->Draw();
}
