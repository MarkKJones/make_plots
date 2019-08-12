void make_chain() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 gStyle->SetPadLeftMargin(0.17);
 // CHain files together
    TString inputroot;
    TChain *fdelta;
    fdelta = new TChain("T");
        static const Int_t nrtot=6;
     Int_t runlist[nrtot]={1815,1816,1819,1820,1821,1822};
    //    static const Int_t nrtot=3;
    //   Int_t runlist[nrtot]={1817,1818,1823};
    //     static const Int_t nrtot=5;
    //	Int_t runlist[nrtot]={1914,1915,1916,1917,1918};
    //   static const Int_t nrtot=8;
    //	 Int_t runlist[nrtot]={1817,1818,1823,1914,1915,1916,1917,1918};
   for (Int_t nr=0;nr<nrtot;nr++) {
      TString inputroot=Form("shms_replay_matrixopt_%d_-1.root",runlist[nr]);
      if (gSystem->FindFile("ROOTfiles",inputroot)) {
		cout << " file = " << inputroot << endl;
		//	TFile fch = new TFile("hist/"+inputroot);
       	fdelta->Add(inputroot);	       
      } else {
	cout << " No file = " << inputroot << endl;
      }
    }
    fdelta->Merge(Form("ROOTfiles/shms_replay_matrixopt_%d_%d_-1.root",runlist[0],runlist[nrtot-1]));
}
