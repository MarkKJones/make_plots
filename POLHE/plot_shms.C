// supposed to be at 11.7 degree and for one-pass (2.185 GeV)
#include <fmt/format.h>

std::string data_fmt = "/home/cdaq/polhe3/production/hallc_replay/ROOTfiles/shms_replay_production_default_{:05d}_-1.root";

void plot_shms(int run = 2555)
{
    gStyle->SetPalette(kBird, 0, 0.5);
    gStyle->SetHistFillStyle(3001);
    // gStyle->SetOptStat(10);

    TFile *f = new TFile(fmt::format(data_fmt, run).c_str());
    TTree *T = static_cast<TTree*>(f->Get("T"));

    auto c = new TCanvas("he3_dis", "DIS Events", 1600, 900);
    c->Divide(2, 2);
    
    std::string cut_win = fmt::format("H.react.z < {:.1f} && H.react.z > {:.1f}", z_max, z_min);
    auto combine = [](const std::vector<std::string> &all_cuts) {
        if (all_cuts.empty()) { return std::string(); }
        std::string res = all_cuts[0];
        for (size_t i = 1; i < all_cuts.size(); ++i) {
            res += "&&" + all_cuts[i];
        }
        return res;
    };

    
    std::string cut_desc = fmt::format(" (Window Cut from {:.1f} to {:.1f}, Elastic Cut)", z_min, z_max);

    c->cd(1);
    auto hcal_cer = new TH2F("hcal_cer", "Cal E/p vs Cerenkov; NG cerenkov; E/p", 200, 0, 60, 75, 0, 1.5);
    T->Draw("P.cal.etottracknorm: >> hwcal_cer", "", "");
    hcal_cer->Draw("colz");
    
 
}

