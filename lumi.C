{

  TCanvas *can1 = new TCanvas("can1", "HMS Cer", 800, 800);
  can1->Divide(1, 2);
  can1->cd(1);
  // gPad->SetStatOpt(1111111);
  T->Draw("H.cer.npeSum","H.cer.npeSum>2.0&&H.cer.npeSum<20.");
  can1->cd(2);
  // gPad->SetStatOpt(1111111);
  T->Draw("H.rb.raster.fr_ya","H.cer.npeSum>2.0&&H.cer.npeSum<20.");
}
