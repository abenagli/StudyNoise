{

//   TFile* fileDAMC = TFile::Open("outPlots_DAMC.root","read");
//   TFile* fileMCMC = TFile::Open("outPlots_MCMC.root","read");

  TFile* fileDAMC = TFile::Open("outPlots_DAMC.root","read");
  TFile* fileMCMC = TFile::Open("outPlots_MCMC.root","read");

  TFile* fileAX = TFile::Open("alexey.root", "read");

  TFile* fileNG = TFile::Open("nuGunGraphs_MC-NuGun-runAB-noPU.root", "read");
  TGraph* gNuGun_recHitE_RMS_vsIRing = (TGraph*)( fileNG->Get("g_recHitE_RMS_vsIRing") );
  gNuGun_recHitE_RMS_vsIRing->SetLineColor(kBlue);
  gNuGun_recHitE_RMS_vsIRing->SetLineWidth(2);


  TGraphErrors* g_recHitE_sigma_vsAbsEta_DA = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigma_vsAbsEta_DA");
  g_recHitE_sigma_vsAbsEta_DA->SetMarkerColor(kBlack);
  g_recHitE_sigma_vsAbsEta_DA->SetMarkerStyle(20);
  TGraphErrors* g_recHitE_sigmaAt0PU_vsAbsEta_DA = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigmaAt0PU_vsAbsEta_DA");
  g_recHitE_sigmaAt0PU_vsAbsEta_DA->SetMarkerColor(kBlack);
  g_recHitE_sigmaAt0PU_vsAbsEta_DA->SetMarkerStyle(28);
  TGraphErrors* g_recHitE_RMSAt0PU_vsIRing_DA = (TGraphErrors*)fileDAMC->Get("g_recHitE_RMSAt0PU_vsIRing_DA");
  g_recHitE_RMSAt0PU_vsIRing_DA->SetMarkerColor(kBlack);
  g_recHitE_RMSAt0PU_vsIRing_DA->SetMarkerStyle(28);
  TGraphErrors* g_recHitE_sigma_vsAbsEtans_DA = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigma_vsAbsEtans_DA");
  g_recHitE_sigma_vsAbsEtans_DA->SetMarkerColor(kBlack);
  g_recHitE_sigma_vsAbsEtans_DA->SetMarkerStyle(3);
  TGraphErrors* g_recHitE_sigmaAt70PU_vsAbsEtans_DA = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigmaAt70PU_vsAbsEtans_DA");
  g_recHitE_sigmaAt70PU_vsAbsEtans_DA->SetMarkerColor(kRed);
  g_recHitE_sigmaAt70PU_vsAbsEtans_DA->SetMarkerStyle(21);
  g_recHitE_sigmaAt70PU_vsAbsEtans_DA->SetMarkerSize(0.7);
  TGraphErrors* g_recHitE_sigmaAt140PU_vsAbsEtans_DA = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigmaAt140PU_vsAbsEtans_DA");
  g_recHitE_sigmaAt140PU_vsAbsEtans_DA->SetMarkerColor(kBlue);
  g_recHitE_sigmaAt140PU_vsAbsEtans_DA->SetMarkerStyle(21);
  g_recHitE_sigmaAt140PU_vsAbsEtans_DA->SetMarkerSize(0.7);

  TGraphErrors* g_recHitE_sigma_vsAbsEta_MC = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigma_vsAbsEta_MC");
  g_recHitE_sigma_vsAbsEta_MC->SetMarkerColor(kGreen+2);
  g_recHitE_sigma_vsAbsEta_MC->SetMarkerStyle(20);
  TGraphErrors* g_recHitE_sigmaAt0PU_vsAbsEta_MC = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigmaAt0PU_vsAbsEta_MC");
  g_recHitE_sigmaAt0PU_vsAbsEta_MC->SetMarkerColor(kGreen+2);
  g_recHitE_sigmaAt0PU_vsAbsEta_MC->SetMarkerStyle(28);
  TGraphErrors* g_recHitE_sigma_vsAbsEtans_MC = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigma_vsAbsEtans_MC");
  g_recHitE_sigma_vsAbsEtans_MC->SetMarkerColor(kGreen+2);
  g_recHitE_sigma_vsAbsEtans_MC->SetMarkerStyle(3);
  TGraphErrors* g_recHitE_sigmaAt70PU_vsAbsEtans_MC = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigmaAt70PU_vsAbsEtans_MC");
  g_recHitE_sigmaAt70PU_vsAbsEtans_MC->SetMarkerColor(kRed);
  g_recHitE_sigmaAt70PU_vsAbsEtans_MC->SetMarkerStyle(25);
  g_recHitE_sigmaAt70PU_vsAbsEtans_MC->SetMarkerSize(0.7);
  TGraphErrors* g_recHitE_sigmaAt140PU_vsAbsEtans_MC = (TGraphErrors*)fileDAMC->Get("g_recHitE_sigmaAt140PU_vsAbsEtans_MC");
  g_recHitE_sigmaAt140PU_vsAbsEtans_MC->SetMarkerColor(kBlue);
  g_recHitE_sigmaAt140PU_vsAbsEtans_MC->SetMarkerStyle(25);
  g_recHitE_sigmaAt140PU_vsAbsEtans_MC->SetMarkerSize(0.7);

  TGraphErrors* g_recHitE_sigma_vsAbsEta_PU0_MC = (TGraphErrors*)fileMCMC->Get("g_recHitE_sigma_vsAbsEta_PU0_MC");
  g_recHitE_sigma_vsAbsEta_PU0_MC->SetMarkerColor(kBlack);
  g_recHitE_sigma_vsAbsEta_PU0_MC->SetMarkerStyle(20);
  TGraphErrors* g_recHitE_sigmaAt0PU_vsAbsEta_070140MC = (TGraphErrors*)fileMCMC->Get("g_recHitE_sigmaAt0PU_vsAbsEta_MC");
  g_recHitE_sigmaAt0PU_vsAbsEta_070140MC->SetMarkerColor(kBlack);
  g_recHitE_sigmaAt0PU_vsAbsEta_070140MC->SetMarkerStyle(20);
  TGraphErrors* g_recHitE_RMSAt0PU_vsIRing_MC = (TGraphErrors*)fileMCMC->Get("g_recHitE_RMSAt0PU_vsIRing_MC");
  g_recHitE_RMSAt0PU_vsIRing_MC->SetMarkerColor(kBlack);
  g_recHitE_RMSAt0PU_vsIRing_MC->SetMarkerStyle(20);

  TGraphErrors* g_recHitE_sigma_vsAbsEta_PU70_MC = (TGraphErrors*)fileMCMC->Get("g_recHitE_sigma_vsAbsEta_PU70_MC");
  g_recHitE_sigma_vsAbsEta_PU70_MC->SetMarkerColor(kRed);
  g_recHitE_sigma_vsAbsEta_PU70_MC->SetMarkerStyle(20);
  TGraphErrors* g_recHitE_sigma_vsAbsEta_PU70ns_MC = (TGraphErrors*)fileMCMC->Get("g_recHitE_sigma_vsAbsEta_PU70ns_MC");
  g_recHitE_sigma_vsAbsEta_PU70ns_MC->SetMarkerColor(kRed);
  g_recHitE_sigma_vsAbsEta_PU70ns_MC->SetMarkerStyle(3);

  TGraphErrors* g_recHitE_sigma_vsAbsEta_PU140_MC = (TGraphErrors*)fileMCMC->Get("g_recHitE_sigma_vsAbsEta_PU140_MC");
  g_recHitE_sigma_vsAbsEta_PU140_MC->SetMarkerColor(kBlue);
  g_recHitE_sigma_vsAbsEta_PU140_MC->SetMarkerStyle(20);

  TGraphErrors* g_recHitE_sigma_vsAbsEta_PU140ns_MC = (TGraphErrors*)fileMCMC->Get("g_recHitE_sigma_vsAbsEta_PU140ns_MC");
  g_recHitE_sigma_vsAbsEta_PU140ns_MC->SetMarkerColor(kBlue);
  g_recHitE_sigma_vsAbsEta_PU140ns_MC->SetMarkerStyle(3);


  TGraphErrors* ax_PU20 = (TGraphErrors*)fileAX->Get("g_20");
  ax_PU20->SetMarkerColor(kBlack);
  ax_PU20->SetLineColor(kBlack);
  ax_PU20->SetLineWidth(2);

  TGraphErrors* ax_PU20_25ns = new TGraphErrors();
  ax_PU20_25ns->SetMarkerColor(kBlack);
  ax_PU20_25ns->SetLineColor(kBlack);
  ax_PU20_25ns->SetLineWidth(2);
  ax_PU20_25ns->SetLineStyle(2);

  ax_PU20_25ns->SetPoint(0, 0.0722, 0.00869);
  ax_PU20_25ns->SetPoint(1, 0.217, 0.00908);
  ax_PU20_25ns->SetPoint(2, 0.361, 0.00804);
  ax_PU20_25ns->SetPoint(3, 0.505, 0.00803);
  ax_PU20_25ns->SetPoint(4, 0.65, 0.00834);
  ax_PU20_25ns->SetPoint(5, 0.794, 0.0114);
  ax_PU20_25ns->SetPoint(6, 0.939, 0.0123);
  ax_PU20_25ns->SetPoint(7, 1.08, 0.0128);
  ax_PU20_25ns->SetPoint(8, 1.23, 0.0161);
  ax_PU20_25ns->SetPoint(9, 1.37, 0.0179);
  ax_PU20_25ns->SetPoint(10, 1.64, 0.0249);
  ax_PU20_25ns->SetPoint(11, 1.78, 0.0406);
  ax_PU20_25ns->SetPoint(12, 1.92, 0.0577);
  ax_PU20_25ns->SetPoint(13, 2.07, 0.0888);
  ax_PU20_25ns->SetPoint(14, 2.21, 0.148);
  ax_PU20_25ns->SetPoint(15, 2.35, 0.211);
  ax_PU20_25ns->SetPoint(16, 2.5, 0.308);
  ax_PU20_25ns->SetPoint(17, 2.64, 0.559);
  ax_PU20_25ns->SetPoint(18, 2.78, 0.792);
  ax_PU20_25ns->SetPoint(19, 2.93, 1.17);



  TGraphErrors* ax_PU70 = new TGraphErrors();
  ax_PU70->SetMarkerColor(kRed);
  ax_PU70->SetLineColor(kRed);
  ax_PU70->SetLineWidth(2);
  ax_PU70->SetMarkerStyle(20);

  ax_PU70->SetPoint(0,0.07220000029,0.03840000182);
  ax_PU70->SetPoint(1,0.2169999927,0.0335999988);
  ax_PU70->SetPoint(2,0.3610000014,0.04060000181);
  ax_PU70->SetPoint(3,0.5049999952,0.0447999984);
  ax_PU70->SetPoint(4,0.6499999762,0.03959999979);
  ax_PU70->SetPoint(5,0.7940000296,0.03889999911);
  ax_PU70->SetPoint(6,0.9390000105,0.04179999977);
  ax_PU70->SetPoint(7,1.080000043,0.05449999869);
  ax_PU70->SetPoint(8,1.230000019,0.05860000104);
  ax_PU70->SetPoint(9,1.370000005,0.06419999897);
  ax_PU70->SetPoint(10,1.639999986,0.09700000286);
  ax_PU70->SetPoint(11,1.779999971,0.1220000014);
  ax_PU70->SetPoint(12,1.919999957,0.1790000051);
  ax_PU70->SetPoint(13,2.069999933,0.2619999945);
  ax_PU70->SetPoint(14,2.210000038,0.3849999905);
  ax_PU70->SetPoint(15,2.349999905,0.5500000119);
  ax_PU70->SetPoint(16,2.5,0.7979999781);
  ax_PU70->SetPoint(17,2.640000105,1.460000038);
  ax_PU70->SetPoint(18,2.779999971,2.039999962);
  ax_PU70->SetPoint(19,2.930000067,2.710000038);


  TGraphErrors* ax_PU140 = new TGraphErrors();
  ax_PU140->SetMarkerColor(kBlue);
  ax_PU140->SetLineColor(kBlue);
  ax_PU140->SetLineWidth(2);
  ax_PU140->SetMarkerStyle(20);

  ax_PU140->SetPoint(0,0.07220000029,0.07190000266);
  ax_PU140->SetPoint(1,0.2169999927,0.07959999889);
  ax_PU140->SetPoint(2,0.3610000014,0.06960000098);
  ax_PU140->SetPoint(3,0.5049999952,0.07209999859);
  ax_PU140->SetPoint(4,0.6499999762,0.07079999894);
  ax_PU140->SetPoint(5,0.7940000296,0.08380000293);
  ax_PU140->SetPoint(6,0.9390000105,0.09229999781);
  ax_PU140->SetPoint(7,1.080000043,0.09390000254);
  ax_PU140->SetPoint(8,1.230000019,0.1030000001);
  ax_PU140->SetPoint(9,1.370000005,0.1080000028);
  ax_PU140->SetPoint(10,1.639999986,0.1490000039);
  ax_PU140->SetPoint(11,1.779999971,0.2329999954);
  ax_PU140->SetPoint(12,1.919999957,0.3120000064);
  ax_PU140->SetPoint(13,2.069999933,0.4269999862);
  ax_PU140->SetPoint(14,2.210000038,0.69599998);
  ax_PU140->SetPoint(15,2.349999905,0.9520000219);
  ax_PU140->SetPoint(16,2.5,1.320000052);
  ax_PU140->SetPoint(17,2.640000105,2.130000114);
  ax_PU140->SetPoint(18,2.779999971,3.279999971);
  ax_PU140->SetPoint(19,2.930000067,4.329999924);

  //// DRAW
  TLegend* leg = new TLegend(0.15,0.70,0.45,0.95,NULL,"brNDC");
  leg->SetTextSize(0.035);
  leg->SetLineWidth(0);
  leg->SetLineColor(kWhite);
  leg->SetTextFont(42);
  leg->SetFillColor(kWhite);
  //  leg->SetFillStyle(1001);
  //  leg->SetBorderSize(1);





  ////c_AxZeeDa
  TCanvas* c_AxZeeDa = new TCanvas();
  c_AxZeeDa->SetLogy();
  c_AxZeeDa->SetGridx();
  c_AxZeeDa->SetGridy();

  ax_PU140->Draw("AL");
  g_recHitE_sigma_vsAbsEtans_DA->Draw("p, same");
  g_recHitE_sigma_vsAbsEta_PU70ns_MC->Draw("p,same");
  g_recHitE_sigma_vsAbsEta_PU140ns_MC->Draw("p,same");
  g_recHitE_sigmaAt70PU_vsAbsEtans_DA->Draw("p,same");
  g_recHitE_sigmaAt70PU_vsAbsEtans_MC->Draw("p,same");
  g_recHitE_sigmaAt140PU_vsAbsEtans_DA->Draw("p,same");
  g_recHitE_sigmaAt140PU_vsAbsEtans_MC->Draw("p,same");

  ax_PU70->Draw("l,same");
  ax_PU20->Draw("l,same");

  ax_PU140->GetXaxis()->SetTitle("|#eta|");
  ax_PU140->GetYaxis()->SetTitle("effective #sigma(E_{recHit}) [GeV]");
  ax_PU140->GetXaxis()->SetRangeUser(0., 3.);
  ax_PU140->SetMinimum(0.001);
  ax_PU140->SetMaximum(10.);

  leg->AddEntry(g_recHitE_sigma_vsAbsEta_PU140ns_MC, "MC PU140 (PU only)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEta_PU70ns_MC, "MC PU70 (PU only)", "p");
  leg->AddEntry(ax_PU140, "MC PU140 (alexey)", "l");
  leg->AddEntry(ax_PU70, "MC PU70 (alexey)", "l");
  leg->AddEntry(ax_PU20, "MC PU20 (alexey)", "l");
  leg->AddEntry(g_recHitE_sigma_vsAbsEtans_DA, "DA RunD (PU only)", "p");
  leg->AddEntry(g_recHitE_sigmaAt70PU_vsAbsEtans_DA, "DA atPU70 (PU only)", "p");
  leg->AddEntry(g_recHitE_sigmaAt70PU_vsAbsEtans_MC, "MC atPU70 (PU only)", "p");
  leg->AddEntry(g_recHitE_sigmaAt140PU_vsAbsEtans_DA, "DA atPU140 (PU only)", "p");
  leg->AddEntry(g_recHitE_sigmaAt140PU_vsAbsEtans_MC, "MC atPU140 (PU only)", "p");

  leg->Draw("same");
  c_AxZeeDa->Print("PLOTS/c_AxZeeDa.png", "png");

  leg->Clear();


  ////c_ZeeDa_noise
  TCanvas* c_ZeeDa_noise = new TCanvas();
  c_ZeeDa_noise->SetLogy();
  c_ZeeDa_noise->SetGridx();
  c_ZeeDa_noise->SetGridy();

  g_recHitE_sigma_vsAbsEta_PU0_MC->Draw("AP");
  g_recHitE_sigma_vsAbsEta_PU70_MC->Draw("p,same");
  g_recHitE_sigma_vsAbsEta_PU140_MC->Draw("p,same");
  g_recHitE_sigma_vsAbsEta_MC->Draw("p,same");
  g_recHitE_sigma_vsAbsEta_PU70ns_MC->Draw("p,same");
  g_recHitE_sigma_vsAbsEta_PU140ns_MC->Draw("p,same");
  g_recHitE_sigma_vsAbsEtans_MC->Draw("p,same");
  //  g_recHitE_sigma_vsAbsEtans_DA->Draw("p,same");

  g_recHitE_sigma_vsAbsEta_PU0_MC->GetXaxis()->SetTitle("|#eta|");
  g_recHitE_sigma_vsAbsEta_PU0_MC->GetYaxis()->SetTitle("effective #sigma(E_{recHit}) [GeV]");
  g_recHitE_sigma_vsAbsEta_PU0_MC->GetXaxis()->SetRangeUser(0., 3.);
  g_recHitE_sigma_vsAbsEta_PU0_MC->SetMinimum(0.01);
  g_recHitE_sigma_vsAbsEta_PU0_MC->SetMaximum(3.);

  leg->AddEntry(g_recHitE_sigma_vsAbsEta_PU140_MC, "MC PU140 (noise 0ageing)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEta_PU140ns_MC, "MC PU140 (PU only)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEta_PU70_MC, "MC PU70 (noise 0ageing)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEta_PU70ns_MC, "MC PU70 (PU only)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEta_MC, "MC PUrunD (noise runD)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEtans_MC, "MC PUrunD (PU only)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEta_PU0_MC, "MC PU0 (noise 0ageing)", "p");
  leg->Draw("same");
  c_ZeeDa_noise->Print("PLOTS/c_ZeeDa_noise.png","png");


  leg->Clear();


  ////c_At0PU
  TCanvas* c_At0PU = new TCanvas();
  c_At0PU->SetLogy();
  c_At0PU->SetGridx();
  c_At0PU->SetGridy();

  g_recHitE_sigma_vsAbsEta_MC->Draw("AP,same");
  g_recHitE_sigma_vsAbsEta_DA->Draw("p,same");
  g_recHitE_sigmaAt0PU_vsAbsEta_MC->Draw("p,same");
  g_recHitE_sigmaAt0PU_vsAbsEta_DA->Draw("p,same");

  g_recHitE_sigma_vsAbsEta_MC->GetXaxis()->SetTitle("|#eta|");
  g_recHitE_sigma_vsAbsEta_MC->GetYaxis()->SetTitle("effective #sigma(E_{recHit}) [GeV]");
  g_recHitE_sigma_vsAbsEta_MC->GetXaxis()->SetRangeUser(0., 3.);
  g_recHitE_sigma_vsAbsEta_MC->SetMinimum(0.04);
  g_recHitE_sigma_vsAbsEta_MC->SetMaximum(0.5);

  leg->AddEntry(g_recHitE_sigmaAt0PU_vsAbsEta_DA, "DA at0PU (noise runD)", "p");
  leg->AddEntry(g_recHitE_sigmaAt0PU_vsAbsEta_MC, "MC at0PU (noise runD)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEta_DA, "DA PUrunD (noise runD)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEta_MC, "MC PUrunD (noise runD)", "p");
  leg->Draw("same");
  c_At0PU->Print("PLOTS/c_At0PU.png","png");

  leg->Clear();


//   leg->SetNColumns(2);
//   leg->SetColumnSeparation(1.);
//   std::cout << " separation = " << leg->GetEntrySeparation() << std::endl;
  leg->SetEntrySeparation(0.01);
  leg->SetTextSize(0.025);
  ////c_AxZeeDa
  TCanvas* c_AxZeeDaOK = new TCanvas();
  c_AxZeeDaOK->SetLogy();
  c_AxZeeDaOK->SetGridx();
  c_AxZeeDaOK->SetGridy();

  g_recHitE_sigma_vsAbsEtans_DA->SetMarkerStyle(20);
  g_recHitE_sigma_vsAbsEtans_MC->SetMarkerStyle(20);
  g_recHitE_sigma_vsAbsEta_PU70ns_MC->SetMarkerStyle(20);
  g_recHitE_sigma_vsAbsEta_PU140ns_MC->SetMarkerStyle(20);

  ax_PU140->Draw("AL");
  g_recHitE_sigma_vsAbsEtans_DA->Draw("p, same");
  g_recHitE_sigma_vsAbsEtans_MC->Draw("p, same");
  g_recHitE_sigma_vsAbsEta_PU70ns_MC->Draw("p,same");
  g_recHitE_sigma_vsAbsEta_PU140ns_MC->Draw("p,same");
  ax_PU70->Draw("l,same");
  ax_PU20->Draw("l,same");
  ax_PU20_25ns->Draw("l,same");

  ax_PU140->GetXaxis()->SetTitle("|#eta|");
  ax_PU140->GetYaxis()->SetTitle("effective #sigma(E_{recHit}) [GeV]");
  ax_PU140->GetXaxis()->SetRangeUser(0., 3.);
  ax_PU140->SetMinimum(0.001);
  ax_PU140->SetMaximum(10.);

  leg->AddEntry(g_recHitE_sigma_vsAbsEta_PU140ns_MC, "MC PU140 (PU only)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEta_PU70ns_MC, "MC PU70 (PU only)", "p");
  leg->AddEntry(ax_PU140, "MC PU140 (alexey)", "l");
  leg->AddEntry(ax_PU70, "MC PU70 (alexey)", "l");
  leg->AddEntry(ax_PU20, "MC PU20 BX 50ns (alexey)", "l");
  leg->AddEntry(ax_PU20_25ns, "MC PU20 BX 25ns (alexey)", "l");
  leg->AddEntry(g_recHitE_sigma_vsAbsEtans_DA, "DA RunD (PU only)", "p");
  leg->AddEntry(g_recHitE_sigma_vsAbsEtans_MC, "MC RunD (PU only)", "p");

  leg->Draw("same");
  c_AxZeeDaOK->Print("PLOTS/c_AxZeeDaOK.png", "png");

  leg->Clear();



  TGraphErrors*  gEB_recHitE_sigma_vsPUrunD_MC = (TGraphErrors*)fileDAMC->Get("gEB_recHitE_sigma_vsNPU_MC");
  TGraphErrors*  gEB_recHitE_sigma_vsPUother_MC = (TGraphErrors*)fileMCMC->Get("gEB_recHitE_sigma_vsNPU_MC");
  TGraphErrors*  gEE_recHitE_sigma_vsPUrunD_MC = (TGraphErrors*)fileDAMC->Get("gEE_recHitE_sigma_vsNPU_MC");
  TGraphErrors*  gEE_recHitE_sigma_vsPUother_MC = (TGraphErrors*)fileMCMC->Get("gEE_recHitE_sigma_vsNPU_MC");


  gEB_recHitE_sigma_vsPUrunD_MC->SetMarkerColor(kGreen+2);
  gEB_recHitE_sigma_vsPUrunD_MC->SetMarkerStyle(20);
  gEE_recHitE_sigma_vsPUrunD_MC->SetMarkerColor(kGreen+2);
  gEE_recHitE_sigma_vsPUrunD_MC->SetMarkerStyle(3);

  gEB_recHitE_sigma_vsPUother_MC->SetMarkerColor(kRed+2);
  gEB_recHitE_sigma_vsPUother_MC->SetMarkerStyle(20);
  gEE_recHitE_sigma_vsPUother_MC->SetMarkerColor(kRed+2);
  gEE_recHitE_sigma_vsPUother_MC->SetMarkerStyle(3);

  TCanvas* cTrends_sigma = new TCanvas();
  cTrends_sigma->SetGridx();
  cTrends_sigma->SetGridy();

  gEB_recHitE_sigma_vsPUrunD_MC->Draw("ap");
  gEE_recHitE_sigma_vsPUrunD_MC->Draw("p, same");
  gEB_recHitE_sigma_vsPUother_MC->Draw("p, same");
  gEE_recHitE_sigma_vsPUother_MC->Draw("p, same");

  gEB_recHitE_sigma_vsPUrunD_MC->GetXaxis()->SetTitle("<nPU>");
  gEB_recHitE_sigma_vsPUrunD_MC->GetYaxis()->SetTitle("effective #sigma [GeV]");
  gEB_recHitE_sigma_vsPUrunD_MC->GetXaxis()->SetRangeUser(5., 150.);
  gEB_recHitE_sigma_vsPUrunD_MC->SetMinimum(0.001);
  gEB_recHitE_sigma_vsPUrunD_MC->SetMaximum(1.);

  leg->AddEntry(gEB_recHitE_sigma_vsPUrunD_MC, "MC runD EB", "p");
  leg->AddEntry(gEE_recHitE_sigma_vsPUrunD_MC, "MC runD EE", "p");
  leg->AddEntry(gEB_recHitE_sigma_vsPUother_MC, "MC EB", "p");
  leg->AddEntry(gEE_recHitE_sigma_vsPUother_MC, "MC EE", "p");

  leg->Draw("same");
  cTrends_sigma->Print("PLOTS/cTrends_sigma.png", "png");



  TGraphErrors*  gEB_recHitE_RMS_vsPUrunD_MC = (TGraphErrors*)fileDAMC->Get("gEB_recHitE_RMS_vsNPU_MC");
  TGraphErrors*  gEB_recHitE_RMS_vsPUother_MC = (TGraphErrors*)fileMCMC->Get("gEB_recHitE_RMS_vsNPU_MC");
  TGraphErrors*  gEE_recHitE_RMS_vsPUrunD_MC = (TGraphErrors*)fileDAMC->Get("gEE_recHitE_RMS_vsNPU_MC");
  TGraphErrors*  gEE_recHitE_RMS_vsPUother_MC = (TGraphErrors*)fileMCMC->Get("gEE_recHitE_RMS_vsNPU_MC");


  gEB_recHitE_RMS_vsPUrunD_MC->SetMarkerColor(kGreen+2);
  gEB_recHitE_RMS_vsPUrunD_MC->SetMarkerStyle(20);
  gEE_recHitE_RMS_vsPUrunD_MC->SetMarkerColor(kGreen+2);
  gEE_recHitE_RMS_vsPUrunD_MC->SetMarkerStyle(3);

  gEB_recHitE_RMS_vsPUother_MC->SetMarkerColor(kRed+2);
  gEB_recHitE_RMS_vsPUother_MC->SetMarkerStyle(20);
  gEE_recHitE_RMS_vsPUother_MC->SetMarkerColor(kRed+2);
  gEE_recHitE_RMS_vsPUother_MC->SetMarkerStyle(3);

  TCanvas* cTrends_RMS = new TCanvas();
  cTrends_RMS->SetGridx();
  cTrends_RMS->SetGridy();

  gEB_recHitE_RMS_vsPUrunD_MC->Draw("ap");
  gEE_recHitE_RMS_vsPUrunD_MC->Draw("p, same");
  gEB_recHitE_RMS_vsPUother_MC->Draw("p, same");
  gEE_recHitE_RMS_vsPUother_MC->Draw("p, same");

  gEB_recHitE_RMS_vsPUrunD_MC->GetXaxis()->SetTitle("<nPU>");
  gEB_recHitE_RMS_vsPUrunD_MC->GetYaxis()->SetTitle("RMS [GeV]");
  gEB_recHitE_RMS_vsPUrunD_MC->GetXaxis()->SetRangeUser(5., 150.);
  gEB_recHitE_RMS_vsPUrunD_MC->SetMinimum(0.001);
  gEB_recHitE_RMS_vsPUrunD_MC->SetMaximum(1.);

//   leg->AddEntry(gEB_recHitE_RMS_vsPUrunD_MC, "MC runD EB", "p");
//   leg->AddEntry(gEE_recHitE_RMS_vsPUrunD_MC, "MC runD EE", "p");
//   leg->AddEntry(gEB_recHitE_RMS_vsPUother_MC, "MC EB", "p");
//   leg->AddEntry(gEE_recHitE_RMS_vsPUother_MC, "MC EE", "p");

  leg->Draw("same");
  cTrends_RMS->Print("PLOTS/cTrends_RMS.png", "png");

}


