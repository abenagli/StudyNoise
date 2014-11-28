// g++  -o drawNoisePlots `root-config --cflags --glibs` drawNoisePlots.cpp

#include "setTDRStyle.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TArrow.h"
#include "TLegend.h"

// std::string label  = "Upgrade_NoAging_PU0-140";
// std::string type = "MCMC";

std::string label  = "runABCD_DA-MC_added";
std::string type = "DAMC";

std::string label2 = "MC-NuGun-runAB-noPU";


int rebin = 120;

void SetHistoStyle(TH1F* h, const std::string& dataset, const std::string& region, const std::string& type, TLegend** leg, const bool& log = false);
void SetGraphStyle(TGraphErrors* g, const std::string& dataset, const std::string& region, const std::string& type, const std::string& type2, TLegend** leg);
float GetEffectiveSigma(TH1F* h);
void GetEffectiveSigma(TH1F* h, double& E1, double& E2);
void FillHistoRMS(TH1F* histo, float rms, float mean);

int main()
{
  setTDRStyle();
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  std::string baseDir = std::string(getenv("STUDYNOISE"));
  TFile* f2 = TFile::Open((baseDir+"/output/nuGunGraphs/nuGunGraphs_"+label2+".root").c_str(),"READ");
  TGraph* gNuGun_recHitE_RMS_vsIRing = (TGraph*)( f2->Get("g_recHitE_RMS_vsIRing") );
  gNuGun_recHitE_RMS_vsIRing -> SetLineColor(kBlue);
  gNuGun_recHitE_RMS_vsIRing -> SetLineWidth(2);
  
  TFile* ax = TFile::Open("alexey.root","read");
  TGraphErrors* ax_PU70 = (TGraphErrors*)ax->Get("g_70");
  TGraphErrors* ax_PU140 = (TGraphErrors*)ax->Get("g_140");

  TFile* f  = TFile::Open((baseDir+"/output/studyNoise_" +label +".root").c_str(),"READ");
  f->cd();

  TFile* fAB  = TFile::Open((baseDir+"/output/studyNoise_runAB_DA-MC.root").c_str(), "READ");
  fAB->cd();

  TFile* fC  = TFile::Open((baseDir+"/output/studyNoise_runC_DA-MC.root").c_str(), "READ");
  fC->cd();

  TFile* fD  = TFile::Open((baseDir+"/output/studyNoise_runD_DA-MC.root").c_str(), "READ");
  //  TFile* fD  = TFile::Open((baseDir+"/output/ROOT_OK/studyNoise_runD_DA-MC.root").c_str(), "READ");
  fD->cd();


  std::string folderName = baseDir+"/output/noisePlots__" + label + "/";
  system(("mkdir "+folderName).c_str());  
  
  TLegend* leg = NULL;
  
  
  std::vector<std::string> datasets;
  datasets.push_back("DA");
  datasets.push_back("MC");
  int nDatasets = int(datasets.size());
  
  std::vector<std::string> regions;
  regions.push_back("EB");
  regions.push_back("EE");
  int nRegions = int(regions.size());
  
  
  
  //-------------------------
  // get occupancy histograms
  
  TH1F* h_occupancy_vsNPU   = (TH1F*)( f->Get("MC/EB/h_occupancy_vsNPU_EB_MC") );
  TH1F* h_occupancy_vsIRing = (TH1F*)( f->Get("MC/EB/h_occupancy_vsIRing_EB_MC") );
  TH1F* h_occupancy_vsEta = (TH1F*)( f->Get("MC/EB/h_occupancy_vsEta_EB_MC") );
  TH1F* h_occupancy_vsAbsEta = (TH1F*)( f->Get("MC/EB/h_occupancy_vsAbsEta_EB_MC") );
  
  
  //-------------------------
  // get inclusive histograms
  
  TH1F** hEB_recHitE         = new TH1F*[datasets.size()];
  TH1F** hEB_recHitPedSubADC = new TH1F*[datasets.size()];
  TH1F** hEE_recHitE         = new TH1F*[datasets.size()];
  TH1F** hEE_recHitPedSubADC = new TH1F*[datasets.size()];
  
  
  // ------------
  // build graphs
    
  TGraphErrors** gEB_recHitE_RMS_vsNPU   = new TGraphErrors*[datasets.size()];
  TGraphErrors** gEB_recHitE_sigma_vsNPU = new TGraphErrors*[datasets.size()];
  TGraphErrors** gEE_recHitE_RMS_vsNPU   = new TGraphErrors*[datasets.size()];
  TGraphErrors** gEE_recHitE_sigma_vsNPU = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsIRing   = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsIRing_PU0   = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsIRing_PU70  = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsIRing_PU140 = new TGraphErrors*[datasets.size()];

  TGraphErrors** g_recHitE_sigma_vsIRing = new TGraphErrors*[datasets.size()];  
  TGraphErrors** g_recHitE_RMS_vsEta     = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigma_vsEta   = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsAbsEta   = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigma_vsAbsEta = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigma_vsAbsEtans = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigmaAt70PU_vsAbsEtans = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigmaAt140PU_vsAbsEtans = new TGraphErrors*[datasets.size()];

  TGraphErrors** g_recHitE_RMSAt0PU_vsIRing   = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigmaAt0PU_vsIRing = new TGraphErrors*[datasets.size()];  

  TGraphErrors** g_recHitE_RMSAt0PU_vsAbsEta   = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigmaAt0PU_vsAbsEta = new TGraphErrors*[datasets.size()];  
  TGraphErrors** g_recHitE_sigmaAt70PU_vsAbsEta = new TGraphErrors*[datasets.size()];  
  TGraphErrors** g_recHitE_sigmaAt140PU_vsAbsEta = new TGraphErrors*[datasets.size()];  

  TGraphErrors** g_recHitE_RMS_vsEta_PU0       = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsEta_PU70      = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsEta_PU140     = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsAbsEta_PU0    = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsAbsEta_PU70   = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsAbsEta_PU140  = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsAbsEta_PU70ns   = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_RMS_vsAbsEta_PU140ns  = new TGraphErrors*[datasets.size()];

  TGraphErrors** g_recHitE_sigma_vsEta_PU0      = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigma_vsEta_PU70     = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigma_vsEta_PU140    = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigma_vsAbsEta_PU0   = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigma_vsAbsEta_PU70  = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigma_vsAbsEta_PU140 = new TGraphErrors*[datasets.size()];

  TGraphErrors** g_recHitE_sigma_vsAbsEta_PU70ns  = new TGraphErrors*[datasets.size()];
  TGraphErrors** g_recHitE_sigma_vsAbsEta_PU140ns = new TGraphErrors*[datasets.size()];

  
  
  
  //-------------------
  // loop over datasets
  
  for(unsigned int datasetIt = 0; datasetIt < datasets.size(); ++datasetIt)
  {
    std::string dataset = datasets.at(datasetIt);
    
    
    char histoName[250];
    
    sprintf(histoName,"%s/EB/h_recHitE_EB_%s",dataset.c_str(),dataset.c_str());
    hEB_recHitE[datasetIt] = (TH1F*)( f->Get(histoName) );
    hEB_recHitE[datasetIt] -> Rebin(rebin);
    
    sprintf(histoName,"%s/EB/h_recHitPedSubADC_EB_%s",dataset.c_str(),dataset.c_str());
    hEB_recHitPedSubADC[datasetIt] = (TH1F*)( f->Get(histoName) );
    hEB_recHitPedSubADC[datasetIt] -> Rebin(2*rebin);
    
    sprintf(histoName,"%s/EE/h_recHitE_EE_%s",dataset.c_str(),dataset.c_str());
    hEE_recHitE[datasetIt] = (TH1F*)( f->Get(histoName) );
    hEE_recHitE[datasetIt] -> Rebin(rebin);
    
    sprintf(histoName,"%s/EE/h_recHitPedSubADC_EE_%s",dataset.c_str(),dataset.c_str());
    hEE_recHitPedSubADC[datasetIt] = (TH1F*)( f->Get(histoName) );
    hEE_recHitPedSubADC[datasetIt] -> Rebin(2*rebin);
    
    
    gEB_recHitE_RMS_vsNPU[datasetIt]   = new TGraphErrors();
    gEB_recHitE_sigma_vsNPU[datasetIt] = new TGraphErrors();
    gEE_recHitE_RMS_vsNPU[datasetIt]   = new TGraphErrors();
    gEE_recHitE_sigma_vsNPU[datasetIt] = new TGraphErrors();
    g_recHitE_RMS_vsIRing[datasetIt]   = new TGraphErrors();
    g_recHitE_RMS_vsIRing_PU0[datasetIt]   = new TGraphErrors();
    g_recHitE_RMS_vsIRing_PU70[datasetIt]   = new TGraphErrors();
    g_recHitE_RMS_vsIRing_PU140[datasetIt]   = new TGraphErrors();

    g_recHitE_sigma_vsIRing[datasetIt] = new TGraphErrors();
    g_recHitE_RMSAt0PU_vsIRing[datasetIt]   = new TGraphErrors();
    g_recHitE_sigmaAt0PU_vsIRing[datasetIt] = new TGraphErrors();    

    g_recHitE_RMSAt0PU_vsAbsEta[datasetIt]   = new TGraphErrors();
    g_recHitE_sigmaAt0PU_vsAbsEta[datasetIt] = new TGraphErrors();    
    g_recHitE_sigmaAt70PU_vsAbsEta[datasetIt] = new TGraphErrors();    
    g_recHitE_sigmaAt140PU_vsAbsEta[datasetIt] = new TGraphErrors();    

    g_recHitE_RMS_vsEta[datasetIt]      = new TGraphErrors();
    g_recHitE_sigma_vsEta[datasetIt]    = new TGraphErrors();
    g_recHitE_RMS_vsAbsEta[datasetIt]   = new TGraphErrors();
    g_recHitE_sigma_vsAbsEta[datasetIt] = new TGraphErrors();
    g_recHitE_sigma_vsAbsEtans[datasetIt] = new TGraphErrors();
    g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] = new TGraphErrors();
    g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] = new TGraphErrors();

    g_recHitE_RMS_vsEta_PU0[datasetIt]     = new TGraphErrors();
    g_recHitE_RMS_vsEta_PU70[datasetIt]    = new TGraphErrors();
    g_recHitE_RMS_vsEta_PU140[datasetIt]   = new TGraphErrors();
    g_recHitE_RMS_vsAbsEta_PU0[datasetIt]     = new TGraphErrors();
    g_recHitE_RMS_vsAbsEta_PU70[datasetIt]    = new TGraphErrors();
    g_recHitE_RMS_vsAbsEta_PU140[datasetIt]   = new TGraphErrors();
    g_recHitE_RMS_vsAbsEta_PU70ns[datasetIt]    = new TGraphErrors();
    g_recHitE_RMS_vsAbsEta_PU140ns[datasetIt]   = new TGraphErrors();

    g_recHitE_sigma_vsEta_PU0[datasetIt]   = new TGraphErrors();
    g_recHitE_sigma_vsEta_PU70[datasetIt]  = new TGraphErrors();
    g_recHitE_sigma_vsEta_PU140[datasetIt] = new TGraphErrors();
    g_recHitE_sigma_vsAbsEta_PU0[datasetIt]   = new TGraphErrors();
    g_recHitE_sigma_vsAbsEta_PU70[datasetIt]  = new TGraphErrors();
    g_recHitE_sigma_vsAbsEta_PU140[datasetIt] = new TGraphErrors();
    g_recHitE_sigma_vsAbsEta_PU70ns[datasetIt]  = new TGraphErrors();
    g_recHitE_sigma_vsAbsEta_PU140ns[datasetIt] = new TGraphErrors();

    
    for(unsigned int regionIt = 0; regionIt < regions.size(); ++regionIt)
    {
      std::string region = regions.at(regionIt);
      
      
      char canvasName[250];
      sprintf(canvasName,"c_recHitE_vsNPU_%s_%s",region.c_str(),dataset.c_str());
      TCanvas* c_temp1 = new TCanvas(canvasName,canvasName);
      
      char fileName[250];
      sprintf(fileName,"%s/output/noisePlots__%s/recHitE_vsNPU_%s_%s.pdf",baseDir.c_str(),label.c_str(),region.c_str(),dataset.c_str());
      std::string fileNameString1(fileName);
      c_temp1 -> Print((fileNameString1+"[").c_str(),"pdf");
      
      sprintf(canvasName,"c_recHitE_vsNPU_%s_%s_log",region.c_str(),dataset.c_str());
      TCanvas* c_temp2 = new TCanvas(canvasName,canvasName);
            
      sprintf(fileName,"%s/output/noisePlots__%s/recHitE_vsNPU_%s_%s_log.pdf",baseDir.c_str(),label.c_str(),region.c_str(),dataset.c_str());
      std::string fileNameString2(fileName);
      c_temp2 -> Print((fileNameString2+"[").c_str(),"pdf");      
      
      
      int point = 0;
      for(int bin = 1; bin <= h_occupancy_vsNPU->GetNbinsX(); ++bin)
      {
        float binCenter  = h_occupancy_vsNPU -> GetBinCenter(bin);
        float binLowEdge = h_occupancy_vsNPU -> GetBinLowEdge(bin);
        float binHigEdge = h_occupancy_vsNPU -> GetBinLowEdge(bin) + h_occupancy_vsNPU->GetBinWidth(bin);
        
        sprintf(histoName,"%s/%s/plots_vs_nPU/h_nPUDistr_nPU%2.1f-%2.1f_%s_%s",dataset.c_str(),region.c_str(),
                                                                               binLowEdge,binHigEdge,
                                                                               region.c_str(),dataset.c_str());

	std::cout << "histoName = " << std::string(histoName) << std::endl;
        TH1F* histo = (TH1F*)( f->Get(histoName) );
	//        TH1F* histoRMS;
        TH1F* histoRMS_AB;
        TH1F* histoRMS_C;
        TH1F* histoRMS_D;
        float binMean    = histo -> GetMean();
        float binRMS     = histo -> GetRMS();
        float binMeanErr = histo -> GetMeanError();
	std::cout << " nPU = " << binLowEdge << " " << binHigEdge << " binMean = " << binMean << " " << binRMS << " " << binMeanErr << std::endl;
        if( type == "DAMC" && (binMean < 6. || binMean > 30.) ) continue;
        
        
        sprintf(histoName,"%s/%s/plots_vs_nPU/h_recHitE_nPU%2.1f-%2.1f_%s_%s",dataset.c_str(),region.c_str(),
                                                                              binLowEdge,binHigEdge,
                                                                              region.c_str(),dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
	std::cout << "histoName = " << std::string(histoName) << std::endl;
        if( histo->GetEntries() < 30 ) continue;
        histo -> Rebin(rebin);
        
        sprintf(histoName,"%s/%s/plots_vs_nPU/h_recHitE_sigmanPU%2.1f-%2.1f_%s_%s",dataset.c_str(),region.c_str(),
		                                                                    binLowEdge,binHigEdge,
                                                                        	    region.c_str(),dataset.c_str());

	std::cout << "histoName = " << std::string(histoName) << std::endl;
	//histoRMS = (TH1F*)( f->Get(histoName) );
	histoRMS_AB = (TH1F*)( fAB->Get(histoName) );
	histoRMS_C = (TH1F*)( fC->Get(histoName) );
	histoRMS_D = (TH1F*)( fD->Get(histoName) );
	// 	FillHistoRMS(histoRMS, histo->GetRMS(), histo->GetMean() );
	// 	std::cout << " histoRMS = " << histo->GetRMS() << " effective = " << GetEffectiveSigma(histo) << std::endl;

	double E1 = (histoRMS_AB->GetBinContent(2) + histoRMS_C->GetBinContent(2) + histoRMS_D->GetBinContent(2))/3.; 
	double E2 = (histoRMS_AB->GetBinContent(3) + histoRMS_C->GetBinContent(3) + histoRMS_D->GetBinContent(3))/3.; 
	double Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.; 
	double EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
				pow(Eabcd, 2) ); 
	double R1 = histo->GetMean()-histo->GetRMS();
	double R2 = histo->GetMean()+histo->GetRMS();
	//	GetEffectiveSigma(histo, E1, E2);

// 	std::cout << " E1 = " << E1 << " E2 = " << E2 << " E2-E1 = " << E2-E1 << " sig = " << histoRMS->GetBinContent(1) << std::endl;
// 	std::cout << " R1 = " << R1 << " R2 = " << R2 << " rms = " << histo->GetRMS() << std::endl;
	TArrow* line1 = new TArrow(E1,0.000,E1,0.500); 
	line1->SetLineColor(kBlue); 	line1->SetLineStyle(2); line1->SetLineWidth(2);
	TArrow* line2 = new TArrow(E2,0.000,E2,0.500);
	line2->SetLineColor(kBlue); 	line2->SetLineStyle(2); line2->SetLineWidth(2);

	TArrow* line1r = new TArrow(R1,0.000,R1,0.500); 
	line1r->SetLineColor(kYellow); 	line1r->SetLineWidth(2);
	TArrow* line2r = new TArrow(R2,0.000,R2,0.500);
	line2r->SetLineColor(kYellow); 	line2r->SetLineWidth(2);
	

        
        TF1* f_gaus;
        if( region == "EB" ) f_gaus = new TF1("f_gaus","[0]*exp(-1.*(x-[1])*(x-[1])/2/[2]/[2])",-0.10,0.10);
        if( region == "EE" ) f_gaus = new TF1("f_gaus","[0]*exp(-1.*(x-[1])*(x-[1])/2/[2]/[2])",-0.25,0.25);
        f_gaus -> SetParameters(histo->GetMaximum(),0.,0.05);
        histo -> Fit("f_gaus","QENRS+");
        
        sprintf(canvasName,"c_recHitE_nPU%2.1f-%2.1f_%s_%s",binLowEdge,binHigEdge,region.c_str(),dataset.c_str());
        TCanvas* c_temp = new TCanvas(canvasName,canvasName);
        SetHistoStyle(histo,dataset,region,"recHitE",&leg);
        histo -> Draw();
        f_gaus -> Draw("same");
	//	histoRMS->Draw("hist, same");
	line1 -> Draw("same"); 
	line2 -> Draw("same");     
	line1r -> Draw("same"); 
	line2r -> Draw("same");     
        c_temp -> Print(fileNameString1.c_str(),"pdf");
        delete c_temp;
        
        sprintf(canvasName,"c_recHitE_nPU%2.1f-%2.1f_%s_%s_log",binLowEdge,binHigEdge,region.c_str(),dataset.c_str());
        c_temp = new TCanvas(canvasName,canvasName);
        c_temp -> SetLogy();
        SetHistoStyle(histo,dataset,region,"recHitE",&leg,true);
        histo -> Draw();
        f_gaus -> Draw("same");
	//	histoRMS->Draw("hist, same");
	line1 -> Draw("same"); 
	line2 -> Draw("same");     
	line1r -> Draw("same"); 
	line2r -> Draw("same");     
        c_temp -> Print(fileNameString2.c_str(),"pdf");
        delete c_temp;        
        
        if( region == "EB" )
        {
          gEB_recHitE_RMS_vsNPU[datasetIt] -> SetPoint(point,binMean,histo->GetRMS());
          gEB_recHitE_RMS_vsNPU[datasetIt] -> SetPointError(point,binRMS,histo->GetRMSError());

	  //	  gEB_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(point,binMean,GetEffectiveSigma(histo));
// 	  gEB_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(point,binMean,histoRMS->GetBinContent(1));
// 	  gEB_recHitE_sigma_vsNPU[datasetIt] -> SetPointError(point,binRMS,histo->GetRMSError());
	  gEB_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(point,binMean,Eabcd);
	  gEB_recHitE_sigma_vsNPU[datasetIt] -> SetPointError(point,binRMS,EabcdRMS);
        }
        if( region == "EE" )
        {
          gEE_recHitE_RMS_vsNPU[datasetIt] -> SetPoint(point,binMean,histo->GetRMS());
          gEE_recHitE_RMS_vsNPU[datasetIt] -> SetPointError(point,binRMS,histo->GetRMSError());
          
	  //          gEE_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(point,binMean,GetEffectiveSigma(histo));
// 	  gEE_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(point,binMean,histoRMS->GetBinContent(1));
//           gEE_recHitE_sigma_vsNPU[datasetIt] -> SetPointError(point,binRMS,histo->GetRMSError());
	  gEE_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(point,binMean,Eabcd);
          gEE_recHitE_sigma_vsNPU[datasetIt] -> SetPointError(point,binRMS,EabcdRMS);
        }
        delete f_gaus;
        ++point;
      }
      
      c_temp1 -> Print((fileNameString1+"]").c_str(),"pdf");
      delete c_temp1;
      
      c_temp2 -> Print((fileNameString2+"]").c_str(),"pdf");
      delete c_temp2;
                  
      gEB_recHitE_RMS_vsNPU[datasetIt] -> SetPoint(gEB_recHitE_RMS_vsNPU[datasetIt]->GetN(),-1.,-1.);
      gEB_recHitE_RMS_vsNPU[datasetIt] -> SetPoint(gEB_recHitE_RMS_vsNPU[datasetIt]->GetN(),200.,-1.);
      
      gEB_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(gEB_recHitE_sigma_vsNPU[datasetIt]->GetN(),-1.,-1.);
      gEB_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(gEB_recHitE_sigma_vsNPU[datasetIt]->GetN(),200.,-1.);
      
      gEE_recHitE_RMS_vsNPU[datasetIt] -> SetPoint(gEE_recHitE_RMS_vsNPU[datasetIt]->GetN(),-1.,-1.);
      gEE_recHitE_RMS_vsNPU[datasetIt] -> SetPoint(gEE_recHitE_RMS_vsNPU[datasetIt]->GetN(),200.,-1.);
      
      gEE_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(gEE_recHitE_sigma_vsNPU[datasetIt]->GetN(),-1.,-1.);
      gEE_recHitE_sigma_vsNPU[datasetIt] -> SetPoint(gEE_recHitE_sigma_vsNPU[datasetIt]->GetN(),200.,-1.);
    }
  
    //////////////////////////////
   std::cout << " >>> h_occupancy_vsIRing iRing " << std::endl;
    
   int  point = 0;
    for(int iRingBin = 1; iRingBin <= h_occupancy_vsIRing->GetNbinsX(); ++iRingBin)
    {  
      float iRingBinCenter  = h_occupancy_vsIRing -> GetBinCenter(iRingBin);
      float iRingBinLowEdge = h_occupancy_vsIRing -> GetBinLowEdge(iRingBin);
      float iRingBinHigEdge = h_occupancy_vsIRing -> GetBinLowEdge(iRingBin) + h_occupancy_vsIRing->GetBinWidth(iRingBin);
      
      if( iRingBinCenter == 0 ) continue;
      
            
      TGraphErrors* g_temp = new TGraphErrors();
      int point2 = 0;
      for(int nPUBin = 1; nPUBin <= h_occupancy_vsNPU->GetNbinsX(); ++nPUBin)
      {
        float nPUBinLowEdge = h_occupancy_vsNPU -> GetBinLowEdge(nPUBin);
        float nPUBinHigEdge = h_occupancy_vsNPU -> GetBinLowEdge(nPUBin) + h_occupancy_vsNPU->GetBinWidth(nPUBin);      
        
        TH1F* histo;
        char histoName[250];
        
        if( fabs(iRingBinCenter) < 86 )
          sprintf(histoName,"%s/EB/plots_vs_nPU_and_iRing/h_nPUDistr_nPU%2.1f-%2.1f__iRing%2.1f-%2.1f_EB_%s",dataset.c_str(),
                                                                                                             nPUBinLowEdge,nPUBinHigEdge,
                                                                                                             iRingBinLowEdge,iRingBinHigEdge,
                                                                                                             dataset.c_str());        
        else
          sprintf(histoName,"%s/EE/plots_vs_nPU_and_iRing/h_nPUDistr_nPU%2.1f-%2.1f__iRing%2.1f-%2.1f_EE_%s",dataset.c_str(),
                                                                                                             nPUBinLowEdge,nPUBinHigEdge,
                                                                                                             iRingBinLowEdge,iRingBinHigEdge,
                                                                                                             dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
        float binMean    = histo -> GetMean();
        float binMeanErr = histo -> GetMeanError();
        if( type == "DAMC" && (binMean < 6. || binMean > 30.) ) continue;        
        
        if( fabs(iRingBinCenter) < 86 )
          sprintf(histoName,"%s/EB/plots_vs_nPU_and_iRing/h_recHitE_nPU%2.1f-%2.1f__iRing%2.1f-%2.1f_EB_%s",dataset.c_str(),
                                                                                                            nPUBinLowEdge,nPUBinHigEdge,
                                                                                                            iRingBinLowEdge,iRingBinHigEdge,
                                                                                                            dataset.c_str());        
        else
          sprintf(histoName,"%s/EE/plots_vs_nPU_and_iRing/h_recHitE_nPU%2.1f-%2.1f__iRing%2.1f-%2.1f_EE_%s",dataset.c_str(),
                                                                                                            nPUBinLowEdge,nPUBinHigEdge,
                                                                                                            iRingBinLowEdge,iRingBinHigEdge,
                                                                                                            dataset.c_str());        
        histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
	histo -> Rebin(rebin);
        
        g_temp -> SetPoint(point2,binMean,histo->GetRMS());
        g_temp -> SetPointError(point2,binMeanErr,histo->GetRMSError());
        ++point2;
      }
      
      TF1* f_temp = new TF1("f_temp","[0]+[1]*x",-0.5,199.5);
      if(type == "DAMC") f_temp->SetRange(5., 30.);
      g_temp -> Fit("f_temp","QNRS");
      if(type == "DAMC") f_temp->SetRange(5., 100.);
      g_temp -> Fit("f_temp","QNRS");
      
      g_recHitE_RMSAt0PU_vsIRing[datasetIt] -> SetPoint(point,iRingBinCenter,f_temp->GetParameter(0));
      g_recHitE_RMSAt0PU_vsIRing[datasetIt] -> SetPointError(point,0.,f_temp->GetParError(0));
      
      delete f_temp;
      delete g_temp;
      ++point;
    }
  

    ///////
    std::cout << " >>> h_occupancy_vsAbsEta absEta " << std::endl;
    
    char cName[250];
    sprintf(cName,"c_recHitE_sigma_vsNPU_%s",dataset.c_str());
    TCanvas* c_t1 = new TCanvas(cName,cName);

    char fName[250];
    sprintf(fName,"%s/output/noisePlots__%s/recHitE_sigma_vsNPU_%s.pdf",baseDir.c_str(),label.c_str(),dataset.c_str());
    std::string fNameString1(fName);
    c_t1 -> Print((fNameString1+"[").c_str(),"pdf");


    point = 0;
    for(int iRingBin = 1; iRingBin <= h_occupancy_vsAbsEta->GetNbinsX(); ++iRingBin)
    {  
      float iRingBinCenter  = h_occupancy_vsAbsEta -> GetBinCenter(iRingBin);
      float iRingBinLowEdge = h_occupancy_vsAbsEta -> GetBinLowEdge(iRingBin);
      float iRingBinHigEdge = h_occupancy_vsAbsEta -> GetBinLowEdge(iRingBin) + h_occupancy_vsAbsEta->GetBinWidth(iRingBin);
      
      if( iRingBinCenter == 0 ) continue;
      
            
      TGraphErrors* g_temp = new TGraphErrors();
      int point2 = 0;
      std::cout << " >>> h_occupancy_vsNPU->GetNbinsX() = " << h_occupancy_vsNPU->GetNbinsX() << std::endl;
      for(int nPUBin = 1; nPUBin <= h_occupancy_vsNPU->GetNbinsX(); ++nPUBin)
      {
        float nPUBinLowEdge = h_occupancy_vsNPU -> GetBinLowEdge(nPUBin);
        float nPUBinHigEdge = h_occupancy_vsNPU -> GetBinLowEdge(nPUBin) + h_occupancy_vsNPU->GetBinWidth(nPUBin);      
        

        TH1F* histo;
	//	TH1F* histoRMS;
	TH1F* histoRMS_AB;
	TH1F* histoRMS_C;
	TH1F* histoRMS_D;
        char histoName[250];
        char histoNameRMS[250];
        
	std::cout << ">>>>>>>> histo definiti " << std::endl;

        if( fabs(iRingBinCenter) < 1.4442 )
          sprintf(histoName,"%s/EB/plots_vs_nPU_and_absEta/h_nPUDistr_nPU%2.1f-%2.1f__absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),
                                                                                                             nPUBinLowEdge,nPUBinHigEdge,
                                                                                                             iRingBinLowEdge,iRingBinHigEdge,
                                                                                                             dataset.c_str());        
        else
          sprintf(histoName,"%s/EE/plots_vs_nPU_and_absEta/h_nPUDistr_nPU%2.1f-%2.1f__absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),
                                                                                                             nPUBinLowEdge,nPUBinHigEdge,
                                                                                                             iRingBinLowEdge,iRingBinHigEdge,
                                                                                                             dataset.c_str());
	std::cout << ">>>>>>>> f->GetName() = " << f->GetName() << std::endl;
	std::cout << ">>>>>>>> histo->GetName() = " << std::string(histoName) << std::endl;

        histo = (TH1F*)( f->Get(histoName) );
        float binMean    = histo -> GetMean();
        float binMeanErr = histo -> GetMeanError();
        if(type == "DAMC" &&  (binMean < 6. || binMean > 30.) ) continue;        
        
        if( fabs(iRingBinCenter) < 1.4442 ){
          sprintf(histoName,"%s/EB/plots_vs_nPU_and_absEta/h_recHitE_nPU%2.1f-%2.1f__absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),
                                                                                                            nPUBinLowEdge,nPUBinHigEdge,
                                                                                                            iRingBinLowEdge,iRingBinHigEdge,
                                                                                                            dataset.c_str());        

          sprintf(histoNameRMS,"%s/EB/plots_vs_nPU_and_absEta/h_recHitE_sigma_nPU%2.1f-%2.1f__absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),
                                                                                                            nPUBinLowEdge,nPUBinHigEdge,
                                                                                                            iRingBinLowEdge,iRingBinHigEdge,
                                                                                                            dataset.c_str());        
	}
        else{
          sprintf(histoName,"%s/EE/plots_vs_nPU_and_absEta/h_recHitE_nPU%2.1f-%2.1f__absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),
                                                                                                            nPUBinLowEdge,nPUBinHigEdge,
                                                                                                            iRingBinLowEdge,iRingBinHigEdge,
                                                                                                            dataset.c_str());        
          sprintf(histoNameRMS,"%s/EE/plots_vs_nPU_and_absEta/h_recHitE_sigma_nPU%2.1f-%2.1f__absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),
                                                                                                            nPUBinLowEdge,nPUBinHigEdge,
                                                                                                            iRingBinLowEdge,iRingBinHigEdge,
                                                                                                            dataset.c_str());        
	}

        histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
	histo -> Rebin(rebin);

	std::cout << ">>>>>>>> prima di histoRMS " << std::endl;
	//	histoRMS = (TH1F*)( f->Get(histoNameRMS) );
	histoRMS_AB = (TH1F*)( fAB->Get(histoNameRMS) );
	histoRMS_C = (TH1F*)( fC->Get(histoNameRMS) );
	histoRMS_D = (TH1F*)( fD->Get(histoNameRMS) );

	std::cout << ">>>>>>>> histo presi " << std::endl;
	double Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.;
        double EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
                                pow(Eabcd, 2) );

	std::cout << " >>> Eabcd " << Eabcd << std::endl;
	std::cout << " >>> EabcdRMS " << EabcdRMS << std::endl;
        g_temp -> SetPoint(point2,binMean,Eabcd);
        g_temp -> SetPointError(point2,binMeanErr,EabcdRMS);

//         g_temp -> SetPoint(point2,binMean,histoRMS->GetBinContent(1));
//         g_temp -> SetPointError(point2,binMeanErr,histo->GetRMSError());
        ++point2;
      }
      
      TF1* f_temp = new TF1("f_temp","[0]+[1]*x",-0.5,199.5);
      if(type == "DAMC") f_temp->SetRange(5., 30.);
      g_temp -> Fit("f_temp","QNRS");
      if(type == "DAMC") f_temp->SetRange(5., 100.);
      g_temp -> Fit("f_temp","QNRS");


      sprintf(cName,"c_sigma_vsPU_absEta%1.4f-%1.4f_%s",iRingBinLowEdge, iRingBinHigEdge, dataset.c_str());
      TCanvas* c_temp = new TCanvas(cName,cName);
      c_temp->SetGridx();
      c_temp->SetGridy();
      g_temp->Draw("ap");
      f_temp->Draw("same");
      if(type == "DAMC") g_temp->GetXaxis()->SetRangeUser(0., 50.);
      else g_temp->GetXaxis()->SetRangeUser(0., 190.);
      g_temp->SetMaximum(0.5);
      g_temp->SetMinimum(0.);
      g_temp->GetXaxis()->SetTitle(Form("|#eta| %1.4f-%1.4f", iRingBinLowEdge, iRingBinHigEdge));

      double y0 = f_temp->GetParameter(0);
      double y0E = f_temp->GetParError(0);
      double m = f_temp->GetParameter(1);
      double mE = f_temp->GetParError(1);
      ///////////
      TPaveText* text = new TPaveText(0.20,0.85,0.50,0.90,"brNDC");
      text -> SetTextSize(0.04);
      text -> SetLineWidth(0);
      text -> SetLineColor(kWhite);
      text -> SetTextFont(42);
      text -> SetTextColor(kBlack);
      text -> SetFillColor(kWhite);
      text -> SetFillStyle(1001);
      text -> SetBorderSize(1);
      char par0Text[50]; sprintf(par0Text,"y_{0} = %1.2e #pm %1.0e",y0, y0E);
      text -> AddText(par0Text);
      text -> Draw("same");

      text = new TPaveText(0.50,0.85,0.80,0.90,"brNDC");
      text -> SetTextSize(0.04);
      text -> SetLineWidth(0);
      text -> SetLineColor(kWhite);
      text -> SetTextFont(42);
      text -> SetTextColor(kBlack);
      text -> SetFillColor(kWhite);
      text -> SetFillStyle(1001);
      text -> SetBorderSize(1);
      char par1Text[50]; sprintf(par1Text,"m = %1.2e #pm %1.0e",m, mE);
      text -> AddText(par1Text);
      text -> Draw("same");
      ///////////

      std::cout << ">>>>>>>> histo disegnati " << std::endl;
      c_temp -> Print(fNameString1.c_str(),"pdf");
      delete c_temp;


      g_recHitE_sigmaAt0PU_vsAbsEta[datasetIt] -> SetPoint(point,iRingBinCenter,y0);
      g_recHitE_sigmaAt0PU_vsAbsEta[datasetIt] -> SetPointError(point,0.,y0E);

      g_recHitE_sigmaAt70PU_vsAbsEta[datasetIt] -> SetPoint(point,iRingBinCenter,70.*m + y0);
      g_recHitE_sigmaAt70PU_vsAbsEta[datasetIt] -> SetPointError(point,0.,70*mE + y0E);

      g_recHitE_sigmaAt140PU_vsAbsEta[datasetIt] -> SetPoint(point,iRingBinCenter,140.*m + y0);
      g_recHitE_sigmaAt140PU_vsAbsEta[datasetIt] -> SetPointError(point,0.,140*mE + y0E);
      
      delete f_temp;
      delete g_temp;
      ++point;
    }
    c_t1 -> Print((fNameString1+"]").c_str(),"pdf");
    delete c_t1;
 
    //////////////////////////////
    std::cout << " >>> h_occupancy_vsIRing bin " << std::endl;

    point = 0;
    for(int bin = 1; bin <= h_occupancy_vsIRing->GetNbinsX(); ++bin)
    {  
      float binCenter  = h_occupancy_vsIRing -> GetBinCenter(bin);
      float binLowEdge = h_occupancy_vsIRing -> GetBinLowEdge(bin);
      float binHigEdge = h_occupancy_vsIRing -> GetBinLowEdge(bin) + h_occupancy_vsIRing->GetBinWidth(bin);
      
      if( binCenter == 0 ) continue;
    
      
      if( fabs(binCenter) < 86 )
      {
        char histoName[50];
        sprintf(histoName,"%s/EB/plots_vs_iRing/h_recHitE_iRing%2.1f-%2.1f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        TH1F* histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
	histo -> Rebin(rebin);

	//         TF1* f_gaus = new TF1("f_gaus","[0]*exp(-1.*(x-[1])*(x-[1])/2/[2]/[2])",-0.10,0.10);
	//         f_gaus -> SetParameters(histo->GetMaximum(),0.,0.05);
        
        g_recHitE_RMS_vsIRing[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsIRing[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        

	//         g_recHitE_sigma_vsIRing[datasetIt] -> SetPoint(point,binCenter,f_gaus->GetParameter(2));
	//         g_recHitE_sigma_vsIRing[datasetIt] -> SetPointError(point,0.,f_gaus->GetParError(2));
        g_recHitE_sigma_vsIRing[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
        g_recHitE_sigma_vsIRing[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());

	delete histo;

	if(type == "MCMC"){
	//PU0        
sprintf(histoName,"%s/EB/plots_vs_nPU_and_iRing/h_recHitE_nPU-1.0-1.0__iRing%2.1f-%2.1f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
	//	histo -> Rebin(rebin);

        g_recHitE_RMS_vsIRing_PU0[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsIRing_PU0[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
	//        delete f_gaus;

	delete histo;

	//PU70        
sprintf(histoName,"%s/EB/plots_vs_nPU_and_iRing/h_recHitE_nPU69.0-71.0__iRing%2.1f-%2.1f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
	//	histo -> Rebin(rebin);

        g_recHitE_RMS_vsIRing_PU70[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsIRing_PU70[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
	//	delete f_gaus;

	delete histo;

	//PU140        
sprintf(histoName,"%s/EB/plots_vs_nPU_and_iRing/h_recHitE_nPU139.0-141.0__iRing%2.1f-%2.1f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
// 	histo -> Rebin(rebin);
// 	histo -> Rebin(2);

        g_recHitE_RMS_vsIRing_PU140[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsIRing_PU140[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
	//        delete f_gaus;

	delete histo;
	}
        ++point;
      }
      else
      {
        char histoName[50];
        sprintf(histoName,"%s/EE/plots_vs_iRing/h_recHitE_iRing%2.1f-%2.1f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        TH1F* histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
        histo -> Rebin(rebin);
        
	//         TF1* f_gaus = new TF1("f_gaus","[0]*exp(-1.*(x-[1])*(x-[1])/2/[2]/[2])",-0.25,0.25);
	//         f_gaus -> SetParameters(histo->GetMaximum(),0.,0.05);
	//         histo -> Fit("f_gaus","QENRS");
              
        g_recHitE_RMS_vsIRing[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsIRing[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	//         g_recHitE_sigma_vsIRing[datasetIt] -> SetPoint(point,binCenter,f_gaus->GetParameter(2));
	//         g_recHitE_sigma_vsIRing[datasetIt] -> SetPointError(point,0.,f_gaus->GetParError(2));
        g_recHitE_sigma_vsIRing[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
	g_recHitE_sigma_vsIRing[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	//        delete f_gaus;
	delete histo;

	if(type == "MCMC"){
	//PU0        
sprintf(histoName,"%s/EE/plots_vs_nPU_and_iRing/h_recHitE_nPU-1.0-1.0__iRing%2.1f-%2.1f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
	//        histo -> Rebin(rebin);

        g_recHitE_RMS_vsIRing_PU0[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsIRing_PU0[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
	//        delete f_gaus;

	delete histo;

	//PU70        
sprintf(histoName,"%s/EE/plots_vs_nPU_and_iRing/h_recHitE_nPU69.0-71.0__iRing%2.1f-%2.1f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
	//	histo -> Rebin(rebin);

        g_recHitE_RMS_vsIRing_PU70[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsIRing_PU70[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
	//        delete f_gaus;

	delete histo;

	//PU140        
sprintf(histoName,"%s/EE/plots_vs_nPU_and_iRing/h_recHitE_nPU139.0-141.0__iRing%2.1f-%2.1f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
//      histo -> Rebin(rebin);
// 	histo -> Rebin(2);

        g_recHitE_RMS_vsIRing_PU140[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsIRing_PU140[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
	//        delete f_gaus;
	//	std::cout << " iRing " << binLowEdge << " " << binHigEdge << "rms = " << histo->GetRMS() << std::endl;
	delete histo;
	}
        ++point;      
      }
    }
  

    std::cout << " >>> h_occupancy_vsEta " << std::endl;

    point = 0;
    for(int bin = 1; bin <= h_occupancy_vsEta->GetNbinsX(); ++bin)
    {  
      float binCenter  = h_occupancy_vsEta -> GetBinCenter(bin);
      float binLowEdge = h_occupancy_vsEta -> GetBinLowEdge(bin);
      float binHigEdge = h_occupancy_vsEta -> GetBinLowEdge(bin) + h_occupancy_vsEta->GetBinWidth(bin);
      
      if( binCenter == 0 ) continue;
      
      if( fabs(binCenter) < 1.4442 )
      {
        char histoName[50];
        sprintf(histoName,"%s/EB/plots_vs_eta/h_recHitE_eta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        TH1F* histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
	histo -> Rebin(rebin);
        
        g_recHitE_RMS_vsEta[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsEta[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
        g_recHitE_sigma_vsEta[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
	g_recHitE_sigma_vsEta[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	delete histo;

	if(type == "MCMC"){
	//PU0
        sprintf(histoName,"%s/EB/plots_vs_nPU_and_eta/h_recHitE_nPU-1.0-1.0__eta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
	if( histo->GetEntries() < 30 ) continue;
	histo -> Rebin(rebin);
        
        g_recHitE_RMS_vsEta_PU0[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsEta_PU0[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
        g_recHitE_sigma_vsEta_PU0[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
	g_recHitE_sigma_vsEta_PU0[datasetIt] -> SetPointError(point,0., histo->GetRMSError());
        
	delete histo;

	//PU70
        sprintf(histoName,"%s/EB/plots_vs_nPU_and_eta/h_recHitE_nPU69.0-71.0__eta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
        histo -> Rebin(rebin);
        
        g_recHitE_RMS_vsEta_PU70[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsEta_PU70[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
        g_recHitE_sigma_vsEta_PU70[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
        g_recHitE_sigma_vsEta_PU70[datasetIt] -> SetPointError(point,0., histo->GetRMSError());
        
	delete histo;

	//PU140
    sprintf(histoName,"%s/EB/plots_vs_nPU_and_eta/h_recHitE_nPU139.0-141.0__eta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
        histo -> Rebin(rebin);
	histo -> Rebin(2);

        g_recHitE_RMS_vsEta_PU140[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsEta_PU140[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
        g_recHitE_sigma_vsEta_PU140[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
	g_recHitE_sigma_vsEta_PU140[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	delete histo;
	}
        ++point;
      }
      else
      {
        char histoName[50];
        sprintf(histoName,"%s/EE/plots_vs_eta/h_recHitE_eta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        TH1F* histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
	histo -> Rebin(rebin);
        
        g_recHitE_RMS_vsEta[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsEta[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
        g_recHitE_sigma_vsEta[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
        g_recHitE_sigma_vsEta[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	delete histo;

	if(type == "MCMC"){
	//PU0
        sprintf(histoName,"%s/EE/plots_vs_nPU_and_eta/h_recHitE_nPU-1.0-1.0__eta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
	histo -> Rebin(rebin);
        
        g_recHitE_RMS_vsEta_PU0[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsEta_PU0[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
        g_recHitE_sigma_vsEta_PU0[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
	g_recHitE_sigma_vsEta_PU0[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	delete histo;

	//PU70
        sprintf(histoName,"%s/EE/plots_vs_nPU_and_eta/h_recHitE_nPU69.0-71.0__eta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
        histo -> Rebin(rebin);
        
        g_recHitE_RMS_vsEta_PU70[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsEta_PU70[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
        g_recHitE_sigma_vsEta_PU70[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
	g_recHitE_sigma_vsEta_PU70[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	delete histo;

	//PU140
    sprintf(histoName,"%s/EE/plots_vs_nPU_and_eta/h_recHitE_nPU139.0-141.0__eta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
        if( histo->GetEntries() < 30 ) continue;
        histo -> Rebin(rebin);
	histo -> Rebin(2);

        g_recHitE_RMS_vsEta_PU140[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsEta_PU140[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
        g_recHitE_sigma_vsEta_PU140[datasetIt] -> SetPoint(point,binCenter,GetEffectiveSigma(histo));
	g_recHitE_sigma_vsEta_PU140[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());

	delete histo;
	}

        ++point;      
      }
    }
  

    std::cout << " >>> h_occupancy_vsAbsEta " << std::endl;

    point = 0;
    for(int bin = 1; bin <= h_occupancy_vsAbsEta->GetNbinsX(); ++bin)
    {  
      float binCenter  = h_occupancy_vsAbsEta -> GetBinCenter(bin);
      float binLowEdge = h_occupancy_vsAbsEta -> GetBinLowEdge(bin);
      float binHigEdge = h_occupancy_vsAbsEta -> GetBinLowEdge(bin) + h_occupancy_vsAbsEta->GetBinWidth(bin);
      
      if( binCenter == 0 ) continue;
      
      if( fabs(binCenter) < 1.4442 )
      {
        char histoName[50];
        sprintf(histoName,"%s/EB/plots_vs_absEta/h_recHitE_absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        TH1F* histo = (TH1F*)( f->Get(histoName) );
	sprintf(histoName,"%s/EB/plots_vs_absEta/h_recHitE_sigma_absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
//         TH1F* histoRMS = (TH1F*)( f->Get(histoName) );
	TH1F* histoRMS_AB = (TH1F*)( fAB->Get(histoName) );
	TH1F* histoRMS_C = (TH1F*)( fC->Get(histoName) );
	TH1F* histoRMS_D = (TH1F*)( fD->Get(histoName) );
	std::cout << ">>>>>>>> EB prima di getBins " << std::endl;
	double Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.;
	double EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
				pow(Eabcd, 2) );
	std::cout << ">>>>>>>> EB dopo getBins " << std::endl;
        if( histo->GetEntries() < 30 ) continue;
        histo -> Rebin(rebin);
        
        g_recHitE_RMS_vsAbsEta[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsAbsEta[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	float sigmaPUContrib = Eabcd; //histoRMS->GetBinContent(1);
	float sigmaPUContribE = EabcdRMS; //histo->GetRMSError();
	float sigmaPUContrib2 = pow(sigmaPUContrib, 2);
        g_recHitE_sigma_vsAbsEta[datasetIt] -> SetPoint(point,binCenter, sigmaPUContrib);
        g_recHitE_sigma_vsAbsEta[datasetIt] -> SetPointError(point,0.,sigmaPUContribE);
        
	double x, y;
	float sigmaAt0PUPUContrib = g_recHitE_sigmaAt0PU_vsAbsEta[datasetIt]->GetPoint(point, x, y);
	double ey = g_recHitE_sigmaAt0PU_vsAbsEta[datasetIt]->GetErrorY(point);
	if(y < sigmaPUContrib) g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, sqrt(sigmaPUContrib2 - y*y));
	else g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, 0. );
	double numE = pow(sigmaPUContrib * sigmaPUContribE, 2) + pow(y*ey, 2);
	double denE = sigmaPUContrib2 - y*y;
	//	g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPointError(point,0., sqrt(numE/denE));
	//	g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPointError(point,0., 0.);
	g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPointError(point,0., sigmaPUContribE);
	std::cout << ">>>>>>>> EB error su sigma = " << sigmaPUContribE << " at PU0 = " << ey << std::endl;

	double xp, yp;
	float sigmaAt70PUPUContrib = g_recHitE_sigmaAt70PU_vsAbsEta[datasetIt]->GetPoint(point, xp, yp);
	double eyp = g_recHitE_sigmaAt70PU_vsAbsEta[datasetIt]->GetErrorY(point);
	if(y < yp) g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, sqrt(yp*yp - y*y));
	else g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, 0. );
	double numEp = pow(yp*eyp, 2) + pow(y*ey, 2);
	double denEp = yp*yp - y*y;
	//        g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,sqrt(numEp/denEp));
	//        g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,0.);
        g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,eyp);


	float sigmaAt140PUPUContrib = g_recHitE_sigmaAt140PU_vsAbsEta[datasetIt]->GetPoint(point, xp, yp);
	eyp = g_recHitE_sigmaAt140PU_vsAbsEta[datasetIt]->GetErrorY(point);
	if(y < yp) g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, sqrt(yp*yp - y*y));
	else g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, 0. );
	numEp = pow(yp*eyp, 2) + pow(y*ey, 2);
	denEp = yp*yp - y*y;
	//        g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,sqrt(numEp/denEp));
	//        g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,0.);
        g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,eyp);


	delete histo;
	delete histoRMS_AB;
	delete histoRMS_C;
	delete histoRMS_D;

	if(type == "MCMC"){
	//PU0
  sprintf(histoName,"%s/EB/plots_vs_nPU_and_absEta/h_recHitE_nPU-1.0-1.0__absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
sprintf(histoName,"%s/EB/plots_vs_nPU_and_absEta/h_recHitE_sigma_nPU-1.0-1.0__absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
//        TH1F* histoRMS = (TH1F*)( f->Get(histoName) );
          histoRMS_AB = (TH1F*)( fAB->Get(histoName) );
	  histoRMS_C = (TH1F*)( fC->Get(histoName) );
	  histoRMS_D = (TH1F*)( fD->Get(histoName) );
	  Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.;
	  EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
				 pow(Eabcd, 2) );
	//        if( histo->GetEntries() < 30 ) continue;
	//        histo -> Rebin(rebin);

	float rmsPU0contrib = histo->GetRMS();
	float rmsPU0contrib2 = pow(rmsPU0contrib,2);      
        g_recHitE_RMS_vsAbsEta_PU0[datasetIt] -> SetPoint(point,binCenter,rmsPU0contrib);
        g_recHitE_RMS_vsAbsEta_PU0[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	float PU0contrib = Eabcd; //histoRMS->GetBinContent(1);
	float PU0ContribE = EabcdRMS; //histo->GetRMSError();
	float PU0contrib2 = pow(PU0contrib,2);
        g_recHitE_sigma_vsAbsEta_PU0[datasetIt]->SetPoint(point,binCenter,PU0contrib);
        g_recHitE_sigma_vsAbsEta_PU0[datasetIt]->SetPointError(point,0., PU0ContribE);
	//	std::cout << " >>> PU 0 - |eta| << " << binLowEdge << " " << binHigEdge << " " << PU0contrib << std::endl;
        
	delete histo;
	//	delete histoRMS;
	delete histoRMS_AB;
	delete histoRMS_C;
	delete histoRMS_D;
	//PU70
  sprintf(histoName,"%s/EB/plots_vs_nPU_and_absEta/h_recHitE_nPU69.0-71.0__absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
  sprintf(histoName,"%s/EB/plots_vs_nPU_and_absEta/h_recHitE_sigma_nPU69.0-71.0__absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
  //        histoRMS = (TH1F*)( f->Get(histoName) );
        histoRMS_AB = (TH1F*)( fAB->Get(histoName) );
	histoRMS_C = (TH1F*)( fC->Get(histoName) );
	histoRMS_D = (TH1F*)( fD->Get(histoName) );
	Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.;
	EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
			 pow(Eabcd, 2) );

	//	if( histo->GetEntries() < 30 ) continue;
	//        histo -> Rebin(rebin);
        
	float rmsPU70contrib = histo->GetRMS();
	float rmsPU70contrib2 = pow(rmsPU70contrib,2);      
        g_recHitE_RMS_vsAbsEta_PU70[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsAbsEta_PU70[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        g_recHitE_RMS_vsAbsEta_PU70ns[datasetIt] -> SetPoint(point,binCenter, sqrt(rmsPU70contrib2 - rmsPU0contrib2) );
        g_recHitE_RMS_vsAbsEta_PU70ns[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	float PU70contrib = Eabcd; //histoRMS->GetBinContent(1);
	float PU70ContribE = EabcdRMS; //histo->GetRMSError();
	float PU70contrib2 = pow(PU70contrib,2);
        g_recHitE_sigma_vsAbsEta_PU70[datasetIt] -> SetPoint(point,binCenter,PU70contrib);
	g_recHitE_sigma_vsAbsEta_PU70[datasetIt] -> SetPointError(point,0.,PU70ContribE);
        g_recHitE_sigma_vsAbsEta_PU70ns[datasetIt]->SetPoint(point,binCenter,sqrt(PU70contrib2 - PU0contrib2));
	numE = pow(PU70contrib * PU70ContribE, 2) + pow(PU0contrib * PU0ContribE, 2);
	denE = PU70ContribE - PU0contrib2;
	//	g_recHitE_sigma_vsAbsEta_PU70ns[datasetIt] -> SetPointError(point,0., sqrt(numE/denE));
	//	g_recHitE_sigma_vsAbsEta_PU70ns[datasetIt] -> SetPointError(point,0., 0.);
	g_recHitE_sigma_vsAbsEta_PU70ns[datasetIt] -> SetPointError(point,0., PU70ContribE);

	//	std::cout << " >>> PU 70 - |eta| << " << binLowEdge << " " << binHigEdge << " " << PU70contrib << std::endl;	
	delete histo;
	//	delete histoRMS;
	delete histoRMS_AB;
	delete histoRMS_C;
	delete histoRMS_D;

	//PU140
sprintf(histoName,"%s/EB/plots_vs_nPU_and_absEta/h_recHitE_nPU139.0-141.0__absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
sprintf(histoName,"%s/EB/plots_vs_nPU_and_absEta/h_recHitE_sigma_nPU139.0-141.0__absEta%1.4f-%1.4f_EB_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
//      histoRMS = (TH1F*)( f->Get(histoName) );
        histoRMS_AB = (TH1F*)( fAB->Get(histoName) );
	histoRMS_C = (TH1F*)( fC->Get(histoName) );
	histoRMS_D = (TH1F*)( fD->Get(histoName) );
	Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.;
	EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
			 pow(Eabcd, 2) );

	//	if( histo->GetEntries() < 30 ) continue;
	//        histo -> Rebin(rebin);
	histo -> Rebin(2);

	float rmsPU140contrib = histo->GetRMS();
	float rmsPU140contrib2 = pow(rmsPU140contrib,2);      
        g_recHitE_RMS_vsAbsEta_PU140[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
	g_recHitE_RMS_vsAbsEta_PU140[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        g_recHitE_RMS_vsAbsEta_PU140ns[datasetIt]->SetPoint(point,binCenter,sqrt(rmsPU140contrib2 - rmsPU0contrib2));
        g_recHitE_RMS_vsAbsEta_PU140ns[datasetIt]->SetPointError(point,0.,histo->GetRMSError());
        
	float PU140contrib = Eabcd; //histoRMS->GetBinContent(1);
	float PU140ContribE = EabcdRMS; //histo->GetRMSError();
	float PU140contrib2 = pow(PU140contrib,2);

        g_recHitE_sigma_vsAbsEta_PU140[datasetIt] -> SetPoint(point,binCenter,PU140contrib);
        g_recHitE_sigma_vsAbsEta_PU140[datasetIt] -> SetPointError(point,0.,PU140ContribE);
        g_recHitE_sigma_vsAbsEta_PU140ns[datasetIt]->SetPoint(point,binCenter,sqrt(PU140contrib2 - PU0contrib2));
	numE = pow(PU140contrib * PU140ContribE, 2) + pow(PU0contrib * PU0ContribE, 2);
        denE = PU140ContribE - PU0contrib2;
	//        g_recHitE_sigma_vsAbsEta_PU140ns[datasetIt] -> SetPointError(point,0., sqrt(numE/denE));
	//        g_recHitE_sigma_vsAbsEta_PU140ns[datasetIt] -> SetPointError(point,0., 0.);
        g_recHitE_sigma_vsAbsEta_PU140ns[datasetIt] -> SetPointError(point,0., PU140ContribE);
        
	//	std::cout << " >>> PU 140 - |eta| << " << binLowEdge << " " << binHigEdge << " " << PU140contrib << std::endl;	
	delete histo;
	}
        ++point;
      }
      else
      {
        char histoName[50];
        sprintf(histoName,"%s/EE/plots_vs_absEta/h_recHitE_absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        TH1F* histo = (TH1F*)( f->Get(histoName) );
	sprintf(histoName,"%s/EE/plots_vs_absEta/h_recHitE_sigma_absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
//         TH1F* histoRMS = (TH1F*)( f->Get(histoName) );

	TH1F* histoRMS_AB = (TH1F*)( fAB->Get(histoName) );
	TH1F* histoRMS_C = (TH1F*)( fC->Get(histoName) );
        TH1F* histoRMS_D = (TH1F*)( fD->Get(histoName) );
	double Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.;
	double EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
                                pow(Eabcd, 2) );

	if( histo->GetEntries() < 30 ) continue;
	histo -> Rebin(rebin);
        
              
        g_recHitE_RMS_vsAbsEta[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsAbsEta[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        
	//	std::cout << " >>> all - |eta| << " << binLowEdge << " " << binHigEdge << " " << GetEffectiveSigma(histo) << std::endl;	
	float sigmaPUContrib = Eabcd; //histoRMS->GetBinContent(1);
	float sigmaPUContribE = EabcdRMS; histo->GetRMSError();
	float sigmaPUContrib2 = pow(sigmaPUContrib, 2);
        g_recHitE_sigma_vsAbsEta[datasetIt] -> SetPoint(point,binCenter,sigmaPUContrib);
	g_recHitE_sigma_vsAbsEta[datasetIt] -> SetPointError(point,0.,sigmaPUContribE);


	double x, y;
	float sigmaAt0PUPUContrib = g_recHitE_sigmaAt0PU_vsAbsEta[datasetIt]->GetPoint(point, x, y);
	double ey = g_recHitE_sigmaAt0PU_vsAbsEta[datasetIt]->GetErrorY(point);
	if(y < sigmaPUContrib) g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, sqrt(sigmaPUContrib2 - y*y));
	else g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, 0. );
	double numE = pow(sigmaPUContrib * sigmaPUContribE, 2) + pow(y*ey, 2);
	double denE = sigmaPUContrib2 - y*y;
	//	g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPointError(point,0., sqrt(numE/denE));
	//	g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPointError(point,0., 0.);
	g_recHitE_sigma_vsAbsEtans[datasetIt] -> SetPointError(point,0., sigmaPUContribE);
	std::cout << ">>>>>>>> EE error su sigma = " << sigmaPUContribE << " at PU0 = " << ey << std::endl;

	double xp, yp;
	float sigmaAt70PUPUContrib = g_recHitE_sigmaAt70PU_vsAbsEta[datasetIt]->GetPoint(point, xp, yp);
	double eyp = g_recHitE_sigmaAt70PU_vsAbsEta[datasetIt]->GetErrorY(point);
	if(y < yp) g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, sqrt(yp*yp - y*y));
	else g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, 0. );
	double numEp = pow(yp*eyp, 2) + pow(y*ey, 2);
	double denEp = yp*yp - y*y;
	//	g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,sqrt(numEp/denEp));
	//	g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,0.);
	g_recHitE_sigmaAt70PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,eyp);

	float sigmaAt140PUPUContrib = g_recHitE_sigmaAt140PU_vsAbsEta[datasetIt]->GetPoint(point, xp, yp);
	eyp = g_recHitE_sigmaAt140PU_vsAbsEta[datasetIt]->GetErrorY(point);
	if(y < yp) g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, sqrt(yp*yp - y*y));
	else g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPoint(point,binCenter, 0. );
	numEp = pow(yp*eyp, 2) + pow(y*ey, 2);
	denEp = yp*yp - y*y;
	//	g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,sqrt(numEp/denEp));
	//	g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,0.);
	g_recHitE_sigmaAt140PU_vsAbsEtans[datasetIt] -> SetPointError(point,0.,eyp);
        
	delete histo;
	delete histoRMS_AB;
	delete histoRMS_C;
	delete histoRMS_D;

	if(type == "MCMC"){
	//PU0
  sprintf(histoName,"%s/EE/plots_vs_nPU_and_absEta/h_recHitE_nPU-1.0-1.0__absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
  sprintf(histoName,"%s/EE/plots_vs_nPU_and_absEta/h_recHitE_sigma_nPU-1.0-1.0__absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
  //        TH1F* histoRMS = (TH1F*)( f->Get(histoName) );
        histoRMS_AB = (TH1F*)( fAB->Get(histoName) );
	histoRMS_C = (TH1F*)( fC->Get(histoName) );
	histoRMS_D = (TH1F*)( fD->Get(histoName) );
	Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.;
	EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
			 pow(Eabcd, 2) );

	//	if( histo->GetEntries() < 30 ) continue;
	//        histo -> Rebin(rebin);
        

	float rmsPU0contrib = histo->GetRMS();
	float rmsPU0contrib2 = pow(rmsPU0contrib,2);      
        g_recHitE_RMS_vsAbsEta_PU0[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsAbsEta_PU0[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        

	float PU0contrib = Eabcd; //histoRMS->GetBinContent(1);
	float PU0ContribE = EabcdRMS; //histo->GetRMSError();
        float PU0contrib2 = pow(PU0contrib,2);

        g_recHitE_sigma_vsAbsEta_PU0[datasetIt] -> SetPoint(point,binCenter,PU0contrib);
        g_recHitE_sigma_vsAbsEta_PU0[datasetIt] -> SetPointError(point,0.,PU0ContribE);
	//	std::cout << " >>> PU 0 - |eta| << " << binLowEdge << " " << binHigEdge << " " << PU0contrib << std::endl;	

	delete histo;
	//	delete histoRMS;
	delete histoRMS_AB;
	delete histoRMS_C;
	delete histoRMS_D;

	//PU70
  sprintf(histoName,"%s/EE/plots_vs_nPU_and_absEta/h_recHitE_nPU69.0-71.0__absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
  sprintf(histoName,"%s/EE/plots_vs_nPU_and_absEta/h_recHitE_sigma_nPU69.0-71.0__absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
  //        histoRMS = (TH1F*)( f->Get(histoName) );
        histoRMS_AB = (TH1F*)( fAB->Get(histoName) );
	histoRMS_C = (TH1F*)( fC->Get(histoName) );
	histoRMS_D = (TH1F*)( fD->Get(histoName) );
	Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.;
	EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
			 pow(Eabcd, 2) );

	//        if( histo->GetEntries() < 30 ) continue;
	//	histo -> Rebin(rebin);
        
	float rmsPU70contrib = histo->GetRMS();
	float rmsPU70contrib2 = pow(rmsPU70contrib,2);            
        g_recHitE_RMS_vsAbsEta_PU70[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsAbsEta_PU70[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        g_recHitE_RMS_vsAbsEta_PU70ns[datasetIt] -> SetPoint(point,binCenter, sqrt(rmsPU70contrib2 - rmsPU0contrib2));
        g_recHitE_RMS_vsAbsEta_PU70ns[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        

	float PU70contrib = Eabcd; //histoRMS->GetBinContent(1);
	float PU70ContribE = EabcdRMS; //histo->GetRMSError();
	float PU70contrib2 = pow(PU70contrib,2);
        g_recHitE_sigma_vsAbsEta_PU70[datasetIt] -> SetPoint(point,binCenter,PU70contrib);
        g_recHitE_sigma_vsAbsEta_PU70[datasetIt] -> SetPointError(point,0.,PU70ContribE);
        g_recHitE_sigma_vsAbsEta_PU70ns[datasetIt]->SetPoint(point,binCenter,sqrt(PU70contrib2 - PU0contrib2));
	numE = pow(PU70contrib * PU70ContribE, 2) + pow(PU0contrib * PU0ContribE, 2);
        denE = PU70ContribE - PU0contrib2;
	//        g_recHitE_sigma_vsAbsEta_PU70ns[datasetIt] -> SetPointError(point,0., sqrt(numE/denE));
	//        g_recHitE_sigma_vsAbsEta_PU70ns[datasetIt] -> SetPointError(point,0., 0.);
        g_recHitE_sigma_vsAbsEta_PU70ns[datasetIt] -> SetPointError(point,0., PU70ContribE);

	//	std::cout << " >>> PU 70 - |eta| << " << binLowEdge << " " << binHigEdge << " " << PU70contrib << std::endl;	
	delete histo;
	//	delete histoRMS;
	delete histoRMS_AB;
	delete histoRMS_C;
	delete histoRMS_D;

	//PU140
sprintf(histoName,"%s/EE/plots_vs_nPU_and_absEta/h_recHitE_nPU139.0-141.0__absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
        histo = (TH1F*)( f->Get(histoName) );
sprintf(histoName,"%s/EE/plots_vs_nPU_and_absEta/h_recHitE_sigma_nPU139.0-141.0__absEta%1.4f-%1.4f_EE_%s",dataset.c_str(),binLowEdge,binHigEdge,dataset.c_str());
//        histoRMS = (TH1F*)( f->Get(histoName) );
        histoRMS_AB = (TH1F*)( fAB->Get(histoName) );
	histoRMS_C = (TH1F*)( fC->Get(histoName) );
	histoRMS_D = (TH1F*)( fD->Get(histoName) );
	Eabcd = (histoRMS_AB->GetBinContent(1) + histoRMS_C->GetBinContent(1) + histoRMS_D->GetBinContent(1))/3.;
	EabcdRMS = sqrt( (pow(histoRMS_AB->GetBinContent(1), 2) + pow(histoRMS_C->GetBinContent(1), 2) + pow(histoRMS_D->GetBinContent(1), 2))/3. -
			 pow(Eabcd, 2) );

	//	if( histo->GetEntries() < 30 ) continue;
	//        histo -> Rebin(rebin);
	histo -> Rebin(2);

	float rmsPU140contrib = histo->GetRMS();
        float rmsPU140contrib2 = pow(rmsPU140contrib,2);
        g_recHitE_RMS_vsAbsEta_PU140[datasetIt] -> SetPoint(point,binCenter,histo->GetRMS());
        g_recHitE_RMS_vsAbsEta_PU140[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());
        g_recHitE_RMS_vsAbsEta_PU140ns[datasetIt] -> SetPoint(point,binCenter,sqrt(rmsPU140contrib2 - rmsPU0contrib2));
        g_recHitE_RMS_vsAbsEta_PU140ns[datasetIt] -> SetPointError(point,0.,histo->GetRMSError());


	float PU140contrib = Eabcd; //histoRMS->GetBinContent(1);
	float PU140ContribE = EabcdRMS; //histo->GetRMSError();
	float PU140contrib2 = pow(PU140contrib,2);        
        g_recHitE_sigma_vsAbsEta_PU140[datasetIt] -> SetPoint(point,binCenter,PU140contrib);
        g_recHitE_sigma_vsAbsEta_PU140[datasetIt] -> SetPointError(point,0.,PU140ContribE);
        g_recHitE_sigma_vsAbsEta_PU140ns[datasetIt]->SetPoint(point,binCenter,sqrt(PU140contrib2 - PU0contrib2));
	numE = pow(PU140contrib * PU140ContribE, 2) + pow(PU0contrib * PU0ContribE, 2);
        denE = PU140ContribE - PU0contrib2;
	//        g_recHitE_sigma_vsAbsEta_PU140ns[datasetIt] -> SetPointError(point,0., sqrt(numE/denE));
	//        g_recHitE_sigma_vsAbsEta_PU140ns[datasetIt] -> SetPointError(point,0., 0.);
        g_recHitE_sigma_vsAbsEta_PU140ns[datasetIt] -> SetPointError(point,0., PU140ContribE);

	//	std::cout << " >>> PU 140 - |eta| << " << binLowEdge << " " << binHigEdge << " " << PU140contrib << std::endl;	
        //histo->GetBinWidth(1);
	delete histo;
	}
        ++point;      
      }
    }

  }
  
  
  std::cout << " >>>>>>>>>>>>>> DRAW " << std::endl;
  
  char canvasName[50];
  
//   std::string folderName = baseDir+"/output/noisePlots__" + label + "/";
//   system(("mkdir "+folderName).c_str());
  
  leg = NULL;
  
  
  //--------------------------
  // draw inclusive histograms
  sprintf(canvasName,"cEB_recHitE");
  TCanvas* c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  hEB_recHitE[0] -> Draw("P");
  hEB_recHitE[1] -> Draw("hist,same");
  SetHistoStyle(hEB_recHitE[0],"DA","EB","recHitE",&leg);
  SetHistoStyle(hEB_recHitE[1],"MC","EB","recHitE",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"hEB_recHitE.png").c_str(),"png");
  delete c;
  leg = NULL;
  
  sprintf(canvasName,"cEB_recHitE_log");
  TCanvas* clog = new TCanvas(canvasName,canvasName);
  clog -> SetLogy();
  clog -> SetGridx();
  clog -> SetGridy();
  hEB_recHitE[0] -> Draw("P");
  hEB_recHitE[1] -> Draw("hist,same");  
  SetHistoStyle(hEB_recHitE[0],"DA","EB","recHitE",&leg,true);
  SetHistoStyle(hEB_recHitE[1],"MC","EB","recHitE",&leg,true);
  leg -> Draw("same");
  clog -> Print((folderName+"hEB_recHitE_log.png").c_str(),"png");
  delete clog;
  leg = NULL;
  
  
  sprintf(canvasName,"cEB_recHitPedSubADC");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  hEB_recHitPedSubADC[0] -> Draw("P");
  hEB_recHitPedSubADC[1] -> Draw("hist,same");
  SetHistoStyle(hEB_recHitPedSubADC[0],"DA","EB","recHitPedSubADC",&leg);
  SetHistoStyle(hEB_recHitPedSubADC[1],"MC","EB","recHitPedSubADC",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"hEB_recHitPedSubADC.png").c_str(),"png");
  delete c;
  leg = NULL;
  
  sprintf(canvasName,"cEB_recHitPedSubADC_log");
  clog = new TCanvas(canvasName,canvasName);
  clog -> SetLogy();
  clog -> SetGridx();
  clog -> SetGridy();
  hEB_recHitPedSubADC[0] -> Draw("P");
  hEB_recHitPedSubADC[1] -> Draw("hist,same");  
  SetHistoStyle(hEB_recHitPedSubADC[0],"DA","EB","recHitPedSubADC",&leg,true);
  SetHistoStyle(hEB_recHitPedSubADC[1],"MC","EB","recHitPedSubADC",&leg,true);
  leg -> Draw("same");
  clog -> Print((folderName+"hEB_recHitPedSubADC_log.png").c_str(),"png");
  delete clog;
  leg = NULL;
  
  
  sprintf(canvasName,"cEE_recHitE");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  hEE_recHitE[0] -> Draw("P");
  hEE_recHitE[1] -> Draw("hist,same");
  SetHistoStyle(hEE_recHitE[0],"DA","EE","recHitE",&leg);
  SetHistoStyle(hEE_recHitE[1],"MC","EE","recHitE",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"hEE_recHitE.png").c_str(),"png");
  delete c;
  leg = NULL;
  
  sprintf(canvasName,"cEE_recHitE_log");
  clog = new TCanvas(canvasName,canvasName);
  clog -> SetLogy();
  clog -> SetGridx();
  clog -> SetGridy();
  hEE_recHitE[0] -> Draw("P");
  hEE_recHitE[1] -> Draw("hist,same");  
  SetHistoStyle(hEE_recHitE[0],"DA","EE","recHitE",&leg,true);
  SetHistoStyle(hEE_recHitE[1],"MC","EE","recHitE",&leg,true);
  leg -> Draw("same");
  clog -> Print((folderName+"hEE_recHitE_log.png").c_str(),"png");
  delete clog;
  leg = NULL;
  
  
  sprintf(canvasName,"cEE_recHitPedSubADC");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  hEE_recHitPedSubADC[0] -> Draw("P");
  hEE_recHitPedSubADC[1] -> Draw("hist,same");
  SetHistoStyle(hEE_recHitPedSubADC[0],"DA","EE","recHitPedSubADC",&leg);
  SetHistoStyle(hEE_recHitPedSubADC[1],"MC","EE","recHitPedSubADC",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"hEE_recHitPedSubADC.png").c_str(),"png");
  delete c;
  leg = NULL;
  
  sprintf(canvasName,"cEE_recHitPedSubADC_log");
  clog = new TCanvas(canvasName,canvasName);
  clog -> SetLogy();
  clog -> SetGridx();
  clog -> SetGridy();
  hEE_recHitPedSubADC[0] -> Draw("P");
  hEE_recHitPedSubADC[1] -> Draw("hist,same");  
  SetHistoStyle(hEE_recHitPedSubADC[0],"DA","EE","recHitPedSubADC",&leg,true);
  SetHistoStyle(hEE_recHitPedSubADC[1],"MC","EE","recHitPedSubADC",&leg,true);
  leg -> Draw("same");
  clog -> Print((folderName+"hEE_recHitPedSubADC_log.png").c_str(),"png");
  delete clog;
  leg = NULL;
  
  
  std::cout << " >>>>>>>>>>>>>> TFILE " << std::endl;

  TFile outPlots(("outPlots_"+type+".root").c_str(),"recreate");
  outPlots.cd();  

  //----------------------
  // draw exclusive graphs
    
  sprintf(canvasName,"cEB_recHitE_RMS_vsNPU");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  gEB_recHitE_RMS_vsNPU[0] -> Draw("AP");
  if(type == "DAMC") gEB_recHitE_RMS_vsNPU[1] -> Draw("P,same");
  SetGraphStyle(gEB_recHitE_RMS_vsNPU[0],"DA","EB","vsNPU","recHitE_RMS",&leg);
  if(type == "DAMC") SetGraphStyle(gEB_recHitE_RMS_vsNPU[1],"MC","EB","vsNPU","recHitE_RMS",&leg);
  //  if(type == "DAMC") leg -> Draw("same");
  if(type == "DAMC") {
    gEB_recHitE_RMS_vsNPU[0]->GetXaxis()->SetRangeUser(0.,50.);
    gEB_recHitE_RMS_vsNPU[0]->SetMaximum(0.5);
  }
  c -> Print((folderName+"gEB_recHitE_RMS_vsNPU.png").c_str(),"png");
  delete c;
  leg = NULL;
  gEB_recHitE_RMS_vsNPU[0]->Write("gEB_recHitE_RMS_vsNPU_DA");
  gEB_recHitE_RMS_vsNPU[1]->Write("gEB_recHitE_RMS_vsNPU_MC");

  sprintf(canvasName,"cEE_recHitE_RMS_vsNPU");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  gEE_recHitE_RMS_vsNPU[0] -> Draw("AP");
  if(type == "DAMC") gEE_recHitE_RMS_vsNPU[1] -> Draw("P, same");
  SetGraphStyle(gEE_recHitE_RMS_vsNPU[0],"DA","EE","vsNPU","recHitE_RMS",&leg);
  if(type == "DAMC") SetGraphStyle(gEE_recHitE_RMS_vsNPU[1],"MC","EE","vsNPU","recHitE_RMS",&leg);
  if(type == "DAMC") {
    gEE_recHitE_RMS_vsNPU[0]->GetXaxis()->SetRangeUser(0.,50.);
    gEE_recHitE_RMS_vsNPU[0]->SetMaximum(0.5);
  }
  //  if(type == "DAMC") leg -> Draw("same");
  c -> Print((folderName+"gEE_recHitE_RMS_vsNPU.png").c_str(),"png");
  delete c;
  leg = NULL;
  gEE_recHitE_RMS_vsNPU[0]->Write("gEE_recHitE_RMS_vsNPU_DA");
  gEE_recHitE_RMS_vsNPU[1]->Write("gEE_recHitE_RMS_vsNPU_MC");  
  
  sprintf(canvasName,"cEB_recHitE_sigma_vsNPU");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  gEB_recHitE_sigma_vsNPU[0] -> Draw("AP");
  if(type == "DAMC") gEB_recHitE_sigma_vsNPU[1] -> Draw("P,same");
  SetGraphStyle(gEB_recHitE_sigma_vsNPU[0],"DA","EB","vsNPU","recHitE_sigma",&leg);
  if(type == "DAMC") SetGraphStyle(gEB_recHitE_sigma_vsNPU[1],"MC","EB","vsNPU","recHitE_sigma",&leg);
  if(type == "DAMC") {
    gEB_recHitE_sigma_vsNPU[0]->GetXaxis()->SetRangeUser(0.,50.);
    gEB_recHitE_sigma_vsNPU[0]->SetMaximum(0.5);
  }
  //  if(type == "DAMC") leg -> Draw("same");
  c -> Print((folderName+"gEB_recHitE_sigma_vsNPU.png").c_str(),"png");
  delete c;
  leg = NULL;
  gEB_recHitE_sigma_vsNPU[0]->Write("gEB_recHitE_sigma_vsNPU_DA");
  gEB_recHitE_sigma_vsNPU[1]->Write("gEB_recHitE_sigma_vsNPU_MC");
  
  sprintf(canvasName,"cEE_recHitE_sigma_vsNPU");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  gEE_recHitE_sigma_vsNPU[0] -> Draw("AP");
  if(type == "DAMC") gEE_recHitE_sigma_vsNPU[1] -> Draw("P,same");
  SetGraphStyle(gEE_recHitE_sigma_vsNPU[0],"DA","EE","vsNPU","recHitE_sigma",&leg);
  if(type == "DAMC") SetGraphStyle(gEE_recHitE_sigma_vsNPU[1],"MC","EE","vsNPU","recHitE_sigma",&leg);
  if(type == "DAMC"){
    gEE_recHitE_sigma_vsNPU[0]->GetXaxis()->SetRangeUser(0.,50.);
    gEE_recHitE_sigma_vsNPU[0]->SetMaximum(1.2);
  }
  //  if(type == "DAMC") leg -> Draw("same");
  c -> Print((folderName+"gEE_recHitE_sigma_vsNPU.png").c_str(),"png");
  delete c;
  leg = NULL;
  gEE_recHitE_sigma_vsNPU[0]->Write("gEE_recHitE_sigma_vsNPU_DA");
  gEE_recHitE_sigma_vsNPU[1]->Write("gEE_recHitE_sigma_vsNPU_MC");
    
  sprintf(canvasName,"c_recHitE_RMS_vsIRing");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
 c -> SetGridy();
  g_recHitE_RMS_vsIRing_PU0[1] -> Draw("AP");
  g_recHitE_RMS_vsIRing_PU70[1] -> Draw("P,same");
  g_recHitE_RMS_vsIRing_PU140[1] -> Draw("P,same");
  if(type == "DAMC"){
    g_recHitE_RMS_vsIRing[0] -> Draw("P,same");
    g_recHitE_RMS_vsIRing[1] -> Draw("P,same");
    SetGraphStyle(g_recHitE_RMS_vsIRing[0],"DA","EBEE","vsIRing","recHitE_RMS",&leg);
    SetGraphStyle(g_recHitE_RMS_vsIRing[1],"MC","EBEE","vsIRing","recHitE_RMS",&leg);
  }
  if(type == "MCMC"){
    SetGraphStyle(g_recHitE_RMS_vsIRing_PU0[1],"MC_PU0","EBEE","vsIRing","recHitE_RMS",&leg);
    SetGraphStyle(g_recHitE_RMS_vsIRing_PU70[1],"MC_PU70","EBEE","vsIRing","recHitE_RMS",&leg);
    SetGraphStyle(g_recHitE_RMS_vsIRing_PU140[1],"MC_PU140","EBEE","vsIRing","recHitE_RMS",&leg);
  }
  gNuGun_recHitE_RMS_vsIRing -> Draw("L,same");
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_RMS_vsIRing.png").c_str(),"png");
  delete c;
  leg = NULL;
  if(type == "DAMC"){
    g_recHitE_RMS_vsIRing[0]->Write("g_recHitE_RMS_vsIRing_DA"); 
    g_recHitE_RMS_vsIRing[1]->Write("g_recHitE_RMS_vsIRing_MC"); 
  }
  g_recHitE_RMS_vsIRing_PU0[1]->Write("g_recHitE_RMS_vsIRing_PU0_MC");
  g_recHitE_RMS_vsIRing_PU70[1]->Write("g_recHitE_RMS_vsIRing_PU70_MC");
  g_recHitE_RMS_vsIRing_PU140[1]->Write("g_recHitE_RMS_vsIRing_PU140_MC");


  sprintf(canvasName,"c_recHitE_RMS_vsIRing_log");
  c = new TCanvas(canvasName,canvasName);
  c->SetLogy();
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_RMS_vsIRing_PU0[1] -> Draw("AP");
  g_recHitE_RMS_vsIRing_PU70[1] -> Draw("P,same");
  g_recHitE_RMS_vsIRing_PU140[1] -> Draw("P,same");
  if(type == "DAMC"){
    g_recHitE_RMS_vsIRing[0] -> Draw("P,same");
    g_recHitE_RMS_vsIRing[1] -> Draw("P,same");
    SetGraphStyle(g_recHitE_RMS_vsIRing[0],"DA","EBEE","vsIRing","recHitE_RMS",&leg);
    SetGraphStyle(g_recHitE_RMS_vsIRing[1],"MC","EBEE","vsIRing","recHitE_RMS",&leg);
  }
  if(type == "MCMC"){
    SetGraphStyle(g_recHitE_RMS_vsIRing_PU0[1],"MC_PU0","EBEE","vsIRing","recHitE_RMS",&leg);
    SetGraphStyle(g_recHitE_RMS_vsIRing_PU70[1],"MC_PU70","EBEE","vsIRing","recHitE_RMS",&leg);
    SetGraphStyle(g_recHitE_RMS_vsIRing_PU140[1],"MC_PU140","EBEE","vsIRing","recHitE_RMS",&leg);
  }
  gNuGun_recHitE_RMS_vsIRing -> Draw("L,same");
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_RMS_vsIRing_log.png").c_str(),"png");
  delete c;
  leg = NULL;


  sprintf(canvasName,"c_recHitE_sigma_vsAbsEta");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_sigma_vsAbsEta_PU0[1] -> Draw("AP");
  g_recHitE_sigma_vsAbsEta_PU70[1] -> Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU140[1] -> Draw("P,same");
  if(type == "DAMC"){
    g_recHitE_sigma_vsAbsEta[0] -> Draw("P,same");
    g_recHitE_sigma_vsAbsEta[1] -> Draw("P,same");
    SetGraphStyle(g_recHitE_sigma_vsAbsEta[0],"DA","EBEE","vsAbsEta","recHitE_sigma",&leg);
    SetGraphStyle(g_recHitE_sigma_vsAbsEta[1],"MC","EBEE","vsAbsEta","recHitE_sigma",&leg);
  }
  if(type == "MCMC"){
    SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU0[1],"MC_PU0","EBEE","vsAbsEta","recHitE_sigma",&leg);
    SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU70[1],"MC_PU70","EBEE","vsAbsEta","recHitE_sigma",&leg);
    SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU140[1],"MC_PU140","EBEE","vsAbsEta","recHitE_sigma",&leg);
  }
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigma_vsAbsEta.png").c_str(),"png");
  delete c;
  leg = NULL;
  if(type == "DAMC"){
    g_recHitE_sigma_vsAbsEta[0]->Write("g_recHitE_sigma_vsAbsEta_DA"); 
    g_recHitE_sigma_vsAbsEta[1]->Write("g_recHitE_sigma_vsAbsEta_MC"); 
  }
  g_recHitE_sigma_vsAbsEta_PU0[1]->Write("g_recHitE_sigma_vsAbsEta_PU0_MC");
  g_recHitE_sigma_vsAbsEta_PU70[1]->Write("g_recHitE_sigma_vsAbsEta_PU70_MC");
  g_recHitE_sigma_vsAbsEta_PU140[1]->Write("g_recHitE_sigma_vsAbsEta_PU140_MC");
  g_recHitE_sigma_vsAbsEtans[0]->Write("g_recHitE_sigma_vsAbsEtans_DA");
  g_recHitE_sigma_vsAbsEtans[1]->Write("g_recHitE_sigma_vsAbsEtans_MC");
  g_recHitE_sigmaAt70PU_vsAbsEtans[0]->Write("g_recHitE_sigmaAt70PU_vsAbsEtans_DA");
  g_recHitE_sigmaAt70PU_vsAbsEtans[1]->Write("g_recHitE_sigmaAt70PU_vsAbsEtans_MC");
  g_recHitE_sigmaAt140PU_vsAbsEtans[0]->Write("g_recHitE_sigmaAt140PU_vsAbsEtans_DA");
  g_recHitE_sigmaAt140PU_vsAbsEtans[1]->Write("g_recHitE_sigmaAt140PU_vsAbsEtans_MC");

  
  sprintf(canvasName,"c_recHitE_RMSAt0PU_vsIRing");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_RMSAt0PU_vsIRing[0] -> Draw("AP");
  g_recHitE_RMSAt0PU_vsIRing[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_RMSAt0PU_vsIRing[0],"DA","EBEE","vsIRing","recHitE_RMSAt0PU",&leg);
  SetGraphStyle(g_recHitE_RMSAt0PU_vsIRing[1],"MC","EBEE","vsIRing","recHitE_RMSAt0PU",&leg);
  gNuGun_recHitE_RMS_vsIRing -> Draw("L,same");
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_RMSAt0PU_vsIRing.png").c_str(),"png");
  delete c;
  leg = NULL;  
  g_recHitE_RMSAt0PU_vsIRing[0]->Write("g_recHitE_RMSAt0PU_vsIRing_DA");
  g_recHitE_RMSAt0PU_vsIRing[1]->Write("g_recHitE_RMSAt0PU_vsIRing_MC");


  
  sprintf(canvasName,"c_recHitE_sigmaAt0PU_vsAbsEta");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_sigmaAt0PU_vsAbsEta[0] -> Draw("AP");
  g_recHitE_sigmaAt0PU_vsAbsEta[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_sigmaAt0PU_vsAbsEta[0],"DA","EBEE","vsAbsEta","recHitE_sigmaAt0PU",&leg);
  SetGraphStyle(g_recHitE_sigmaAt0PU_vsAbsEta[1],"MC","EBEE","vsAbsEta","recHitE_sigmaAt0PU",&leg);
  gNuGun_recHitE_RMS_vsIRing -> Draw("L,same");
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigmaAt0PU_vsAbsEta.png").c_str(),"png");
  delete c;
  leg = NULL;  
  g_recHitE_sigmaAt0PU_vsAbsEta[0]->Write("g_recHitE_sigmaAt0PU_vsAbsEta_DA");
  g_recHitE_sigmaAt0PU_vsAbsEta[1]->Write("g_recHitE_sigmaAt0PU_vsAbsEta_MC");
  g_recHitE_sigmaAt70PU_vsAbsEta[0]->Write("g_recHitE_sigmaAt70PU_vsAbsEta_DA");
  g_recHitE_sigmaAt70PU_vsAbsEta[1]->Write("g_recHitE_sigmaAt70PU_vsAbsEta_MC");
  g_recHitE_sigmaAt140PU_vsAbsEta[0]->Write("g_recHitE_sigmaAt140PU_vsAbsEta_DA");
  g_recHitE_sigmaAt140PU_vsAbsEta[1]->Write("g_recHitE_sigmaAt140PU_vsAbsEta_MC");

  //rms
  sprintf(canvasName,"c_recHitE_RMS_vsEta");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_RMS_vsEta[0] -> Draw("AP");
  g_recHitE_RMS_vsEta[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_RMS_vsEta[0],"DA","EBEE","vsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsEta[1],"MC","EBEE","vsEta","recHitE_RMS",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_RMS_vsEta.png").c_str(),"png");
  delete c;
  leg = NULL;
  g_recHitE_RMS_vsEta[0]->Write("g_recHitE_RMS_vsEta_DA");
  g_recHitE_RMS_vsEta[1]->Write("g_recHitE_RMS_vsEta_MC");


  sprintf(canvasName,"c_recHitE_RMS_vsAbsEta");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_RMS_vsAbsEta[0] -> Draw("AP");
  g_recHitE_RMS_vsAbsEta[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_RMS_vsAbsEta[0],"DA","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta[1],"MC","EBEE","vsAbsEta","recHitE_RMS",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_RMS_vsAbsEta.png").c_str(),"png");
  delete c;
  leg = NULL;
  g_recHitE_RMS_vsAbsEta[0]->Write("g_recHitE_RMS_vsAbsEta_DA");
  g_recHitE_RMS_vsAbsEta[1]->Write("g_recHitE_RMS_vsAbsEta_MC");

  sprintf(canvasName,"c_recHitE_RMS_vsEta_vsPU");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_RMS_vsEta_PU0[1] -> Draw("AP");
  g_recHitE_RMS_vsEta_PU70[1] -> Draw("P,same");
  g_recHitE_RMS_vsEta_PU140[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_RMS_vsEta_PU0[1],"MC_PU0","EBEE","vsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsEta_PU70[1],"MC_PU70","EBEE","vsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsEta_PU140[1],"MC_PU140","EBEE","vsEta","recHitE_RMS",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_RMS_vsEta_vsPU.png").c_str(),"png");
  delete c;
  leg = NULL;
  g_recHitE_RMS_vsEta_PU0[1]->Write("g_recHitE_RMS_vsEta_PU0_MC");
  g_recHitE_RMS_vsEta_PU70[1]->Write("g_recHitE_RMS_vsEta_PU70_MC");
  g_recHitE_RMS_vsEta_PU140[1]->Write("g_recHitE_RMS_vsEta_PU140_MC");

  sprintf(canvasName,"c_recHitE_RMS_vsAbsEta_vsPU");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_RMS_vsAbsEta_PU0[1] -> Draw("AP");
  g_recHitE_RMS_vsAbsEta_PU70[1] -> Draw("P,same");
  g_recHitE_RMS_vsAbsEta_PU140[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU0[1],"MC_PU0","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU70[1],"MC_PU70","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU140[1],"MC_PU140","EBEE","vsAbsEta","recHitE_RMS",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"c_recHitE_RMS_vsAbsEta_vsPU.png").c_str(),"png");
  delete c;
  leg = NULL;
  g_recHitE_RMS_vsAbsEta_PU0[1]->Write("g_recHitE_RMS_vsAbsEta_PU0_MC");
  g_recHitE_RMS_vsAbsEta_PU70[1]->Write("g_recHitE_RMS_vsAbsEta_PU70_MC");
  g_recHitE_RMS_vsAbsEta_PU140[1]->Write("g_recHitE_RMS_vsAbsEta_PU140_MC");

  /// sigma
  sprintf(canvasName,"c_recHitE_sigma_vsEta");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_sigma_vsEta[0] -> Draw("AP");
  g_recHitE_sigma_vsEta[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_sigma_vsEta[0],"DA","EBEE","vsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsEta[1],"MC","EBEE","vsEta","recHitE_sigma",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigma_vsEta.png").c_str(),"png");
  delete c;
  leg = NULL;
  g_recHitE_sigma_vsEta[0]->Write("g_recHitE_sigma_vsEta_DA");
  g_recHitE_sigma_vsEta[1]->Write("g_recHitE_sigma_vsEta_MC");

  sprintf(canvasName,"c_recHitE_sigma_vsEta_log");
  c = new TCanvas(canvasName,canvasName);
  c->SetLogy();
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_sigma_vsEta[0] -> Draw("AP");
  g_recHitE_sigma_vsEta[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_sigma_vsEta[0],"DA","EBEE","vsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsEta[1],"MC","EBEE","vsEta","recHitE_sigma",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigma_vsEta_log.png").c_str(),"png");
  delete c;
  leg = NULL;

  sprintf(canvasName,"c_recHitE_sigma_vsAbsEta");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_sigma_vsAbsEta[0] -> Draw("AP");
  g_recHitE_sigma_vsAbsEta[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_sigma_vsAbsEta[0],"DA","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta[1],"MC","EBEE","vsAbsEta","recHitE_sigma",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigma_vsAbsEta.png").c_str(),"png");
  delete c;
  leg = NULL;
  g_recHitE_sigma_vsAbsEta[0]->Write("g_recHitE_sigma_vsAbsEta_DA");
  g_recHitE_sigma_vsAbsEta[1]->Write("g_recHitE_sigma_vsAbsEta_MC");


  sprintf(canvasName,"c_recHitE_sigma_vsEta_vsPU");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_sigma_vsEta_PU0[1] -> Draw("AP");
  g_recHitE_sigma_vsEta_PU70[1] -> Draw("P,same");
  g_recHitE_sigma_vsEta_PU140[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_sigma_vsEta_PU0[1],"MC_PU0","EBEE","vsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsEta_PU70[1],"MC_PU70","EBEE","vsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsEta_PU140[1],"MC_PU140","EBEE","vsEta","recHitE_sigma",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigma_vsEta_vsPU.png").c_str(),"png");
  delete c;
  leg = NULL;
  g_recHitE_sigma_vsEta_PU0[1]->Write("g_recHitE_sigma_vsEta_PU0_MC");
  g_recHitE_sigma_vsEta_PU70[1]->Write("g_recHitE_sigma_vsEta_PU70_MC");
  g_recHitE_sigma_vsEta_PU140[1]->Write("g_recHitE_sigma_vsEta_PU140_MC");


  sprintf(canvasName,"c_recHitE_RMS_vsAbsEta_vsPUns");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_RMS_vsAbsEta_PU0[1] -> Draw("AP");
  g_recHitE_RMS_vsAbsEta_PU70[1] -> Draw("P,same");
  g_recHitE_RMS_vsAbsEta_PU140[1] -> Draw("P,same");
  g_recHitE_RMS_vsAbsEta_PU70ns[1] -> Draw("P,same");
  g_recHitE_RMS_vsAbsEta_PU140ns[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU0[1],"MC_PU0","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU70[1],"MC_PU70","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU140[1],"MC_PU140","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU70ns[1],"MC_PU70ns","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU140ns[1],"MC_PU140ns","EBEE","vsAbsEta","recHitE_RMS",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"c_recHitE_RMS_vsAbsEta_vsPUns.png").c_str(),"png");
  delete c;
  leg = NULL;
  g_recHitE_RMS_vsAbsEta_PU0[1]->Write("g_recHitE_RMS_vsAbsEta_PU0_MC");
  g_recHitE_RMS_vsAbsEta_PU70[1]->Write("g_recHitE_RMS_vsAbsEta_PU70_MC");
  g_recHitE_RMS_vsAbsEta_PU140[1]->Write("g_recHitE_RMS_vsAbsEta_PU140_MC");
  g_recHitE_RMS_vsAbsEta_PU70ns[1]->Write("g_recHitE_RMS_vsAbsEta_PU70ns_MC");
  g_recHitE_RMS_vsAbsEta_PU140ns[1]->Write("g_recHitE_RMS_vsAbsEta_PU140ns_MC");


  sprintf(canvasName,"c_recHitE_RMS_vsAbsEta_vsPUns_log");
  c = new TCanvas(canvasName,canvasName);
  c->SetLogy();
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_RMS_vsAbsEta_PU0[1] -> Draw("AP");
  g_recHitE_RMS_vsAbsEta_PU70[1] -> Draw("P,same");
  g_recHitE_RMS_vsAbsEta_PU140[1] -> Draw("P,same");
  g_recHitE_RMS_vsAbsEta_PU70ns[1] -> Draw("P,same");
  g_recHitE_RMS_vsAbsEta_PU140ns[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU0[1],"MC_PU0","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU70[1],"MC_PU70","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU140[1],"MC_PU140","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU70ns[1],"MC_PU70ns","EBEE","vsAbsEta","recHitE_RMS",&leg);
  SetGraphStyle(g_recHitE_RMS_vsAbsEta_PU140ns[1],"MC_PU140ns","EBEE","vsAbsEta","recHitE_RMS",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"c_recHitE_RMS_vsAbsEta_vsPUns_log.png").c_str(),"png");
  delete c;
  leg = NULL;


  /////////////////
  sprintf(canvasName,"c_recHitE_sigma_vsAbsEta_vsPU");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_sigma_vsAbsEta_PU0[1] -> Draw("AP");
  g_recHitE_sigma_vsAbsEta_PU70[1] -> Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU140[1] -> Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU70ns[1] -> Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU140ns[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU0[1],"MC_PU0","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU70[1],"MC_PU70","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU140[1],"MC_PU140","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU70ns[1],"MC_PU70ns","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU140ns[1],"MC_PU140ns","EBEE","vsAbsEta","recHitE_sigma",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigma_vsAbsEta_vsPU.png").c_str(),"png");
  delete c;
  leg = NULL;
  g_recHitE_sigma_vsAbsEta_PU0[1]->Write("g_recHitE_sigma_vsAbsEta_PU0_MC");
  g_recHitE_sigma_vsAbsEta_PU70[1]->Write("g_recHitE_sigma_vsAbsEta_PU70_MC");
  g_recHitE_sigma_vsAbsEta_PU140[1]->Write("g_recHitE_sigma_vsAbsEta_PU140_MC");
  g_recHitE_sigma_vsAbsEta_PU70ns[1]->Write("g_recHitE_sigma_vsAbsEta_PU70ns_MC");
  g_recHitE_sigma_vsAbsEta_PU140ns[1]->Write("g_recHitE_sigma_vsAbsEta_PU140ns_MC");


  sprintf(canvasName,"c_recHitE_sigma_vsAbsEta_vsPU_log");
  c = new TCanvas(canvasName,canvasName);
  c->SetLogy();
  c -> SetGridx();
  c -> SetGridy();
  g_recHitE_sigma_vsAbsEta_PU0[1] -> Draw("AP");
  g_recHitE_sigma_vsAbsEta_PU70[1] -> Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU140[1] -> Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU70ns[1] -> Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU140ns[1] -> Draw("P,same");
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU0[1],"MC_PU0","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU70[1],"MC_PU70","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU140[1],"MC_PU140","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU70ns[1],"MC_PU70ns","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU140ns[1],"MC_PU140ns","EBEE","vsAbsEta","recHitE_sigma",&leg);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigma_vsAbsEta_vsPU_log.png").c_str(),"png");
  delete c;
  leg = NULL;


  sprintf(canvasName,"c_recHitE_sigma_vsAbsEta_vsPU_comparison_log");
  c = new TCanvas(canvasName,canvasName);
  c->SetLogy();
  c -> SetGridx();
  c -> SetGridy();
  ax_PU70->Draw("AL");
  ax_PU140->Draw("L,same");
  g_recHitE_sigma_vsAbsEta_PU70ns[1] -> Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU140ns[1] -> Draw("P,same");
  SetGraphStyle(ax_PU70,"MC_PU70_ax","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(ax_PU140,"MC_PU140_ax","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU70ns[1],"MC_PU70nsc","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU140ns[1],"MC_PU140nsc","EBEE","vsAbsEta","recHitE_sigma",&leg);
  ax_PU70->SetMaximum(10.);
  ax_PU140->SetMaximum(10.);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigma_vsAbsEta_vsPU_comparison_log.png").c_str(),"png");
  delete c;
  leg = NULL;

  sprintf(canvasName,"c_recHitE_sigma_vsAbsEta_vsPU_comparison");
  c = new TCanvas(canvasName,canvasName);
  c -> SetGridx();
  c -> SetGridy();
  ax_PU70->Draw("AP");
  ax_PU140->Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU70ns[1] -> Draw("P,same");
  g_recHitE_sigma_vsAbsEta_PU140ns[1] -> Draw("P,same");
  SetGraphStyle(ax_PU70,"MC_PU70_ax","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(ax_PU140,"MC_PU140_ax","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU70ns[1],"MC_PU70nsc","EBEE","vsAbsEta","recHitE_sigma",&leg);
  SetGraphStyle(g_recHitE_sigma_vsAbsEta_PU140ns[1],"MC_PU140nsc","EBEE","vsAbsEta","recHitE_sigma",&leg);
  ax_PU70->SetMaximum(10.);
  ax_PU140->SetMaximum(10.);
  leg -> Draw("same");
  c -> Print((folderName+"g_recHitE_sigma_vsAbsEta_vsPU_comparison.png").c_str(),"png");
  delete c;
  leg = NULL;

  outPlots.Close();
}



void SetHistoStyle(TH1F* h, const std::string& dataset, const std::string& region, const std::string& type, TLegend** leg, const bool& log)
{
  if( (*leg) == NULL )
  {
    (*leg) = new TLegend(0.72,0.14,0.82,0.25,NULL,"brNDC");
    (*leg)->SetTextSize(0.04);
    (*leg)->SetLineWidth(0);
    (*leg)->SetLineColor(kWhite);
    (*leg)->SetTextFont(42);
    (*leg)->SetFillColor(kWhite);
    (*leg)->SetFillStyle(1001);
    (*leg)->SetBorderSize(1);
  }
  
  if( dataset == "DA" )
  {
    (*leg) -> AddEntry(h,"data","P");
    
    h -> SetMarkerStyle(20);
    h -> SetMarkerColor(kBlack);
    h -> SetMarkerSize(0.7);
    if(!log) h -> SetMinimum(0.);
    if(log)  h -> SetMinimum(0.00000001);
    h -> SetMaximum(1.1*h->GetMaximum());
    
    TPaveText* text = new TPaveText(0.50,0.85,0.80,0.90,"brNDC");
    text -> SetTextSize(0.04);
    text -> SetLineWidth(0);
    text -> SetLineColor(kWhite);
    text -> SetTextFont(42);
    text -> SetTextColor(kBlack);    
    text -> SetFillColor(kWhite);
    text -> SetFillStyle(1001);
    text -> SetBorderSize(1);
    char RMSText[50]; sprintf(RMSText,"RMS = %1.2e #pm %1.0e",h->GetRMS(),h->GetRMSError());
    text -> AddText(RMSText);
    text -> Draw();
  }
  if( dataset == "MC" )
  {
    (*leg) -> AddEntry(h,"MC","F");
    
    h -> SetMarkerStyle(20);
    h -> SetMarkerColor(kRed);
    h -> SetMarkerSize(0.7);
    if(!log) h -> SetMinimum(0.);
    if(log)  h -> SetMinimum(0.00000001);
    h -> SetMaximum(1.1*h->GetMaximum());
    h -> SetLineColor(kRed);
    h -> SetFillColor(kRed);
    h -> SetFillStyle(3003);
    
    TPaveText* text = new TPaveText(0.50,0.80,0.80,0.85,"brNDC");
    text -> SetTextSize(0.04);
    text -> SetLineWidth(0);
    text -> SetLineColor(kWhite);
    text -> SetTextFont(42);
    text -> SetTextColor(kRed);
    text -> SetFillColor(kWhite);
    text -> SetFillStyle(1001);
    text -> SetBorderSize(1);
    char RMSText[50]; sprintf(RMSText,"RMS = %1.2e #pm %1.0e",h->GetRMS(),h->GetRMSError());
    text -> AddText(RMSText);
    text -> Draw();    
  }
  
  if( type == "recHitE" )
  {
    if( region == "EB" ) h -> GetXaxis() -> SetRangeUser(-0.5,1.0);
    if( region == "EE" ) h -> GetXaxis() -> SetRangeUser(-1.0,2.0);
    h -> GetXaxis() -> SetTitle("E_{recHit} [GeV]");
    
    TF1* f_gaus;
    if( region == "EB" ) f_gaus = new TF1("f_gaus","[0]*exp(-1.*(x-[1])*(x-[1])/2/[2]/[2])",-0.08,0.08);
    if( region == "EE" ) f_gaus = new TF1("f_gaus","[0]*exp(-1.*(x-[1])*(x-[1])/2/[2]/[2])",-0.20,0.20);
    f_gaus -> SetParameters(h->GetMaximum(),0.,0.05);
    f_gaus -> SetNpx(10000);
    if( dataset == "DA" ) f_gaus -> SetLineColor(kBlack);
    if( dataset == "MC" ) f_gaus -> SetLineColor(kRed);
    h -> Fit("f_gaus","QENRS+");
    f_gaus -> Draw("same");
    
    if( dataset == "DA" )
    {
      TPaveText* text = new TPaveText(0.50,0.70,0.80,0.75,"brNDC");
      text -> SetTextSize(0.04);
      text -> SetLineWidth(0);
      text -> SetLineColor(kWhite);
      text -> SetTextFont(42);
      text -> SetTextColor(kBlack);    
      text -> SetFillColor(kWhite);
      text -> SetFillStyle(1001);
      text -> SetBorderSize(1);
      char RMSText[50]; sprintf(RMSText,"#sigma = %1.2e #pm %1.0e",f_gaus->GetParameter(2),f_gaus->GetParError(2));
      text -> AddText(RMSText);
      text -> Draw();    
    }
    if( dataset == "MC" )
    {
      TPaveText* text = new TPaveText(0.50,0.65,0.80,0.70,"brNDC");
      text -> SetTextSize(0.04);
      text -> SetLineWidth(0);
      text -> SetLineColor(kWhite);
      text -> SetTextFont(42);
      text -> SetTextColor(kRed);    
      text -> SetFillColor(kWhite);
      text -> SetFillStyle(1001);
      text -> SetBorderSize(1);
      char RMSText[50]; sprintf(RMSText,"#sigma = %1.2e #pm %1.0e",f_gaus->GetParameter(2),f_gaus->GetParError(2));
      text -> AddText(RMSText);
      text -> Draw();    
    }
  }
  if( type == "recHitPedSubADC" )
  {
    if( region == "EB" ) h -> GetXaxis() -> SetRangeUser(-25.,50.);
    if( region == "EE" ) h -> GetXaxis() -> SetRangeUser(-50.,100.);
    h -> GetXaxis() -> SetTitle("sample-pedestal [ADC]");
  }
  
  h -> GetYaxis() -> SetTitle("event fraction");
}



void SetGraphStyle(TGraphErrors* g, const std::string& dataset, const std::string& region, const std::string& type, const std::string& type2, TLegend** leg)
{
  if( (*leg) == NULL )
  {
    (*leg) = new TLegend(0.25,0.75,0.45,0.95,NULL,"brNDC");
    (*leg)->SetTextSize(0.04);
    (*leg)->SetLineWidth(0);
    (*leg)->SetLineColor(kWhite);
    (*leg)->SetTextFont(42);
    (*leg)->SetFillColor(kWhite);
    (*leg)->SetFillStyle(1001);
    (*leg)->SetBorderSize(1);
  }
  
  if( dataset == "DA" )
  {
    (*leg) -> AddEntry(g,"data","P");
    
    g -> SetLineColor(kBlack);
    g -> SetMarkerColor(kBlack);
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.);
  }
  if( dataset == "MC" )
  {
    (*leg) -> AddEntry(g,"MC","P");
    
    g -> SetLineColor(kRed);
    g -> SetMarkerColor(kRed);
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.);
  }
  if( dataset == "MC_PU0" )
  {
    (*leg) -> AddEntry(g,"MC PU0","P");
    
    g -> SetLineColor(kBlack);
    g -> SetMarkerColor(kBlack);
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.);
  }
  if( dataset == "MC_PU70" )
  {
    (*leg) -> AddEntry(g,"MC PU70 + noise (0 ageing)","P");
    
    g -> SetLineColor(kRed);
    g -> SetMarkerColor(kRed);
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.);
  }
  if( dataset == "MC_PU140" )
  {
    (*leg) -> AddEntry(g,"MC PU140 + noise (0 ageing)","P");
    
    g -> SetLineColor(kBlue);
    g -> SetMarkerColor(kBlue);
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.);
  }
  if( dataset == "MC_PU70ns" )
  {
    (*leg) -> AddEntry(g,"MC PU70 only","P");
    
    g -> SetLineColor(kRed);
    g -> SetMarkerColor(kRed);
    g -> SetMarkerStyle(25);
    g -> SetMarkerSize(1.);
  }
  if( dataset == "MC_PU140ns" )
  {
    (*leg) -> AddEntry(g,"MC PU140 only","P");
    
    g -> SetLineColor(kBlue);
    g -> SetMarkerColor(kBlue);
    g -> SetMarkerStyle(25);
    g -> SetMarkerSize(1.);
  }

  if( dataset == "MC_PU70nsc" )
  {
    (*leg) -> AddEntry(g,"MC PU70 only","P");
    
    g -> SetLineColor(kRed);
    g -> SetMarkerColor(kRed);
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.);
  }
  if( dataset == "MC_PU140nsc" )
  {
    (*leg) -> AddEntry(g,"MC PU140 only","P");
    
    g -> SetLineColor(kBlue);
    g -> SetMarkerColor(kBlue);
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.);
  }

  if( dataset == "MC_PU70_ax" )
  {
    (*leg) -> AddEntry(g,"MC PU70 (alexey)","L");
    
    g -> SetLineColor(kRed);
    g -> SetLineWidth(2);
    g -> SetMarkerColor(kRed);
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.);
  }
  if( dataset == "MC_PU140_ax" )
  {
    (*leg) -> AddEntry(g,"MC PU140 (alexey)","L");
    
    g -> SetLineColor(kBlue);
    g -> SetLineWidth(2);
    g -> SetMarkerColor(kBlue);
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.);
  }

  
  if( type == "vsNPU" )
  {
    if( dataset == "DA" )
    {
      g -> GetXaxis() -> SetRangeUser(0.,200.);
      if(type == "DAMC")       g -> GetXaxis() -> SetRangeUser(0.,50.);
      g -> GetXaxis() -> SetTitle("#LT N_{PU} #GT");
      if( region == "EB" )
      {
        g -> SetMinimum(0.000);
        g -> SetMaximum(0.40);
      }
      if( region == "EE" )
      {
        g -> SetMinimum(0.000);
        g -> SetMaximum(1.500);
      }
      
      TF1* f_prefit_DA = new TF1("f_prefit_DA","[0]+[1]*x",-0.5,199.5);
      if(type == "DAMC") f_prefit_DA->SetRange(5., 30.);
      g -> Fit("f_prefit_DA","QNRS");
      if(type == "DAMC") f_prefit_DA->SetRange(5., 100.);
      g -> Fit("f_prefit_DA","QNRS");
            
      for(int point = 0; point < g->GetN(); ++point)
      {
        double ey = g -> GetErrorY(point);
        //g -> SetPointError(point,g->GetErrorX(point),ey*sqrt(f_prefit_DA->GetChisquare()/f_prefit_DA->GetNDF()));
      }
            
      TF1* f_pol1_DA = new TF1("f_pol1_DA","[0]+[1]*x",-0.5,199.5);
      if(type == "DAMC") f_pol1_DA->SetRange(5., 30.);
      g -> Fit("f_pol1_DA","QNRS");
      if(type == "DAMC") f_pol1_DA->SetRange(5., 100.);
      g -> Fit("f_pol1_DA","QNRS");
      f_pol1_DA -> SetLineColor(kBlack);
      f_pol1_DA -> Draw("same");
      
      TPaveText* text = new TPaveText(0.20,0.85,0.50,0.90,"brNDC");
      text -> SetTextSize(0.04);
      text -> SetLineWidth(0);
      text -> SetLineColor(kWhite);
      text -> SetTextFont(42);
      text -> SetTextColor(kBlack);    
      text -> SetFillColor(kWhite);
      text -> SetFillStyle(1001);
      text -> SetBorderSize(1);
      char par0Text[50]; sprintf(par0Text,"y_{0} = %1.2e #pm %1.0e",f_pol1_DA->GetParameter(0),f_pol1_DA->GetParError(0));
      text -> AddText(par0Text);
      text -> Draw();
      
      text = new TPaveText(0.50,0.85,0.80,0.90,"brNDC");
      text -> SetTextSize(0.04);
      text -> SetLineWidth(0);
      text -> SetLineColor(kWhite);
      text -> SetTextFont(42);
      text -> SetTextColor(kBlack);    
      text -> SetFillColor(kWhite);
      text -> SetFillStyle(1001);
      text -> SetBorderSize(1);
      char par1Text[50]; sprintf(par1Text,"m = %1.2e #pm %1.0e",f_pol1_DA->GetParameter(1),f_pol1_DA->GetParError(1));
      text -> AddText(par1Text);
      text -> Draw();
    }
    
    if( dataset == "MC" )
    {
      TF1* f_prefit_MC = new TF1("f_prefit_MC","[0]+[1]*x",-0.5,199.5);
      if(type == "DAMC") f_prefit_MC->SetRange(5., 30.);
      g -> Fit("f_prefit_MC","QNRS");
      if(type == "DAMC") f_prefit_MC->SetRange(5., 100.);
      g -> Fit("f_prefit_MC","QNRS");

            
      for(int point = 0; point < g->GetN(); ++point)
      {
        double ey = g -> GetErrorY(point);
        //g -> SetPointError(point,g->GetErrorX(point),ey*sqrt(f_prefit_MC->GetChisquare()/f_prefit_MC->GetNDF()));
      }
            
      TF1* f_pol1_MC = new TF1("f_pol1_MC","[0]+[1]*x",-0.5,199.5);
      if(type == "DAMC") f_pol1_MC->SetRange(5., 30.);
      g -> Fit("f_pol1_MC","QNRS");
      if(type == "DAMC") f_pol1_MC->SetRange(5., 100.);
      g -> Fit("f_pol1_MC","QNRS");
      f_pol1_MC -> SetLineColor(kRed);
      f_pol1_MC -> Draw("same");
      
      TPaveText* text = new TPaveText(0.20,0.80,0.50,0.85,"brNDC");
      text -> SetTextSize(0.04);
      text -> SetLineWidth(0);
      text -> SetLineColor(kWhite);
      text -> SetTextFont(42);
      text -> SetTextColor(kRed);    
      text -> SetFillColor(kWhite);
      text -> SetFillStyle(1001);
      text -> SetBorderSize(1);
      char par0Text[50]; sprintf(par0Text,"y_{0} = %1.2e #pm %1.0e",f_pol1_MC->GetParameter(0),f_pol1_MC->GetParError(0));
      text -> AddText(par0Text);
      text -> Draw();
      
      text = new TPaveText(0.50,0.80,0.80,0.85,"brNDC");
      text -> SetTextSize(0.04);
      text -> SetLineWidth(0);
      text -> SetLineColor(kWhite);
      text -> SetTextFont(42);
      text -> SetTextColor(kRed);    
      text -> SetFillColor(kWhite);
      text -> SetFillStyle(1001);
      text -> SetBorderSize(1);
      char par1Text[50]; sprintf(par1Text,"m = %1.2e #pm %1.0e",f_pol1_MC->GetParameter(1),f_pol1_MC->GetParError(1));
      text -> AddText(par1Text);
      text -> Draw();        
    }
  }
  if( type == "vsIRing" )
  {
    g -> GetXaxis() -> SetRangeUser(-150.,150.);
    g -> GetXaxis() -> SetTitle("i_{ring}");
    
    if( region == "EBEE" )
    {
      g -> SetMinimum(0.01);
      g -> SetMaximum(2.);
      
      /*
      TArrow* line1 = new TArrow(-85.5,0.000,-85.5,0.500);
      line1 -> Draw("same");
      TArrow* line2 = new TArrow(+85.5,0.000,+85.5,0.500);
      line2 -> Draw("same");
      */
    }
  }
  if( type == "vsEta" )
  {
    g -> GetXaxis() -> SetRangeUser(-150.,150.);
    g -> GetXaxis() -> SetTitle("#eta");
    
    if( region == "EBEE" )
    {
      g -> SetMinimum(0.01);
      //g -> SetMaximum(0.6);
      g -> SetMaximum(2.);
      /*      
      TArrow* line1 = new TArrow(-1.4442,0.000,-1.4442,0.500);
      line1 -> Draw("same");
      TArrow* line2 = new TArrow(+1.4442,0.000,+1.4442,0.500);
      line2 -> Draw("same");
      */
    }
  }
  if( type == "vsAbsEta" )
  {
    g -> GetXaxis() -> SetRangeUser(-150.,150.);
    g -> GetXaxis() -> SetTitle("|#eta|");
    
    if( region == "EBEE" )
    {
      g -> SetMinimum(0.01);
      //      g -> SetMaximum(0.6);
      g -> SetMaximum(2.);
      
      /*
      TArrow* line1 = new TArrow(-1.4442,0.000,-1.4442,0.500);
      line1 -> Draw("same");
      TArrow* line2 = new TArrow(+1.4442,0.000,+1.4442,0.500);
      line2 -> Draw("same");
      */
    }
  }

  
  if( type2 == "recHitE_RMS" )
  {
    g -> GetYaxis() -> SetTitle("RMS(E_{recHit}) [GeV]");
  }
  if( type2 == "recHitE_RMSAt0PU" )
  {
    g -> GetYaxis() -> SetTitle("RMS(E_{recHit}) @ #LT N_{PU} #GT = 0 [GeV]");
  }
  if( type2 == "recHitE_sigma" )
  {
    g -> GetYaxis() -> SetTitle("effective #sigma(E_{recHit}) [GeV]");
  }
  if( type2 == "recHitPedSubADC_RMS" )
  {
    g -> GetYaxis() -> SetTitle("RMS(sample-pedestal) [ADC]");
  }
  if( type2 == "recHitPedSubADC_sigma" )
  {
    g -> GetYaxis() -> SetTitle("#sigma(sample-pedestal) [ADC]");
  }
}




float GetEffectiveSigma(TH1F* h){

  double TotEvents = h->Integral(1, h->GetNbinsX()-1);
  float LocEvents = 0.;
  int binI = h->FindBin(h->GetMean());
  int binF = h->GetNbinsX()-1;
  bool keepGoing = false;

  for(int jBin=binI; jBin>0; --jBin){
    LocEvents = 0.;
    keepGoing = false;
    for(int iBin=jBin; iBin<binF; ++iBin){
      LocEvents += h->GetBinContent(iBin);
      if(LocEvents/TotEvents >= 0.68) {
	if(iBin-jBin < binF-binI) {
	  binF = iBin;
	  binI = jBin;
	  keepGoing = true;
	}
        break;
      }
      if(iBin == binF-1 && binF == h->GetNbinsX()-1){
	keepGoing = true;
	--binI;
      }
      if(iBin == binF-1 && binF != h->GetNbinsX()-1) break;
    }
    if(keepGoing == false) break;
  }

  float sigma = (h->GetBinCenter(binF) - h->GetBinCenter(binI))/2.;
  //  std::cout << " >>> sigma" << sigma << std::endl;
  return sigma;
}


void GetEffectiveSigma(TH1F* h, double& E1, double& E2){

  double TotEvents = h->Integral(1, h->GetNbinsX()-1);
  float LocEvents = 0.;
  int binI = h->FindBin(h->GetMean());
  int binF = h->GetNbinsX()-1;
  bool keepGoing = false;

  for(int jBin=binI; jBin>0; --jBin){
    LocEvents = 0.;
    keepGoing = false;
    for(int iBin=jBin; iBin<binF; ++iBin){
      LocEvents += h->GetBinContent(iBin);
      if(LocEvents/TotEvents >= 0.68) {
	if(iBin-jBin < binF-binI) {
	  binF = iBin;
	  binI = jBin;
	  keepGoing = true;
	}
        break;
      }
      if(iBin == binF-1 && binF == h->GetNbinsX()-1){
	keepGoing = true;
	--binI;
      }
      if(iBin == binF-1 && binF != h->GetNbinsX()-1) break;
    }
    if(keepGoing == false) break;
  }

  E1 = h->GetBinCenter(binI);
  E2 = h->GetBinCenter(binF);
  float sigma = (h->GetBinCenter(binF) - h->GetBinCenter(binI))/2.;
  //  std::cout << " >>> sigma" << sigma << std::endl;
  return;
}



void FillHistoRMS(TH1F* histo, float rms, float mean){
  std::cout << " (mean - rms)/2. = " << (mean - rms)/2. << std::endl;
  std::cout << " (mean + rms)/2. = " << (mean + rms)/2. << std::endl;
  for(int iB=1; iB<histo->GetNbinsX(); ++iB) {
    if(iB < histo->FindBin(mean - rms) ) histo->SetBinContent(iB, 0.);
    else if(iB > histo->FindBin(mean + rms) ) histo->SetBinContent(iB, 0.);
  }
  std::cout << " histo->FindBin((mean - rms)/2.) = " << histo->FindBin((mean - rms)/2.) << std::endl;
  histo->SetFillColor(kYellow);
  histo->SetFillStyle(3001);
}

