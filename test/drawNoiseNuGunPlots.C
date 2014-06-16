#include "setTDRStyle.cc"
#include "setTDRStyle.h"

std::string label = "MC-NuGun-noPU";

int rebin = 50;



void drawNoiseNuGunPlots()
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
    
  TFile* f = TFile::Open(("studyNoiseNuGun_"+label+".root").c_str(),"READ");
  f -> cd();
  
  
  //-------------------------
  // get occupancy histograms
  
  TH1F* h_occupancy_vsIRing = (TH1F*)( f->Get("h_occupancy_vsIRing") );
  
  
  // ------------
  // build graphs
    
  TGraphErrors* g_recHitADC_RMS_vsIRing       = new TGraphErrors();
  TGraphErrors* g_recHitPedSubADC_RMS_vsIRing = new TGraphErrors();  
  TGraphErrors* g_recHitE_RMS_vsIRing         = new TGraphErrors();  
  
  int point = 0;
  for(int bin = 1; bin <= h_occupancy_vsIRing->GetNbinsX(); ++bin)
  {  
    float binCenter  = h_occupancy_vsIRing -> GetBinCenter(bin);
    float binLowEdge = h_occupancy_vsIRing -> GetBinLowEdge(bin);
    float binHigEdge = h_occupancy_vsIRing -> GetBinLowEdge(bin) + h_occupancy_vsIRing->GetBinWidth(bin);
    
    char histoName[50];
    
    sprintf(histoName,"plots_vs_iRing/h_recHitADC_iRing%2.1f-%2.1f",binLowEdge,binHigEdge);
    TH1F* histo = (TH1F*)( f->Get(histoName) );
    if( histo->GetEntries() < 30 ) continue;
    g_recHitADC_RMS_vsIRing -> SetPoint(point,binCenter,histo->GetRMS());
    g_recHitADC_RMS_vsIRing -> SetPointError(point,0.,histo->GetRMSError());
    
    sprintf(histoName,"plots_vs_iRing/h_recHitPedSubADC_iRing%2.1f-%2.1f",binLowEdge,binHigEdge);
    TH1F* histo = (TH1F*)( f->Get(histoName) );
    if( histo->GetEntries() < 30 ) continue;
    g_recHitPedSubADC_RMS_vsIRing -> SetPoint(point,binCenter,histo->GetRMS());
    g_recHitPedSubADC_RMS_vsIRing -> SetPointError(point,0.,histo->GetRMSError());
        
    sprintf(histoName,"plots_vs_iRing/h_recHitE_iRing%2.1f-%2.1f",binLowEdge,binHigEdge);
    TH1F* histo = (TH1F*)( f->Get(histoName) );
    if( histo->GetEntries() < 30 ) continue;
    g_recHitE_RMS_vsIRing -> SetPoint(point,binCenter,histo->GetRMS());
    g_recHitE_RMS_vsIRing -> SetPointError(point,0.,histo->GetRMSError());    
    
    ++point;
  }
  
  
  TFile* outFile = TFile::Open(("nuGunGraphs_"+label+".root").c_str(),"RECREATE");
  
  g_recHitADC_RMS_vsIRing       -> Write("g_recHitADC_RMS_vsIRing");
  g_recHitPedSubADC_RMS_vsIRing -> Write("g_recHitPedSubADC_RMS_vsIRing");
  g_recHitE_RMS_vsIRing         -> Write("g_recHitE_RMS_vsIRing");
}