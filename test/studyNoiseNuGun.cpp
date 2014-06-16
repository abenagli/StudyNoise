#include "setTDRStyle.h"
#include "ConfigParser.h"
#include "ParserUtils.h"
#include "PUreweightingUtils.h"
#include "avgPUList.h"
#include "TEndcapRings.h"

#include <iostream>
#include <fstream>
#include <map>
#include <cmath>

#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"



bool FillChain(TChain& chain, const std::string& inputFileList);

void WriteNormalized(TH1F* histo);



int nBins_iRing = 249;
float iRingMin = -124.5;
float iRingMax = 124.5;

int nBins_recHitADC = 300;
float recHitADCMin = 100.5;
float recHitADCMax = 400.5;

int nBins_recHitPedSubADC = 96000;
float recHitPedSubADCMin = -50.;
float recHitPedSubADCMax = 100.;

int nBins_recHitE = 192000;
float recHitEMin = -4.;
float recHitEMax = +8.;

int pedMaxSample = 10;






int main(int argc, char** argv)
{
  // Input parameters
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>> studyNoiseNuGun::usage: " << argv[0] << " configFileName" << std::endl;
    return 1;
  }
  
  // Parse the config file
  parseConfigFile(argv[1]) ;
  
  std::string inputListMC = gConfigParser -> readStringOption("Input::inputListMC");
  std::string inputTreeMC = gConfigParser -> readStringOption("Input::inputTreeMC");
  
  std::string EERingsFile   = gConfigParser -> readStringOption("Input::EERingsFile");
  
  std::string label = gConfigParser -> readStringOption("Options::label");
  
  std::string outputDir = gConfigParser -> readStringOption("Output::outputDir");
  
  
  
  
  
  
  // open files and fill the tree chain
  TChain* chain = new TChain(inputTreeMC.c_str());
  FillChain(*chain,inputListMC);
  
  // define outfile
  system(("mkdir "+outputDir).c_str());
  std::string outfileName = outputDir + "/studyNoiseNuGun_" + label + ".root";
  TFile* outFile = new TFile(outfileName.c_str(),"RECREATE");
  
  
  
  
  
  
  // EE geometry
  TEndcapRings* EERings = new TEndcapRings(EERingsFile);
  
  
  
  
  
  
  // define histograms
  TH1F* h_occupancy_vsIRing;
  
  std::map<std::string,TH1F*> h_recHitADC_vsIRing;
  std::map<std::string,TH1F*> h_recHitPedSubADC_vsIRing;
  std::map<std::string,TH1F*> h_recHitE_vsIRing;
  
  std::string histoName = "h_occupancy_vsIRing";
  h_occupancy_vsIRing = new TH1F(histoName.c_str(),"",nBins_iRing,iRingMin,iRingMax);
  h_occupancy_vsIRing -> Sumw2();
  
  
  for(int bin = 1; bin <= nBins_iRing; ++bin)
  {
    float binLowEdge = h_occupancy_vsIRing -> GetBinLowEdge(bin);
    float binHigEdge = h_occupancy_vsIRing -> GetBinLowEdge(bin) + h_occupancy_vsIRing->GetBinWidth(bin);
    
    char iRingChar[50];
    sprintf(iRingChar,"iRing%2.1f-%2.1f",binLowEdge,binHigEdge);
    std::string iRingString(iRingChar);
    
    histoName = "h_recHitADC_" + iRingString;
    h_recHitADC_vsIRing[iRingString] = new TH1F(histoName.c_str(),"",nBins_recHitADC,recHitADCMin,recHitADCMax);
    h_recHitADC_vsIRing[iRingString] -> Sumw2();
    
    histoName = "h_recHitPedSubADC_" + iRingString;
    h_recHitPedSubADC_vsIRing[iRingString] = new TH1F(histoName.c_str(),"",nBins_recHitPedSubADC,recHitPedSubADCMin,recHitPedSubADCMax);
    h_recHitPedSubADC_vsIRing[iRingString] -> Sumw2();
    
    histoName = "h_recHitE_" + iRingString;
    h_recHitE_vsIRing[iRingString] = new TH1F(histoName.c_str(),"",nBins_recHitE,recHitEMin,recHitEMax);
    h_recHitE_vsIRing[iRingString] -> Sumw2();
  }
  
  
  
  
  
  
  // define variables;
  std::vector<float>* EBRecHit_ADC = new std::vector<float>;
  std::vector<float>* EERecHit_ADC = new std::vector<float>;  
  std::vector<float>* EBRecHit_E = new std::vector<float>;
  std::vector<float>* EERecHit_E = new std::vector<float>;
  std::vector<float>* EBRecHit_flag = new std::vector<float>;
  std::vector<float>* EERecHit_flag = new std::vector<float>;
  std::vector<int>* EBRecHit_ieta = new std::vector<int>;
  std::vector<int>* EERecHit_ieta = new std::vector<int>;
  std::vector<int>* EBRecHit_iphi = new std::vector<int>;
  std::vector<int>* EERecHit_iphi = new std::vector<int>;
  std::vector<int>* EBRecHit_iz = new std::vector<int>;
  std::vector<int>* EERecHit_iz = new std::vector<int>;
  
  chain -> SetBranchStatus("EBRecHit_samples",1);  chain -> SetBranchAddress("EBRecHit_samples",  &EBRecHit_ADC);
  chain -> SetBranchStatus("EERecHit_samples",1);  chain -> SetBranchAddress("EERecHit_samples",  &EERecHit_ADC);
  chain -> SetBranchStatus("EBRecHit_E",1);        chain -> SetBranchAddress("EBRecHit_E",        &EBRecHit_E);
  chain -> SetBranchStatus("EERecHit_E",1);        chain -> SetBranchAddress("EERecHit_E",        &EERecHit_E);
  chain -> SetBranchStatus("EBRecHit_flag",1);     chain -> SetBranchAddress("EBRecHit_flag",     &EBRecHit_flag);
  chain -> SetBranchStatus("EERecHit_flag",1);     chain -> SetBranchAddress("EERecHit_flag",     &EERecHit_flag);
  chain -> SetBranchStatus("EBRecHit_ietaORix",1); chain -> SetBranchAddress("EBRecHit_ietaORix", &EBRecHit_ieta);
  chain -> SetBranchStatus("EERecHit_ietaORix",1); chain -> SetBranchAddress("EERecHit_ietaORix", &EERecHit_ieta);
  chain -> SetBranchStatus("EBRecHit_iphiORiy",1); chain -> SetBranchAddress("EBRecHit_iphiORiy", &EBRecHit_iphi);
  chain -> SetBranchStatus("EERecHit_iphiORiy",1); chain -> SetBranchAddress("EERecHit_iphiORiy", &EERecHit_iphi);
  chain -> SetBranchStatus("EBRecHit_zside",1);    chain -> SetBranchAddress("EBRecHit_zside",    &EBRecHit_iz);
  chain -> SetBranchStatus("EERecHit_zside",1);    chain -> SetBranchAddress("EERecHit_zside",    &EERecHit_iz);
  
  
  
  // loop over events
  std::cout << ">>> start of loop: " << chain->GetEntries() << " entries" << std::endl;
  for(int entry = 0; entry < chain->GetEntries(); ++entry)
  {
    //if( entry >= 10 ) break;
    if( entry%1 == 0 ) std::cout << ">>>>>> reading entry " << entry << " / " << chain->GetEntries() << "\r" << std::flush;
    chain -> GetEntry(entry);
    
    
    // fill EBRecHits
    for(unsigned int rhIt = 0; rhIt < EBRecHit_E->size(); ++rhIt)
    {
      float recHit_E    = EBRecHit_E   ->at(rhIt);
      float recHit_ieta = EBRecHit_ieta->at(rhIt);
      float recHit_iRing = recHit_ieta;
      
      if( EBRecHit_flag->at(rhIt) != 0 ) continue;
      
      
      // fill occupancy histograms
      int iRingBin = h_occupancy_vsIRing -> Fill(recHit_iRing);
      if( iRingBin < 1 || iRingBin > nBins_iRing ) continue;
      
      float binLowEdge = h_occupancy_vsIRing -> GetBinLowEdge(iRingBin);
      float binHigEdge = h_occupancy_vsIRing -> GetBinLowEdge(iRingBin) + h_occupancy_vsIRing->GetBinWidth(iRingBin);
      
      char iRingChar[50];
      sprintf(iRingChar,"iRing%2.1f-%2.1f",binLowEdge,binHigEdge);
      std::string iRingString(iRingChar);
      
      
      // fill noise histograms
      h_recHitE_vsIRing[iRingString] -> Fill(recHit_E);
      
      float pedestal = 0.;
      for(int sampleIt = 0; sampleIt < 10; ++sampleIt)
      {
        float recHit_ADC = EBRecHit_ADC->at(rhIt*10+sampleIt);
        h_recHitADC_vsIRing[iRingString] -> Fill(recHit_ADC);
        if( sampleIt < pedMaxSample ) pedestal += recHit_ADC;
      }
      pedestal /= pedMaxSample;
      
      for(int sampleIt = 0; sampleIt < 10; ++sampleIt)
      {
        float recHit_ADC = EBRecHit_ADC->at(rhIt*10+sampleIt);
        h_recHitPedSubADC_vsIRing[iRingString] -> Fill(recHit_ADC-pedestal);
      }
    }
    // end of fill EBRecHits
    
    
    // fill EERecHits
    for(unsigned int rhIt = 0; rhIt < EERecHit_E->size(); ++rhIt)
    {
      float recHit_E     = EERecHit_E   ->at(rhIt);
      float recHit_ieta  = EERecHit_ieta->at(rhIt);
      float recHit_iphi  = EERecHit_iphi->at(rhIt);
      float recHit_iz    = EERecHit_iz->at(rhIt);
      float recHit_iRing = recHit_iz * ( 86 + EERings->GetEndcapRing(recHit_ieta,recHit_iphi,recHit_iz) );
      
      if( EERecHit_flag->at(rhIt) != 0 ) continue;
      
      
      // fill occupancy histograms
      int iRingBin = h_occupancy_vsIRing -> Fill(recHit_iRing);
      if( iRingBin < 1 || iRingBin > nBins_iRing ) continue;
      
      float binLowEdge = h_occupancy_vsIRing -> GetBinLowEdge(iRingBin);
      float binHigEdge = h_occupancy_vsIRing -> GetBinLowEdge(iRingBin) + h_occupancy_vsIRing->GetBinWidth(iRingBin);
      
      char iRingChar[50];
      sprintf(iRingChar,"iRing%2.1f-%2.1f",binLowEdge,binHigEdge);
      std::string iRingString(iRingChar);
      
      
      // fill noise histograms
      h_recHitE_vsIRing[iRingString] -> Fill(recHit_E);
      
      float pedestal = 0.;
      for(int sampleIt = 0; sampleIt < 10; ++sampleIt)
      {
        float recHit_ADC = EERecHit_ADC->at(rhIt*10+sampleIt);
        h_recHitADC_vsIRing[iRingString] -> Fill(recHit_ADC);
        if( sampleIt < pedMaxSample ) pedestal += recHit_ADC;
      }
      pedestal /= pedMaxSample;
      
      for(int sampleIt = 0; sampleIt < 10; ++sampleIt)
      {
        float recHit_ADC = EERecHit_ADC->at(rhIt*10+sampleIt);
        h_recHitPedSubADC_vsIRing[iRingString] -> Fill(recHit_ADC-pedestal);
      }
    }
    // end of fill EERecHits
  }
  
  std::cout << std::endl;
  std::cout << ">>> end of loop" << std::endl;
  
  
  
  outFile -> cd();
  
  h_occupancy_vsIRing -> Write();
  
  TDirectory* c1 = outFile -> mkdir("plots_vs_iRing");
  c1 -> cd();
  for(int iRingBin = 1; iRingBin <= nBins_iRing; ++iRingBin)
  {
    float binLowEdge = h_occupancy_vsIRing -> GetBinLowEdge(iRingBin);
    float binHigEdge = h_occupancy_vsIRing -> GetBinLowEdge(iRingBin) + h_occupancy_vsIRing->GetBinWidth(iRingBin);
    
    char iRingChar[50];
    sprintf(iRingChar,"iRing%2.1f-%2.1f",binLowEdge,binHigEdge);
    std::string iRingString(iRingChar);
    WriteNormalized( h_recHitADC_vsIRing[iRingString] );
    WriteNormalized( h_recHitPedSubADC_vsIRing[iRingString] );
    WriteNormalized( h_recHitE_vsIRing[iRingString] );
  }
      
  outFile -> cd("");
  
  outFile -> Close();
  std::cout << ">>> outfile " << outfileName << " created" << std::endl;
  
  
  return 0;
}






bool FillChain(TChain& chain, const std::string& inputFileList)
{
  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;

  if(!inFile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return false;
  }

  while(1)
  {
    inFile >> buffer;
    if(!inFile.good()) break;
    if( buffer.at(0) == '#' ) continue;
    chain.Add(buffer.c_str());
    std::cout << ">>> ntupleUtils::FillChain - treeName = " << chain.GetName() << " from file " << buffer << std::endl;
  }

  return true;
}






void ReadAxis(const std::string& fileName, std::vector<float>& axis)
{
  std::ifstream inFile(fileName.c_str(),std::ios::in);
  float buffer;
  while(1)
  {
    inFile >> buffer;
    if(!inFile.good()) break;
    
    axis.push_back(buffer);
  }
}



void MySetAxisBins(TH1F* h, const std::vector<float>& axis)
{
  double* new_bins = new double[axis.size()];
  
  for(unsigned int it = 0; it < axis.size(); ++it)
    new_bins[it] = axis.at(it);
  
  h -> SetBins(axis.size()-1,new_bins);
}



void WriteNormalized(TH1F* histo)
{
  if( histo->Integral() != 0.)
    histo -> Scale(1./histo->Integral());
  histo -> Write();
}



void FillHisto(TH1F* histo,
              int seedIeta, int seedIphi, 
              std::vector<float> recHitMatrix_E, std::vector<int> recHitMatrix_ieta, std::vector<int> recHitMatrix_iphi)
{
  //std::cout << ">>>>>> seed: (" << seedIeta << "," << seedIphi << ")" << std::endl;
  
  for(unsigned int it = 0; it < recHitMatrix_E.size(); ++it)
  {
    //std::cout << ">>>>>>>>> recHit: (" << recHitMatrix_ieta.at(it) << "," << recHitMatrix_iphi.at(it) << ")" << std::endl;
    if( abs(seedIeta - recHitMatrix_ieta.at(it)) > 5 ) histo -> Fill( recHitMatrix_E.at(it) );
  }
}


bool CheckSeedDistance(const int& recHit_ieta, const int& recHit_iphi,
                       const int& seed_ieta, const int& seed_iphi, const int& seed_iz)
{
  if( recHit_ieta == -40 && recHit_iphi == 121 && seed_iz == 0  ) return false;
  if( recHit_ieta == -14 && recHit_iphi == 333 && seed_iz == 0  ) return false;
  if( recHit_ieta ==  86 && recHit_iphi ==  47 && seed_iz == -1 ) return false;
  if( recHit_ieta ==  87 && recHit_iphi ==  28 && seed_iz == -1 ) return false;
  if( recHit_ieta ==  63 && recHit_iphi ==  49 && seed_iz == -1 ) return false;
  if( recHit_ieta ==  89 && recHit_iphi ==  29 && seed_iz == +1 ) return false;
  if( recHit_ieta ==  36 && recHit_iphi ==  36 && seed_iz == +1 ) return false;
  if( recHit_ieta ==  55 && recHit_iphi ==  65 && seed_iz == +1 ) return false;
  if( recHit_ieta ==  39 && recHit_iphi ==  53 && seed_iz == +1 ) return false;
  
  bool farRecHit = false;
  
  if( seed_iz == 0 )
  {
    if( fabs(seed_ieta - recHit_ieta) > 5 ) farRecHit = true;
  }
  else
  {
    if( fabs(sqrt(recHit_ieta*recHit_ieta+recHit_iphi*recHit_iphi) - sqrt(seed_ieta*seed_ieta+seed_iphi*seed_iphi)) > 5 ) farRecHit = true;
  }
  
  return farRecHit;
}
