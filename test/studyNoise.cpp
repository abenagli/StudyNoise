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

void ReadAxis(const std::string& fileName, std::vector<float>& axis);
void MySetAxisBins(TH1F* h, const std::vector<float>& axis);

void WriteNormalized(TH1F* histo);

void FillHisto(TH1F* histo,
              int seedIeta, int seedIphi, 
              std::vector<float> recHitMatrix_E, std::vector<int> recHitMatrix_ieta, std::vector<int> recHitMatrix_iphi);

bool CheckSeedDistance(const int& recHit_ieta, const int& recHit_iphi,
                       const int& seed_ieta, const int& seed_iphi, const int& seed_iz);



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
    std::cerr << ">>> studyNoise::usage: " << argv[0] << " configFileName" << std::endl;
    return 1;
  }
  
  // Parse the config file
  parseConfigFile(argv[1]) ;
  
  std::string inputListDA = gConfigParser -> readStringOption("Input::inputListDA");
  std::string inputListMC = gConfigParser -> readStringOption("Input::inputListMC");
  std::string inputTreeDA = gConfigParser -> readStringOption("Input::inputTreeDA");
  std::string inputTreeMC = gConfigParser -> readStringOption("Input::inputTreeMC");
  
  std::string PUDistrFile   = gConfigParser -> readStringOption("Input::PUDistrFile");
  std::string avgPUListFile = gConfigParser -> readStringOption("Input::avgPUListFile");
  std::string EERingsFile   = gConfigParser -> readStringOption("Input::EERingsFile");
  
  std::string nPUAxisFile    = gConfigParser -> readStringOption("Input::nPUAxisFile");
  std::string iRingAxisFile  = gConfigParser -> readStringOption("Input::iRingAxisFile");
  std::string etaAxisFile    = gConfigParser -> readStringOption("Input::etaAxisFile");
  std::string absEtaAxisFile = gConfigParser -> readStringOption("Input::absEtaAxisFile");
  
  std::string label = gConfigParser -> readStringOption("Options::label");
  std::string type = gConfigParser -> readStringOption("Options::type");
  int PUReweighting = gConfigParser -> readIntOption("Options::PUReweighting");
  
  std::string outputDir = gConfigParser -> readStringOption("Output::outputDir");
  
  
  
  
  
  
  // define datasets
  std::vector<std::string> datasets;
  datasets.push_back("DA");
  datasets.push_back("MC");
  int nDatasets = int(datasets.size());
  
  // define regions
  std::vector<std::string> regions;
  regions.push_back("EB");
  regions.push_back("EE");
  int nRegions = int(regions.size());
  
  // open files and fill the tree chain
  std::map<std::string,std::string> inputFileList;
  inputFileList["DA"] = inputListDA;
  inputFileList["MC"] = inputListMC;
  
  std::map<std::string,TChain*> chain;
  chain["DA"] = new TChain(inputTreeDA.c_str());
  chain["MC"] = new TChain(inputTreeMC.c_str());
  
  for(int datasetIt = 0; datasetIt < nDatasets; ++datasetIt)
    FillChain(*chain[datasets.at(datasetIt)],inputFileList[datasets.at(datasetIt)]);
  
  // define outfile
  system(("mkdir "+outputDir).c_str());
  std::string outfileName = outputDir + "/studyNoise_" + label + ".root";
  TFile* outFile = new TFile(outfileName.c_str(),"RECREATE");
  
  
  
  
  
  
  // PU reweighting
  std::map<float,float> PUWeights = *(ComputePUweights(chain["MC"],PUDistrFile,false));
  
  // define list of avg PU
  avgPUList* myAvgPUList = new avgPUList(avgPUListFile);
  
  // EE geometry
  TEndcapRings* EERings = new TEndcapRings(EERingsFile);
  
  // read axes
  std::vector<float> nPU_axis; ReadAxis(nPUAxisFile,nPU_axis);
  int nBins_nPU = int(nPU_axis.size()-1);
  
  std::vector<float> iRing_axis; ReadAxis(iRingAxisFile,iRing_axis);
  int nBins_iRing = int(iRing_axis.size()-1);
  
  std::vector<float> eta_axis; ReadAxis(etaAxisFile,eta_axis);
  int nBins_eta = int(eta_axis.size()-1);
  
  std::vector<float> absEta_axis; ReadAxis(absEtaAxisFile,absEta_axis);
  int nBins_absEta = int(absEta_axis.size()-1);
  
  
  
  
  
  // define histograms
  std::map<std::string,std::map<std::string,TH1F*> > h_nVtx;
  std::map<std::string,std::map<std::string,TH1F*> > h_occupancy_vsNPU;
  std::map<std::string,std::map<std::string,TH1F*> > h_occupancy_vsIRing;
  std::map<std::string,std::map<std::string,TH1F*> > h_occupancy_vsEta;
  std::map<std::string,std::map<std::string,TH1F*> > h_occupancy_vsAbsEta;
  
  std::map<std::string,std::map<std::string,TProfile2D*> > p_recHitADC;
  
  std::map<std::string,std::map<std::string,TH1F*> > h_recHitADC;
  std::map<std::string,std::map<std::string,TH1F*> > h_recHitPedSubADC;
  std::map<std::string,std::map<std::string,TH1F*> > h_recHitE;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_nPUDistr_vsNPU;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_nPUDistr_vsNPUIRing;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_nPUDistr_vsNPUEta;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_nPUDistr_vsNPUAbsEta;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitPedSubADC_vsNPU;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitPedSubADC_vsIRing;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitPedSubADC_vsEta;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitPedSubADC_vsAbsEta;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitPedSubADC_vsNPUIRing;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitPedSubADC_vsNPUEta;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitPedSubADC_vsNPUAbsEta;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitE_vsNPU;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitE_vsIRing;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitE_vsEta;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitE_vsAbsEta;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitE_vsNPUIRing;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitE_vsNPUEta;
  std::map<std::string,std::map<std::string,std::map<std::string,TH1F*> > > h_recHitE_vsNPUAbsEta;
  
  std::cout << ">>> define histograms" << std::endl;
  for(int datasetIt = 0; datasetIt < nDatasets; ++datasetIt)
  {
    for(int regionIt = 0; regionIt < nRegions; ++regionIt)
    {
      std::string dataset = datasets.at(datasetIt);
      std::string region = regions.at(regionIt);
      std::cout << ">>>>>> dataset: " << dataset << "   region: " << region << std::endl;
      
      std::string histoName = "h_nVtx_" + region + "_" + dataset;
      (h_nVtx[dataset])[region] = new TH1F(histoName.c_str(),"",200,-0.5,199.5);
      (h_nVtx[dataset])[region] -> Sumw2();
      
      histoName = "h_occupancy_vsNPU_" + region + "_" + dataset;
      (h_occupancy_vsNPU[dataset])[region] = new TH1F(histoName.c_str(),"",nBins_nPU,0.,200.);
      (h_occupancy_vsNPU[dataset])[region] -> Sumw2();
      MySetAxisBins((h_occupancy_vsNPU[dataset])[region],nPU_axis);
      
      histoName = "h_occupancy_vsIRing_" + region + "_" + dataset;
      (h_occupancy_vsIRing[dataset])[region] = new TH1F(histoName.c_str(),"",nBins_iRing,1.,126.);
      (h_occupancy_vsIRing[dataset])[region] -> Sumw2();
      MySetAxisBins((h_occupancy_vsIRing[dataset])[region],iRing_axis);
      
      histoName = "h_occupancy_vsEta_" + region + "_" + dataset;
      (h_occupancy_vsEta[dataset])[region] = new TH1F(histoName.c_str(),"",nBins_eta,-2.5,2.5);
      (h_occupancy_vsEta[dataset])[region] -> Sumw2();
      MySetAxisBins((h_occupancy_vsEta[dataset])[region],eta_axis);
      
      histoName = "h_occupancy_vsAbsEta_" + region + "_" + dataset;
      (h_occupancy_vsAbsEta[dataset])[region] = new TH1F(histoName.c_str(),"",nBins_absEta,0,2.5);
      (h_occupancy_vsAbsEta[dataset])[region] -> Sumw2();
      MySetAxisBins((h_occupancy_vsAbsEta[dataset])[region],absEta_axis);
      
      
      histoName = "p_recHitADC_" + region + "_" + dataset;
      if( region == "EB" ) (p_recHitADC[dataset])[region] = new TProfile2D(histoName.c_str(),"",171,-85.,86.,360,1.,361.);
      if( region == "EE" ) (p_recHitADC[dataset])[region] = new TProfile2D(histoName.c_str(),"",201,-100.,101.,100,1,101.);
      (p_recHitADC[dataset])[region] -> Sumw2();
      
      
      histoName = "h_recHitADC_" + region + "_" + dataset;
      (h_recHitADC[dataset])[region] = new TH1F(histoName.c_str(),"",nBins_recHitADC,recHitADCMin,recHitADCMax);
      (h_recHitADC[dataset])[region] -> Sumw2();
      
      histoName = "h_recHitPedSubADC_" + region + "_" + dataset;
      (h_recHitPedSubADC[dataset])[region] = new TH1F(histoName.c_str(),"",nBins_recHitPedSubADC,recHitPedSubADCMin,recHitPedSubADCMax);
      (h_recHitPedSubADC[dataset])[region] -> Sumw2();
      
      histoName = "h_recHitE_" + region + "_" + dataset;
      (h_recHitE[dataset])[region] = new TH1F(histoName.c_str(),"",nBins_recHitE,recHitEMin,recHitEMax);
      (h_recHitE[dataset])[region] -> Sumw2();
      
      std::cout << ">>>>>> global histos defined" << std::endl;
      
      
      for(int nPUIt = 0; nPUIt < nBins_nPU; ++nPUIt)
      {
        char nPUChar[50];
        sprintf(nPUChar,"nPU%2.1f-%2.1f",nPU_axis.at(nPUIt),nPU_axis.at(nPUIt+1));
        std::string nPUString(nPUChar);
        
        histoName = "h_nPUDistr_" + nPUString + "_" + region + "_" + dataset;
        ((h_nPUDistr_vsNPU[dataset])[region])[nPUString] = new TH1F(histoName.c_str(),"",500,0.,500.);
        ((h_nPUDistr_vsNPU[dataset])[region])[nPUString] -> Sumw2();
                
        histoName = "h_recHitPedSubADC_" + nPUString + "_" + region + "_" + dataset;
        ((h_recHitPedSubADC_vsNPU[dataset])[region])[nPUString] = new TH1F(histoName.c_str(),"",nBins_recHitPedSubADC,recHitPedSubADCMin,recHitPedSubADCMax);
        ((h_recHitPedSubADC_vsNPU[dataset])[region])[nPUString] -> Sumw2();
        
        histoName = "h_recHitE_" + nPUString + "_" + region + "_" + dataset;
        ((h_recHitE_vsNPU[dataset])[region])[nPUString] = new TH1F(histoName.c_str(),"",nBins_recHitE,recHitEMin,recHitEMax);
        ((h_recHitE_vsNPU[dataset])[region])[nPUString] -> Sumw2();
      }
      std::cout << ">>>>>> histos vs nPU defined" << std::endl;
      
      for(int iRingIt = 0; iRingIt < nBins_iRing; ++iRingIt)
      {
        char iRingChar[50];
        sprintf(iRingChar,"iRing%2.1f-%2.1f",iRing_axis.at(iRingIt),iRing_axis.at(iRingIt+1));
        std::string iRingString(iRingChar);
        
        histoName = "h_recHitPedSubADC_" + iRingString + "_" + region + "_" + dataset;
        ((h_recHitPedSubADC_vsIRing[dataset])[region])[iRingString] = new TH1F(histoName.c_str(),"",nBins_recHitPedSubADC,recHitPedSubADCMin,recHitPedSubADCMax);
        ((h_recHitPedSubADC_vsIRing[dataset])[region])[iRingString] -> Sumw2();
        
        histoName = "h_recHitE_" + iRingString + "_" + region + "_" + dataset;
        ((h_recHitE_vsIRing[dataset])[region])[iRingString] = new TH1F(histoName.c_str(),"",nBins_recHitE,recHitEMin,recHitEMax);
        ((h_recHitE_vsIRing[dataset])[region])[iRingString] -> Sumw2();
      }
      std::cout << ">>>>>> histos vs iRing defined" << std::endl;
      
      for(int etaIt = 0; etaIt < nBins_eta; ++etaIt)
      {
        char etaChar[50];
        sprintf(etaChar,"eta%1.4f-%1.4f",eta_axis.at(etaIt),eta_axis.at(etaIt+1));
        std::string etaString(etaChar);
        
        histoName = "h_recHitPedSubADC_" + etaString + "_" + region + "_" + dataset;
        ((h_recHitPedSubADC_vsEta[dataset])[region])[etaString] = new TH1F(histoName.c_str(),"",nBins_recHitPedSubADC,recHitPedSubADCMin,recHitPedSubADCMax);
        ((h_recHitPedSubADC_vsEta[dataset])[region])[etaString] -> Sumw2();
        
        histoName = "h_recHitE_" + etaString + "_" + region + "_" + dataset;
        ((h_recHitE_vsEta[dataset])[region])[etaString] = new TH1F(histoName.c_str(),"",nBins_recHitE,recHitEMin,recHitEMax);
        ((h_recHitE_vsEta[dataset])[region])[etaString] -> Sumw2();
      }
      std::cout << ">>>>>> histos vs eta defined" << std::endl;
      
      for(int absEtaIt = 0; absEtaIt < nBins_absEta; ++absEtaIt)
      {
        char absEtaChar[50];
        sprintf(absEtaChar,"absEta%1.4f-%1.4f",absEta_axis.at(absEtaIt),absEta_axis.at(absEtaIt+1));
        std::string absEtaString(absEtaChar);
        
        histoName = "h_recHitPedSubADC_" + absEtaString + "_" + region + "_" + dataset;
        ((h_recHitPedSubADC_vsAbsEta[dataset])[region])[absEtaString] = new TH1F(histoName.c_str(),"",nBins_recHitPedSubADC,recHitPedSubADCMin,recHitPedSubADCMax);
        ((h_recHitPedSubADC_vsAbsEta[dataset])[region])[absEtaString] -> Sumw2();
        
        histoName = "h_recHitE_" + absEtaString + "_" + region + "_" + dataset;
        ((h_recHitE_vsAbsEta[dataset])[region])[absEtaString] = new TH1F(histoName.c_str(),"",nBins_recHitE,recHitEMin,recHitEMax);
        ((h_recHitE_vsAbsEta[dataset])[region])[absEtaString] -> Sumw2();
      }
      std::cout << ">>>>>> histos vs absEta defined" << std::endl;
      
      for(int nPUIt = 0; nPUIt < nBins_nPU; ++nPUIt)
        for(int iRingIt = 0; iRingIt < nBins_iRing; ++iRingIt)
        {
          char nPUIRingChar[50];
          sprintf(nPUIRingChar,"nPU%2.1f-%2.1f__iRing%2.1f-%2.1f",nPU_axis.at(nPUIt),nPU_axis.at(nPUIt+1),iRing_axis.at(iRingIt),iRing_axis.at(iRingIt+1));        
          std::string nPUIRingString(nPUIRingChar);
          
          histoName = "h_nPUDistr_" + nPUIRingString + "_" + region + "_" + dataset;
          ((h_nPUDistr_vsNPUIRing[dataset])[region])[nPUIRingString] = new TH1F(histoName.c_str(),"",500,0.,500.);
          ((h_nPUDistr_vsNPUIRing[dataset])[region])[nPUIRingString] -> Sumw2();
          
          histoName = "h_recHitPedSubADC_" + nPUIRingString + "_" + region + "_" + dataset;
          ((h_recHitPedSubADC_vsNPUIRing[dataset])[region])[nPUIRingString] = new TH1F(histoName.c_str(),"",nBins_recHitPedSubADC,recHitPedSubADCMin,recHitPedSubADCMax);
          ((h_recHitPedSubADC_vsNPUIRing[dataset])[region])[nPUIRingString] -> Sumw2();
          
          histoName = "h_recHitE_" + nPUIRingString + "_" + region + "_" + dataset;
          ((h_recHitE_vsNPUIRing[dataset])[region])[nPUIRingString] = new TH1F(histoName.c_str(),"",nBins_recHitE,recHitEMin,recHitEMax);
          ((h_recHitE_vsNPUIRing[dataset])[region])[nPUIRingString] -> Sumw2();
        }
      std::cout << ">>>>>> histos vs nPU-iRing defined" << std::endl;
      
      for(int nPUIt = 0; nPUIt < nBins_nPU; ++nPUIt)
        for(int etaIt = 0; etaIt < nBins_eta; ++etaIt)
        {
          char nPUEtaChar[50];
          sprintf(nPUEtaChar,"nPU%2.1f-%2.1f__eta%1.4f-%1.4f",nPU_axis.at(nPUIt),nPU_axis.at(nPUIt+1),eta_axis.at(etaIt),eta_axis.at(etaIt+1));        
          std::string nPUEtaString(nPUEtaChar);
          
          histoName = "h_nPUDistr_" + nPUEtaString + "_" + region + "_" + dataset;
          ((h_nPUDistr_vsNPUEta[dataset])[region])[nPUEtaString] = new TH1F(histoName.c_str(),"",500,0.,500.);
          ((h_nPUDistr_vsNPUEta[dataset])[region])[nPUEtaString] -> Sumw2();
          
          histoName = "h_recHitPedSubADC_" + nPUEtaString + "_" + region + "_" + dataset;
          ((h_recHitPedSubADC_vsNPUEta[dataset])[region])[nPUEtaString] = new TH1F(histoName.c_str(),"",nBins_recHitPedSubADC,recHitPedSubADCMin,recHitPedSubADCMax);
          ((h_recHitPedSubADC_vsNPUEta[dataset])[region])[nPUEtaString] -> Sumw2();
          
          histoName = "h_recHitE_" + nPUEtaString + "_" + region + "_" + dataset;
          ((h_recHitE_vsNPUEta[dataset])[region])[nPUEtaString] = new TH1F(histoName.c_str(),"",nBins_recHitE,recHitEMin,recHitEMax);
          ((h_recHitE_vsNPUEta[dataset])[region])[nPUEtaString] -> Sumw2();
        }
      std::cout << ">>>>>> histos vs nPU-eta defined" << std::endl;
      
      for(int nPUIt = 0; nPUIt < nBins_nPU; ++nPUIt)
        for(int absEtaIt = 0; absEtaIt < nBins_absEta; ++absEtaIt)
        {
          char nPUAbsEtaChar[50];
          sprintf(nPUAbsEtaChar,"nPU%2.1f-%2.1f__absEta%1.4f-%1.4f",nPU_axis.at(nPUIt),nPU_axis.at(nPUIt+1),absEta_axis.at(absEtaIt),absEta_axis.at(absEtaIt+1));        
          std::string nPUAbsEtaString(nPUAbsEtaChar);
          
          histoName = "h_nPUDistr_" + nPUAbsEtaString + "_" + region + "_" + dataset;
          ((h_nPUDistr_vsNPUAbsEta[dataset])[region])[nPUAbsEtaString] = new TH1F(histoName.c_str(),"",500,0.,500.);
          ((h_nPUDistr_vsNPUAbsEta[dataset])[region])[nPUAbsEtaString] -> Sumw2();
          
          histoName = "h_recHitPedSubADC_" + nPUAbsEtaString + "_" + region + "_" + dataset;
          ((h_recHitPedSubADC_vsNPUAbsEta[dataset])[region])[nPUAbsEtaString] = new TH1F(histoName.c_str(),"",nBins_recHitPedSubADC,recHitPedSubADCMin,recHitPedSubADCMax);
          ((h_recHitPedSubADC_vsNPUAbsEta[dataset])[region])[nPUAbsEtaString] -> Sumw2();
          
          histoName = "h_recHitE_" + nPUAbsEtaString + "_" + region + "_" + dataset;
          ((h_recHitE_vsNPUAbsEta[dataset])[region])[nPUAbsEtaString] = new TH1F(histoName.c_str(),"",nBins_recHitE,recHitEMin,recHitEMax);
          ((h_recHitE_vsNPUAbsEta[dataset])[region])[nPUAbsEtaString] -> Sumw2();
        }
      std::cout << ">>>>>> histos vs nPU-absEta defined" << std::endl;
    }
  }
  
  
  
  
  
  
  // define variables;
  int isZ,runId,lumiId;
  int nVtx;
  float nPU, mee;
  int ele1_isEB, ele2_isEB;
  float ele1_scERaw, ele2_scERaw;
  float ele1_e3x3, ele2_e3x3;
  int ele1_seedIeta, ele2_seedIeta;
  int ele1_seedIphi, ele2_seedIphi;
  int ele1_seedIx, ele2_seedIx;
  int ele1_seedIy, ele2_seedIy;
  int ele1_seedIz, ele2_seedIz;

  std::vector<float>* ele1_recHit_ADC = new std::vector<float>;
  std::vector<float>* ele2_recHit_ADC = new std::vector<float>;  
  std::vector<float>* ele1_recHit_E = new std::vector<float>;
  std::vector<float>* ele2_recHit_E = new std::vector<float>;
  std::vector<int>* ele1_recHit_flag = new std::vector<int>;
  std::vector<int>* ele2_recHit_flag = new std::vector<int>;
  std::vector<int>* ele1_recHit_ieta = new std::vector<int>;
  std::vector<int>* ele2_recHit_ieta = new std::vector<int>;
  std::vector<int>* ele1_recHit_iphi = new std::vector<int>;
  std::vector<int>* ele2_recHit_iphi = new std::vector<int>;
  std::vector<int>* ele1_recHit_iz   = new std::vector<int>;
  std::vector<int>* ele2_recHit_iz   = new std::vector<int>;
  
  for(int datasetIt = 0; datasetIt < nDatasets; ++datasetIt)
  {
    std::string dataset = datasets.at(datasetIt);
    std::cout << ">>> define branches for dataset " << dataset << std::endl;
    
    chain[dataset] -> SetBranchStatus("*",0);                                                                 
    chain[dataset] -> SetBranchStatus("isZ",1);                     chain[dataset] -> SetBranchAddress("isZ",                    &isZ);
    chain[dataset] -> SetBranchStatus("runId",1);                   chain[dataset] -> SetBranchAddress("runId",                  &runId);
    chain[dataset] -> SetBranchStatus("lumiId",1);                  chain[dataset] -> SetBranchAddress("lumiId",                 &lumiId);
    chain[dataset] -> SetBranchStatus("PV_n",1);                    chain[dataset] -> SetBranchAddress("PV_n",                   &nVtx);
    chain[dataset] -> SetBranchStatus("ele1ele2_scM_regression",1); chain[dataset] -> SetBranchAddress("ele1ele2_scM_regression",&mee);
    chain[dataset] -> SetBranchStatus("ele1_isEB",1);               chain[dataset] -> SetBranchAddress("ele1_isEB",              &ele1_isEB);
    chain[dataset] -> SetBranchStatus("ele2_isEB",1);               chain[dataset] -> SetBranchAddress("ele2_isEB",              &ele2_isEB);
    chain[dataset] -> SetBranchStatus("ele1_scERaw",1);             chain[dataset] -> SetBranchAddress("ele1_scERaw",            &ele1_scERaw);
    chain[dataset] -> SetBranchStatus("ele2_scERaw",1);             chain[dataset] -> SetBranchAddress("ele2_scERaw",            &ele2_scERaw);
    chain[dataset] -> SetBranchStatus("ele1_e3x3",1);               chain[dataset] -> SetBranchAddress("ele1_e3x3",              &ele1_e3x3);
    chain[dataset] -> SetBranchStatus("ele2_e3x3",1);               chain[dataset] -> SetBranchAddress("ele2_e3x3",              &ele2_e3x3);
    chain[dataset] -> SetBranchStatus("ele1_seedIeta",1);           chain[dataset] -> SetBranchAddress("ele1_seedIeta",          &ele1_seedIeta);
    chain[dataset] -> SetBranchStatus("ele2_seedIeta",1);           chain[dataset] -> SetBranchAddress("ele2_seedIeta",          &ele2_seedIeta);
    chain[dataset] -> SetBranchStatus("ele1_seedIphi",1);           chain[dataset] -> SetBranchAddress("ele1_seedIphi",          &ele1_seedIphi);
    chain[dataset] -> SetBranchStatus("ele2_seedIphi",1);           chain[dataset] -> SetBranchAddress("ele2_seedIphi",          &ele2_seedIphi);
    chain[dataset] -> SetBranchStatus("ele1_seedIx",1);             chain[dataset] -> SetBranchAddress("ele1_seedIx",            &ele1_seedIx);
    chain[dataset] -> SetBranchStatus("ele2_seedIx",1);             chain[dataset] -> SetBranchAddress("ele2_seedIx",            &ele2_seedIx);
    chain[dataset] -> SetBranchStatus("ele1_seedIy",1);             chain[dataset] -> SetBranchAddress("ele1_seedIy",            &ele1_seedIy);
    chain[dataset] -> SetBranchStatus("ele2_seedIy",1);             chain[dataset] -> SetBranchAddress("ele2_seedIy",            &ele2_seedIy);
    chain[dataset] -> SetBranchStatus("ele1_seedZside",1);          chain[dataset] -> SetBranchAddress("ele1_seedZside",         &ele1_seedIz);
    chain[dataset] -> SetBranchStatus("ele2_seedZside",1);          chain[dataset] -> SetBranchAddress("ele2_seedZside",         &ele2_seedIz);
    chain[dataset] -> SetBranchStatus("ele1_recHitMatrix_samples",1);  chain[dataset] -> SetBranchAddress("ele1_recHitMatrix_samples",  &ele1_recHit_ADC);
    chain[dataset] -> SetBranchStatus("ele2_recHitMatrix_samples",1);  chain[dataset] -> SetBranchAddress("ele2_recHitMatrix_samples",  &ele2_recHit_ADC);
    chain[dataset] -> SetBranchStatus("ele1_recHitMatrix_E",1);        chain[dataset] -> SetBranchAddress("ele1_recHitMatrix_E",        &ele1_recHit_E);
    chain[dataset] -> SetBranchStatus("ele2_recHitMatrix_E",1);        chain[dataset] -> SetBranchAddress("ele2_recHitMatrix_E",        &ele2_recHit_E);
    chain[dataset] -> SetBranchStatus("ele1_recHitMatrix_flag",1);     chain[dataset] -> SetBranchAddress("ele1_recHitMatrix_flag",     &ele1_recHit_flag);
    chain[dataset] -> SetBranchStatus("ele2_recHitMatrix_flag",1);     chain[dataset] -> SetBranchAddress("ele2_recHitMatrix_flag",     &ele2_recHit_flag);
    chain[dataset] -> SetBranchStatus("ele1_recHitMatrix_ietaORix",1); chain[dataset] -> SetBranchAddress("ele1_recHitMatrix_ietaORix", &ele1_recHit_ieta);
    chain[dataset] -> SetBranchStatus("ele2_recHitMatrix_ietaORix",1); chain[dataset] -> SetBranchAddress("ele2_recHitMatrix_ietaORix", &ele2_recHit_ieta);
    chain[dataset] -> SetBranchStatus("ele1_recHitMatrix_iphiORiy",1); chain[dataset] -> SetBranchAddress("ele1_recHitMatrix_iphiORiy", &ele1_recHit_iphi);
    chain[dataset] -> SetBranchStatus("ele2_recHitMatrix_iphiORiy",1); chain[dataset] -> SetBranchAddress("ele2_recHitMatrix_iphiORiy", &ele2_recHit_iphi);  
    chain[dataset] -> SetBranchStatus("ele1_recHitMatrix_zside",1);    chain[dataset] -> SetBranchAddress("ele1_recHitMatrix_zside",    &ele1_recHit_iz);
    chain[dataset] -> SetBranchStatus("ele2_recHitMatrix_zside",1);    chain[dataset] -> SetBranchAddress("ele2_recHitMatrix_zside",    &ele2_recHit_iz);  
    
    if( (dataset == "MC" && type == "DAMC") || ( type == "MCMC"))
    {
      chain[dataset] -> SetBranchStatus("PUit_TrueNumInteractions",1); chain[dataset] -> SetBranchAddress("PUit_TrueNumInteractions",&nPU);
    }
    
    
    
    // loop over events
    std::cout << ">>> start of loop for dataset " << dataset << ": " << chain[dataset]->GetEntries() << " entries" << std::endl;
    for(int entry = 0; entry < chain[dataset]->GetEntries(); ++entry)
    {
      //if( entry >= 1000 ) break;
      if( entry%10000 == 0 ) std::cout << ">>>>>> reading entry " << entry << " / " << chain[dataset]->GetEntries() << "\r" << std::flush;
      chain[dataset] -> GetEntry(entry);
      
      
      // selections
      if(type == "DAMC" && dataset == "DA" ) nPU = myAvgPUList -> GetAvgPU(runId,lumiId); // da comentare per MC - MC
      if(type == "DAMC" && dataset == "DA" && nPU == 0.) continue;                        // da comentare per MC - MC
      if( isZ == 0 ) continue;
      if(type == "DAMC" && fabs(mee-91.18) > 3. ) continue;
      
      float PUWeight = 1.;
      if(type == "DAMC" && PUReweighting == 1 && dataset == "MC" ) PUWeight *= PUWeights[int(nPU+0.5)];
      
      
      // fill ele1
      std::string ele1_region = "";
      
      if( ele1_e3x3/ele1_scERaw > 0.94 )
      {
        ele1_region = "EB";
        
        if( ele1_isEB != 1 )
        {
          ele1_region = "EE";
          ele1_seedIeta = ele1_seedIx;
          ele1_seedIphi = ele1_seedIy;
        }
        
        
        (h_nVtx[dataset])[ele1_region] -> Fill(nVtx,PUWeight);
        
        
        for(unsigned int rhIt = 0; rhIt < ele1_recHit_E->size(); ++rhIt)
        {
          float recHit_E     = ele1_recHit_E   ->at(rhIt);
          float recHit_ieta  = ele1_recHit_ieta->at(rhIt);
          float recHit_iphi  = ele1_recHit_iphi->at(rhIt);
          float recHit_iRing = recHit_ieta;
          if( ele1_region == "EE" ) recHit_iRing = ele1_seedIz * (86 + EERings->GetEndcapRing(recHit_ieta,recHit_iphi,ele1_seedIz));
          float recHit_eta = GetEtaFromIRing(recHit_iRing);
          float recHit_absEta = fabs(recHit_eta);
          
          bool farRecHit = CheckSeedDistance(recHit_ieta,recHit_iphi,ele1_seedIeta,ele1_seedIphi,ele1_seedIz);
          
	  if( ele1_recHit_flag->at(rhIt) != 3000 ) continue;
	  if( farRecHit != true )                  continue;
          
          
          int nPUBin = (h_occupancy_vsNPU[dataset])[ele1_region] -> Fill(nPU);
          if( nPUBin < 1 || nPUBin > nBins_nPU ) continue;
          char nPUChar[50];
          sprintf(nPUChar,"nPU%2.1f-%2.1f",nPU_axis.at(nPUBin-1),nPU_axis.at(nPUBin));
          std::string nPUString(nPUChar);
         
           int iRingBin = (h_occupancy_vsIRing[dataset])[ele1_region] -> Fill(recHit_iRing,PUWeight);
          if( iRingBin < 1 || iRingBin > nBins_iRing ) continue;
          char iRingChar[50];
          sprintf(iRingChar,"iRing%2.1f-%2.1f",iRing_axis.at(iRingBin-1),iRing_axis.at(iRingBin));
          std::string iRingString(iRingChar);
          
          int etaBin = (h_occupancy_vsEta[dataset])[ele1_region] -> Fill(recHit_eta,PUWeight);
          if( etaBin < 1 || etaBin > nBins_eta ) continue;
          char etaChar[50];
          sprintf(etaChar,"eta%1.4f-%1.4f",eta_axis.at(etaBin-1),eta_axis.at(etaBin));
          std::string etaString(etaChar);
          
          int absEtaBin = (h_occupancy_vsAbsEta[dataset])[ele1_region] -> Fill(recHit_absEta,PUWeight);
          if( absEtaBin < 1 || absEtaBin > nBins_absEta ) continue;
          char absEtaChar[50];
          sprintf(absEtaChar,"absEta%1.4f-%1.4f",absEta_axis.at(absEtaBin-1),absEta_axis.at(absEtaBin));
          std::string absEtaString(absEtaChar);
          
          char nPUIRingChar[50];
          sprintf(nPUIRingChar,"nPU%2.1f-%2.1f__iRing%2.1f-%2.1f",nPU_axis.at(nPUBin-1),nPU_axis.at(nPUBin),iRing_axis.at(iRingBin-1),iRing_axis.at(iRingBin));
          std::string nPUIRingString(nPUIRingChar);
          
          char nPUEtaChar[50];
          sprintf(nPUEtaChar,"nPU%2.1f-%2.1f__eta%1.4f-%1.4f",nPU_axis.at(nPUBin-1),nPU_axis.at(nPUBin),eta_axis.at(etaBin-1),eta_axis.at(etaBin));
          std::string nPUEtaString(nPUEtaChar);
          
          char nPUAbsEtaChar[50];
          sprintf(nPUAbsEtaChar,"nPU%2.1f-%2.1f__absEta%1.4f-%1.4f",nPU_axis.at(nPUBin-1),nPU_axis.at(nPUBin),absEta_axis.at(absEtaBin-1),absEta_axis.at(absEtaBin));
          std::string nPUAbsEtaString(nPUAbsEtaChar);
          
          ((h_nPUDistr_vsNPU[dataset])      [ele1_region])[nPUString]       -> Fill(nPU);
          ((h_nPUDistr_vsNPUIRing[dataset]) [ele1_region])[nPUIRingString]  -> Fill(nPU);
          ((h_nPUDistr_vsNPUEta[dataset])   [ele1_region])[nPUEtaString]    -> Fill(nPU);
          ((h_nPUDistr_vsNPUAbsEta[dataset])[ele1_region])[nPUAbsEtaString] -> Fill(nPU);
          
          (h_recHitE[dataset])             [ele1_region]                   -> Fill(recHit_E,PUWeight);
          ((h_recHitE_vsNPU[dataset])      [ele1_region])[nPUString]       -> Fill(recHit_E);
          ((h_recHitE_vsIRing[dataset])    [ele1_region])[iRingString]     -> Fill(recHit_E,PUWeight);
          ((h_recHitE_vsEta[dataset])      [ele1_region])[etaString]       -> Fill(recHit_E,PUWeight);
          ((h_recHitE_vsAbsEta[dataset])   [ele1_region])[absEtaString]    -> Fill(recHit_E,PUWeight);
          ((h_recHitE_vsNPUIRing[dataset]) [ele1_region])[nPUIRingString]  -> Fill(recHit_E);
          ((h_recHitE_vsNPUEta[dataset])   [ele1_region])[nPUEtaString]    -> Fill(recHit_E);
          ((h_recHitE_vsNPUAbsEta[dataset])[ele1_region])[nPUAbsEtaString] -> Fill(recHit_E);
          
          float pedestal = 0.;
          for(int sampleIt = 0; sampleIt < 10; ++sampleIt)
          {
            float recHit_ADC = ele1_recHit_ADC->at(rhIt*10+sampleIt);
            float shift = ele1_seedIz == 0 ? 1. : ele1_seedIz;
            (p_recHitADC[dataset])[ele1_region] -> Fill(shift*recHit_ieta,recHit_iphi,recHit_ADC);
            (h_recHitADC[dataset])[ele1_region] -> Fill(recHit_ADC,PUWeight);
            if( sampleIt < pedMaxSample ) pedestal += recHit_ADC;
          }
          pedestal /= pedMaxSample;
          
          for(int sampleIt = 0; sampleIt < 10; ++sampleIt)
          {
            float recHit_ADC = ele1_recHit_ADC->at(rhIt*10+sampleIt);
            (h_recHitPedSubADC[dataset])             [ele1_region]                   -> Fill(recHit_ADC-pedestal,PUWeight);
            ((h_recHitPedSubADC_vsNPU[dataset])      [ele1_region])[nPUString]       -> Fill(recHit_ADC-pedestal);
            ((h_recHitPedSubADC_vsIRing[dataset])    [ele1_region])[iRingString]     -> Fill(recHit_ADC-pedestal,PUWeight);
            ((h_recHitPedSubADC_vsEta[dataset])      [ele1_region])[etaString]       -> Fill(recHit_ADC-pedestal,PUWeight);
            ((h_recHitPedSubADC_vsAbsEta[dataset])   [ele1_region])[absEtaString]    -> Fill(recHit_ADC-pedestal,PUWeight);
            ((h_recHitPedSubADC_vsNPUIRing[dataset]) [ele1_region])[nPUIRingString]  -> Fill(recHit_ADC-pedestal);
            ((h_recHitPedSubADC_vsNPUEta[dataset])   [ele1_region])[nPUEtaString]    -> Fill(recHit_ADC-pedestal);
            ((h_recHitPedSubADC_vsNPUAbsEta[dataset])[ele1_region])[nPUAbsEtaString] -> Fill(recHit_ADC-pedestal);
          }
        }
      }    
      // end of fill ele1
      
      
      // fill ele2
      std::string ele2_region = "";
      
      if( ele2_e3x3/ele2_scERaw > 0.94 )
      {
        ele2_region = "EB";
        
        if( ele2_isEB != 1 )
        {
          ele2_region = "EE";
          ele2_seedIeta = ele2_seedIx;
          ele2_seedIphi = ele2_seedIy;
        }
        
        
        (h_nVtx[dataset])[ele2_region] -> Fill(nVtx,PUWeight);
        
        
        for(unsigned int rhIt = 0; rhIt < ele2_recHit_E->size(); ++rhIt)
        {
          float recHit_E     = ele2_recHit_E   ->at(rhIt);
          float recHit_ieta  = ele2_recHit_ieta->at(rhIt);
          float recHit_iphi  = ele2_recHit_iphi->at(rhIt);
          float recHit_iRing = recHit_ieta;
          if( ele2_region == "EE" ) recHit_iRing = ele2_seedIz * (86 + EERings->GetEndcapRing(recHit_ieta,recHit_iphi,ele2_seedIz));
          float recHit_eta = GetEtaFromIRing(recHit_iRing);
          float recHit_absEta = fabs(recHit_eta);
          
          bool farRecHit = CheckSeedDistance(recHit_ieta,recHit_iphi,ele2_seedIeta,ele2_seedIphi,ele2_seedIz);
          
          if( ele2_recHit_flag->at(rhIt) != 3000 ) continue;
          if( farRecHit != true )                  continue;
          
          
          int nPUBin = (h_occupancy_vsNPU[dataset])[ele2_region] -> Fill(nPU);
          if( nPUBin < 1 || nPUBin > nBins_nPU ) continue;
          char nPUChar[50];
          sprintf(nPUChar,"nPU%2.1f-%2.1f",nPU_axis.at(nPUBin-1),nPU_axis.at(nPUBin));
          std::string nPUString(nPUChar);
         
           int iRingBin = (h_occupancy_vsIRing[dataset])[ele2_region] -> Fill(recHit_iRing,PUWeight);
          if( iRingBin < 1 || iRingBin > nBins_iRing ) continue;
          char iRingChar[50];
          sprintf(iRingChar,"iRing%2.1f-%2.1f",iRing_axis.at(iRingBin-1),iRing_axis.at(iRingBin));
          std::string iRingString(iRingChar);
          
          int etaBin = (h_occupancy_vsEta[dataset])[ele2_region] -> Fill(recHit_eta,PUWeight);
          if( etaBin < 1 || etaBin > nBins_eta ) continue;
          char etaChar[50];
          sprintf(etaChar,"eta%1.4f-%1.4f",eta_axis.at(etaBin-1),eta_axis.at(etaBin));
          std::string etaString(etaChar);
          
          int absEtaBin = (h_occupancy_vsAbsEta[dataset])[ele2_region] -> Fill(recHit_absEta,PUWeight);
          if( absEtaBin < 1 || absEtaBin > nBins_absEta ) continue;
          char absEtaChar[50];
          sprintf(absEtaChar,"absEta%1.4f-%1.4f",absEta_axis.at(absEtaBin-1),absEta_axis.at(absEtaBin));
          std::string absEtaString(absEtaChar);
          
          char nPUIRingChar[50];
          sprintf(nPUIRingChar,"nPU%2.1f-%2.1f__iRing%2.1f-%2.1f",nPU_axis.at(nPUBin-1),nPU_axis.at(nPUBin),iRing_axis.at(iRingBin-1),iRing_axis.at(iRingBin));
          std::string nPUIRingString(nPUIRingChar);
          
          char nPUEtaChar[50];
          sprintf(nPUEtaChar,"nPU%2.1f-%2.1f__eta%1.4f-%1.4f",nPU_axis.at(nPUBin-1),nPU_axis.at(nPUBin),eta_axis.at(etaBin-1),eta_axis.at(etaBin));
          std::string nPUEtaString(nPUEtaChar);
          
          char nPUAbsEtaChar[50];
          sprintf(nPUAbsEtaChar,"nPU%2.1f-%2.1f__absEta%1.4f-%1.4f",nPU_axis.at(nPUBin-1),nPU_axis.at(nPUBin),absEta_axis.at(absEtaBin-1),absEta_axis.at(absEtaBin));
          std::string nPUAbsEtaString(nPUAbsEtaChar);
          
          ((h_nPUDistr_vsNPU[dataset])      [ele2_region])[nPUString]       -> Fill(nPU);
          ((h_nPUDistr_vsNPUIRing[dataset]) [ele2_region])[nPUIRingString]  -> Fill(nPU);
          ((h_nPUDistr_vsNPUEta[dataset])   [ele2_region])[nPUEtaString]    -> Fill(nPU);
          ((h_nPUDistr_vsNPUAbsEta[dataset])[ele2_region])[nPUAbsEtaString] -> Fill(nPU);
          
          (h_recHitE[dataset])             [ele2_region]                   -> Fill(recHit_E,PUWeight);
          ((h_recHitE_vsNPU[dataset])      [ele2_region])[nPUString]       -> Fill(recHit_E);
          ((h_recHitE_vsIRing[dataset])    [ele2_region])[iRingString]     -> Fill(recHit_E,PUWeight);
          ((h_recHitE_vsEta[dataset])      [ele2_region])[etaString]       -> Fill(recHit_E,PUWeight);
          ((h_recHitE_vsAbsEta[dataset])   [ele2_region])[absEtaString]    -> Fill(recHit_E,PUWeight);
          ((h_recHitE_vsNPUIRing[dataset]) [ele2_region])[nPUIRingString]  -> Fill(recHit_E);
          ((h_recHitE_vsNPUEta[dataset])   [ele2_region])[nPUEtaString]    -> Fill(recHit_E);
          ((h_recHitE_vsNPUAbsEta[dataset])[ele2_region])[nPUAbsEtaString] -> Fill(recHit_E);
          
          float pedestal = 0.;
          for(int sampleIt = 0; sampleIt < 10; ++sampleIt)
          {
            float recHit_ADC = ele2_recHit_ADC->at(rhIt*10+sampleIt);
            float shift = ele2_seedIz == 0 ? 1. : ele2_seedIz;
            (p_recHitADC[dataset])[ele2_region] -> Fill(shift*recHit_ieta,recHit_iphi,recHit_ADC);
            (h_recHitADC[dataset])[ele2_region] -> Fill(recHit_ADC,PUWeight);
            if( sampleIt < pedMaxSample ) pedestal += recHit_ADC;
          }
          pedestal /= pedMaxSample;
          
          for(int sampleIt = 0; sampleIt < 10; ++sampleIt)
          {
            float recHit_ADC = ele2_recHit_ADC->at(rhIt*10+sampleIt);
            (h_recHitPedSubADC[dataset])             [ele2_region]                   -> Fill(recHit_ADC-pedestal,PUWeight);
            ((h_recHitPedSubADC_vsNPU[dataset])      [ele2_region])[nPUString]       -> Fill(recHit_ADC-pedestal);
            ((h_recHitPedSubADC_vsIRing[dataset])    [ele2_region])[iRingString]     -> Fill(recHit_ADC-pedestal,PUWeight);
            ((h_recHitPedSubADC_vsEta[dataset])      [ele2_region])[etaString]       -> Fill(recHit_ADC-pedestal,PUWeight);
            ((h_recHitPedSubADC_vsAbsEta[dataset])   [ele2_region])[absEtaString]    -> Fill(recHit_ADC-pedestal,PUWeight);
            ((h_recHitPedSubADC_vsNPUIRing[dataset]) [ele2_region])[nPUIRingString]  -> Fill(recHit_ADC-pedestal);
            ((h_recHitPedSubADC_vsNPUEta[dataset])   [ele2_region])[nPUEtaString]    -> Fill(recHit_ADC-pedestal);
            ((h_recHitPedSubADC_vsNPUAbsEta[dataset])[ele2_region])[nPUAbsEtaString] -> Fill(recHit_ADC-pedestal);
          }
        }
      }    
      // end of fill ele2
    }
    std::cout << std::endl;
    std::cout << ">>> end of loop for dataset " << dataset << std::endl;
  }
  
  
  
  outFile -> cd();
  
  for(int datasetIt = 0; datasetIt < nDatasets; ++datasetIt)
  {
    std::string dataset = datasets.at(datasetIt);
    TDirectory* c1 = outFile -> mkdir(dataset.c_str());
    c1 -> cd();
    
    for(int regionIt = 0; regionIt < nRegions; ++regionIt)
    {
      std::string region = regions.at(regionIt);  
      TDirectory* c2 = c1 -> mkdir(region.c_str());
      c2 -> cd();
      
      WriteNormalized( (h_nVtx[dataset])[region] );
      
      (h_occupancy_vsNPU[dataset])[region] -> Write();
      (h_occupancy_vsIRing[dataset])[region] -> Write();
      (h_occupancy_vsEta[dataset])[region] -> Write();
      (h_occupancy_vsAbsEta[dataset])[region] -> Write();
      
      (p_recHitADC[dataset])[region] -> Write();
      WriteNormalized( (h_recHitADC[dataset])[region] );
      WriteNormalized( (h_recHitPedSubADC[dataset])[region] );
      WriteNormalized( (h_recHitE[dataset])[region] );
      
      TDirectory* c3 = c2 -> mkdir("plots_vs_nPU");
      c3 -> cd();
      for(int nPUIt = 0; nPUIt < nBins_nPU; ++nPUIt)
      {
        char nPUChar[50];
        sprintf(nPUChar,"nPU%2.1f-%2.1f",nPU_axis.at(nPUIt),nPU_axis.at(nPUIt+1));
        std::string nPUString(nPUChar);
        WriteNormalized( ((h_nPUDistr_vsNPU[dataset])[region])[nPUString] );
        WriteNormalized( ((h_recHitE_vsNPU[dataset])[region])[nPUString] );
        WriteNormalized( ((h_recHitPedSubADC_vsNPU[dataset])[region])[nPUString] );
      }
      
      c3 = c2 -> mkdir("plots_vs_iRing");
      c3 -> cd();
      for(int iRingIt = 0; iRingIt < nBins_iRing; ++iRingIt)
      {
        char iRingChar[50];
        sprintf(iRingChar,"iRing%2.1f-%2.1f",iRing_axis.at(iRingIt),iRing_axis.at(iRingIt+1));
        std::string iRingString(iRingChar);
        WriteNormalized( ((h_recHitPedSubADC_vsIRing[dataset])[region])[iRingString] );
        WriteNormalized( ((h_recHitE_vsIRing[dataset])[region])[iRingString] );
      }
      
      c3 = c2 -> mkdir("plots_vs_eta");
      c3 -> cd();
      for(int etaIt = 0; etaIt < nBins_eta; ++etaIt)
      {
        char etaChar[50];
        sprintf(etaChar,"eta%1.4f-%1.4f",eta_axis.at(etaIt),eta_axis.at(etaIt+1));
        std::string etaString(etaChar);
        WriteNormalized( ((h_recHitPedSubADC_vsEta[dataset])[region])[etaString] );
        WriteNormalized( ((h_recHitE_vsEta[dataset])[region])[etaString] );
      }
      
      c3 = c2 -> mkdir("plots_vs_absEta");
      c3 -> cd();
      for(int absEtaIt = 0; absEtaIt < nBins_absEta; ++absEtaIt)
      {
        char absEtaChar[50];
        sprintf(absEtaChar,"absEta%1.4f-%1.4f",absEta_axis.at(absEtaIt),absEta_axis.at(absEtaIt+1));
        std::string absEtaString(absEtaChar);
        WriteNormalized( ((h_recHitPedSubADC_vsAbsEta[dataset])[region])[absEtaString] );
        WriteNormalized( ((h_recHitE_vsAbsEta[dataset])[region])[absEtaString] );
      }
      
      c3 = c2 -> mkdir("plots_vs_nPU_and_iRing");
      c3 -> cd();
      for(int nPUIt = 0; nPUIt < nBins_nPU; ++nPUIt)
        for(int iRingIt = 0; iRingIt < nBins_iRing; ++iRingIt)
        {
          char nPUIRingChar[50];
          sprintf(nPUIRingChar,"nPU%2.1f-%2.1f__iRing%2.1f-%2.1f",nPU_axis.at(nPUIt),nPU_axis.at(nPUIt+1),iRing_axis.at(iRingIt),iRing_axis.at(iRingIt+1));
          std::string nPUIRingString(nPUIRingChar);
          WriteNormalized( ((h_nPUDistr_vsNPUIRing[dataset])[region])[nPUIRingString] );
          WriteNormalized( ((h_recHitPedSubADC_vsNPUIRing[dataset])[region])[nPUIRingString] );
          WriteNormalized( ((h_recHitE_vsNPUIRing[dataset])[region])[nPUIRingString] );
        }
      
      c3 = c2 -> mkdir("plots_vs_nPU_and_eta");
      c3 -> cd();
      for(int nPUIt = 0; nPUIt < nBins_nPU; ++nPUIt)
        for(int etaIt = 0; etaIt < nBins_eta; ++etaIt)
        {
          char nPUEtaChar[50];
          sprintf(nPUEtaChar,"nPU%2.1f-%2.1f__eta%1.4f-%1.4f",nPU_axis.at(nPUIt),nPU_axis.at(nPUIt+1),eta_axis.at(etaIt),eta_axis.at(etaIt+1));
          std::string nPUEtaString(nPUEtaChar);
          WriteNormalized( ((h_nPUDistr_vsNPUEta[dataset])[region])[nPUEtaString] );
          WriteNormalized( ((h_recHitPedSubADC_vsNPUEta[dataset])[region])[nPUEtaString] );
          WriteNormalized( ((h_recHitE_vsNPUEta[dataset])[region])[nPUEtaString] );
        }
      
      c3 = c2 -> mkdir("plots_vs_nPU_and_absEta");
      c3 -> cd();
      for(int nPUIt = 0; nPUIt < nBins_nPU; ++nPUIt)
        for(int absEtaIt = 0; absEtaIt < nBins_absEta; ++absEtaIt)
        {
          char nPUAbsEtaChar[50];
          sprintf(nPUAbsEtaChar,"nPU%2.1f-%2.1f__absEta%1.4f-%1.4f",nPU_axis.at(nPUIt),nPU_axis.at(nPUIt+1),absEta_axis.at(absEtaIt),absEta_axis.at(absEtaIt+1));
          std::string nPUAbsEtaString(nPUAbsEtaChar);
          WriteNormalized( ((h_nPUDistr_vsNPUAbsEta[dataset])[region])[nPUAbsEtaString] );
          WriteNormalized( ((h_recHitPedSubADC_vsNPUAbsEta[dataset])[region])[nPUAbsEtaString] );
          WriteNormalized( ((h_recHitE_vsNPUAbsEta[dataset])[region])[nPUAbsEtaString] );
        }
            
      outFile -> cd("");
    }
  }
  
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

//  LocalWords:  for
