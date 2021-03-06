#ifndef PUreweightingUtils_h
#define PUreweightingUtils_h

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TRandom.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TPaveLabel.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

void makePUDistributionMC(std::string, std::string, std::string, std::string);

void ComputePUweights(std::string, std::string, std::string, std::string, std::string);
std::map<float,float>* ComputePUweights(TTree* t_mc, const std::string& PUFileName_da, const bool& verbosity = false);

class TPileupReweighting
{
 private:
  float w[100];

 public:
  TPileupReweighting(std::string,std::string); 
  ~TPileupReweighting();
  
  double GetWeight(int);
  //ClassDef(TPileupReweighting,1); 
};

#endif
