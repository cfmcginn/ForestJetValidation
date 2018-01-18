#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TNamed.h"

#include "include/returnRootFileContentsList.h"

int plotPFEta(const std::string inFileName)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> nameList = returnRootFileContentsList(inFile_p, "TNamed", "");
  std::vector<std::string> histList = returnRootFileContentsList(inFile_p, "TH1F", "pfFracVEta");
  std::vector<TNamed*> names;

  unsigned int pos = 0;
  while(pos < nameList.size()){
    bool isGood = true;
    for(unsigned int i = pos+1; i < nameList.size(); ++i){
      if(nameList.at(i).size() == nameList.at(pos).size() && nameList.at(i).find(nameList.at(pos)) != std::string::npos){nameList.erase(nameList.begin()+i); isGood = false; break;}
    }
    
    if(isGood) ++pos;
  }

  Int_t nCentBins=-1;
  std::vector<Int_t> centBinsLow;
  std::vector<Int_t> centBinsHi;
  std::vector<std::string> jetAlgos;

  std::cout << "Names..." << std::endl;
  for(unsigned int i = 0; i < nameList.size(); ++i){
    std::cout << " " << i << "/" << nameList.size() << ": " << nameList.at(i) << std::endl;
    names.push_back((TNamed*)inFile_p->Get(nameList.at(i).c_str()));
    std::cout << "  " << names.at(names.size()-1)->GetTitle() << std::endl;

    if(nameList.at(i).find("nCentBins") != std::string::npos) nCentBins = std::stoi(names.at(names.size()-1)->GetTitle());
    else if(nameList.at(i).find("centBinsLow") != std::string::npos){
      std::string tempStr = std::string(names.at(names.size()-1)->GetTitle());   
      while(tempStr.find(",") != std::string::npos){
	centBinsLow.push_back(std::stoi(tempStr.substr(0,tempStr.find(","))));
	tempStr.replace(0, tempStr.find(",")+1, "");
      }
      if(tempStr.size() != 0) centBinsLow.push_back(std::stoi(tempStr));
    }
    else if(nameList.at(i).find("centBinsHi") != std::string::npos){
      std::string tempStr = std::string(names.at(names.size()-1)->GetTitle());   
      while(tempStr.find(",") != std::string::npos){
	centBinsHi.push_back(std::stoi(tempStr.substr(0,tempStr.find(","))));
	tempStr.replace(0, tempStr.find(",")+1, "");
      }
      if(tempStr.size() != 0) centBinsHi.push_back(std::stoi(tempStr));
    }
    else if(nameList.at(i).find("jetAlgos") != std::string::npos){
      std::string tempStr = std::string(names.at(names.size()-1)->GetTitle());
      while(tempStr.find(",") != std::string::npos){
        jetAlgos.push_back(tempStr.substr(0,tempStr.find(",")));
        tempStr.replace(0, tempStr.find(",")+1, "");
      }
      if(tempStr.size() != 0) jetAlgos.push_back(tempStr);
    }
  }

  std::cout << "nCentBins: " << nCentBins << std::endl;

  std::cout << "centBinsLow: ";
  for(unsigned int i = 0; i < centBinsLow.size(); ++i){std::cout << centBinsLow.at(i) << ",";}
  std::cout << std::endl;

  std::cout << "centBinsHi: ";
  for(unsigned int i = 0; i < centBinsHi.size(); ++i){std::cout << centBinsHi.at(i) << ",";}
  std::cout << std::endl;

  std::cout << "jetAlgos: ";
  for(unsigned int i = 0; i < jetAlgos.size(); ++i){std::cout << jetAlgos.at(i) << ",";}
  std::cout << std::endl;

 


  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage ./plotPFEta.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotPFEta(argv[1]);
  return retVal;
}
