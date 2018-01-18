#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"

#include "include/returnRootFileContentsList.h"
#include "include/kirchnerPalette.h"

int plotPFEta(const std::string inFileName)
{
  kirchnerPalette col;

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

  pos=0;
  while(pos < histList.size()){
    bool isGood = true;
    
    if(histList.at(pos).find("_Eta") != std::string::npos || histList.at(pos).find("_Neg") != std::string::npos){
      histList.erase(histList.begin()+pos);
      isGood = false;
    }
    else{
      for(unsigned int i = pos+1; i < histList.size(); ++i){
	if(histList.at(i).size() == histList.at(pos).size() && histList.at(i).find(histList.at(pos)) != std::string::npos){histList.erase(histList.begin()+i); isGood = false; break;}
      }
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

  std::cout << "Hists..." << std::endl;
  for(unsigned int i = 0; i < histList.size(); ++i){
    std::cout << " " << i << "/" << histList.size() << ": " << histList.at(i) << std::endl;
  }

  const Int_t nCentHist = nCentBins;
  const Int_t nPFID = 7;
  TH1F* hists_p[nCentHist][nPFID];

  for(unsigned int i = 0; i < histList.size(); ++i){
    bool isJet = false;
    for(unsigned int jI = 0; jI < jetAlgos.size(); ++jI){
      if(histList.at(i).find(jetAlgos.at(jI)) != std::string::npos){isJet = true; break;}
    }
    if(isJet) continue;

    Int_t centPos = -1;
    Int_t pfIdPos = -1;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr1 = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      if(histList.at(i).find(centStr1) != std::string::npos){centPos = cI; break;}
    }

    for(Int_t pI = 0; pI < nPFID; ++pI){
      std::string pfStr = "PFID" + std::to_string(pI + 1);
      if(histList.at(i).find(pfStr) != std::string::npos){pfIdPos = pI; break;}
    }

    hists_p[centPos][pfIdPos] = (TH1F*)inFile_p->Get(histList.at(i).c_str());
    hists_p[centPos][pfIdPos]->Sumw2();
  }

  TDatime* date = new TDatime();
  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    TCanvas* canv_p = new TCanvas("canv_c", "canv_c", 500, 500);

    canv_p->SetTopMargin(canv_p->GetLeftMargin()/2.);
    canv_p->SetRightMargin(canv_p->GetLeftMargin()/2.);
    canv_p->SetBottomMargin(canv_p->GetLeftMargin());
    
    const std::string centStr1 = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string centStr2 =  std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
    std::cout << centStr1 << std::endl;
    for(Int_t pI = 0; pI < nPFID; ++pI){
      std::cout << " " << hists_p[cI][pI]->GetName() << std::endl;

      hists_p[cI][pI]->SetMaximum(1.1);
      hists_p[cI][pI]->SetMinimum(0.0);
      hists_p[cI][pI]->SetFillColor(col.getColor(pI));
      hists_p[cI][pI]->SetLineColor(1);

      gStyle->SetOptStat(0);
      if(pI != 0){ 
	hists_p[cI][pI]->Add(hists_p[cI][pI-1]);
	hists_p[cI][pI]->DrawCopy("SAME HIST E1");
      }
      else{
	hists_p[cI][pI]->GetYaxis()->SetTitle(("PF Fraction (" + centStr2 + ")").c_str());
	hists_p[cI][pI]->DrawCopy("HIST E1");
      }
    }

    for(Int_t pI = nPFID-1; pI >= 0; --pI){
      hists_p[cI][pI]->DrawCopy("SAME HIST E1");
      hists_p[cI][pI]->DrawCopy("SAME E1");
    }

    gPad->RedrawAxis();
    std::string keyStr = "PF Key: ";
    for(Int_t pI = 0; pI < nPFID; ++pI){
      keyStr = keyStr + "#color[" + std::to_string(col.getColor(pI)) + "]{" + std::to_string(pI+1) + "}, ";
    }
    keyStr.replace(keyStr.size()-2,2,"");
    label_p->DrawLatex(.15, .96, keyStr.c_str());

    canv_p->SaveAs(("pdfDir/pfFracVEta_" + centStr1 + "_" + std::to_string(date->GetDate()) + ".pdf").c_str());
    delete canv_p;
  }
  delete date;
  delete label_p;

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
