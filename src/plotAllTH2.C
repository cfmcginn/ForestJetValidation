#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"

#include <iostream>
#include <string>
#include <vector>

#include "include/returnRootFileContentsList.h"

const std::string th2fStr = "TH2F";

std::string getStringExt(const std::string inFileName)
{
  std::string retString = "";
  if(inFileName.find("PAMB") != std::string::npos) retString = "PAMB";
  else if(inFileName.find("EPOSMB_PbP") != std::string::npos) retString = "EPOSMB_PbP";
  else if(inFileName.find("EPOSMB_PPb") != std::string::npos) retString = "EPOSMB_PPb";

  return retString;
}

int plotAllTH2F(const std::string inFileName)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  std::vector<std::string> inTH2FNames = returnRootFileContentsList(inFile_p, "TH2F");
    
  const Int_t nTH2F = (Int_t)inTH2FNames.size();

  TH2F* inTH2F_p[nTH2F];
  TCanvas* inTH2F_c[nTH2F];

  for(Int_t iter = 0; iter < nTH2F; iter++){
    inTH2F_p[iter] = (TH2F*)inFile_p->Get(inTH2FNames.at(iter).c_str());

    inTH2F_c[iter] = new TCanvas(inTH2FNames.at(iter).c_str(), inTH2FNames.at(iter).c_str(), 1000, 1000);
    gPad->SetLeftMargin(gPad->GetLeftMargin()*1.25);
    gPad->SetRightMargin(gPad->GetLeftMargin());
    gPad->SetBottomMargin(gPad->GetLeftMargin());
    gPad->SetTopMargin(.01);

    TH2F* dummyHist_p = new TH2F("dummyHist_h", ";;", 10, inTH2F_p[iter]->GetXaxis()->GetBinLowEdge(1), inTH2F_p[iter]->GetXaxis()->GetBinLowEdge(inTH2F_p[iter]->GetXaxis()->GetNbins()+1), 10, inTH2F_p[iter]->GetYaxis()->GetBinLowEdge(1), inTH2F_p[iter]->GetYaxis()->GetBinLowEdge(inTH2F_p[iter]->GetYaxis()->GetNbins()+1));
    //    dummyHist_p->SetMinimum(0);
    dummyHist_p->SetMinimum(inTH2F_p[iter]->GetMinimum()/5.);
    if(dummyHist_p->GetMinimum() == 0) dummyHist_p->SetMinimum(.5);
    dummyHist_p->SetMaximum(inTH2F_p[iter]->GetMaximum() + TMath::Sqrt(inTH2F_p[iter]->GetMaximum()));

    dummyHist_p->SetXTitle(inTH2F_p[iter]->GetXaxis()->GetTitle());
    dummyHist_p->SetYTitle(inTH2F_p[iter]->GetYaxis()->GetTitle());

    dummyHist_p->GetXaxis()->CenterTitle();
    dummyHist_p->GetYaxis()->CenterTitle();

    dummyHist_p->GetXaxis()->SetNdivisions(505);

    gStyle->SetOptStat(0);
    dummyHist_p->DrawCopy();

    inTH2F_p[iter]->DrawCopy("COLZ SAME");

    gPad->SetLogz();
    /*
    TLine* line_p = new TLine();
    line_p->SetLineStyle(2);
    line_p->DrawLine(30, dummyHist_p->GetYaxis()->GetBinLowEdge(1), 30, dummyHist_p->GetYaxis()->GetBinLowEdge(dummyHist_p->GetYaxis()->GetNbins()+1));
    */
    gPad->Modified();

    TDatime* date = new TDatime();

    std::string tempOutName = inTH2FNames.at(iter).substr(0, inTH2FNames.at(iter).size()-1);

    while(tempOutName.find("/") != std::string::npos){
      tempOutName.replace(0, tempOutName.find("/")+1, "");
    }

    std::string extensionStr = getStringExt(inFileName);

    std::string saveNamePDF = "pdfDir/" + tempOutName + std::to_string(date->GetDate());
    if(extensionStr.size() != 0)  saveNamePDF = saveNamePDF + "_" + extensionStr;
    saveNamePDF = saveNamePDF + ".pdf";

    inTH2F_c[iter]->SaveAs(saveNamePDF.c_str());

    delete date;
    delete dummyHist_p;
    delete inTH2F_c[iter];
  }

  inTH2FNames.clear();

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./plotAllTH2.exe <inFileName>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += plotAllTH2F(argv[1]);
  return retVal;
}
