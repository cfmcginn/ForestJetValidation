#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNamed.h"
#include "TRandom3.h"

#include "include/returnRootFileContentsList.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"

int forestJetValidation(const std::string inFileName, std::string outFileName = "")
{
  if(outFileName.size() == 0){
    outFileName = inFileName;
    while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  }
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");

  TDatime* date = new TDatime();
  outFileName = outFileName + "_JetValidation_" + std::to_string(date->GetDate()) + ".root";
  delete date;

  std::vector<std::string> fileList; 
  if(inFileName.find(".root") != std::string::npos) fileList.push_back(inFileName);
  else if(inFileName.find(".txt") != std::string::npos){
    std::ifstream file(inFileName.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      while(tempStr.substr(0,1).find(" ") != std::string::npos) tempStr.replace(0,1,"");
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") == std::string::npos) continue;

      fileList.push_back(tempStr);
    }

    file.close();
  }

  if(fileList.size() == 0){
    std::cout << "Given input \'" << inFileName << "\' doesn't produce valid input. return 1" << std::endl;
    return 1;
  }

  UInt_t runMin = 1000000000;
  UInt_t runMax = 0;
  UInt_t lumiMin = 1000000000;
  UInt_t lumiMax = 0;
  UInt_t evtMin = 1000000000;
  UInt_t evtMax = 0;

  bool hasPFTreeGlobal = true;
  bool hasHITreeGlobal = true;
  bool hasJetTreeGlobal = true;

  std::vector<std::string> jetTreeListGlobal;
  std::vector<std::string> jetTreeAlgGlobal;
  
  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");
    std::vector<std::string> treeList = returnRootFileContentsList(inFile_p, "TTree", "");
    std::vector<std::string> jetTreeList;
    std::vector<std::string> jetTreeAlg;

    unsigned int pos = 0;
    while(treeList.size() > pos){
      bool isUnique = true;
      for(unsigned int pI = pos+1; pI < treeList.size(); ++pI){
	if(treeList.at(pos).size() == treeList.at(pI).size() && treeList.at(pos).find(treeList.at(pI)) != std::string::npos){
	  treeList.erase(treeList.begin()+pI);
	  isUnique = false;
	  break;
	}
      }
      if(isUnique) ++pos;
    }
    
    bool hasPFTree = false;
    bool hasHITree = false;
    bool hasJetTree = false;
    
    if(fI == 0) std::cout << "Trees..." << std::endl;
    for(unsigned int pI = 0; pI < treeList.size(); ++pI){
      if(fI == 0) std::cout << " Tree " << pI << "/" << treeList.size() << ": " << treeList.at(pI) << std::endl;
      if(treeList.at(pI).find("JetAnalyzer") != std::string::npos){
	std::string tempJet = treeList.at(pI);
	jetTreeList.push_back(tempJet);
	tempJet.replace(tempJet.find("JetAnalyzer"), tempJet.size()-tempJet.find("JetAnalyzer"), "");
	jetTreeAlg.push_back(tempJet);
	hasJetTree = true;
      }
      else if(treeList.at(pI).find("pfcandAnalyzer") != std::string::npos) hasPFTree = true;
      else if(treeList.at(pI).find("hiEvtAnalyzer") != std::string::npos){
	hasHITree = true;
	TTree* temp = (TTree*)inFile_p->Get(treeList.at(pI).c_str());
	runMin = TMath::Min((Int_t)runMin, (Int_t)temp->GetMinimum("run"));
	runMax = TMath::Max((Int_t)runMax, (Int_t)temp->GetMaximum("run"));
	
	lumiMin = TMath::Min((Int_t)lumiMin, (Int_t)temp->GetMinimum("lumi"));
	lumiMax = TMath::Max((Int_t)lumiMax, (Int_t)temp->GetMaximum("lumi"));
	
	evtMin = TMath::Min((Int_t)evtMin, (Int_t)temp->GetMinimum("evt"));
	evtMax = TMath::Max((Int_t)evtMax, (Int_t)temp->GetMaximum("evt"));
      }
    }
    
    inFile_p->Close();
    delete inFile_p;

    hasPFTreeGlobal = hasPFTree;
    hasHITreeGlobal = hasHITree;
    hasJetTreeGlobal = hasJetTree;

    if(!hasHITree) break;
    if(!hasJetTree) break;

    if(jetTreeListGlobal.size() == 0){
      jetTreeListGlobal = jetTreeList;
      jetTreeAlgGlobal = jetTreeAlg;
    }
  }


  if(!hasHITreeGlobal) std::cout << "No hiEvtAnalyzer/HiTree in input \'" << inFileName << "\'. Required. return 1." << std::endl;
  if(!hasJetTreeGlobal) std::cout << "No *JetAnalyzer/t in input \'" << inFileName << "\'. Required. return 1." << std::endl;

  if(!hasHITreeGlobal || !hasJetTreeGlobal) return 1;

  TRandom3* randGen_p = new TRandom3(0);
  const Double_t rcParam = 0.4;
  const Int_t nRC = 11;

  const Int_t nJetTrees = jetTreeListGlobal.size();

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  const Int_t nCentBins = 5;
  const Int_t centBinsLow[nCentBins] = { 0,  0, 10, 30,  50};
  const Int_t centBinsHi[nCentBins] = {100, 10, 30, 50, 100};

  const Int_t nJtPtBins = 4;
  const Float_t jtPtBinsLow[nJtPtBins] = {30., 30., 50.,  80.};
  const Float_t jtPtBinsHi[nJtPtBins] = {200., 50., 80., 200.};

  const Int_t nJtEtaBins = 7;
  const Float_t jtEtaBinsLow[nJtEtaBins] = {-2.0, -2.0, -1.6, -1.0, 0.0, 1.0, 1.6};
  const Float_t jtEtaBinsHi[nJtEtaBins] = {  2.0, -1.6, -1.0,  0.0, 1.0, 1.6, 2.0};

  const Int_t nJtAbsEtaBins = 4;
  const Float_t jtAbsEtaBinsLow[nJtAbsEtaBins] = {0.0, 0.0, 1.0, 1.6};
  const Float_t jtAbsEtaBinsHi[nJtAbsEtaBins] = { 2.0, 1.0, 1.6, 2.0};

  TH1F* run_h = 0;
  TH1F* lumi_h = 0;
  TH1F* evt_h = 0;
  runMin==runMax? run_h = new TH1F("run_h", ";Run;Counts", 3, runMin-1, runMax+1) : run_h = new TH1F("run_h", ";Run;Counts", TMath::Min(200, Int_t(runMax - runMin)), runMin, runMax);
  lumiMin==lumiMax? lumi_h = new TH1F("lumi_h", ";Lumi;Counts", 3, lumiMin-1, lumiMax+1) : lumi_h = new TH1F("lumi_h", ";Lumi;Counts", TMath::Min(200, Int_t(lumiMax - lumiMin)), lumiMin, lumiMax);
  evtMin==evtMax? evt_h = new TH1F("evt_h", ";Evt;Counts", 3, evtMin-1, evtMax+1) : evt_h = new TH1F("evt_h", ";Evt;Counts", TMath::Min(200, Int_t(evtMax - evtMin)), evtMin, evtMax);

  TH1F* hiBin_h = new TH1F("hiBin_h", ";hiBin;Counts", 200, -0.5, 199.5);
  TH1F* vz_h = new TH1F("vz_h", ";vz;Counts", 81, -20, 20);

  centerTitles({run_h, lumi_h, evt_h, hiBin_h, vz_h});

  TH1F* jtPt_h[nJetTrees][nCentBins][nJtAbsEtaBins];
  TH1F* jtRawPt_h[nJetTrees][nCentBins][nJtAbsEtaBins];
  TH1F* jtPhi_h[nJetTrees][nCentBins][nJtPtBins];
  TH1F* jtEta_h[nJetTrees][nCentBins][nJtPtBins];

  TH2F* jtPtOverRawPt_h[nJetTrees][nCentBins];

  TH1F* jtCHF_h[nJetTrees][nCentBins];
  TH1F* jtNHF_h[nJetTrees][nCentBins];
  TH1F* jtCEF_h[nJetTrees][nCentBins];
  TH1F* jtNEF_h[nJetTrees][nCentBins];
  TH1F* jtMUF_h[nJetTrees][nCentBins];
  TH1F* jtTF_h[nJetTrees][nCentBins];

  TH2F* randomConeSumVsHiBin_h[nJtEtaBins];
  TH2F* randomConeRhoVsHiBin_h[nJtEtaBins];

  const Int_t nPFID = 7;
  const Int_t nPFEtaBins = 50;
  const Float_t pfEtaLow = -5.;
  const Float_t pfEtaHi = 5.;
  Double_t pfEtaBins[nPFEtaBins+1];
  getLinBins(pfEtaLow, pfEtaHi, nPFEtaBins, pfEtaBins);

  TH1F* pfFracVEta_h[nCentBins][nPFID];
  TH1F* pfFracVEta_Bins_h[nCentBins][nPFID][nPFEtaBins];

  TH1F* pfFracVEta_InJet_h[nJetTrees][nCentBins][nPFID];
  TH1F* pfFracVEta_Bins_InJet_h[nJetTrees][nCentBins][nPFID][nPFEtaBins];


  for(Int_t jI = 0; jI < nJetTrees; ++jI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr1 = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      const std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";

      for(Int_t eI = 0; eI < nJtAbsEtaBins; ++eI){
	std::string absEtaStr1 = "AbsEta" + prettyString(jtAbsEtaBinsLow[eI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[eI], 1, true);
	std::string absEtaStr2 = prettyString(jtAbsEtaBinsLow[eI], 1, true) + " |#eta_{Jet}| " + prettyString(jtAbsEtaBinsHi[eI], 1, false);

	jtPt_h[jI][cI][eI] = new TH1F(("jtPt_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_" + absEtaStr1 + "_h").c_str(), (";Jet p_{T} (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ", " + absEtaStr2 + ")").c_str(), 50, 30, 130);
	jtRawPt_h[jI][cI][eI] = new TH1F(("rawPt_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_" + absEtaStr1 + "_h").c_str(), (";Jet Raw p_{T} (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ", " + absEtaStr2 + ")").c_str(), 50, 30, 130);
	centerTitles({jtPt_h[jI][cI][eI], jtRawPt_h[jI][cI][eI]});
      }

      for(Int_t pI = 0; pI < nJtPtBins; ++pI){
	std::string ptStr1 = "Pt" + prettyString(jtPtBinsLow[pI], 1, true) + "to" + prettyString(jtPtBinsHi[pI], 1, true);
	std::string ptStr2 = prettyString(jtPtBinsLow[pI], 1, true) + " |#eta_{Jet}| " + prettyString(jtPtBinsHi[pI], 1, false);

	jtPhi_h[jI][cI][pI] = new TH1F(("jtPhi_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_" + ptStr1 + "_h").c_str(), (";Jet #phi (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ", " + ptStr2 + ")").c_str(), 50, -TMath::Pi(), TMath::Pi());
	jtEta_h[jI][cI][pI] = new TH1F(("jtEta_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_" + ptStr1 + "_h").c_str(), (";Jet #eta (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ", " + ptStr2 + ")").c_str(), 50, -2., 2.);
	centerTitles({jtPhi_h[jI][cI][pI], jtEta_h[jI][cI][pI]});
      }

      jtPtOverRawPt_h[jI][cI] = new TH2F(("jtPtOverRawPt_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1  +"_h").c_str(), (";Jet p_{T};Jet p_{T}/Raw p_{T} (" + centStr2 + ")").c_str(), 50, 30, 130, 200, 1., 1.25);
      
      jtCHF_h[jI][cI] = new TH1F(("jtCHF_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_h").c_str(), (";Jet CHF (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ")").c_str(), 55, 0, 1.1);
      jtNHF_h[jI][cI] = new TH1F(("jtNHF_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_h").c_str(), (";Jet NHF (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ")").c_str(), 55, 0, 1.1);
      jtCEF_h[jI][cI] = new TH1F(("jtCEF_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_h").c_str(), (";Jet CEF (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ")").c_str(), 55, 0, 1.1);
      jtNEF_h[jI][cI] = new TH1F(("jtNEF_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_h").c_str(), (";Jet NEF (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ")").c_str(), 55, 0, 1.1);
      jtMUF_h[jI][cI] = new TH1F(("jtMUF_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_h").c_str(), (";Jet MUF (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ")").c_str(), 55, 0, 1.1);
      jtTF_h[jI][cI] = new TH1F(("jtTF_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_h").c_str(), (";Jet TF (" + jetTreeAlgGlobal.at(jI) + ");Counts (" + centStr2 + ")").c_str(), 55, 0, 1.1);
      
      centerTitles({jtPtOverRawPt_h[jI][cI], jtCHF_h[jI][cI], jtNHF_h[jI][cI], jtCEF_h[jI][cI], jtNEF_h[jI][cI], jtMUF_h[jI][cI], jtTF_h[jI][cI]});
    }
  }


  for(Int_t jI = 0; jI < nJtEtaBins; ++jI){
    randomConeSumVsHiBin_h[jI] = 0;
    randomConeRhoVsHiBin_h[jI] = 0;

    if(hasPFTreeGlobal){
      std::string jetEtaStr = "JetEta" + prettyString(jtEtaBinsLow[jI],2,true) + "to" + prettyString(jtEtaBinsHi[jI],2,true);
      while(jetEtaStr.find("-") != std::string::npos){jetEtaStr.replace(jetEtaStr.find("-"), 1, "Neg");}

      randomConeSumVsHiBin_h[jI] = new TH2F(("randomConeSumVsHiBin_" + jetEtaStr + "_h").c_str(), (";hiBin;RC_{Sum} (" + prettyString(jtEtaBinsLow[jI],2,false) + " < #eta_{Cone} < " + prettyString(jtEtaBinsHi[jI],2,false) + ")").c_str(), 200, -0.5, 199.5, 200, 0, 200);

      randomConeRhoVsHiBin_h[jI] = new TH2F(("randomConeRhoVsHiBin_" + jetEtaStr + "_h").c_str(), (";hiBin;#rho (" + prettyString(jtEtaBinsLow[jI],2,false) + " < #eta_{Cone} < " + prettyString(jtEtaBinsHi[jI],2,false) + ")").c_str(), 200, -0.5, 199.5, 200, 0, 400);

      centerTitles({randomConeSumVsHiBin_h[jI], randomConeRhoVsHiBin_h[jI]});
    }
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){ 
    const std::string centStr1 = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";

    for(Int_t pI = 0; pI < nPFID; ++pI){

      pfFracVEta_h[cI][pI] = 0;

      for(Int_t jI = 0; jI < nJetTrees; ++jI){
	pfFracVEta_InJet_h[jI][cI][pI] = 0;
      }

      if(hasPFTreeGlobal){
	pfFracVEta_h[cI][pI] = new TH1F(("pfFracVEta_" + centStr1 + "_PFID" + std::to_string(pI + 1) + "_h").c_str(), (";#eta;Fraction ID=" + std::to_string(pI+1) + " (" + centStr2 + ")").c_str(), nPFEtaBins, pfEtaBins);
	centerTitles({pfFracVEta_h[cI][pI]});

	for(Int_t jI = 0; jI < nJetTrees; ++jI){
	  pfFracVEta_InJet_h[jI][cI][pI] = new TH1F(("pfFracVEta_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_PFID" + std::to_string(pI + 1) + "_InJet_h").c_str(), (";#eta_{" + jetTreeAlgGlobal.at(jI) + "};Fraction ID=" + std::to_string(pI+1) + " (" + centStr2 + ")").c_str(), nPFEtaBins, pfEtaBins);
	  centerTitles({pfFracVEta_InJet_h[jI][cI][pI]});
	}
      }

      for(Int_t eI = 0; eI < nPFEtaBins; ++eI){
	std::string etaStr1 = "Eta" + prettyString(pfEtaBins[eI],1,true) + "to" + prettyString(pfEtaBins[eI+1],1,true);
	std::string etaStr2 = prettyString(pfEtaBins[eI],1,false) + " < #eta < " + prettyString(pfEtaBins[eI+1],1,false);

	while(etaStr1.find("-") != std::string::npos){etaStr1.replace(etaStr1.find("-"), 1, "Neg");}

	pfFracVEta_Bins_h[cI][pI][eI] = 0;

	for(Int_t jI = 0; jI < nJetTrees; ++jI){
	  pfFracVEta_Bins_InJet_h[jI][cI][pI][eI] = 0;
	}

	if(hasPFTreeGlobal){
	  pfFracVEta_Bins_h[cI][pI][eI] = new TH1F(("pfFracVEta_" + centStr1 + "_PFID" + std::to_string(pI + 1) + "_" + etaStr1 + "_h").c_str(), (";Fraction ID=" + std::to_string(pI+1) + " (" + centStr2 + ", " + etaStr2 + ");Counts").c_str(), 101, -0.005, 1.005);

	  for(Int_t jI = 0; jI < nJetTrees; ++jI){
	    pfFracVEta_Bins_InJet_h[jI][cI][pI][eI] = new TH1F(("pfFracVEta_" + jetTreeAlgGlobal.at(jI) + "_" + centStr1 + "_PFID" + std::to_string(pI + 1) + "_" + etaStr1 + "_InJet_h").c_str(), (";Fraction ID=" + std::to_string(pI+1) + " (" + centStr2 + ", " + etaStr2 + ");Counts").c_str(), 101, -0.005, 1.005);
	  }
	}
      }
    }
  }

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");
    
    TTree* pfTree_p=0;
    if(hasPFTreeGlobal) pfTree_p = (TTree*)inFile_p->Get("pfcandAnalyzer/pfTree");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* jetTrees_p[nJetTrees];
    for(Int_t jI = 0; jI < nJetTrees; ++jI){
      jetTrees_p[jI] = (TTree*)inFile_p->Get(jetTreeListGlobal.at(jI).c_str());
    }

    
    UInt_t run_,lumi_;
    ULong64_t evt_;
    Int_t hiBin_;
    Float_t vz_;
    
    std::vector<float>* pfPt_p=0;
    std::vector<float>* pfPhi_p=0;
    std::vector<float>* pfEta_p=0;
    std::vector<int>* pfId_p=0;
    
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("vz", 1);

    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    
    if(hasPFTreeGlobal){
      pfTree_p->SetBranchStatus("*", 0);
      pfTree_p->SetBranchStatus("pfPt", 1);
      pfTree_p->SetBranchStatus("pfPhi", 1);
      pfTree_p->SetBranchStatus("pfEta", 1);
      pfTree_p->SetBranchStatus("pfId", 1);
      
      pfTree_p->SetBranchAddress("pfPt", &pfPt_p);
      pfTree_p->SetBranchAddress("pfPhi", &pfPhi_p);
      pfTree_p->SetBranchAddress("pfEta", &pfEta_p);
      pfTree_p->SetBranchAddress("pfId", &pfId_p);
    }
    
    const Int_t nJetMax = 500;
    Int_t nref_[nJetTrees];
    Float_t jtpt_[nJetTrees][nJetMax];
    Float_t rawpt_[nJetTrees][nJetMax];
    Float_t jtphi_[nJetTrees][nJetMax];
    Float_t jteta_[nJetTrees][nJetMax];
    
    //PU specific vars;
    Float_t chargedSum_[nJetTrees][nJetMax];
    Float_t photonSum_[nJetTrees][nJetMax];
    Float_t neutralSum_[nJetTrees][nJetMax];
    Float_t eSum_[nJetTrees][nJetMax];
    Float_t muSum_[nJetTrees][nJetMax];
    
    //CS specific vars;  
    Float_t jtPfCHF_[nJetTrees][nJetMax];
    Float_t jtPfNHF_[nJetTrees][nJetMax];
    Float_t jtPfCEF_[nJetTrees][nJetMax];
    Float_t jtPfNEF_[nJetTrees][nJetMax];
    Float_t jtPfMUF_[nJetTrees][nJetMax];

    for(Int_t jI = 0; jI < nJetTrees; ++jI){
      jetTrees_p[jI]->SetBranchStatus("*", 0);
      jetTrees_p[jI]->SetBranchStatus("nref", 1);
      jetTrees_p[jI]->SetBranchStatus("jtpt", 1);
      jetTrees_p[jI]->SetBranchStatus("rawpt", 1);
      jetTrees_p[jI]->SetBranchStatus("jtphi", 1);
      jetTrees_p[jI]->SetBranchStatus("jteta", 1);
      if(jetTreeAlgGlobal.at(jI).find("akPu") != std::string::npos){
	jetTrees_p[jI]->SetBranchStatus("chargedSum", 1);
	jetTrees_p[jI]->SetBranchStatus("neutralSum", 1);
	jetTrees_p[jI]->SetBranchStatus("photonSum", 1);
	jetTrees_p[jI]->SetBranchStatus("eSum", 1);
	jetTrees_p[jI]->SetBranchStatus("muSum", 1);
      }
      else if(jetTreeAlgGlobal.at(jI).find("akCs") != std::string::npos){
	jetTrees_p[jI]->SetBranchStatus("jtPfCHF", 1);
	jetTrees_p[jI]->SetBranchStatus("jtPfNHF", 1);
	jetTrees_p[jI]->SetBranchStatus("jtPfCEF", 1);
	jetTrees_p[jI]->SetBranchStatus("jtPfNEF", 1);
	jetTrees_p[jI]->SetBranchStatus("jtPfMUF", 1);
      }
      
      jetTrees_p[jI]->SetBranchAddress("nref", (&nref_[jI]));
      jetTrees_p[jI]->SetBranchAddress("jtpt", jtpt_[jI]);
      jetTrees_p[jI]->SetBranchAddress("rawpt", rawpt_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtphi", jtphi_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jteta", jteta_[jI]);
      if(jetTreeAlgGlobal.at(jI).find("akPu") != std::string::npos){
	jetTrees_p[jI]->SetBranchAddress("chargedSum", chargedSum_[jI]);
	jetTrees_p[jI]->SetBranchAddress("neutralSum", neutralSum_[jI]);
	jetTrees_p[jI]->SetBranchAddress("photonSum", photonSum_[jI]);
	jetTrees_p[jI]->SetBranchAddress("eSum", eSum_[jI]);
	jetTrees_p[jI]->SetBranchAddress("muSum", muSum_[jI]);
      }
      else if(jetTreeAlgGlobal.at(jI).find("akCs") != std::string::npos){
	jetTrees_p[jI]->SetBranchAddress("jtPfCHF", jtPfCHF_[jI]);
	jetTrees_p[jI]->SetBranchAddress("jtPfNHF", jtPfNHF_[jI]);
	jetTrees_p[jI]->SetBranchAddress("jtPfCEF", jtPfCEF_[jI]);
	jetTrees_p[jI]->SetBranchAddress("jtPfNEF", jtPfNEF_[jI]);
	jetTrees_p[jI]->SetBranchAddress("jtPfMUF", jtPfMUF_[jI]);
      }
    }
    
    const Int_t nEntries = hiTree_p->GetEntries();
    
    std::cout << "Processing file " << fI << "/" << fileList.size() << ", \'" << fileList.at(fI) << "\'..." << std::endl;

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
      
      hiTree_p->GetEntry(entry);
      if(hasPFTreeGlobal) pfTree_p->GetEntry(entry);
      for(Int_t jI = 0; jI < nJetTrees; ++jI){
	jetTrees_p[jI]->GetEntry(entry);
      }
      
      if(TMath::Abs(vz_) > 15) continue;
      if(run_ == 304899 && lumi_ <= 32) continue;

      run_h->Fill(run_);
      lumi_h->Fill(lumi_);
      evt_h->Fill(evt_);
      vz_h->Fill(vz_);
      hiBin_h->Fill(hiBin_);
      
      std::vector<Int_t> centPos;
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHi[cI]){
	  centPos.push_back(cI);
	}
      }
 
      std::vector< std::vector<Float_t> > jetEtas;
      std::vector< std::vector<Float_t> > jetPhis;
     
      for(Int_t jI = 0; jI < nJetTrees; ++jI){
	std::vector<Float_t> jetEtasTemp;
	std::vector<Float_t> jetPhisTemp;

	for(Int_t sI = 0; sI < nref_[jI]; ++sI){
	  if(jtpt_[jI][sI] < 30.) continue;
	  
	  jetEtasTemp.push_back(jteta_[jI][sI]);
	  jetPhisTemp.push_back(jtphi_[jI][sI]);

	  if(TMath::Abs(jteta_[jI][sI]) > 2.) continue;

	  std::vector<Int_t> absEtaPos;
	  std::vector<Int_t> ptPos;

	  
	  for(Int_t eI = 0; eI < nJtAbsEtaBins; ++eI){
	    if(TMath::Abs(jteta_[jI][sI]) >= jtAbsEtaBinsLow[eI] && TMath::Abs(jteta_[jI][sI]) < jtAbsEtaBinsHi[eI]){
	      absEtaPos.push_back(eI);
	    }
	  }
	  
	  for(Int_t eI = 0; eI < nJtPtBins; ++eI){
	    if(jtpt_[jI][sI] >= jtPtBinsLow[eI] && jtpt_[jI][sI] < jtPtBinsHi[eI]){
	      ptPos.push_back(eI);
	    }
	  }
	  
	  for(unsigned int cI = 0; cI < centPos.size(); ++cI){

	    for(unsigned int eI = 0; eI < absEtaPos.size(); ++eI){
	      jtPt_h[jI][centPos.at(cI)][absEtaPos.at(eI)]->Fill(jtpt_[jI][sI]);
	      jtRawPt_h[jI][centPos.at(cI)][absEtaPos.at(eI)]->Fill(rawpt_[jI][sI]);
	    }
	    
	    for(unsigned int eI = 0; eI < ptPos.size(); ++eI){
	      jtPhi_h[jI][centPos.at(cI)][ptPos.at(eI)]->Fill(jtphi_[jI][sI]);
	      jtEta_h[jI][centPos.at(cI)][ptPos.at(eI)]->Fill(jteta_[jI][sI]);
	    }
	    
	    jtPtOverRawPt_h[jI][centPos.at(cI)]->Fill(jtpt_[jI][sI], jtpt_[jI][sI]/rawpt_[jI][sI]);
	    
	    if(jetTreeAlgGlobal.at(jI).find("akPu") != std::string::npos){
	      Float_t total = chargedSum_[jI][sI] + photonSum_[jI][sI] + neutralSum_[jI][sI] + eSum_[jI][sI] + muSum_[jI][sI];
	   
	      jtCHF_h[jI][centPos.at(cI)]->Fill(chargedSum_[jI][sI]/total);
	      jtNHF_h[jI][centPos.at(cI)]->Fill(neutralSum_[jI][sI]/total);
	      jtCEF_h[jI][centPos.at(cI)]->Fill(eSum_[jI][sI]/total);
	      jtNEF_h[jI][centPos.at(cI)]->Fill(photonSum_[jI][sI]/total);
	      jtMUF_h[jI][centPos.at(cI)]->Fill(muSum_[jI][sI]/total);
	      jtTF_h[jI][centPos.at(cI)]->Fill(rawpt_[jI][sI]/total);
	    }
	    else if(jetTreeAlgGlobal.at(jI).find("akCs") != std::string::npos){
	      jtCHF_h[jI][centPos.at(cI)]->Fill(jtPfCHF_[jI][sI]);
	      jtNHF_h[jI][centPos.at(cI)]->Fill(jtPfNHF_[jI][sI]);
	      jtCEF_h[jI][centPos.at(cI)]->Fill(jtPfCEF_[jI][sI]);
	      jtNEF_h[jI][centPos.at(cI)]->Fill(jtPfNEF_[jI][sI]);
	      jtMUF_h[jI][centPos.at(cI)]->Fill(jtPfMUF_[jI][sI]);
	      jtTF_h[jI][centPos.at(cI)]->Fill(jtPfCHF_[jI][sI] + jtPfNHF_[jI][sI] + jtPfCEF_[jI][sI] + jtPfNEF_[jI][sI] + jtPfMUF_[jI][sI]);
	    }
	  }
	}

	jetEtas.push_back(jetEtasTemp);
	jetPhis.push_back(jetPhisTemp);
      }

      
      if(hasPFTreeGlobal){
	Float_t pfEtaSum[nPFID+1][nPFEtaBins];
	Float_t pfEtaSumJet[nJetTrees][nPFID+1][nPFEtaBins];

	for(Int_t iI = 0; iI < nPFID+1; ++iI){
	  for(Int_t eI = 0; eI < nPFEtaBins; ++eI){
	    pfEtaSum[iI][eI] = 0.0;

	    for(Int_t jI = 0; jI < nJetTrees; ++jI){
	      pfEtaSumJet[jI][iI][eI] = 0.0;
	    }
	  }
	}

	std::vector<Float_t> etaRC;
	std::vector<Float_t> phiRC;
	std::vector<Float_t> sumRC;
	
	while((Int_t)etaRC.size() < nRC){
	  Float_t etaVal = randGen_p->Uniform(jtEtaBinsLow[0], jtEtaBinsHi[nJtEtaBins-1]);
	  Float_t phiVal = randGen_p->Uniform(-TMath::Pi(), TMath::Pi());
	  
	  bool isGood = true;
	  for(unsigned int rI = 0; rI < etaRC.size(); ++rI){
	    if(getDR(etaVal, phiVal, etaRC.at(rI), phiRC.at(rI)) < rcParam){
	      isGood = false;
	      break;
	    }
	  }
	  
	  if(!isGood) continue;
	  
	  etaRC.push_back(etaVal);
	  phiRC.push_back(phiVal);
	  sumRC.push_back(0.0);
	}
	
	for(unsigned int pI = 0; pI < pfPt_p->size(); ++pI){
	  for(Int_t eI = 0; eI < nPFEtaBins; ++eI){
	    if(pfEta_p->at(pI) >= pfEtaBins[eI] && pfEta_p->at(pI) < pfEtaBins[eI+1]){
	      pfEtaSum[0][eI] += pfPt_p->at(pI);
	      pfEtaSum[pfId_p->at(pI)][eI] += pfPt_p->at(pI);

	      for(unsigned int jI = 0; jI < jetEtas.size(); ++jI){
		for(unsigned int sI = 0; sI < jetEtas.at(jI).size(); ++sI){
		  if(getDR(jetEtas.at(jI).at(sI), jetPhis.at(jI).at(sI), pfEta_p->at(pI), pfPhi_p->at(pI)) < rcParam){
		    pfEtaSumJet[jI][0][eI] += pfPt_p->at(pI);
		    pfEtaSumJet[jI][pfId_p->at(pI)][eI] += pfPt_p->at(pI);
		    break;
		  }
		}
	      }

	      break;
	    }
	  }

	  for(unsigned int rI = 0; rI < etaRC.size(); ++rI){
	    if(getDR(etaRC.at(rI), phiRC.at(rI), pfEta_p->at(pI), pfPhi_p->at(pI)) < rcParam){
	      sumRC.at(rI) += pfPt_p->at(pI);
	    }
	  }
	}

	
	for(unsigned int cI = 0; cI < centPos.size(); ++cI){
	  for(Int_t iI = 0; iI < nPFID; ++iI){
	    for(Int_t eI = 0; eI < nPFEtaBins; ++eI){
	      if(pfEtaSum[0][eI] > .01) pfFracVEta_Bins_h[centPos.at(cI)][iI][eI]->Fill(pfEtaSum[iI + 1][eI]/pfEtaSum[0][eI]);
	      for(Int_t jI = 0; jI < nJetTrees; ++jI){
		if(pfEtaSumJet[jI][0][eI] > 0.01) pfFracVEta_Bins_InJet_h[jI][centPos.at(cI)][iI][eI]->Fill(pfEtaSumJet[jI][iI + 1][eI]/pfEtaSumJet[jI][0][eI]);
	      }
	    }
	  }	  
	}

	unsigned int rcPos = 0;
	while(rcPos < sumRC.size()){
	  bool isGood = true;
	  
	  for(unsigned int rI = rcPos + 1; rI < sumRC.size(); ++rI){
	    if(sumRC.at(rcPos) < sumRC.at(rI)){
	      Float_t tempEta = etaRC.at(rI);
	      Float_t tempPhi = phiRC.at(rI);
	      Float_t tempSum = sumRC.at(rI);
	      
	      etaRC.at(rI) =  etaRC.at(rcPos);
	      phiRC.at(rI) =  phiRC.at(rcPos);
	      sumRC.at(rI) =  sumRC.at(rcPos);
	      
	      etaRC.at(rcPos) = tempEta;
	      phiRC.at(rcPos) = tempPhi;
	      sumRC.at(rcPos) = tempSum;
	      
	      isGood = false;
	      break;
	    }
	  }
	  
	  if(isGood) rcPos++;
	}

	std::vector<Int_t> rcEtaPos;
	
	for(Int_t jI = 0; jI < nJtEtaBins; ++jI){
	  if(etaRC.at(etaRC.size()/2) >= jtEtaBinsLow[jI] && etaRC.at(etaRC.size()/2) < jtEtaBinsHi[jI]) rcEtaPos.push_back(jI);
	}
	
	if(hiBin_ < 20 && sumRC.at(sumRC.size()/2) < 30){
	  std::cout << "hiBin, sumRC, etaRC, phiRC: " << hiBin_ << ", " << sumRC.at(sumRC.size()/2) << ", " << etaRC.at(etaRC.size()/2) << ", " << phiRC.at(phiRC.size()/2) << std::endl;
	  std::cout << " " << run_ << ", " << lumi_ << std::endl;
	}
	
	for(unsigned int rI = 0; rI < rcEtaPos.size(); ++rI){
	  randomConeSumVsHiBin_h[rcEtaPos.at(rI)]->Fill(hiBin_, sumRC.at(sumRC.size()/2));
	  randomConeRhoVsHiBin_h[rcEtaPos.at(rI)]->Fill(hiBin_, sumRC.at(sumRC.size()/2)/(rcParam*rcParam*TMath::Pi()));
	}
      }
    }
  
    inFile_p->Close();
    delete inFile_p;
  }

  outFile_p->cd();

  run_h->Write("", TObject::kOverwrite);
  delete run_h;

  lumi_h->Write("", TObject::kOverwrite);
  delete lumi_h;

  evt_h->Write("", TObject::kOverwrite);
  delete evt_h;

  vz_h->Write("", TObject::kOverwrite);
  delete vz_h;

  hiBin_h->Write("", TObject::kOverwrite);
  delete hiBin_h;

  for(Int_t jI = 0; jI < nJetTrees; ++jI){
    outFile_p->cd();
    outFile_p->mkdir(jetTreeAlgGlobal.at(jI).c_str());
    outFile_p->cd(jetTreeAlgGlobal.at(jI).c_str());

    for(Int_t cI = 0; cI < nCentBins; ++cI){      
      const std::string centStr1 = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      outFile_p->cd();
      outFile_p->mkdir((jetTreeAlgGlobal.at(jI) + "/" + centStr1).c_str());
      outFile_p->cd((jetTreeAlgGlobal.at(jI) + "/" + centStr1).c_str());

      for(Int_t eI = 0; eI < nJtAbsEtaBins; ++eI){
	jtPt_h[jI][cI][eI]->Write("", TObject::kOverwrite);
	jtRawPt_h[jI][cI][eI]->Write("", TObject::kOverwrite);
      }

      for(Int_t eI = 0; eI < nJtPtBins; ++eI){
	jtPhi_h[jI][cI][eI]->Write("", TObject::kOverwrite);
	jtEta_h[jI][cI][eI]->Write("", TObject::kOverwrite);
      }

      jtPtOverRawPt_h[jI][cI]->Write("", TObject::kOverwrite);
      jtCHF_h[jI][cI]->Write("", TObject::kOverwrite);
      jtNHF_h[jI][cI]->Write("", TObject::kOverwrite);
      jtCEF_h[jI][cI]->Write("", TObject::kOverwrite);
      jtNEF_h[jI][cI]->Write("", TObject::kOverwrite);
      jtMUF_h[jI][cI]->Write("", TObject::kOverwrite);
      jtTF_h[jI][cI]->Write("", TObject::kOverwrite);
      
      for(Int_t eI = 0; eI < nJtAbsEtaBins; ++eI){
	delete jtPt_h[jI][cI][eI];
	delete jtRawPt_h[jI][cI][eI];
      }

      for(Int_t eI = 0; eI < nJtPtBins; ++eI){
	delete jtPhi_h[jI][cI][eI];
	delete jtEta_h[jI][cI][eI];
      }

      delete jtPtOverRawPt_h[jI][cI];
      delete jtCHF_h[jI][cI];
      delete jtNHF_h[jI][cI];
      delete jtCEF_h[jI][cI];
      delete jtNEF_h[jI][cI];
      delete jtMUF_h[jI][cI];
      delete jtTF_h[jI][cI];
    }
  }

  if(hasPFTreeGlobal){
    outFile_p->cd();
    outFile_p->mkdir("pfcandDir");
    outFile_p->cd("pfcandDir");

    for(Int_t rI = 0; rI < nJtEtaBins; ++rI){
      randomConeSumVsHiBin_h[rI]->Write("", TObject::kOverwrite);
      randomConeRhoVsHiBin_h[rI]->Write("", TObject::kOverwrite);

      delete randomConeSumVsHiBin_h[rI];
      delete randomConeRhoVsHiBin_h[rI];
    }

    for(Int_t cI = 0; cI < nCentBins; ++cI){ 
      const std::string centStr1 = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      outFile_p->cd();
      outFile_p->mkdir(("pfcandDir/" + centStr1).c_str());
      outFile_p->cd(("pfcandDir/" + centStr1).c_str());

      for(Int_t pI = 0; pI < nPFID; ++pI){
	for(Int_t eI = 0; eI < pfFracVEta_h[cI][pI]->GetNbinsX(); ++eI){
	  pfFracVEta_h[cI][pI]->SetBinContent(eI+1, pfFracVEta_Bins_h[cI][pI][eI]->GetMean());
	  pfFracVEta_h[cI][pI]->SetBinError(eI+1, pfFracVEta_Bins_h[cI][pI][eI]->GetMeanError());
	}

	pfFracVEta_h[cI][pI]->Write("", TObject::kOverwrite);
	delete pfFracVEta_h[cI][pI];
	
	for(Int_t jI = 0; jI < nJetTrees; ++jI){
	  for(Int_t eI = 0; eI < pfFracVEta_InJet_h[jI][cI][pI]->GetNbinsX(); ++eI){
	    pfFracVEta_InJet_h[jI][cI][pI]->SetBinContent(eI+1, pfFracVEta_Bins_InJet_h[jI][cI][pI][eI]->GetMean());
	    pfFracVEta_InJet_h[jI][cI][pI]->SetBinError(eI+1, pfFracVEta_Bins_InJet_h[jI][cI][pI][eI]->GetMeanError());
	  }

	  pfFracVEta_InJet_h[jI][cI][pI]->Write("", TObject::kOverwrite);
	  delete pfFracVEta_InJet_h[jI][cI][pI];
	}
      }
    }

    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr1 = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      outFile_p->cd();
      outFile_p->mkdir(("pfcandDir/" + centStr1 + "/Bins").c_str());
      outFile_p->cd(("pfcandDir/" + centStr1 + "/Bins").c_str());

      for(Int_t pI = 0; pI < nPFID; ++pI){
	for(Int_t eI = 0; eI < nPFEtaBins; ++eI){
	  pfFracVEta_Bins_h[cI][pI][eI]->Write("", TObject::kOverwrite);
	  delete pfFracVEta_Bins_h[cI][pI][eI];
	}

	for(Int_t jI = 0; jI < nJetTrees; ++jI){
	  for(Int_t eI = 0; eI < nPFEtaBins; ++eI){
	    pfFracVEta_Bins_InJet_h[jI][cI][pI][eI]->Write("", TObject::kOverwrite);
	    delete pfFracVEta_Bins_InJet_h[jI][cI][pI][eI];
	  }
	}
      }
    }
  }

  outFile_p->cd();
  outFile_p->mkdir("binDir");
  outFile_p->cd("binDir");

  std::string jetAlgStr = "";
  for(Int_t jI = 0; jI < nJetTrees; ++jI){
    if(jI < nJetTrees-1) jetAlgStr = jetAlgStr + jetTreeAlgGlobal.at(jI) + ", ";
    else jetAlgStr = jetAlgStr + jetTreeAlgGlobal.at(jI);
  }
  
  TNamed jetAlgName("jetAlgos", jetAlgStr.c_str());
  jetAlgName.Write("", TObject::kOverwrite);

  TNamed nCentBinsName("nCentBins", std::to_string(nCentBins).c_str());
  nCentBinsName.Write("", TObject::kOverwrite);

  std::string centBinsLowStr = "";
  std::string centBinsHiStr = "";

  for(unsigned int i = 0; i < nCentBins; ++i){centBinsLowStr = centBinsLowStr + "," + prettyString(centBinsLow[i], 0, false);}
  for(unsigned int i = 0; i < nCentBins; ++i){centBinsHiStr = centBinsHiStr + "," + prettyString(centBinsHi[i], 0, false);}
  if(centBinsLowStr.substr(0,1).find(",") != std::string::npos) centBinsLowStr.replace(0,1,"");
  if(centBinsLowStr.substr(centBinsLowStr.size()-1,1).find(",") != std::string::npos) centBinsLowStr.replace(centBinsLowStr.size()-1,1,"");
  if(centBinsHiStr.substr(0,1).find(",") != std::string::npos) centBinsHiStr.replace(0,1,"");
  if(centBinsHiStr.substr(centBinsHiStr.size()-1,1).find(",") != std::string::npos) centBinsHiStr.replace(centBinsHiStr.size()-1,1,"");
  
  TNamed centBinsLowName("centBinsLow", centBinsLowStr.c_str());
  centBinsLowName.Write("", TObject::kOverwrite);
  TNamed centBinsHiName("centBinsHi", centBinsHiStr.c_str());
  centBinsHiName.Write("", TObject::kOverwrite);


  TNamed nJtPtBinsName("nJtPtBins", std::to_string(nJtPtBins).c_str());
  nJtPtBinsName.Write("", TObject::kOverwrite);

  std::string jtPtBinsLowStr = "";
  std::string jtPtBinsHiStr = "";

  for(unsigned int i = 0; i < nJtPtBins; ++i){jtPtBinsLowStr = jtPtBinsLowStr + "," + prettyString(jtPtBinsLow[i], 1, false);}
  for(unsigned int i = 0; i < nJtPtBins; ++i){jtPtBinsHiStr = jtPtBinsHiStr + "," + prettyString(jtPtBinsHi[i], 1, false);}
  if(jtPtBinsLowStr.substr(0,1).find(",") != std::string::npos) jtPtBinsLowStr.replace(0,1,"");
  if(jtPtBinsLowStr.substr(jtPtBinsLowStr.size()-1,1).find(",") != std::string::npos) jtPtBinsLowStr.replace(jtPtBinsLowStr.size()-1,1,"");
  if(jtPtBinsHiStr.substr(0,1).find(",") != std::string::npos) jtPtBinsHiStr.replace(0,1,"");
  if(jtPtBinsHiStr.substr(jtPtBinsHiStr.size()-1,1).find(",") != std::string::npos) jtPtBinsHiStr.replace(jtPtBinsHiStr.size()-1,1,"");
  
  TNamed jtPtBinsLowName("jtPtBinsLow", jtPtBinsLowStr.c_str());
  jtPtBinsLowName.Write("", TObject::kOverwrite);
  TNamed jtPtBinsHiName("jtPtBinsHi", jtPtBinsHiStr.c_str());
  jtPtBinsHiName.Write("", TObject::kOverwrite);

  TNamed nJtEtaBinsName("nJtEtaBins", std::to_string(nJtEtaBins).c_str());
  nJtEtaBinsName.Write("", TObject::kOverwrite);

  std::string jtEtaBinsLowStr = "";
  std::string jtEtaBinsHiStr = "";

  for(unsigned int i = 0; i < nJtEtaBins; ++i){jtEtaBinsLowStr = jtEtaBinsLowStr + "," + prettyString(jtEtaBinsLow[i], 2, false);}
  for(unsigned int i = 0; i < nJtEtaBins; ++i){jtEtaBinsHiStr = jtEtaBinsHiStr + "," + prettyString(jtEtaBinsHi[i], 2, false);}
  if(jtEtaBinsLowStr.substr(0,1).find(",") != std::string::npos) jtEtaBinsLowStr.replace(0,1,"");
  if(jtEtaBinsLowStr.substr(jtEtaBinsLowStr.size()-1,1).find(",") != std::string::npos) jtEtaBinsLowStr.replace(jtEtaBinsLowStr.size()-1,1,"");
  if(jtEtaBinsHiStr.substr(0,1).find(",") != std::string::npos) jtEtaBinsHiStr.replace(0,1,"");
  if(jtEtaBinsHiStr.substr(jtEtaBinsHiStr.size()-1,1).find(",") != std::string::npos) jtEtaBinsHiStr.replace(jtEtaBinsHiStr.size()-1,1,"");
  
  TNamed jtEtaBinsLowName("jtEtaBinsLow", jtEtaBinsLowStr.c_str());
  jtEtaBinsLowName.Write("", TObject::kOverwrite);
  TNamed jtEtaBinsHiName("jtEtaBinsHi", jtEtaBinsHiStr.c_str());
  jtEtaBinsHiName.Write("", TObject::kOverwrite);

  TNamed nJtAbsEtaBinsName("nJtAbsEtaBins", std::to_string(nJtAbsEtaBins).c_str());
  nJtAbsEtaBinsName.Write("", TObject::kOverwrite);

  std::string jtAbsEtaBinsLowStr = "";
  std::string jtAbsEtaBinsHiStr = "";

  for(unsigned int i = 0; i < nJtAbsEtaBins; ++i){jtAbsEtaBinsLowStr = jtAbsEtaBinsLowStr + "," + prettyString(jtAbsEtaBinsLow[i], 2, false);}
  for(unsigned int i = 0; i < nJtAbsEtaBins; ++i){jtAbsEtaBinsHiStr = jtAbsEtaBinsHiStr + "," + prettyString(jtAbsEtaBinsHi[i], 2, false);}
  if(jtAbsEtaBinsLowStr.substr(0,1).find(",") != std::string::npos) jtAbsEtaBinsLowStr.replace(0,1,"");
  if(jtAbsEtaBinsLowStr.substr(jtAbsEtaBinsLowStr.size()-1,1).find(",") != std::string::npos) jtAbsEtaBinsLowStr.replace(jtAbsEtaBinsLowStr.size()-1,1,"");
  if(jtAbsEtaBinsHiStr.substr(0,1).find(",") != std::string::npos) jtAbsEtaBinsHiStr.replace(0,1,"");
  if(jtAbsEtaBinsHiStr.substr(jtAbsEtaBinsHiStr.size()-1,1).find(",") != std::string::npos) jtAbsEtaBinsHiStr.replace(jtAbsEtaBinsHiStr.size()-1,1,"");

  TNamed jtAbsEtaBinsLowName("jtAbsEtaBinsLow", jtAbsEtaBinsLowStr.c_str());
  jtAbsEtaBinsLowName.Write("", TObject::kOverwrite);
  TNamed jtAbsEtaBinsHiName("jtAbsEtaBinsHi", jtAbsEtaBinsHiStr.c_str());
  jtAbsEtaBinsHiName.Write("", TObject::kOverwrite);

  TNamed nPFEtaBinsName("nPFEtaBins", std::to_string(nPFEtaBins).c_str());
  nPFEtaBinsName.Write("", TObject::kOverwrite);

  std::string pfEtaBinsStr = "";

  for(unsigned int i = 0; i < nPFEtaBins+1; ++i){pfEtaBinsStr = pfEtaBinsStr + "," + prettyString(pfEtaBins[i], 2, false);}
  if(pfEtaBinsStr.substr(0,1).find(",") != std::string::npos) pfEtaBinsStr.replace(0,1,"");
  if(pfEtaBinsStr.substr(pfEtaBinsStr.size()-1,1).find(",") != std::string::npos) pfEtaBinsStr.replace(pfEtaBinsStr.size()-1,1,"");

  TNamed pfEtaBinsName("pfEtaBins", pfEtaBinsStr.c_str());
  pfEtaBinsName.Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./forestJetValidation.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += forestJetValidation(argv[1]);
  else if(argc == 3) retVal += forestJetValidation(argv[1], argv[2]);
  return retVal;
}
