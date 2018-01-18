#ifndef RETURNROOTFILECONTENTSLIST_H
#define RETURNROOTFILECONTENTSLIST_H

#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TDirectory.h"

#include <string>
#include <vector>

const std::string tdirStr = "TDirectoryFile";

std::vector<std::string> returnTDirContentsList(TFile*, const std::string, const std::string, const std::string);

std::vector<std::string> returnRootFileContentsList(TFile *inFile_p, const std::string classFilter = "", const std::string nameFilter = "")
{
  inFile_p->cd();

  std::vector<std::string> returnList;

  TList* fileContents_p = (TList*)inFile_p->GetListOfKeys();
  const Int_t nKeys = (Int_t)fileContents_p->GetEntries();
  
  for(Int_t keyIter = 0; keyIter < nKeys; keyIter++){
    TKey* tempKey = (TKey*)fileContents_p->At(keyIter);

    const std::string name = tempKey->GetName();
    const std::string className = tempKey->GetClassName();

    if(className.size() == tdirStr.size() && className.find(tdirStr) != std::string::npos){
      std::vector<std::string> tempReturnList = returnTDirContentsList(inFile_p, name, classFilter, nameFilter);

      const Int_t nTempKeys = tempReturnList.size();
      for(Int_t tempIter = 0; tempIter < nTempKeys; tempIter++){
	returnList.push_back(tempReturnList.at(tempIter));
      }
    }
  
    if(classFilter.size() != 0)
      if(className.size() != classFilter.size() || className.find(classFilter) == std::string::npos) continue;   
  
    if(nameFilter.size() != 0)
      if(name.find(nameFilter) == std::string::npos) continue;   
    
    returnList.push_back(name);
  }

  return returnList;
}


std::vector<std::string> returnTDirContentsList(TFile* inFile_p, const std::string dirName, const std::string classFilter, const std::string nameFilter)
{
  inFile_p->cd();

  std::vector<std::string> returnList;
  
  TDirectoryFile* dir_p = (TDirectoryFile*)inFile_p->Get(dirName.c_str());

  TList* fileContents_p = (TList*)dir_p->GetListOfKeys();
  const Int_t nKeys = (Int_t)fileContents_p->GetEntries();

  for(Int_t keyIter = 0; keyIter < nKeys; keyIter++){
    TKey* tempKey = (TKey*)fileContents_p->At(keyIter);

    const std::string name = tempKey->GetName();
    const std::string className = tempKey->GetClassName();

    if(className.size() == tdirStr.size() && className.find(tdirStr) != std::string::npos){
      std::vector<std::string> tempReturnList = returnTDirContentsList(inFile_p, dirName + "/" + name, classFilter, nameFilter);
      
      const Int_t nTempKeys = tempReturnList.size();
      for(Int_t tempIter = 0; tempIter < nTempKeys; tempIter++){
        returnList.push_back(tempReturnList.at(tempIter));
      }
    }

    if(classFilter.size() != 0)
      if(className.size() != classFilter.size() || className.find(classFilter) == std::string::npos) continue;
    
    if(nameFilter.size() != 0)
      if(name.find(nameFilter) == std::string::npos && dirName.find(nameFilter) == std::string::npos) continue;


    returnList.push_back(dirName + "/" + name);
  }  

  return returnList;
}


#endif
