#ifndef HISTDEFUTILITY_H
#define HISTDEFUTILITY_H

#include <vector>

#include "TH1.h"

void centerTitles(TH1* hist_p)
{
  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();
  
  return;
}

void centerTitles(std::vector<TH1*> hists_)
{
  for(unsigned int pI = 0; pI < hists_.size(); ++pI){
    centerTitles(hists_.at(pI));
  }
  
  return;
}

#endif
