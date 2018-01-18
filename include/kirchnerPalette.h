#ifndef KIRCHNERPALETTE_H
#define KIRCHNERPALETTE_H

#include <vector>

#include "TColor.h"

#include "include/doGlobalDebug.h"

class kirchnerPalette{
 public:
  kirchnerPalette();
  ~kirchnerPalette(){};

  std::vector<Int_t> kirchColors;

  Int_t getColor(unsigned int colPos);
};

kirchnerPalette::kirchnerPalette()
{
  TColor kirchCol;
  kirchColors.push_back(kirchCol.GetColor(171, 56, 29)); //RED
  kirchColors.push_back(kirchCol.GetColor(79, 120, 182)); //BLUE
  kirchColors.push_back(kirchCol.GetColor(205, 209, 114)); //YELLOW
  kirchColors.push_back(kirchCol.GetColor(200, 93, 41)); //ORANGE
  kirchColors.push_back(kirchCol.GetColor(100, 128, 43)); //GREEN
  kirchColors.push_back(kirchCol.GetColor(208, 128, 139)); //PINK
  kirchColors.push_back(kirchCol.GetColor(120, 179, 173)); //CYAN

  return;
}


Int_t kirchnerPalette::getColor(unsigned int colPos){return kirchColors.at(colPos);}


#endif
