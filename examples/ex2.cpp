#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "LoveNumbers/DeckModel.h"
using Real = double;

int main() {

  auto fileName = std::string("../data/prem.200.no");

  auto model = LoveNumbers::DeckModel<Real>(fileName);

  for (auto [r1, r2] : model.LayerRadii())
    std::cout << r1 << " " << r2 << std::endl;

  for (auto solid : model.LayerFluid())
    std::cout << solid << std::endl;
}