
#include "mfem.hpp"

#include "mfem.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>

#include "LoveNumbers/DeckModel.hpp"
#include "LoveNumbers/RadialModel.hpp"

int main(int argc, char *argv[]) {

  auto model = LoveNumbers::DeckModel("../data/prem.200.no", 0.05);

  model.PrintMesh("prem.mesh");
}
