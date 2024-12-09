
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

#include "LoveNumbers/LoveNumbers.hpp"
// #include "LoveNumbers/LoveNumbers.hpp"
// #include "LoveNumbers/RadialModel.hpp"

int main(int argc, char *argv[]) {

  // auto model = LoveNumbers::DeckModel::FromMaximumElementSize(
  //     "../data/prem.200.no", 0.1);

  auto model =
      LoveNumbers::DeckModel::FromMaximumDegree("../data/prem.200.no", 128);

  // auto model = LoveNumbers::HomogeneousModel::FromIsotropicWaveSpeeds(
  //     6.371e6, 5000, 10000, 8000, 100, 100);

  model.WriteAsDeckModel("prem.200", 1.e4);

  model.SetFiniteElementSpace(2);

  auto &FESpace = model.FiniteElementSpace();

  std::cout << FESpace.GetNDofs() << std::endl;

  //  model.PrintMesh("prem.mesh");
}
