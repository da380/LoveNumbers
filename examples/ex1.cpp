

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numbers>
#include <ranges>
#include <string>
#include <vector>

#include "LoveNumbers/LoveNumbers.hpp"

using namespace LoveNumbers;
using namespace mfem;

int main(int argc, char *argv[]) {

  auto dimensions = Dimensions();
  Int order = 3;
  Real characteristicLengthScale = 0.5;

  auto model = DeckModel(Dimensions(), "../data/prem.200.no");

  auto modelMesh = RadialModelMesh(model, order, characteristicLengthScale);

  auto x = GridFunction(&modelMesh.H1Space());
  auto phi = modelMesh.GravitationalPotentialCoefficient();
  x.ProjectCoefficient(phi);

  modelMesh.WriteDerivative(modelMesh.GravitationalPotentialCoefficient(),
                            "ex1.out", model.PotentialScale());
}
