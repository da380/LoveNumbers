

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
using namespace std::numbers;

int main(int argc, char* argv[]) {
  auto dimensions = Dimensions();
  Int order = 3;
  Real characteristicLengthScale = 0.01;

  auto model = DeckModel(Dimensions(), "../data/prem.200.no");
  //  auto model = DeckModel(Dimensions(), "../data/model1.txt");

  auto modelMesh = RadialModelMesh(model, order, characteristicLengthScale);

  auto x = GridFunction(&modelMesh.H1Space());
  auto phi = modelMesh.GravitationalPotentialCoefficient();
  // x.ProjectCoefficient(phi);

  modelMesh.Write(modelMesh.GravitationalPotentialCoefficient(),
                  "PREMPotential.out", model.PotentialScale());
  modelMesh.WriteDerivative(modelMesh.GravitationalPotentialCoefficient(),
                            "PREMAcceleration.out", model.PotentialScale());
  modelMesh.Write(modelMesh.DensityCoefficient(), "PREMDensity.out",
                  model.DensityScale());

  /*
  auto rho = model.Density();

  auto phi_exact = FunctionCoefficient([&model](const Vector& x) {
    auto rho = model.Density()(0, 1);
    auto G = model.GravitationalConstant();
    auto b = model.SurfaceRadius();
    return -2 * pi * G * rho * b * b * (1 - std::pow(x(0) / b, 2) / 3);
  });

  modelMesh.Write(phi_exact, "ex1e.out", model.PotentialScale());
  modelMesh.WriteDerivative(phi_exact, "ex1ed.out", model.PotentialScale());
  */
}
