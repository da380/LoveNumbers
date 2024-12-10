

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numbers>
#include <ranges>
#include <string>
#include <vector>

#include "LoveNumbers/LoveNumbers.hpp"

int main(int argc, char *argv[]) {

  auto model = LoveNumbers::DeckModel::FromFEMOrderAndMaximumDegree(
      "../data/prem.200.no", 4, 64);

  auto *L2Space = model.L2Space();

  auto b = mfem::LinearForm(L2Space);

  auto kernel = LoveNumbers::RadialCoefficient(
      [](auto r, auto attribute) { return r * r; });
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(kernel));
  b.Assemble();

  auto rho = mfem::GridFunction(L2Space);
  auto rhoCoefficient = model.RhoCoefficient();
  rho.ProjectCoefficient(rhoCoefficient);

  auto factor = 4 * std::numbers::pi_v<LoveNumbers::Real> *
                model.GravitationalConstant() /
                (std::pow(model.SurfaceRadius(), 2));

  auto g = factor * (b(rho));

  std::cout << g * model.AccelerationScale() << std::endl;
}
