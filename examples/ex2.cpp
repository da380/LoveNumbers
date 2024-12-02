#include <fstream>
#include <iostream>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>

#include "LoveNumbers/DeckModel.h"

using Real = double;

int main() {

  auto fileName = std::string("../data/prem.200.no");

  auto model = LoveNumbers::DeckModel<Real>(fileName);

  auto out = std::ofstream("rho.out");

  for (auto [i, radii] : std::ranges::views::enumerate(model.LayerRadii())) {
    auto [r1, r2] = radii;
    auto n = 100;
    auto dr = (r2 - r1) / (n - 1);
    auto [i1, i2] = model.LayerKnots()[i];
    auto f = model.L(i);
    for (auto j = 0; j < n; j++) {
      auto r = r1 + j * dr;
      out << r << " " << f(r) << std::endl;
    }
  }
}