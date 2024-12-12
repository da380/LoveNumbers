

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
  Real characteristicLengthScale = 0.1;
  auto model = DeckModel(Dimensions(), order, characteristicLengthScale,
                         "../data/prem.200.no");

  model.WriteCoefficient(model.GravitationalAccelerationCoefficient(),
                         "ex1.out", model.AccelerationScale());
}
