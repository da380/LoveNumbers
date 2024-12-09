#pragma once

#include "Configure.hpp"
#include "LoveNumbers/Configure.hpp"
#include "LoveNumbers/RadialModel.hpp"
#include "RadialModel.hpp"

namespace LoveNumbers {

class PWHomogeneousModel : RadialModel {
private:
  std::vector<Real> _r, _rho, _A, _C, _F, _L, _N;

public:
};

} // namespace LoveNumbers