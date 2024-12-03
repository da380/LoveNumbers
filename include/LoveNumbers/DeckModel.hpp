#pragma once

#include "LoveNumbers/RadialModel.hpp"
#include "LoveNumbers/Types.hpp"
#include <Interpolation/CubicSpline>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <ranges>
#include <string>
#include <vector>

namespace LoveNumbers {

class DeckModel : public RadialModel {
private:
  // Model values at the radial knots.
  std::vector<Real> _r;      // radius
  std::vector<Real> _rho;    // density
  std::vector<Real> _A;      // Love parameter, A
  std::vector<Real> _C;      // Love parameter, C
  std::vector<Real> _F;      // Love parameter, F
  std::vector<Real> _L;      // Love parameter, L
  std::vector<Real> _N;      // Love parameter, N
  std::vector<Real> _QKappa; // Bulk quality factor
  std::vector<Real> _QMu;    // Shear quality factor

  std::vector<std::pair<Int, Int>> _boundaryIndices;
  std::vector<std::pair<Real, Real>> _boundaryRadii;
  std::vector<bool> _layerSolid;

public:
  DeckModel() = delete;

  DeckModel(const std::string &fileName, Real maximumElementSize);

  Int NumberOfLayers() const override { return _boundaryIndices.size(); }

  std::pair<Real, Real> LayerRadii(Int i) const override {
    return _boundaryRadii[i];
  }

  bool LayerIsSolid(Int i) const override { return _layerSolid[i]; }

  Int NumberOfKnots() const { return _r.size(); }

private:
};

} // namespace LoveNumbers