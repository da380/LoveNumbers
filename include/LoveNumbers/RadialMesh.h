#pragma once

#include "mfem.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <ranges>
#include <string>
#include <vector>

namespace LoveNumbers {

class RadialMesh : public mfem::Mesh {
  using Int = int;
  using Real = mfem::real_t;

public:
  // Remove the default constructor.
  RadialMesh() = delete;

  // Return a range over the layer radii.
  auto Layers() const {
    using namespace std::ranges::views;
    return all(_layerRadii) | adjacent<2>;
  }

private:
  Int _numberOfLayers;
  std::vector<Real> _layerRadii;
};

} // namespace LoveNumbers