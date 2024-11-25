#pragma once

#include "mfem.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

namespace LoveNumbers {

class RadialPlanetaryModel {
  using Int = int;
  using Real = mfem::real_t;

public:
  RadialPlanetaryModel() = default;

  //----------------------------------------------------------------------------//
  //  Constructs the model from the layer radii and a vector of maximum element
  //  sizes for the layers. Structural parameters are not set.
  //----------------------------------------------------------------------------//
  template <std::ranges::range R1, std::ranges::range R2>
  RadialPlanetaryModel(const R1 &layerRadii, const R2 &maximumElementSizes) {

    std::ranges::copy(layerRadii, std::back_inserter(_layerRadii));
    std::ranges::copy(maximumElementSizes,
                      std::back_inserter(_maximumElementSizes));
    assert(NumberOfLayers() == MaximumElementSizes().size());
    BuildMesh();
  }

  //----------------------------------------------------------------------------//
  //  Constructs the model from the layer radii and a maximum element
  //  size for each of the layers. Structural parameters are not set.
  //----------------------------------------------------------------------------//
  template <std::ranges::range R1>
  RadialPlanetaryModel(const R1 &layerRadii, Real maximumElementSize) {

    std::ranges::copy(layerRadii, std::back_inserter(_layerRadii));
    _maximumElementSizes =
        std::vector<Real>(NumberOfLayers(), maximumElementSize);
    BuildMesh();
  }

  auto NumberOfLayers() const { return _layerRadii.size() - 1; }

  auto Layers() const {
    return std::ranges::views::all(_layerRadii) |
           std::ranges::views::adjacent<2>;
  }

  auto MaximumElementSizes() const {
    return std::ranges::views::all(_maximumElementSizes);
  }

  void PrintMesh(const std::string &mesh_file) {
    std::ofstream ofs(mesh_file);
    ofs.precision(8);
    _mesh.Print(ofs);
    ofs.close();
  }

private:
  // Layer information.
  std::vector<Real> _layerRadii;
  std::vector<Real> _maximumElementSizes;

  // MFEM mesh.
  mfem::Mesh _mesh;

  // Density data.
  mfem::PWCoefficient _rho;
  mfem::Array<mfem::FunctionCoefficient *> _rho_functions;

  void BuildMesh();
};

} // namespace LoveNumbers