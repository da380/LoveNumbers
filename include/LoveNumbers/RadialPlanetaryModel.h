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

    // Store the layer data.
    std::ranges::copy(layerRadii, std::back_inserter(_layerRadii));
    std::ranges::copy(maximumElementSizes,
                      std::back_inserter(_maximumElementSizes));
    assert(NumberOfLayers() == MaximumElementSizes().size());

    // Build up the mesh.
    BuildMesh();

    // Allocate space for the FunctionCoefficients.
    _rho_functions.SetSize(NumberOfLayers());

    // Set the function coefficients.
    for (auto i : LayerIndices()) {
      _rho_functions[i] = new mfem::FunctionCoefficient(
          [](mfem::Vector r) -> Real { return 1; });
    }

    _rho = mfem::PWCoefficient(DomainAttributes(), _rho_functions);
  }

  template <std::ranges::range R1>
  RadialPlanetaryModel(const R1 &layerRadii, Real maximumElementSize)
      : RadialPlanetaryModel(
            layerRadii,
            std::vector<Real>(layerRadii.size() - 1, maximumElementSize)) {}

  // The number of layers in the model.
  auto NumberOfLayers() const { return _layerRadii.size() - 1; }

  // A range indexing the layers.
  auto LayerIndices() const {
    return std::ranges::views::iota(Int(0), static_cast<Int>(NumberOfLayers()));
  }

  // A range over the layers returning the lower and upper radii as a std::pair;
  auto Layers() const {
    return std::ranges::views::all(_layerRadii) |
           std::ranges::views::adjacent<2>;
  }

  // A range over the layers returning the specified maximum
  // element size.
  auto &MaximumElementSizes() const { return _maximumElementSizes; }

  // A range over the layers that equals true if the layer is solid
  // and false otherwise.
  auto &IsSolid() const { return _isSolid; }

  // Return a vector of domain attributes for the layers.
  auto &DomainAttributes() const { return _domainAttributes; }

  // Return a reference to the mesh.
  auto &Mesh() { return _mesh; }

  // Return a const reference to the mesh.
  auto &Mesh() const { return _mesh; }

  // Print the radial mesh to the given file.
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
  std::vector<bool> _isSolid;

  // MFEM mesh and linked data.
  mfem::Mesh _mesh;
  mfem::Array<Int> _domainAttributes;

  // Density data.
  mfem::PWCoefficient _rho;
  mfem::Array<mfem::FunctionCoefficient *> _rho_functions;

  void BuildMesh();
};

} // namespace LoveNumbers