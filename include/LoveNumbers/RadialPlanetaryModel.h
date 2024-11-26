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

using Int = int;
using Real = mfem::real_t;

class RadialPlanetaryModel {
public:
  RadialPlanetaryModel() = default;

  //----------------------------------------------------------------------------//
  //  Constructs the model from the layer radii and a vector of maximum element
  //  sizes for the layers. Structural parameters are not set.
  //----------------------------------------------------------------------------//

  RadialPlanetaryModel(const std::vector<Real> &layerRadii,
                       const std::vector<Real> &maximumElementSizes,
                       const std::vector<bool> &isSolid)
      : _layerRadii{layerRadii}, _isSolid{isSolid} {

    // Store the layer data.
    assert(NumberOfLayers() == maximumElementSizes.size());

    // Build up the mesh.
    BuildMesh(maximumElementSizes);
  }

  RadialPlanetaryModel(const std::vector<Real> &layerRadii,
                       Real maximumElementSize)
      : RadialPlanetaryModel(
            layerRadii,
            std::vector<Real>(layerRadii.size() - 1, maximumElementSize),
            std::vector<bool>(layerRadii.size() - 1, true)) {}

  // The number of layers in the model.
  Int NumberOfLayers() const { return _layerRadii.size() - 1; }

  // The number of boundary elements.
  Int NumberOfBoundaryElements() const { return NumberOfLayers() + 1; }

  // A range indexing the layers.
  auto LayerIndices() const {
    return std::ranges::views::iota(Int(0), static_cast<Int>(NumberOfLayers()));
  }

  // A range over the layers returning the lower and upper radii as a std::pair;
  auto Layers() const {
    return std::ranges::views::all(_layerRadii) |
           std::ranges::views::adjacent<2>;
  }

  // A range over the layers that equals true if the layer is solid
  // and false otherwise.
  auto &IsSolid() const { return _isSolid; }

  // Return a vector of domain attributes for the layers.
  auto &DomainAttributes() const { return _domainAttributes; }

  // Return a marker array for solid regions.
  mfem::Array<Int> SolidMarker() const;

  // Return a marker array for fluid regions.
  mfem::Array<Int> FluidMarker() const;

  // Return a marker for the free surface.
  mfem::Array<Int> FreeSurfaceMarker() const;

  // Return a marker for fluid-solid boundaries.
  mfem::Array<Int> FluidSolidMarker() const;

  // Return a marker for solid-fluid boundaries.
  mfem::Array<Int> SolidFluidMarker() const;

  // Return a marker for solid-solid boundaries.
  mfem::Array<Int> SolidSolidMarker() const;

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
  std::vector<bool> _isSolid;

  // MFEM mesh and linked data.
  mfem::Mesh _mesh;
  mfem::Array<Int> _domainAttributes;

  // Density data.
  mfem::PWCoefficient _rho;
  mfem::Array<mfem::FunctionCoefficient *> _rho_functions;

  void BuildMesh(const std::vector<Real> &);
};

} // namespace LoveNumbers