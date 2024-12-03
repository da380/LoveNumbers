#pragma once

#include "mfem.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <memory>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

#include "LoveNumbers/Types.hpp"

namespace LoveNumbers {

class RadialModel {
private:
  // Constants for non-dimensionalisation, with
  // default values set.
  Real _lengthScale = 6371000;
  Real _massScale = 5.972e24;
  Real _timeScale = 3600;

  // MFEM mesh.
  mfem::Mesh _mesh;

  // Contains mesh domain attribute for each layer.
  mfem::Array<Int> _layerAttributes;

  // Finite element collection.
  std::unique_ptr<mfem::FiniteElementCollection> _FECollection;

  // Finite element space.
  std::unique_ptr<mfem::FiniteElementSpace> _FESpace;

public:
  RadialModel() = default;

  RadialModel(Real lengthScale, Real massScale, Real timeScale)
      : _lengthScale{lengthScale}, _massScale{massScale},
        _timeScale{timeScale} {}

  // Return the number of layers.
  virtual Int NumberOfLayers() const = 0;

  // Return the bounding radii of the ith layer.
  virtual std::pair<Real, Real> LayerRadii(Int i) const = 0;

  // Return true if the ith layer is solid.
  virtual bool LayerIsSolid(Int i) const = 0;

  // Returns the number of boundaries.
  Int NumberOfBoundaries() const { return NumberOfLayers() + 1; }

  // Returns a range over the layer Indices.
  auto LayerIndices() const {
    return std::ranges::views::iota(0, NumberOfLayers());
  }

  // Return a vector of domain attributes for the layers.
  auto &LayerAttributes() const { return _layerAttributes; }

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

  // Return a marker for solid-solid boundaries.
  mfem::Array<Int> FluidFluidMarker() const;

  // Return the length-scale
  auto LengthScale() const { return _lengthScale; }

  // Return the mass-scale.
  auto MassScale() const { return _massScale; }

  // Return the time-scale.
  auto TimeScale() const { return _timeScale; }

  // Return the density scale.
  auto DensityScale() const { return MassScale() / std::pow(LengthScale(), 3); }

  // Return the velocity scale.
  auto VelocityScale() const { return LengthScale() / TimeScale(); }

  // Return the acceleration scale.
  auto AccelerationScale() const { return VelocityScale() / TimeScale(); }

  // Return the force scale.
  auto ForceScale() const { return MassScale() * AccelerationScale(); }

  // Return the traction scale.
  auto TractionScale() const {
    return ForceScale() / std::pow(LengthScale(), 2);
  }

  // Build the mesh.
  void BuildMesh(Real maximumElementSize);

  // Print the radial mesh to the given file.
  void PrintMesh(const std::string &mesh_file) {
    std::ofstream ofs(mesh_file);
    ofs.precision(8);
    _mesh.Print(ofs);
    ofs.close();
  }

private:
};

/*

using Int = int;
using Real = mfem::real_t;

class RadialPlanetaryModel : public DeckModel<Real> {
private:

  // Store the linked deck model
  DeckModel<Real> _deckModel;

  // Layer information.
  std::vector<Real> _layerRadii;
  std::vector<bool> _isSolid;

  // MFEM mesh and linked data.
  mfem::Mesh _mesh;
  mfem::Array<Int> _domainAttributes;

  // Finite element data.
  Int _order;

  // Density data.
  mfem::PWCoefficient _rho;
  mfem::Array<mfem::FunctionCoefficient *> _rho_functions;



public:

RadialPlanetaryModel() = default;

RadialPlanetaryModel(const std::string &deckModelFile,
                     Real maximumElementSize);

// The number of layers in the model.
Int NumberOfLayers() const { return _deckModel.LayerRadii().size(); }

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

// Return the finite element collection.
// auto &FECollection() const { return H1_FECollection(order, dim); }




private:
};

*/

} // namespace LoveNumbers