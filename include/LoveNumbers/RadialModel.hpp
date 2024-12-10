#pragma once

#include "Dimensions.hpp"
#include "RadialCoefficient.hpp"
#include "mfem.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <format>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <numbers>
#include <numeric>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

#include "LoveNumbers/Configure.hpp"
#include "LoveNumbers/Dimensions.hpp"

namespace LoveNumbers {

class RadialModel : public Dimensions {
private:
  // MFEM mesh.
  mfem::Mesh _mesh;

  // Contains mesh domain attribute for each layer.
  mfem::Array<Int> _layerAttributes;

  // Finite element data.
  std::unique_ptr<mfem::L2_FECollection> _L2;
  std::unique_ptr<mfem::H1_FECollection> _H1;

  // Finite element space.
  std::unique_ptr<mfem::FiniteElementSpace> _L2Space;
  std::unique_ptr<mfem::FiniteElementSpace> _H1Space;

public:
  // Construct with default dimension scheme.
  RadialModel() = default;

  // Construct with explicit dimensions scheme.
  RadialModel(const Dimensions &dimensions) : Dimensions(dimensions) {}

  // Return the number of layers.
  virtual Int NumberOfLayers() const = 0;

  // Return the bounding radii of the ith layer.
  virtual std::pair<Real, Real> LayerRadii(Int i) const = 0;

  // Return true if the ith layer is solid.
  virtual bool LayerIsSolid(Int i) const = 0;

  // Returns the number of boundaries.
  Int NumberOfBoundaries() const;

  // Returns a range over the layer Indices.
  auto LayerIndices() const {
    return std::ranges::views::iota(0, NumberOfLayers());
  }

  // Return the surface radius of the planet.
  Real SurfaceRadius() const;

  // Return the Jean length for a given degree (in non-dimensional form).
  Real JeanLength(Int degree) const;

  // Return a vector of domain attributes for the layers.
  mfem::Array<Int> LayerAttributes() const;

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

  // Return a marker for all boundaries.
  mfem::Array<Int> AllBoundaryMarker() const;

  // Return a marker for all boundaries but that at the model's centre
  // (which is not a real physical boundary).
  mfem::Array<Int> AllBoundaryMarkerCentreExcluded() const;

  // Build the mesh given a maximum element sizes for each layer (in
  // non-dimensionalised form).
  void BuildMesh(const std::vector<Real> &maximumElementSizes);

  mfem::Mesh &Mesh() { return _mesh; }

  // Build the mesh given a maximum element size (in non-dimensionalised form).
  void BuildMesh(Real maximumElementSize);

  // Set up the finite element space.
  void SetFiniteElementSpaces(Int order);

  // Return pointers to the finite element spaces.
  mfem::FiniteElementSpace *L2Space() const { return _L2Space.get(); }
  mfem::FiniteElementSpace *H1Space() const { return _H1Space.get(); }

  // Print the radial mesh to the given file.
  void PrintMesh(const std::string &mesh_file);

  // Write the model out in deck format.
  void WriteAsDeckModel(const std::string &fileName,
                        const std::array<std::string, 3> &header,
                        Real maximumKnotSpacing);

  // Write the model out in deck format using default header lines.
  void WriteAsDeckModel(const std::string &fileName, Real maximumKnotSpacing);

  // Material parameter functions for override.
  virtual std::function<Real(Real, Int)> Rho() const = 0;
  virtual std::function<Real(Real, Int)> A() const = 0;
  virtual std::function<Real(Real, Int)> C() const = 0;
  virtual std::function<Real(Real, Int)> F() const = 0;
  virtual std::function<Real(Real, Int)> L() const = 0;
  virtual std::function<Real(Real, Int)> N() const = 0;
  virtual std::function<Real(Real, Int)> QKappa() const = 0;
  virtual std::function<Real(Real, Int)> QMu() const = 0;

  // Derived material parameter functions.
  std::function<Real(Real, Int)> VPV() const;
  std::function<Real(Real, Int)> VSV() const;
  std::function<Real(Real, Int)> VPH() const;
  std::function<Real(Real, Int)> VSH() const;
  std::function<Real(Real, Int)> Eta() const;
  std::function<Real(Real, Int)> Kappa() const;
  std::function<Real(Real, Int)> Mu() const;

  // MFEM Coefficients for material parameters.
  auto RhoCoefficient() const { return RadialCoefficient(Rho()); }
  auto ACoefficient() const { return RadialCoefficient(A()); }
  auto CCoefficient() const { return RadialCoefficient(C()); }
  auto FCoefficient() const { return RadialCoefficient(F()); }
  auto LCoefficient() const { return RadialCoefficient(L()); }
  auto NCoefficient() const { return RadialCoefficient(N()); }
  auto QKappaCoefficient() const { return RadialCoefficient(QKappa()); }
  auto QMuCoefficient() const { return RadialCoefficient(QMu()); }
  auto KappaCoefficient() const { return RadialCoefficient(Kappa()); }
  auto MuCoefficient() const { return RadialCoefficient(QMu()); }

private:
};

} // namespace LoveNumbers