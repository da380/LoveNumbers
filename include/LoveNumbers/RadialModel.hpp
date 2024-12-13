#pragma once

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

#include "Coefficients.hpp"
#include "Configure.hpp"
#include "Dimensions.hpp"
#include "mfem.hpp"

namespace LoveNumbers {

class RadialModel : public Dimensions {
private:
  // MFEM mesh.
  mfem::Mesh _mesh;

  // Contains mesh domain attribute for each layer.
  mfem::Array<Int> _layerAttributes;

  // Finite element data.
  Int _order;
  std::unique_ptr<mfem::L2_FECollection> _L2;
  std::unique_ptr<mfem::H1_FECollection> _H1;

  // Finite element space.
  std::unique_ptr<mfem::FiniteElementSpace> _L2Space;
  std::unique_ptr<mfem::FiniteElementSpace> _H1Space;

  // Computed and stored properties of the model
  Real _surfaceGravity;
  Real _momentOfInertiaFactor;

  std::unique_ptr<mfem::GridFunction> _gravitationalPotential;
  std::unique_ptr<mfem::GridFunction> _gravitationalAcceleration;

public:
  // Constructor to be called from derived classes.
  RadialModel(const Dimensions &dimensions, Int order)
      : Dimensions(dimensions), _order{order},
        _L2{std::make_unique<mfem::L2_FECollection>(order - 1, 1)},
        _H1{std::make_unique<mfem::H1_FECollection>(order, 1)} {}

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

  // Return the surface gravitational acceleration.
  Real SurfaceGravity() const { return _surfaceGravity; }

  // Return the moment of inertia factor.
  Real MomentOfInertiaFactor() const { return _momentOfInertiaFactor; }

  // Return a vector of domain attributes for the layers.
  mfem::Array<Int> LayerAttributes() const;

  // Return a marker array for solid regions.
  mfem::Array<Int> SolidMarker() const;

  // Return a marker array for fluid regions.
  mfem::Array<Int> FluidMarker() const;

  // Return a marker for the free surface.
  mfem::Array<Int> SurfaceMarker() const;

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

  // Return a marker just for the central boundary.
  mfem::Array<Int> CentreMarker() const;

  // Return a marker just for the central boundary.
  mfem::Array<Int> CentreAndSurfaceMarker() const;

  // Build the mesh given a maximum element sizes for each layer (in
  // non-dimensionalised form).
  void BuildMesh(Real characteristicLengthScale);

  // Return a reference to the mesh.
  mfem::Mesh &Mesh() { return _mesh; }
  const mfem::Mesh &Mesh() const { return _mesh; }

  // Return pointers to the finite element spaces.
  auto &L2Space() const { return *_L2Space; }
  auto &H1Space() const { return *_H1Space; }

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
  auto RhoCoefficient() const { return RadialModelCoefficient(Rho()); }
  auto ACoefficient() const { return RadialModelCoefficient(A()); }
  auto CCoefficient() const { return RadialModelCoefficient(C()); }
  auto FCoefficient() const { return RadialModelCoefficient(F()); }
  auto LCoefficient() const { return RadialModelCoefficient(L()); }
  auto NCoefficient() const { return RadialModelCoefficient(N()); }
  auto QKappaCoefficient() const { return RadialModelCoefficient(QKappa()); }
  auto QMuCoefficient() const { return RadialModelCoefficient(QMu()); }
  auto KappaCoefficient() const { return RadialModelCoefficient(Kappa()); }
  auto MuCoefficient() const { return RadialModelCoefficient(QMu()); }

  // MFEM Coefficients for computed properties.
  auto GravitationalPotentialCoefficient() const {
    return mfem::GridFunctionCoefficient(_gravitationalPotential.get());
  }

  auto GravitationalAccelerationCoefficient() const {
    return mfem::GridFunctionCoefficient(_gravitationalAcceleration.get());
  }

  // Write a scalar GridFunction to a file using a simple format for plotting.
  void WriteGridFunction(const mfem::GridFunction &f, const std::string &file,
                         Real scale = 1) const;

  // Write a RadialModelCoefficient to a file using a simple format for
  // plotting.
  void WriteCoefficient(mfem::Coefficient &&f, const std::string &file,
                        Real scale = 1) const;

private:
  // Compute the surface gravitational acceleration and moment of inertia
  // factor via the radial integrals.
  void ComputeSurfaceGravityAndMomentOfInertiaFactor();

  // Compute the gravitational potential field through solution of the radial
  // Poisson equation.
  void ComputeGravitationalPotential();
};

} // namespace LoveNumbers