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
#include "LinearForms.hpp"
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
  std::unique_ptr<mfem::GridFunction> _densityDerivative;

public:
  // Constructor to be called from derived classes.
  RadialModel(const Dimensions &dimensions, Int order)
      : Dimensions(dimensions), _order{order},
        _L2{std::make_unique<mfem::L2_FECollection>(order - 1, 1)},
        _H1{std::make_unique<mfem::H1_FECollection>(order, 1)} {}

  // Return the polynomial order.
  auto Order() const { return _order; }

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
  virtual std::function<Real(Real, Int)> Density() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusA() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusC() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusF() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusL() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusN() const = 0;
  virtual std::function<Real(Real, Int)> BulkQualityFactor() const = 0;
  virtual std::function<Real(Real, Int)> ShearQualityFactor() const = 0;

  // Derived material parameter functions.
  std::function<Real(Real, Int)> VerticalPVelocity() const;
  std::function<Real(Real, Int)> VerticalSVelocity() const;
  std::function<Real(Real, Int)> HorizontalPVelocity() const;
  std::function<Real(Real, Int)> HorizontalSVelocity() const;
  std::function<Real(Real, Int)> AnisotropicEtaParameter() const;
  std::function<Real(Real, Int)> BulkModulus() const;
  std::function<Real(Real, Int)> ShearModulus() const;

  // MFEM Coefficients for material parameters.
  auto DensityCoefficient() const { return RadialCoefficient(Density()); }
  auto LoveModulusACoefficient() const {
    return RadialCoefficient(LoveModulusA());
  }
  auto LoveModulusCCoefficient() const {
    return RadialCoefficient(LoveModulusC());
  }
  auto LoveModulusFCoefficient() const {
    return RadialCoefficient(LoveModulusF());
  }
  auto LoveModulusLCoefficient() const {
    return RadialCoefficient(LoveModulusL());
  }
  auto LoveModulusNCoefficient() const {
    return RadialCoefficient(LoveModulusN());
  }
  auto BulkQualityFactorCoefficient() const {
    return RadialCoefficient(BulkQualityFactor());
  }
  auto ShearQualityFactorCoefficient() const {
    return RadialCoefficient(ShearQualityFactor());
  }
  auto BulkModulusCoefficient() const {
    return RadialCoefficient(BulkModulus());
  }
  auto ShearModulusCoefficient() const {
    return RadialCoefficient(ShearModulus());
  }

  // MFEM Coefficients for computed properties.
  auto GravitationalPotentialCoefficient() const {
    return mfem::GridFunctionCoefficient(_gravitationalPotential.get());
  }

  auto GravitationalAccelerationCoefficient() const {
    return mfem::GridFunctionCoefficient(_gravitationalAcceleration.get());
  }

  auto DensitDerivativeCoefficient() const {
    return mfem::GridFunctionCoefficient(_densityDerivative.get());
  }

  // Write a scalar GridFunction to a file using a simple format for plotting.
  void Write(const mfem::GridFunction &f, const std::string &file,
             Real scale = 1) const;

  // Write the derivagtive of a scalar Gridfunction to a file  using a simple
  // format for plotting.
  void WriteDerivative(const mfem::GridFunction &f, const std::string &file,
                       Real scale = 1) const;

  // Write a RadialCoefficient to a file using a simple format for
  // plotting.
  void Write(mfem::Coefficient &f, const std::string &file,
             Real scale = 1) const;

  void Write(mfem::Coefficient &&f, const std::string &file,
             Real scale = 1) const {
    Write(f, file, scale);
  }

  // Write the derivagtive of a RadialCoefficient to a file  using a simple
  // format for plotting.
  void WriteDerivative(mfem::Coefficient &f, const std::string &file,
                       Real scale = 1) const;

  void WriteDerivative(mfem::Coefficient &&f, const std::string &file,
                       Real scale = 1) const {
    WriteDerivative(f, file, scale);
  }

private:
  // Compute the surface gravitational acceleration and moment of inertia
  // factor via the radial integrals.
  void ComputeSurfaceGravityAndMomentOfInertiaFactor();

  // Compute the gravitational potential field through solution of the radial
  // Poisson equation.
  void ComputeGravitationalPotential();

  // Compute radial derivative of the density.
};

} // namespace LoveNumbers