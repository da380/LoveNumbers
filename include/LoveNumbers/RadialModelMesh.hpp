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
#include "RadialModel.hpp"
#include "mfem.hpp"

namespace LoveNumbers {

class RadialModelMesh {
private:
  // const Reference of the underlying radial model.
  const RadialModel &_model;

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
  RadialModelMesh(const RadialModel &model, Int order,
                  Real characteristicLengthScale)
      : _model{model}, _order{order},
        _L2{std::make_unique<mfem::L2_FECollection>(order - 1, 1)},
        _H1{std::make_unique<mfem::H1_FECollection>(order, 1)} {
    BuildMesh(characteristicLengthScale);
  }

  // Return a reference to the radial model
  const auto &GetRadialModel() const { return _model; }

  // Return a reference to the mesh.
  auto &GetMesh() { return _mesh; }
  const auto &GetMesh() const { return _mesh; }

  // Return the polynomial order.
  auto Order() const { return _order; }

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

  // Return pointers to the finite element spaces.
  auto &L2Space() const { return *_L2Space; }
  auto &H1Space() const { return *_H1Space; }

  // Print the radial mesh to the given file.
  void PrintMesh(const std::string &mesh_file);

  // MFEM Coefficients for material parameters.
  auto DensityCoefficient() const {
    return RadialCoefficient(_model.Density());
  }
  auto LoveModulusACoefficient() const {
    return RadialCoefficient(_model.LoveModulusA());
  }
  auto LoveModulusCCoefficient() const {
    return RadialCoefficient(_model.LoveModulusC());
  }
  auto LoveModulusFCoefficient() const {
    return RadialCoefficient(_model.LoveModulusF());
  }
  auto LoveModulusLCoefficient() const {
    return RadialCoefficient(_model.LoveModulusL());
  }
  auto LoveModulusNCoefficient() const {
    return RadialCoefficient(_model.LoveModulusN());
  }
  auto BulkQualityFactorCoefficient() const {
    return RadialCoefficient(_model.BulkQualityFactor());
  }
  auto ShearQualityFactorCoefficient() const {
    return RadialCoefficient(_model.ShearQualityFactor());
  }
  auto BulkModulusCoefficient() const {
    return RadialCoefficient(_model.BulkModulus());
  }
  auto ShearModulusCoefficient() const {
    return RadialCoefficient(_model.ShearModulus());
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