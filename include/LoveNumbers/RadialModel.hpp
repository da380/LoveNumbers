#pragma once

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

namespace LoveNumbers {

class RadialModel {
private:
  // Constants for non-dimensionalisation.
  Real _lengthScale;
  Real _massScale;
  Real _timeScale;

  // MFEM mesh.
  mfem::Mesh _mesh;

  // Contains mesh domain attribute for each layer.
  mfem::Array<Int> _layerAttributes;

  // PWCoefficients for material parameters.
  mfem::PWCoefficient _rhoCoefficient;
  mfem::PWCoefficient _ACoefficient;
  mfem::PWCoefficient _CCoefficient;
  mfem::PWCoefficient _FCoefficient;
  mfem::PWCoefficient _LCoefficient;
  mfem::PWCoefficient _NCoefficient;

  // Function coefficients for each layer
  std::vector<mfem::FunctionCoefficient> _rhoFunctionCoefficients;
  std::vector<mfem::FunctionCoefficient> _AFunctionCoefficients;
  std::vector<mfem::FunctionCoefficient> _CFunctionCoefficients;
  std::vector<mfem::FunctionCoefficient> _FFunctionCoefficients;
  std::vector<mfem::FunctionCoefficient> _LFunctionCoefficients;
  std::vector<mfem::FunctionCoefficient> _NFunctionCoefficients;

  // Finite element data.
  std::unique_ptr<mfem::H1_FECollection> _FECollection;

  // Finite element space.
  std::unique_ptr<mfem::FiniteElementSpace> _FESpace;

public:
  RadialModel() = delete;

  RadialModel(Real lengthScale, Real massScale, Real timeScale)
      : _lengthScale{lengthScale}, _massScale{massScale},
        _timeScale{timeScale} {}

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

  // Return the density in the ith layer. Must be defined in the
  // derived class.
  virtual std::function<Real(Real)> Rho(Int i) const = 0;

  // Return the modulus A in the ith layer. Must be defined in the
  // derived class.
  virtual std::function<Real(Real)> A(Int i) const = 0;

  // Return the modulus C in the ith layer. Must be defined in the
  // derived class.
  virtual std::function<Real(Real)> C(Int i) const = 0;

  // Return the modulus F in the ith layer. Must be defined in the
  // derived class.
  virtual std::function<Real(Real)> F(Int i) const = 0;

  // Return the modulus L in the ith layer. Must be defined in the
  // derived class.
  virtual std::function<Real(Real)> L(Int i) const = 0;

  // Return the modulus N in the ith layer. Must be defined in the
  // derived class.
  virtual std::function<Real(Real)> N(Int i) const = 0;

  // Return the bulk quality factor in the ith layer. Must be defined in the
  // derived class.
  virtual std::function<Real(Real)> QKappa(Int i) const = 0;

  // Return the shear quality factor in the ith layer. Must be defined in the
  // derived class.
  virtual std::function<Real(Real)> QMu(Int i) const = 0;

  // Return vertical P-velocity in the ith layer.
  std::function<Real(Real)> VPV(Int i) const;

  // Return horizontal P-velocity in the ith layer.
  std::function<Real(Real)> VPH(Int i) const;

  // Return vertical S-velocity in the ith layer.
  std::function<Real(Real)> VSV(Int i) const;

  // Return horizontal S-velocity in the ith layer.
  std::function<Real(Real)> VSH(Int i) const;

  // Return transversely isotropic eta parameter in the ith layer.
  std::function<Real(Real)> Eta(Int i) const;

  // Return the effective bulk modulus in the ith layer.
  std::function<Real(Real)> Kappa(Int i) const;

  // Return the effective shear modulus in the ith layer.
  std::function<Real(Real)> Mu(Int i) const;

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

  // Build the mesh given a maximum element size (in non-dimensionalised form).
  void BuildMesh(Real maximumElementSize);

  // Set up the finite element space.
  void SetFiniteElementSpace(Int order);

  // Return a const reference ot the finite element space.
  auto &FiniteElementSpace() const { return *_FESpace; }

  // Print the radial mesh to the given file.
  void PrintMesh(const std::string &mesh_file) {
    std::ofstream ofs(mesh_file);
    ofs.precision(8);
    _mesh.Print(ofs);
    ofs.close();
  }

  // Write the model out in deck format.
  void WriteAsDeckModel(const std::string &fileName,
                        const std::array<std::string, 3> &header,
                        Real maximumKnotSpacing);

  // Write the model out in deck format using default header lines.
  void WriteAsDeckModel(const std::string &fileName, Real maximumKnotSpacing);

private:
};

} // namespace LoveNumbers