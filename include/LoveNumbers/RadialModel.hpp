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

public:
  // Constructor to be called from derived classes.
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

  // Write the model out in deck format.
  void WriteAsDeckModel(const std::string &fileName,
                        const std::array<std::string, 3> &header,
                        Real maximumKnotSpacing) const;

  // Write the model out in deck format using default header lines.
  void WriteAsDeckModel(const std::string &fileName,
                        Real maximumKnotSpacing) const;

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
};

} // namespace LoveNumbers