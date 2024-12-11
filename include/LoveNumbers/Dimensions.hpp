#pragma once

#include "LoveNumbers/Configure.hpp"
#include <cmath>
#include <mutex>

namespace LoveNumbers {

class Dimensions {
private:
  const Real _gravitationalConstant = 6.67430e-11;
  Real _lengthScale = 6.371e6;
  Real _densityScale = 5.514e3;
  Real _timeScale;

public:
  // Construct using default values.
  Dimensions() : _timeScale{TimeScaleFromGravitationalConstant()} {}

  // Construct with user-defined time scale.
  Dimensions(Real timeScale) : _timeScale{timeScale} {}

  // Construction with user-defined length and density scales,
  // with the time scale then set using default methods.
  Dimensions(Real lengthScale, Real densityScale)
      : _lengthScale{lengthScale}, _densityScale{densityScale},
        _timeScale{TimeScaleFromGravitationalConstant()} {}

  // Construct with user-defined length, density, and time scales.
  Dimensions(Real lengthScale, Real densityScale, Real timeScale)
      : _lengthScale{lengthScale}, _densityScale{densityScale},
        _timeScale{timeScale} {}

  // Default copy and move constructors.
  Dimensions(const Dimensions &) = default;
  Dimensions(Dimensions &&) = default;

  // Return the length-scale
  Real LengthScale() const { return _lengthScale; }

  // Return the density scale.
  Real DensityScale() const { return _densityScale; };

  // Return the time-scale.
  Real TimeScale() const { return _timeScale; }

  // Return the gravitational constant.
  Real GravitationalConstant() const {
    return _gravitationalConstant * DensityScale() * std::pow(TimeScale(), 2);
  }

  // Return the mass-scale.
  Real MassScale() const { return _densityScale * std::pow(_lengthScale, 3); }

  // Return the velocity scale.
  Real VelocityScale() const { return LengthScale() / TimeScale(); }

  // Return the acceleration scale.
  Real AccelerationScale() const { return VelocityScale() / TimeScale(); }

  // Return the force scale.
  Real ForceScale() const { return MassScale() * AccelerationScale(); }

  // Return the traction scale.
  Real TractionScale() const {
    return ForceScale() / std::pow(LengthScale(), 2);
  }

  // Return the potential energy scale.
  Real PotentialScale() const { return AccelerationScale() * LengthScale(); }

private:
  constexpr Real TimeScaleFromGravitationalConstant() const {
    return static_cast<Real>(1) /
           std::sqrt(std::numbers::pi_v<Real> * _gravitationalConstant *
                     DensityScale());
  }
};

} // namespace LoveNumbers