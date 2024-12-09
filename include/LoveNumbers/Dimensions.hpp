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

  // Construct with user-defined time-scale.
  Dimensions(Real timeScale) : _timeScale{timeScale} {}

  Dimensions(Real lengthScale, Real densityScale)
      : _lengthScale{lengthScale}, _densityScale{densityScale},
        _timeScale{TimeScaleFromGravitationalConstant()} {}

  Dimensions(Real lengthScale, Real densityScale, Real timeScale)
      : _lengthScale{lengthScale}, _densityScale{densityScale},
        _timeScale{timeScale} {}

  // Return the density scale.
  Real DensityScale() const { return _densityScale; };

  // Return the length-scale
  Real LengthScale() const { return _lengthScale; }

  // Return the time-scale.
  Real TimeScale() const { return _timeScale; }

  // Return the mass-scale.
  Real MassScale() const { return _densityScale * std::pow(_lengthScale, 3); }

  // Return the velocity scale.
  Real Dimensions::VelocityScale() const { return LengthScale() / TimeScale(); }

  // Return the acceleration scale.
  Real Dimensions::AccelerationScale() const {
    return VelocityScale() / TimeScale();
  }

  // Return the force scale.
  Real Dimensions::ForceScale() const {
    return MassScale() * AccelerationScale();
  }

  // Return the traction scale.
  Real Dimensions::TractionScale() const {
    return ForceScale() / std::pow(LengthScale(), 2);
  }

private:
  constexpr Real TimeScaleFromGravitationalConstant() const;
};

} // namespace LoveNumbers