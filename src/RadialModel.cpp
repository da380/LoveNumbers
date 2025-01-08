#include "LoveNumbers/RadialModel.hpp"

namespace LoveNumbers {

Int RadialModel::NumberOfBoundaries() const { return NumberOfLayers() + 1; }

Real RadialModel::JeanLength(Int degree) const {
  return 2 * std::numbers::pi_v<Real> * SurfaceRadius() /
         static_cast<Real>(degree + 1);
}

Real RadialModel::SurfaceRadius() const {
  return LayerRadii(NumberOfLayers() - 1).second;
}

void RadialModel::WriteAsDeckModel(const std::string &fileName,
                                   const std::array<std::string, 3> &header,
                                   Real maximumKnotSpacing) const {

  // Open the file for output.
  auto modelFile = std::ofstream(fileName);

  // Write the header information.
  for (auto line : header) {
    modelFile << line << "\n";
  }

  // Loop over the layers
  for (auto i = 0; i < NumberOfLayers(); i++) {
    auto attribute = i + 1;
    auto [r0, r1] = LayerRadii(i);
    auto n = std::max<Int>(2, std::round((r1 - r0) / maximumKnotSpacing));
    auto dr = (r1 - r0) / (n - 1);
    for (auto j = 0; j < n; j++) {
      auto r = r0 + j * dr;
      modelFile << std::format(
          "{:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f} "
          "{:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f}\n",
          r * LengthScale(), Density()(r, attribute) * DensityScale(),
          VerticalPVelocity()(r, attribute) * VelocityScale(),
          VerticalSVelocity()(r, attribute) * VelocityScale(),
          BulkQualityFactor()(r, attribute), ShearQualityFactor()(r, attribute),
          HorizontalPVelocity()(r, attribute) * VelocityScale(),
          HorizontalSVelocity()(r, attribute) * VelocityScale(),
          AnisotropicEtaParameter()(r, attribute));
    }
  }
}

void RadialModel::WriteAsDeckModel(const std::string &fileName,
                                   Real maximumKnotSpacing) const {
  auto message = std::string("Header line is ignored!");
  auto header = std::array<std::string, 3>{message, message, message};
  WriteAsDeckModel(fileName, header, maximumKnotSpacing);
}

std::function<Real(Real, Int)> RadialModel::VerticalPVelocity() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(LoveModulusC()(r, attribute) / Density()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::VerticalSVelocity() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(LoveModulusL()(r, attribute) / Density()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::HorizontalPVelocity() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(LoveModulusA()(r, attribute) / Density()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::HorizontalSVelocity() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(LoveModulusN()(r, attribute) / Density()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::AnisotropicEtaParameter() const {
  return [this](Real r, Real attribute) {
    return LoveModulusF()(r, attribute) /
           (LoveModulusA()(r, attribute) - 2 * LoveModulusL()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::BulkModulus() const {
  return [this](Real r, Real attribute) {
    constexpr auto ninth = static_cast<Real>(1) / static_cast<Real>(9);
    return ninth *
           (LoveModulusC()(r, attribute) +
            4 * (LoveModulusA()(r, attribute) - LoveModulusN()(r, attribute) +
                 LoveModulusF()(r, attribute)));
  };
}

std::function<Real(Real, Int)> RadialModel::ShearModulus() const {
  return [this](Real r, Real attribute) {
    constexpr auto fifteenth = static_cast<Real>(1) / static_cast<Real>(15);
    return fifteenth *
           (LoveModulusC()(r, attribute) + LoveModulusA()(r, attribute) +
            6 * LoveModulusL()(r, attribute) +
            5 * LoveModulusN()(r, attribute) -
            2 * LoveModulusF()(r, attribute));
  };
}

} // namespace LoveNumbers