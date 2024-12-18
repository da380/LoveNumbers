#include "LoveNumbers/DeckModel.hpp"

namespace LoveNumbers {

Int DeckModel::NumberOfLayers() const { return _boundaryRadius.size() - 1; }

std::pair<Real, Real> DeckModel::LayerRadii(Int i) const {
  return {_boundaryRadius[i], _boundaryRadius[i + 1]};
}

bool DeckModel::LayerIsSolid(Int i) const { return _layerSolid[i]; }

Int DeckModel::NumberOfKnots() const { return _r.size(); }

void DeckModel::ReadModelFile(const std::string &fileName) {
  // Check the file exists.
  if (!std::filesystem::exists(fileName)) {
    throw std::runtime_error("deck file does not exist");
  }

  // Open the file.
  auto modelFile = std::ifstream(fileName);

  // Store the header information
  for (auto i = 0; i < 3; i++) {
    auto line = std::string();
    std::getline(modelFile, line);
    _header.push_back(line);
  }

  // Read in the model data.
  auto line = std::string();
  while (std::getline(modelFile, line)) {

    // Get the model parameters at the knot.
    Real r, rho, vpv, vsv, QKappa, QMu, vph, vsh, eta;
    if (!(std::stringstream(line) >> r >> rho >> vpv >> vsv >> QKappa >> QMu >>
          vph >> vsh >> eta)) {
      break;
    }

    // Covert velocities to modulii.
    auto A = rho * vph * vph;
    auto C = rho * vpv * vpv;
    auto L = rho * vsv * vsv;
    auto N = rho * vsh * vsh;
    auto F = eta * (A - 2 * L);

    // Store the non-dimensionalised values.
    _r.push_back(r / LengthScale());
    _density.push_back(rho / DensityScale());
    _LoveModulusA.push_back(A / TractionScale());
    _LoveModulusC.push_back(C / TractionScale());
    _LoveModulusF.push_back(F / TractionScale());
    _LoveModulusL.push_back(L / TractionScale());
    _LoveModulusN.push_back(N / TractionScale());
    _bulkQualityFactor.push_back(QKappa);
    _shearQualityFactor.push_back(QMu);
  }

  // Work out the layering.
  constexpr auto eps = std::numeric_limits<Real>::epsilon();
  _boundaryRadius.push_back(_r.front());
  _boundaryIndex.push_back(0);
  for (auto i = 1; i < NumberOfKnots() - 1; i++) {
    auto r0 = _r[i - 1];
    auto r1 = _r[i];
    if (std::abs(r1 - r0) < eps * (r0 + r1)) {
      _boundaryRadius.push_back(r0);
      _boundaryIndex.push_back(i);
    }
  }
  _boundaryRadius.push_back(_r.back());
  _boundaryIndex.push_back(NumberOfKnots());

  // Set up the cubic splines.
  for (auto i = 0; i < NumberOfLayers(); i++) {
    auto i0 = _boundaryIndex[i];
    auto i1 = _boundaryIndex[i + 1];
    auto rS = std::next(_r.begin(), i0);
    auto rF = std::next(_r.begin(), i1);
    _densities.push_back(Spline(rS, rF, std::next(_density.begin(), i0)));
    _LoveModuliiA.push_back(
        Spline(rS, rF, std::next(_LoveModulusA.begin(), i0)));
    _LoveModuliiC.push_back(
        Spline(rS, rF, std::next(_LoveModulusC.begin(), i0)));
    _LoveModuliiF.push_back(
        Spline(rS, rF, std::next(_LoveModulusF.begin(), i0)));
    _LoveModuliiL.push_back(
        Spline(rS, rF, std::next(_LoveModulusL.begin(), i0)));
    _LoveModuliiN.push_back(
        Spline(rS, rF, std::next(_LoveModulusN.begin(), i0)));
    _bulkQualityFactors.push_back(
        Spline(rS, rF, std::next(_bulkQualityFactor.begin(), i0)));
    _shearQualityFactors.push_back(
        Spline(rS, rF, std::next(_shearQualityFactor.begin(), i0)));
  }
}

std::function<Real(Real, Int)> DeckModel::Density() const {
  return [this](Real r, Int attribute) { return _densities[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::LoveModulusA() const {
  return
      [this](Real r, Int attribute) { return _LoveModuliiA[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::LoveModulusC() const {
  return
      [this](Real r, Int attribute) { return _LoveModuliiC[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::LoveModulusF() const {
  return
      [this](Real r, Int attribute) { return _LoveModuliiF[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::LoveModulusL() const {
  return
      [this](Real r, Int attribute) { return _LoveModuliiL[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::LoveModulusN() const {
  return
      [this](Real r, Int attribute) { return _LoveModuliiN[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::BulkQualityFactor() const {
  return [this](Real r, Int attribute) {
    return _bulkQualityFactors[attribute - 1](r);
  };
}

std::function<Real(Real, Int)> DeckModel::ShearQualityFactor() const {
  return [this](Real r, Int attribute) {
    return _shearQualityFactors[attribute - 1](r);
  };
}

} // namespace LoveNumbers