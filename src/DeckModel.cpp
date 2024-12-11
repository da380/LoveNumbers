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
    _rho.push_back(rho / DensityScale());
    _A.push_back(A / TractionScale());
    _C.push_back(C / TractionScale());
    _F.push_back(F / TractionScale());
    _L.push_back(L / TractionScale());
    _N.push_back(N / TractionScale());
    _QKappa.push_back(QKappa);
    _QMu.push_back(QMu);
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
    _rhoSplines.push_back(Spline(rS, rF, std::next(_rho.begin(), i0)));
    _ASplines.push_back(Spline(rS, rF, std::next(_A.begin(), i0)));
    _CSplines.push_back(Spline(rS, rF, std::next(_C.begin(), i0)));
    _FSplines.push_back(Spline(rS, rF, std::next(_F.begin(), i0)));
    _LSplines.push_back(Spline(rS, rF, std::next(_L.begin(), i0)));
    _NSplines.push_back(Spline(rS, rF, std::next(_N.begin(), i0)));
    _QKappaSplines.push_back(Spline(rS, rF, std::next(_QKappa.begin(), i0)));
    _QMuSplines.push_back(Spline(rS, rF, std::next(_QMu.begin(), i0)));
  }
}

std::function<Real(Real, Int)> DeckModel::Rho() const {
  return
      [this](Real r, Int attribute) { return _rhoSplines[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::A() const {
  return [this](Real r, Int attribute) { return _ASplines[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::C() const {
  return [this](Real r, Int attribute) { return _CSplines[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::F() const {
  return [this](Real r, Int attribute) { return _FSplines[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::L() const {
  return [this](Real r, Int attribute) { return _LSplines[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::N() const {
  return [this](Real r, Int attribute) { return _NSplines[attribute - 1](r); };
}

std::function<Real(Real, Int)> DeckModel::QKappa() const {
  return [this](Real r, Int attribute) {
    return _QKappaSplines[attribute - 1](r);
  };
}

std::function<Real(Real, Int)> DeckModel::QMu() const {
  return
      [this](Real r, Int attribute) { return _QMuSplines[attribute - 1](r); };
}

} // namespace LoveNumbers