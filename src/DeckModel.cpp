#include "LoveNumbers/DeckModel.hpp"

namespace LoveNumbers {

DeckModel DeckModel::FromMaximumElementSize(const std::string &fileName,
                                            Real maximumElementSize) {
  auto model = DeckModel(fileName);
  model.BuildMesh(maximumElementSize);
  return model;
}

DeckModel DeckModel::FromMaximumDegree(const std::string &fileName,
                                       Int maximumDegree) {
  auto model = DeckModel(fileName);
  auto maximumElementSize =
      model.JeanLength(maximumDegree) / static_cast<Real>(5);
  model.BuildMesh(maximumElementSize);
  return model;
}

Int DeckModel::NumberOfLayers() const { return _boundaryIndices.size(); }

std::pair<Real, Real> DeckModel::LayerRadii(Int i) const {
  return _boundaryRadii[i];
}

bool DeckModel::LayerIsSolid(Int i) const { return _layerSolid[i]; }

Int DeckModel::NumberOfKnots() const { return _r.size(); }

std::function<Real(Real)> DeckModel::Rho(Int i) const {
  return std::function<Real(Real)>(*_rhoSplines[i]);
};

std::function<Real(Real)> DeckModel::A(Int i) const {
  return std::function<Real(Real)>(*_ASplines[i]);
};

std::function<Real(Real)> DeckModel::C(Int i) const {
  return std::function<Real(Real)>(*_CSplines[i]);
};

std::function<Real(Real)> DeckModel::F(Int i) const {
  return std::function<Real(Real)>(*_FSplines[i]);
};

std::function<Real(Real)> DeckModel::L(Int i) const {
  return std::function<Real(Real)>(*_LSplines[i]);
};

std::function<Real(Real)> DeckModel::N(Int i) const {
  return std::function<Real(Real)>(*_NSplines[i]);
};

std::function<Real(Real)> DeckModel::QKappa(Int i) const {
  return std::function<Real(Real)>(*_QKappaSplines[i]);
};

std::function<Real(Real)> DeckModel::QMu(Int i) const {
  return std::function<Real(Real)>(*_QMuSplines[i]);
};

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
  Int iL = 0;
  Real rL = 0;
  for (auto iU = 1; iU < NumberOfKnots(); iU++) {
    auto r0 = _r[iU - 1];
    auto r1 = _r[iU];
    if (std::abs(r1 - r0) < eps * (r0 + r1)) {
      _boundaryIndices.push_back({iL, iU});
      _boundaryRadii.push_back({rL, r0});
      _layerSolid.push_back(
          std::all_of(&_L[iL], &_L[iU], [](auto L) { return L > 0; }));
      iL = iU;
      rL = r0;
    }
  }

  // Set up the cubic splines.
  for (auto i = 0; i < NumberOfLayers(); i++) {
    auto [i0, i1] = _boundaryIndices[i];
    auto rS = std::next(_r.begin(), i0);
    auto rF = std::next(_r.begin(), i1);
    _rhoSplines.push_back(
        std::make_unique<Spline>(rS, rF, std::next(_rho.begin(), i0)));
    _ASplines.push_back(
        std::make_unique<Spline>(rS, rF, std::next(_A.begin(), i0)));
    _CSplines.push_back(
        std::make_unique<Spline>(rS, rF, std::next(_C.begin(), i0)));
    _FSplines.push_back(
        std::make_unique<Spline>(rS, rF, std::next(_F.begin(), i0)));
    _LSplines.push_back(
        std::make_unique<Spline>(rS, rF, std::next(_L.begin(), i0)));
    _NSplines.push_back(
        std::make_unique<Spline>(rS, rF, std::next(_N.begin(), i0)));
    _QKappaSplines.push_back(
        std::make_unique<Spline>(rS, rF, std::next(_QKappa.begin(), i0)));
    _QMuSplines.push_back(
        std::make_unique<Spline>(rS, rF, std::next(_QMu.begin(), i0)));
  }
}

} // namespace LoveNumbers