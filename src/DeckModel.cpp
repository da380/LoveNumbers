#include "LoveNumbers/DeckModel.hpp"
#include <memory>

namespace LoveNumbers {

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

void DeckModel::WriteModelToFile(const std::string &fileName) {

  auto modelFile = std::ofstream(fileName);

  // Write the header information.
  for (auto line : _header) {
    modelFile << line << "\n";
  }

  // Write the values
  for (auto i = 0; i < NumberOfKnots(); i++) {

    auto r = _r[i] * LengthScale();
    auto rho = _rho[i] * DensityScale();
    auto A = _A[i] * TractionScale();
    auto C = _C[i] * TractionScale();
    auto F = _F[i] * TractionScale();
    auto L = _L[i] * TractionScale();
    auto N = _N[i] * TractionScale();
    auto QKappa = _QKappa[i];
    auto QMu = _QMu[i];

    auto vpv = std::sqrt(C / rho);
    auto vsv = std::sqrt(L / rho);
    auto vph = std::sqrt(A / rho);
    auto vsh = std::sqrt(N / rho);
    auto eta = F / (A - 2 * L);

    modelFile << std::format("{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} "
                             "{:10.3f} {:10.3f} {:10.3f} {:10.3f}\n",
                             r, rho, vpv, vsv, QKappa, QMu, vph, vsh, eta);
  }
}

} // namespace LoveNumbers