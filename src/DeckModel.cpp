#include "LoveNumbers/DeckModel.hpp"

namespace LoveNumbers {

DeckModel::DeckModel(const std::string &fileName, Real maximumElementSize) {
  // Check the file exists.
  if (!std::filesystem::exists(fileName)) {
    throw std::runtime_error("deck file does not exist");
  }

  // Open the file.
  auto modelFile = std::ifstream(fileName);

  // Skip the header information.
  for (auto i = 0; i < 3; i++) {
    auto line = std::string();
    std::getline(modelFile, line);
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

  // Build the FEM mesh.
  BuildMesh(maximumElementSize);
}

} // namespace LoveNumbers