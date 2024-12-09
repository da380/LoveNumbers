#include "LoveNumbers/HomogeneousModel.hpp"

namespace LoveNumbers {

HomogeneousModel HomogeneousModel::FromIsotropicWaveSpeeds(
    Real r, Real rho, Real vp, Real vs, Real QKappa, Real QMu, Real lengthScale,
    Real massScale, Real timeScale) {
  auto A = rho * vp * vp;
  auto L = rho * vs * vs;
  auto F = A - 2 * L;
  return HomogeneousModel(r, rho, A, A, F, L, L, QKappa, QMu, lengthScale,
                          massScale, timeScale);
}

HomogeneousModel HomogeneousModel::FromIsotropicWaveSpeeds(Real r, Real rho,
                                                           Real vp, Real vs,
                                                           Real QKappa,
                                                           Real QMu) {
  return HomogeneousModel::FromIsotropicWaveSpeeds(
      r, rho, vp, vs, QKappa, QMu, _LENGTH_SCALE, _MASS_SCALE, _TIME_SCALE);
}

HomogeneousModel
HomogeneousModel::FromIsotropicModulii(Real r, Real rho, Real kappa, Real mu,
                                       Real QKappa, Real QMu, Real lengthScale,
                                       Real massScale, Real timeScale) {
  auto fourThirds = static_cast<Real>(4) / static_cast<Real>(3);
  auto vp = std::sqrt((kappa + fourThirds * mu) / rho);
  auto vs = std::sqrt(mu / rho);
  return HomogeneousModel::FromIsotropicWaveSpeeds(
      r, rho, vp, vs, QKappa, QMu, lengthScale, massScale, timeScale);
}

HomogeneousModel HomogeneousModel::FromIsotropicModulii(Real r, Real rho,
                                                        Real kappa, Real mu,
                                                        Real QKappa, Real QMu) {
  return HomogeneousModel::FromIsotropicModulii(
      r, rho, kappa, mu, QKappa, QMu, _LENGTH_SCALE, _MASS_SCALE, _TIME_SCALE);
}

Int HomogeneousModel::NumberOfLayers() const { return 1; }

std::pair<Real, Real> HomogeneousModel::LayerRadii(Int i) const {
  return {0, _r};
}

bool HomogeneousModel::LayerIsSolid(Int i) const { return _L > 0; }

std::function<Real(Real)> HomogeneousModel::Rho(Int i) const {
  return std::function<Real(Real)>([this](auto r) { return _rho; });
};

std::function<Real(Real)> HomogeneousModel::A(Int i) const {
  return std::function<Real(Real)>([this](auto r) { return _A; });
};

std::function<Real(Real)> HomogeneousModel::C(Int i) const {
  return std::function<Real(Real)>([this](auto r) { return _C; });
};

std::function<Real(Real)> HomogeneousModel::F(Int i) const {
  return std::function<Real(Real)>([this](auto r) { return _F; });
};

std::function<Real(Real)> HomogeneousModel::L(Int i) const {
  return std::function<Real(Real)>([this](auto r) { return _L; });
};

std::function<Real(Real)> HomogeneousModel::N(Int i) const {
  return std::function<Real(Real)>([this](auto r) { return _N; });
};

std::function<Real(Real)> HomogeneousModel::QKappa(Int i) const {
  return std::function<Real(Real)>([this](auto r) { return _QKappa; });
};

std::function<Real(Real)> HomogeneousModel::QMu(Int i) const {
  return std::function<Real(Real)>([this](auto r) { return _QMu; });
};

} // namespace LoveNumbers