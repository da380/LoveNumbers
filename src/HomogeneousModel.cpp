#include "LoveNumbers/HomogeneousModel.hpp"

namespace LoveNumbers {

Int HomogeneousModel::NumberOfLayers() const { return 1; }

std::pair<Real, Real> HomogeneousModel::LayerRadii(Int i) const {
  return {0, _r};
}

bool HomogeneousModel::LayerIsSolid(Int i) const { return _L > 0; }

/*
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
*/

} // namespace LoveNumbers