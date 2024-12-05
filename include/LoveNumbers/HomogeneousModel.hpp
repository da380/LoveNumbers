#pragma once

#include "Configure.hpp"
#include "LoveNumbers/Configure.hpp"
#include "LoveNumbers/RadialModel.hpp"
#include "RadialModel.hpp"

namespace LoveNumbers {

class HomogeneousModel : public RadialModel {
private:
  Real _r;
  Real _rho;
  Real _A;
  Real _C;
  Real _F;
  Real _L;
  Real _N;
  Real _QKappa;
  Real _QMu;

public:
  HomogeneousModel() = delete;

  HomogeneousModel(Real r, Real rho, Real A, Real C, Real F, Real L, Real N,
                   Real QKappa, Real QMu, Real lengthScale, Real massScale,
                   Real timeScale)
      : RadialModel(lengthScale, massScale, timeScale), _r{r / LengthScale()},
        _rho{rho / DensityScale()}, _A{A / TractionScale()},
        _C{C / TractionScale()}, _F{F / TractionScale()},
        _L{L / TractionScale()}, _N{N / TractionScale()}, _QKappa{QKappa},
        _QMu{QMu} {}

  HomogeneousModel(Real r, Real rho, Real A, Real C, Real F, Real L, Real N,
                   Real QKappa, Real QMu)
      : HomogeneousModel(r, rho, A, C, F, L, N, QKappa, QMu, _LENGTH_SCALE,
                         _MASS_SCALE, _TIME_SCALE) {}

  // Return the number of layers. Override of pure virtual function in base
  // class.
  Int NumberOfLayers() const override { return 1; }

  // Return the bounding radii of the ith layer. Override of pure
  // virtual function in base class.
  std::pair<Real, Real> LayerRadii(Int i) const override { return {0, _r}; }

  // Return true if ith layer is solid. Override of pure virtual function in
  // base class.
  bool LayerIsSolid(Int i) const override { return _L > 0; }

  // Return material parameter functions in each layer.
  std::function<Real(Real)> Rho(Int i) const override {
    return std::function<Real(Real)>([this](auto r) { return _rho; });
  };

  std::function<Real(Real)> A(Int i) const override {
    return std::function<Real(Real)>([this](auto r) { return _A; });
  };

  std::function<Real(Real)> C(Int i) const override {
    return std::function<Real(Real)>([this](auto r) { return _C; });
  };

  std::function<Real(Real)> F(Int i) const override {
    return std::function<Real(Real)>([this](auto r) { return _F; });
  };

  std::function<Real(Real)> L(Int i) const override {
    return std::function<Real(Real)>([this](auto r) { return _L; });
  };

  std::function<Real(Real)> N(Int i) const override {
    return std::function<Real(Real)>([this](auto r) { return _N; });
  };

  std::function<Real(Real)> QKappa(Int i) const override {
    return std::function<Real(Real)>([this](auto r) { return _QKappa; });
  };

  std::function<Real(Real)> QMu(Int i) const override {
    return std::function<Real(Real)>([this](auto r) { return _QMu; });
  };
};

} // namespace LoveNumbers