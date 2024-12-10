#pragma once

#include "Configure.hpp"
#include "Dimensions.hpp"
#include "LoveNumbers/Configure.hpp"
#include "LoveNumbers/RadialModel.hpp"
#include "RadialModel.hpp"

#include <cmath>
#include <functional>

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
                   Real QKappa, Real QMu,
                   const Dimensions &dimensions = Dimensions())
      : RadialModel(dimensions), _r{r / LengthScale()},
        _rho{rho / DensityScale()}, _A{A / TractionScale()},
        _C{C / TractionScale()}, _F{F / TractionScale()},
        _L{L / TractionScale()}, _N{N / TractionScale()}, _QKappa{QKappa},
        _QMu{QMu} {}

  static HomogeneousModel
  FromFEMOrderAndLengthScale(Real r, Real rho, Real A, Real C, Real F, Real L,
                             Real N, Real QKappa, Real QMu, Int Order,
                             Real characteristicLength,
                             const Dimensions &dimensions = Dimensions());

  static HomogeneousModel IsotropicFromFEMOrderAndLengthScale(
      Real r, Real rho, Real kappa, Real mu, Real QKappa, Real QMu, Int Order,
      Real characteristicLength, const Dimensions &dimensions = Dimensions());

  // Return the number of layers. Override of pure virtual function in base
  // class.
  Int NumberOfLayers() const override;

  // Return the bounding radii of the ith layer. Override of pure
  // virtual function in base class.
  std::pair<Real, Real> LayerRadii(Int i) const override;

  // Return true if ith layer is solid. Override of pure virtual function in
  // base class.
  bool LayerIsSolid(Int i) const override;

  // Return material parameter functions in each layer.
  /*
  std::function<Real(Real)> Rho(Int i) const override;

  std::function<Real(Real)> A(Int i) const override;

  std::function<Real(Real)> C(Int i) const override;

  std::function<Real(Real)> F(Int i) const override;

  std::function<Real(Real)> L(Int i) const override;

  std::function<Real(Real)> N(Int i) const override;

  std::function<Real(Real)> QKappa(Int i) const override;

  std::function<Real(Real)> QMu(Int i) const override;
*/
};

} // namespace LoveNumbers