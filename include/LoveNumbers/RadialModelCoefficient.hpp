#pragma once

#include "LoveNumbers/Configure.hpp"
#include "mfem.hpp"
#include <complex>
#include <concepts>
#include <iostream>

namespace LoveNumbers {

template <typename Function>
  requires requires() {
    requires std::regular_invocable<Function, Real, Int>;
    requires std::convertible_to<std::invoke_result_t<Function, Real, Int>,
                                 Real>;
  }
class RadialModelCoefficient : public mfem::Coefficient {
private:
  Function _f;

public:
  RadialModelCoefficient() = default;
  RadialModelCoefficient(Function &&f) : _f{f} {}

  Real Eval(mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override {

    // Map integration point to spatial point.
    Real data[3];
    auto x = mfem::Vector(data, 3);
    T.Transform(ip, x);

    // Evaluate the function of radius and the attribute number.
    auto r = x.Norml2();
    auto attribute = T.Attribute;
    return _f(r, attribute);
  }
};

} // namespace LoveNumbers