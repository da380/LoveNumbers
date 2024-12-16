#pragma once

#include "LoveNumbers/Configure.hpp"
#include "mfem.hpp"
#include <cassert>
#include <concepts>
#include <iostream>

namespace LoveNumbers {

// Coefficient class for piece-wise smooth functions on a radial mesh.
template <typename Function>
  requires requires() {
    requires std::regular_invocable<Function, Real, Int>;
    requires std::convertible_to<std::invoke_result_t<Function, Real, Int>,
                                 Real>;
  }
class RadialCoefficient : public mfem::Coefficient {
private:
  Function _f;

public:
  RadialCoefficient() = default;
  RadialCoefficient(Function &&f) : _f{f} {}

  Real Eval(mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override {
    Real data[3];
    auto x = mfem::Vector(data, 3);
    T.Transform(ip, x);
    auto r = x.Norml2();
    auto attribute = T.Attribute;
    return _f(r, attribute);
  }
};

template <typename Function>
  requires requires() {
    requires std::regular_invocable<Function, Real, Int>;
    requires std::convertible_to<std::invoke_result_t<Function, Real, Int>,
                                 Real>;
  }
class RadialDiscontinuityCoefficient : public mfem::Coefficient {
private:
  Function _f;

public:
  RadialDiscontinuityCoefficient() = default;
  RadialDiscontinuityCoefficient(Function &&f) : _f{f} {}

  Real Eval(mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override {
    Real data[3];
    auto attribute = T.Attribute;
    auto x = mfem::Vector(data, 3);
    T.Transform(ip, x);
    auto r = x.Norml2();
    return _f(r, attribute) - _f(r, attribute - 1);
  }
};

} // namespace LoveNumbers