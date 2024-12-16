#pragma once

#include "Coefficients.hpp"
#include "LoveNumbers/Configure.hpp"
#include "mfem.hpp"

namespace LoveNumbers {

/// Class for domain integration $ L(v) := (f, v) $
class TestIntegrator : public mfem::DeltaLFIntegrator {
  mfem::Vector shape;
  mfem::Coefficient &Q;
  int oa, ob;

public:
  /// Constructs a domain integrator with a given Coefficient
  TestIntegrator(mfem::Coefficient &QF, int a = 2, int b = 0)
      // the old default was a = 1, b = 1
      // for simple elliptic problems a = 2, b = -2 is OK
      : mfem::DeltaLFIntegrator(QF), Q(QF), oa(a), ob(b) {}

  /// Constructs a domain integrator with a given Coefficient
  TestIntegrator(mfem::Coefficient &QF, const mfem::IntegrationRule *ir)
      : mfem::DeltaLFIntegrator(QF, ir), Q(QF), oa(1), ob(1) {}

  bool SupportsDevice() const override { return false; }

  /** Given a particular Finite Element and a transformation (Tr)
      computes the element right hand side element vector, elvect. */
  void AssembleRHSElementVect(const mfem::FiniteElement &el,
                              mfem::ElementTransformation &Tr,
                              mfem::Vector &elvect) override;

  void AssembleDeltaElementVect(const mfem::FiniteElement &fe,
                                mfem::ElementTransformation &Trans,
                                mfem::Vector &elvect) override;

  using LinearFormIntegrator::AssembleRHSElementVect;
};

} // namespace LoveNumbers