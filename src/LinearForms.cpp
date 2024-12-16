#include "LoveNumbers/LinearForms.hpp"

namespace LoveNumbers {

void TestIntegrator::AssembleRHSElementVect(const mfem::FiniteElement &el,
                                            mfem::ElementTransformation &Tr,
                                            mfem::Vector &elvect) {

  using namespace mfem;
  int dof = el.GetDof();

  shape.SetSize(dof); // vector of size dof
  elvect.SetSize(dof);
  elvect = 0.0;

  auto radius = mfem::SphericalRadialCoefficient();

  const IntegrationRule *ir = IntRule;
  if (ir == NULL) {
    ir = &IntRules.Get(el.GetGeomType(), oa * el.GetOrder() + ob);
  }

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);

    Tr.SetIntPoint(&ip);
    auto r = radius.Eval(Tr, ip);
    real_t val = Tr.Weight() * Q.Eval(Tr, ip) * r * r;

    el.CalcPhysShape(Tr, shape);

    add(elvect, ip.weight * val, shape, elvect);
  }
}

void TestIntegrator::AssembleDeltaElementVect(
    const mfem::FiniteElement &fe, mfem::ElementTransformation &Trans,
    mfem::Vector &elvect) {
  MFEM_ASSERT(delta != NULL, "coefficient must be DeltaCoefficient");
  elvect.SetSize(fe.GetDof());
  fe.CalcPhysShape(Trans, elvect);
  elvect *= delta->EvalDelta(Trans, Trans.GetIntPoint());
}

} // namespace LoveNumbers