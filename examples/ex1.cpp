

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numbers>
#include <ranges>
#include <string>
#include <vector>

#include "LoveNumbers/LoveNumbers.hpp"

using namespace LoveNumbers;
using namespace mfem;

int main(int argc, char *argv[]) {

  auto dimensions = Dimensions();
  Int order = 10;
  auto characteristicLengthScale = Real{1};
  auto model = DeckModel(Dimensions(), order, characteristicLengthScale,
                         "../data/prem.200.no");

  // Get pointer to finite element space.
  auto &mesh = model.Mesh();
  auto *fes = model.H1Space();

  // Set Dirichlet condition at the centre.
  auto ess_tdof_list = Array<Int>();
  // auto ess_bdr = model.CentreMarker();
  // fes->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  // Set up the linear form.
  auto rhoTimesRadiusSquared =
      RadialModelCoefficient([&](auto r, auto attribute) {
        return model.Rho()(r, attribute) * r * r;
      });
  auto b = LinearForm(fes);
  b.AddDomainIntegrator(new DomainLFIntegrator(rhoTimesRadiusSquared));
  b.Assemble();
  b *= -(4 * std::numbers::pi_v<Real> * model.GravitationalConstant());

  // Set up the bilinear form.
  auto a = BilinearForm(fes);
  auto radiusSquared =
      RadialModelCoefficient([&](auto r, auto attribute) { return r * r; });
  a.AddDomainIntegrator(new DiffusionIntegrator(radiusSquared));
  auto DtN = ConstantCoefficient(model.SurfaceRadius());
  auto surfaceMarker = model.SurfaceMarker();
  a.AddBoundaryIntegrator(new BoundaryMassIntegrator(DtN), surfaceMarker);
  a.Assemble();

  // Set up the solution vector with appropriate boundary values.
  GridFunction x(fes);
  x = 1.0;

  // Set up the linear system.
  OperatorPtr A;
  Vector B, X;
  a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

  // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
  GSSmoother M((SparseMatrix &)(*A));
  PCG(*A, M, B, X, 1, fes->GetNDofs() * 2, 1e-12, 0.0);

  a.RecoverFEMSolution(X, b, x);

  std::ofstream mesh_ofs("refined.mesh");
  mesh_ofs.precision(8);
  mesh.Print(mesh_ofs);
  std::ofstream sol_ofs("sol.gf");
  sol_ofs.precision(8);
  x.Save(sol_ofs);

  auto fout = std::ofstream("ex1.out");

  auto point = Vector();
  auto values = Array<Real>();
  for (auto i = 0; i < mesh.GetNE(); i++) {
    auto *el = fes->GetFE(i);
    auto *eltrans = fes->GetElementTransformation(i);
    auto &ir = el->GetNodes();
    for (auto j = 0; j < ir.GetNPoints(); j++) {
      auto &ip = ir.IntPoint(j);
      eltrans->SetIntPoint(&ip);
      eltrans->Transform(ip, point);
      fout << point[0] * model.LengthScale() << " "
           << x.GetValue(i, ip) * model.PotentialScale() << std::endl;
    }
  }
}
