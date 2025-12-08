#include "LoveNumbers/RadialModelMesh.hpp"

namespace LoveNumbers {

mfem::Array<Int> RadialModelMesh::LayerAttributes() const {
  auto attributes = mfem::Array<Int>(GetRadialModel().NumberOfLayers());
  for (auto i : GetRadialModel().LayerIndices()) {
    attributes[i] = i + 1;
  }
  return attributes;
}

mfem::Array<Int> RadialModelMesh::SolidMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfLayers());
  for (auto i : GetRadialModel().LayerIndices()) {
    marker[i] = GetRadialModel().LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModelMesh::FluidMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfLayers());
  for (auto i : GetRadialModel().LayerIndices()) {
    marker[i] = GetRadialModel().LayerIsSolid(i) ? 0 : 1;
  }
  return marker;
}

mfem::Array<Int> RadialModelMesh::SurfaceMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfBoundaries());
  marker = 0;
  marker[GetRadialModel().NumberOfBoundaries() - 1] = 1;
  return marker;
}

mfem::Array<Int> RadialModelMesh::FluidSolidMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfBoundaries());
  marker = 0;
  for (auto i = 1; i < GetRadialModel().NumberOfLayers(); i++) {
    marker[i] = !GetRadialModel().LayerIsSolid(i - 1) &&
                        GetRadialModel().LayerIsSolid(i)
                    ? 1
                    : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModelMesh::SolidFluidMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfBoundaries());
  marker = 0;
  for (auto i = 1; i < GetRadialModel().NumberOfLayers(); i++) {
    marker[i] = GetRadialModel().LayerIsSolid(i - 1) &&
                        !GetRadialModel().LayerIsSolid(i)
                    ? 1
                    : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModelMesh::SolidSolidMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfBoundaries());
  marker = 0;
  for (auto i = 1; i < GetRadialModel().NumberOfLayers(); i++) {
    marker[i] =
        GetRadialModel().LayerIsSolid(i - 1) && GetRadialModel().LayerIsSolid(i)
            ? 1
            : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModelMesh::FluidFluidMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfBoundaries());
  marker = 0;
  for (auto i = 1; i < GetRadialModel().NumberOfLayers(); i++) {
    marker[i] = !GetRadialModel().LayerIsSolid(i - 1) &&
                        !GetRadialModel().LayerIsSolid(i)
                    ? 1
                    : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModelMesh::AllBoundaryMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfBoundaries());
  marker = 1;
  return marker;
}

mfem::Array<Int> RadialModelMesh::AllBoundaryMarkerCentreExcluded() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfBoundaries());
  marker = 1;
  marker[0] = 0;
  return marker;
}

mfem::Array<Int> RadialModelMesh::CentreMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfBoundaries());
  marker = 0;
  marker[0] = 1;
  return marker;
}

mfem::Array<Int> RadialModelMesh::CentreAndSurfaceMarker() const {
  auto marker = mfem::Array<Int>(GetRadialModel().NumberOfBoundaries());
  marker = 0;
  marker[0] = 1;
  marker[GetRadialModel().NumberOfBoundaries() - 1] = 1;
  return marker;
}

void RadialModelMesh::BuildMesh(Real characteristicLengthScale) {
  // Clear the mesh in case already set up.
  _mesh.Clear();

  // Work out number of elements in each layer.
  auto elementsInLayer = std::vector<Int>();
  for (auto i : GetRadialModel().LayerIndices()) {
    auto [r1, r2] = GetRadialModel().LayerRadii(i);
    elementsInLayer.push_back(
        std::max<int>(1, std::round((r2 - r1) / characteristicLengthScale)));
  }

  // Initiallise the mesh.
  auto numberOfElements =
      std::accumulate(elementsInLayer.begin(), elementsInLayer.end(), 0);
  auto numberOfVertices = numberOfElements + 1;
  auto numberOfBoundaryElements = GetRadialModel().NumberOfLayers() + 1;
  _mesh = mfem::Mesh(1, numberOfVertices, numberOfElements,
                     numberOfBoundaryElements);

  // Loop over the layers building up the mesh.
  auto vertex1 = 0;
  auto layerAttribute = 1;
  auto boundaryAttribute = 1;
  for (auto i : GetRadialModel().LayerIndices()) {
    // Get number of vertices and spacing.
    auto [r1, r2] = GetRadialModel().LayerRadii(i);
    numberOfVertices = elementsInLayer[i] + 1;
    auto dr = (r2 - r1) / static_cast<Real>(numberOfVertices);

    // Add in the vertices.
    for (auto j = layerAttribute > 1 ? 1 : 0; j <= numberOfVertices; j++) {
      auto r = r1 + j * dr;
      _mesh.AddVertex(r);
    }

    // Add in the elements.
    auto vertex2 = vertex1 + numberOfVertices;
    for (auto j = vertex1; j < vertex2; j++) {
      _mesh.AddSegment(j, j + 1, layerAttribute);
    }

    // Label the boundary elements.
    if (layerAttribute == 1) {
      _mesh.AddBdrPoint(vertex1, boundaryAttribute++);
      _mesh.AddBdrPoint(vertex2, boundaryAttribute++);
    } else {
      _mesh.AddBdrPoint(vertex2, boundaryAttribute++);
    }
    vertex1 = vertex2;

    layerAttribute++;
  }

  // Finalise the mesh construction.
  _mesh.FinalizeMesh();

  // Set the finite element spaces.
  _L2Space = std::make_unique<mfem::FiniteElementSpace>(&_mesh, _L2.get());
  _H1Space = std::make_unique<mfem::FiniteElementSpace>(&_mesh, _H1.get());

  // Compute the surface gravity and moment of inertia factors.
  ComputeSurfaceGravityAndMomentOfInertiaFactor();

  // Compute the interior gravitational potential and acceleration.
  ComputeGravitationalPotential();
}

void RadialModelMesh::ComputeSurfaceGravityAndMomentOfInertiaFactor() {
  auto densityCoefficient = DensityCoefficient();
  auto density = mfem::GridFunction(&L2Space());
  density.ProjectCoefficient(densityCoefficient);
  auto surfaceRadiusSquared = std::pow(GetRadialModel().SurfaceRadius(), 2);
  {
    auto kernel = RadialCoefficient(
        [](auto r, auto attribute) { return std::pow(r, 2); });
    auto b = mfem::LinearForm(&L2Space());
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(kernel));
    b.Assemble();
    auto factor = 4 * std::numbers::pi_v<Real> *
                  GetRadialModel().GravitationalConstant() /
                  surfaceRadiusSquared;
    _surfaceGravity = factor * (b(density));
  }

  {
    auto kernel = RadialCoefficient(
        [](auto r, auto attribute) { return std::pow(r, 4); });
    auto b = mfem::LinearForm(&L2Space());
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(kernel));
    b.Assemble();
    auto factor = 8 * std::numbers::pi_v<Real> *
                  GetRadialModel().GravitationalConstant() /
                  (3 * _surfaceGravity * std::pow(surfaceRadiusSquared, 2));
    _momentOfInertiaFactor = factor * (b * density);
  }
}

void RadialModelMesh::ComputeGravitationalPotential() {
  using namespace mfem;

  // Set up the linear form.
  auto densityTimesRadiusSquared =
      RadialCoefficient([this](auto r, auto attribute) {
        const auto factor = -(4 * std::numbers::pi_v<Real> *
                              GetRadialModel().GravitationalConstant());
        return factor * GetRadialModel().Density()(r, attribute) * r * r;
      });
  auto b = LinearForm(&H1Space());
  b.AddDomainIntegrator(new DomainLFIntegrator(densityTimesRadiusSquared));
  b.Assemble();

  // Set up the bilinear form.
  auto a = BilinearForm(&H1Space());
  auto radiusSquared =
      RadialCoefficient([](auto r, auto attribute) { return r * r; });
  a.AddDomainIntegrator(new DiffusionIntegrator(radiusSquared));
  auto DtN = ConstantCoefficient(GetRadialModel().SurfaceRadius());
  // auto radius = RadialCoefficient([](auto r, auto attribute) { return r; });
  auto surfaceMarker = SurfaceMarker();
  a.AddBoundaryIntegrator(new BoundaryMassIntegrator(DtN), surfaceMarker);
  a.Assemble();

  // Set up the solution vector with appropriate boundary values.
  _gravitationalPotential = std::make_unique<GridFunction>(&H1Space());
  *_gravitationalPotential = 0;

  // Set up the linear system.
  OperatorPtr A;
  Vector B, X;
  auto ess_tdof_list = Array<Int>();
  a.FormLinearSystem(ess_tdof_list, *_gravitationalPotential, b, A, X, B);

  // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
  // GSSmoother M((SparseMatrix &)(*A));
  GSSmoother M(dynamic_cast<SparseMatrix &>(*A));
  PCG(*A, M, B, X, 0, H1Space().GetNDofs() * 2, 1e-12, 0.0);
  a.RecoverFEMSolution(X, b, *_gravitationalPotential);

  // Form the gravitational acceleration.
  _gravitationalAcceleration = std::make_unique<GridFunction>(&L2Space());
  _gravitationalPotential->GetDerivative(1, 0, *_gravitationalAcceleration);

  // Form the density derivative.
  _densityDerivative = std::make_unique<GridFunction>(&L2Space());
}

void RadialModelMesh::PrintMesh(const std::string &mesh_file) {
  std::ofstream ofs(mesh_file);
  ofs.precision(8);
  _mesh.Print(ofs);
  ofs.close();
}

void RadialModelMesh::Write(const mfem::GridFunction &f,
                            const std::string &file, Real scale) const {
  using namespace mfem;
  auto fout = std::ofstream(file);
  auto point = Vector(1);
  auto values = Vector(Order() + 1);
  auto pairs = std::vector<std::pair<Real, Real>>(Order() + 1);
  auto *fes = f.FESpace();
  for (auto i = 0; i < GetMesh().GetNE(); i++) {
    auto *el = fes->GetFE(i);
    auto *eltrans = fes->GetElementTransformation(i);
    auto &ir = el->GetNodes();
    f.GetValues(*eltrans, ir, values);
    for (auto j = 0; j < ir.GetNPoints(); j++) {
      auto &ip = ir.IntPoint(j);
      eltrans->SetIntPoint(&ip);
      eltrans->Transform(ip, point);
      pairs[j] = {point[0], values[j]};
    }
    std::ranges::sort(pairs,
                      [](auto p1, auto p2) { return p1.first < p2.first; });
    for (auto [r, v] : pairs) {
      fout << r * GetRadialModel().LengthScale() << " " << v * scale
           << std::endl;
    }
  }
}

void RadialModelMesh::WriteDerivative(const mfem::GridFunction &f,
                                      const std::string &file,
                                      Real scale) const {
  using namespace mfem;
  auto fout = std::ofstream(file);
  auto point = Vector(1);
  auto grad = DenseMatrix(1, Order() + 1);
  auto pairs = std::vector<std::pair<Real, Real>>(Order() + 1);
  auto *fes = f.FESpace();
  for (auto i = 0; i < GetMesh().GetNE(); i++) {
    auto *el = fes->GetFE(i);
    auto *eltrans = fes->GetElementTransformation(i);
    auto &ir = el->GetNodes();
    f.GetGradients(*eltrans, ir, grad);
    for (auto j = 0; j < ir.GetNPoints(); j++) {
      auto &ip = ir.IntPoint(j);
      eltrans->SetIntPoint(&ip);
      eltrans->Transform(ip, point);
      pairs[j] = {point[0], grad(0, j)};
    }
    std::ranges::sort(pairs,
                      [](auto p1, auto p2) { return p1.first < p2.first; });
    for (auto [r, v] : pairs) {
      fout << r * GetRadialModel().LengthScale() << " "
           << v * scale / GetRadialModel().LengthScale() << std::endl;
    }
  }
}

void RadialModelMesh::Write(mfem::Coefficient &f, const std::string &file,
                            Real scale) const {
  // Loop over the mesh storing values.
  auto fout = std::ofstream(file);
  auto &fes = H1Space();
  auto point = mfem::Vector(1);
  auto pairs = std::vector<std::pair<Real, Real>>(Order() + 1);
  for (auto i = 0; i < GetMesh().GetNE(); i++) {
    auto *el = fes.GetFE(i);
    auto *eltrans = fes.GetElementTransformation(i);
    auto &ir = el->GetNodes();
    for (auto j = 0; j < ir.GetNPoints(); j++) {
      auto &ip = ir.IntPoint(j);
      eltrans->SetIntPoint(&ip);
      eltrans->Transform(ip, point);
      auto value = f.Eval(*eltrans, ip);
      pairs[j] = {point[0], value};
    }
    std::ranges::stable_sort(
        pairs, [](auto p1, auto p2) { return p1.first < p2.first; });
    for (auto [r, v] : pairs) {
      fout << r * GetRadialModel().LengthScale() << " " << v * scale
           << std::endl;
    }
  }
}

void RadialModelMesh::WriteDerivative(mfem::Coefficient &f,
                                      const std::string &file,
                                      Real scale) const {
  using namespace mfem;
  auto fout = std::ofstream(file);
  auto &fes = H1Space();
  auto lval = Vector(Order() + 1);
  auto dshape = DenseMatrix(Order() + 1, 1);
  auto gh = Vector(1);
  auto grad = Vector(1);
  auto point = Vector(1);
  auto pairs = std::vector<std::pair<Real, Real>>(Order() + 1);
  for (auto i = 0; i < GetMesh().GetNE(); i++) {
    auto *el = fes.GetFE(i);
    auto *tr = fes.GetElementTransformation(i);
    el->Project(f, *tr, lval);
    auto &ir = el->GetNodes();
    for (auto j = 0; j < ir.GetNPoints(); j++) {
      auto &ip = ir.IntPoint(j);
      tr->SetIntPoint(&ip);
      tr->Transform(ip, point);
      el->CalcDShape(ip, dshape);
      dshape.MultTranspose(lval, gh);
      tr->InverseJacobian().MultTranspose(gh, grad);
      pairs[j] = {point[0], grad[0]};
    }
    // Sort and write the values.
    std::ranges::stable_sort(
        pairs, [](auto p1, auto p2) { return p1.first < p2.first; });
    for (auto [r, v] : pairs) {
      fout << r * GetRadialModel().LengthScale() << " "
           << v * scale / GetRadialModel().LengthScale() << std::endl;
    }
  }
}

}  // namespace LoveNumbers