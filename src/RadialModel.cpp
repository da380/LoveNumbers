#include "LoveNumbers/RadialModel.hpp"

namespace LoveNumbers {

mfem::Array<Int> RadialModel::LayerAttributes() const {
  auto attributes = mfem::Array<Int>(NumberOfLayers());
  for (auto i : LayerIndices()) {
    attributes[i] = i + 1;
  }
  return attributes;
}

mfem::Array<Int> RadialModel::SolidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfLayers());
  for (auto i : LayerIndices()) {
    marker[i] = LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModel::FluidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfLayers());
  for (auto i : LayerIndices()) {
    marker[i] = LayerIsSolid(i) ? 0 : 1;
  }
  return marker;
}

mfem::Array<Int> RadialModel::SurfaceMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  marker = 0;
  marker[NumberOfBoundaries() - 1] = 1;
  return marker;
}

mfem::Array<Int> RadialModel::FluidSolidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  marker = 0;
  for (auto i = 1; i < NumberOfLayers(); i++) {
    marker[i] = !LayerIsSolid(i - 1) && LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModel::SolidFluidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  marker = 0;
  for (auto i = 1; i < NumberOfLayers(); i++) {
    marker[i] = LayerIsSolid(i - 1) && !LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModel::SolidSolidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  marker = 0;
  for (auto i = 1; i < NumberOfLayers(); i++) {
    marker[i] = LayerIsSolid(i - 1) && LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModel::FluidFluidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  marker = 0;
  for (auto i = 1; i < NumberOfLayers(); i++) {
    marker[i] = !LayerIsSolid(i - 1) && !LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModel::AllBoundaryMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  marker = 1;
  return marker;
}

mfem::Array<Int> RadialModel::AllBoundaryMarkerCentreExcluded() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  marker = 1;
  marker[0] = 0;
  return marker;
}

mfem::Array<Int> RadialModel::CentreMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  marker = 0;
  marker[0] = 1;
  return marker;
}

mfem::Array<Int> RadialModel::CentreAndSurfaceMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  marker = 0;
  marker[0] = 1;
  marker[NumberOfBoundaries() - 1] = 1;
  return marker;
}

Int RadialModel::NumberOfBoundaries() const { return NumberOfLayers() + 1; }

Real RadialModel::JeanLength(Int degree) const {
  return 2 * std::numbers::pi_v<Real> * SurfaceRadius() /
         static_cast<Real>(degree + 1);
}

Real RadialModel::SurfaceRadius() const {
  return LayerRadii(NumberOfLayers() - 1).second;
}

void RadialModel::BuildMesh(Real characteristicLengthScale) {

  // Clear the mesh in case already set up.
  _mesh.Clear();

  // Work out number of elements in each layer.
  auto elementsInLayer = std::vector<Int>();
  for (auto i : LayerIndices()) {
    auto [r1, r2] = LayerRadii(i);
    elementsInLayer.push_back(
        std::max<int>(1, std::round((r2 - r1) / characteristicLengthScale)));
  }

  // Initiallise the mesh.
  auto numberOfElements =
      std::accumulate(elementsInLayer.begin(), elementsInLayer.end(), 0);
  auto numberOfVertices = numberOfElements + 1;
  auto numberOfBoundaryElements = NumberOfLayers() + 1;
  _mesh = mfem::Mesh(1, numberOfVertices, numberOfElements,
                     numberOfBoundaryElements);

  // Loop over the layers building up the mesh.
  auto vertex1 = 0;
  auto layerAttribute = 1;
  auto boundaryAttribute = 1;
  for (auto i : LayerIndices()) {

    // Get number of vertices and spacing.
    auto [r1, r2] = LayerRadii(i);
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

  // Compute the derivative of the density.
}

void RadialModel::ComputeSurfaceGravityAndMomentOfInertiaFactor() {
  auto densityCoefficient = DensityCoefficient();
  auto density = mfem::GridFunction(&L2Space());
  density.ProjectCoefficient(densityCoefficient);
  {
    auto kernel = RadialCoefficient(
        [](auto r, auto attribute) { return std::pow(r, 2); });
    auto b = mfem::LinearForm(&L2Space());
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(kernel));
    b.Assemble();
    auto factor = 4 * std::numbers::pi_v<Real> * GravitationalConstant() /
                  std::pow(SurfaceRadius(), 2);
    _surfaceGravity = factor * (b(density));
  }

  {
    auto kernel = RadialCoefficient(
        [](auto r, auto attribute) { return std::pow(r, 4); });
    auto b = mfem::LinearForm(&L2Space());
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(kernel));
    b.Assemble();
    auto factor = 8 * std::numbers::pi_v<Real> * GravitationalConstant() /
                  (3 * _surfaceGravity * std::pow(SurfaceRadius(), 4));
    _momentOfInertiaFactor = factor * (b * density);
  }
}

void RadialModel::ComputeGravitationalPotential() {

  using namespace mfem;

  // Set up the linear form.
  auto densityTimesRadiusSquared =
      RadialCoefficient([this](auto r, auto attribute) {
        const auto factor =
            -(4 * std::numbers::pi_v<Real> * GravitationalConstant());
        return factor * Density()(r, attribute) * r * r;
      });
  auto b = LinearForm(&H1Space());
  b.AddDomainIntegrator(new DomainLFIntegrator(densityTimesRadiusSquared));
  b.Assemble();

  // Set up the bilinear form.
  auto a = BilinearForm(&H1Space());
  auto radiusSquared =
      RadialCoefficient([](auto r, auto attribute) { return r * r; });
  a.AddDomainIntegrator(new DiffusionIntegrator(radiusSquared));
  auto DtN = ConstantCoefficient(SurfaceRadius());
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

void RadialModel::WriteAsDeckModel(const std::string &fileName,
                                   const std::array<std::string, 3> &header,
                                   Real maximumKnotSpacing) {

  // Open the file for output.
  auto modelFile = std::ofstream(fileName);

  // Write the header information.
  for (auto line : header) {
    modelFile << line << "\n";
  }

  // Loop over the layers
  for (auto i = 0; i < NumberOfLayers(); i++) {
    auto attribute = i + 1;
    auto [r0, r1] = LayerRadii(i);
    auto n = std::max<Int>(2, std::round((r1 - r0) / maximumKnotSpacing));
    auto dr = (r1 - r0) / (n - 1);
    for (auto j = 0; j < n; j++) {
      auto r = r0 + j * dr;
      modelFile << std::format(
          "{:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f} "
          "{:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f}\n",
          r * LengthScale(), Density()(r, attribute) * DensityScale(),
          VerticalPVelocity()(r, attribute) * VelocityScale(),
          VerticalSVelocity()(r, attribute) * VelocityScale(),
          BulkQualityFactor()(r, attribute), ShearQualityFactor()(r, attribute),
          HorizontalPVelocity()(r, attribute) * VelocityScale(),
          HorizontalSVelocity()(r, attribute) * VelocityScale(),
          AnisotropicEtaParameter()(r, attribute));
    }
  }
}

void RadialModel::WriteAsDeckModel(const std::string &fileName,
                                   Real maximumKnotSpacing) {
  auto message = std::string("Header line is ignored!");
  auto header = std::array<std::string, 3>{message, message, message};
  WriteAsDeckModel(fileName, header, maximumKnotSpacing);
}

void RadialModel::PrintMesh(const std::string &mesh_file) {
  std::ofstream ofs(mesh_file);
  ofs.precision(8);
  _mesh.Print(ofs);
  ofs.close();
}

std::function<Real(Real, Int)> RadialModel::VerticalPVelocity() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(LoveModulusC()(r, attribute) / Density()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::VerticalSVelocity() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(LoveModulusL()(r, attribute) / Density()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::HorizontalPVelocity() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(LoveModulusA()(r, attribute) / Density()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::HorizontalSVelocity() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(LoveModulusN()(r, attribute) / Density()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::AnisotropicEtaParameter() const {
  return [this](Real r, Real attribute) {
    return LoveModulusF()(r, attribute) /
           (LoveModulusA()(r, attribute) - 2 * LoveModulusL()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::BulkModulus() const {
  return [this](Real r, Real attribute) {
    constexpr auto ninth = static_cast<Real>(1) / static_cast<Real>(9);
    return ninth *
           (LoveModulusC()(r, attribute) +
            4 * (LoveModulusA()(r, attribute) - LoveModulusN()(r, attribute) +
                 LoveModulusF()(r, attribute)));
  };
}

std::function<Real(Real, Int)> RadialModel::ShearModulus() const {
  return [this](Real r, Real attribute) {
    constexpr auto fifteenth = static_cast<Real>(1) / static_cast<Real>(15);
    return fifteenth *
           (LoveModulusC()(r, attribute) + LoveModulusA()(r, attribute) +
            6 * LoveModulusL()(r, attribute) +
            5 * LoveModulusN()(r, attribute) -
            2 * LoveModulusF()(r, attribute));
  };
}

void RadialModel::Write(const mfem::GridFunction &f, const std::string &file,
                        Real scale) const {

  using namespace mfem;
  auto fout = std::ofstream(file);
  auto point = Vector(1);
  auto values = Vector(Order() + 1);
  auto pairs = std::vector<std::pair<Real, Real>>(Order() + 1);
  auto *fes = f.FESpace();
  for (auto i = 0; i < Mesh().GetNE(); i++) {
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
      fout << r * LengthScale() << " " << v * scale << std::endl;
    }
  }
}

void RadialModel::WriteDerivative(const mfem::GridFunction &f,
                                  const std::string &file, Real scale) const {

  using namespace mfem;
  auto fout = std::ofstream(file);
  auto point = Vector(1);
  auto grad = DenseMatrix(1, Order() + 1);
  auto pairs = std::vector<std::pair<Real, Real>>(Order() + 1);
  auto *fes = f.FESpace();
  for (auto i = 0; i < Mesh().GetNE(); i++) {
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
      fout << r * LengthScale() << " " << v * scale / LengthScale()
           << std::endl;
    }
  }
}

void RadialModel::Write(mfem::Coefficient &f, const std::string &file,
                        Real scale) const {

  // Loop over the mesh storing values.
  auto fout = std::ofstream(file);
  auto &fes = H1Space();
  auto point = mfem::Vector(1);
  auto pairs = std::vector<std::pair<Real, Real>>(Order() + 1);
  for (auto i = 0; i < Mesh().GetNE(); i++) {
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
      fout << r * LengthScale() << " " << v * scale << std::endl;
    }
  }
}

void RadialModel::WriteDerivative(mfem::Coefficient &f, const std::string &file,
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
  for (auto i = 0; i < Mesh().GetNE(); i++) {
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
      fout << r * LengthScale() << " " << v * scale / LengthScale()
           << std::endl;
    }
  }
}

} // namespace LoveNumbers