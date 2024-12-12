#include "LoveNumbers/RadialModel.hpp"
#include <cstddef>
#include <memory>

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
}

void RadialModel::ComputeSurfaceGravityAndMomentOfInertiaFactor() {
  auto rhoCoefficient = RhoCoefficient();
  auto rho = mfem::GridFunction(L2Space());
  rho.ProjectCoefficient(rhoCoefficient);
  {
    auto kernel = RadialModelCoefficient(
        [](auto r, auto attribute) { return std::pow(r, 2); });
    auto b = mfem::LinearForm(L2Space());
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(kernel));
    b.Assemble();
    auto factor = 4 * std::numbers::pi_v<Real> * GravitationalConstant() /
                  std::pow(SurfaceRadius(), 2);
    _surfaceGravity = factor * (b * rho);
  }

  {
    auto kernel = RadialModelCoefficient(
        [](auto r, auto attribute) { return std::pow(r, 4); });
    auto b = mfem::LinearForm(L2Space());
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(kernel));
    b.Assemble();
    auto factor = 8 * std::numbers::pi_v<Real> * GravitationalConstant() /
                  (3 * _surfaceGravity * std::pow(SurfaceRadius(), 4));
    _momentOfInertiaFactor = factor * (b * rho);
  }
}

void RadialModel::ComputeGravitationalPotential() {

  using namespace mfem;

  // Set up the linear form.
  auto rhoTimesRadiusSquared = RadialModelCoefficient(
      [this](auto r, auto attribute) { return Rho()(r, attribute) * r * r; });
  auto b = LinearForm(H1Space());
  b.AddDomainIntegrator(new DomainLFIntegrator(rhoTimesRadiusSquared));
  b.Assemble();
  b *= -(4 * std::numbers::pi_v<Real> * GravitationalConstant());

  // Set up the bilinear form.
  auto a = BilinearForm(H1Space());
  auto radiusSquared =
      RadialModelCoefficient([](auto r, auto attribute) { return r * r; });
  a.AddDomainIntegrator(new DiffusionIntegrator(radiusSquared));
  auto DtN = ConstantCoefficient(SurfaceRadius());
  auto surfaceMarker = SurfaceMarker();
  a.AddBoundaryIntegrator(new BoundaryMassIntegrator(DtN), surfaceMarker);
  a.Assemble();

  // Set up the solution vector with appropriate boundary values.
  _gravitationalPotential = std::make_unique<GridFunction>(H1Space());
  *_gravitationalPotential = 0;

  // Set up the linear system.
  OperatorPtr A;
  Vector B, X;
  auto ess_tdof_list = Array<Int>();
  a.FormLinearSystem(ess_tdof_list, *_gravitationalPotential, b, A, X, B);

  // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
  GSSmoother M((SparseMatrix &)(*A));
  PCG(*A, M, B, X, 0, H1Space()->GetNDofs() * 2, 1e-12, 0.0);
  a.RecoverFEMSolution(X, b, *_gravitationalPotential);

  // Form the gravitational acceleration.
  _gravitationalAcceleration = std::make_unique<GridFunction>(L2Space());
  _gravitationalPotential->GetDerivative(1, 0, *_gravitationalAcceleration);
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
          r * LengthScale(), Rho()(r, attribute) * DensityScale(),
          VPV()(r, attribute) * VelocityScale(),
          VSV()(r, attribute) * VelocityScale(), QKappa()(r, attribute),
          QMu()(r, attribute), VPH()(r, attribute) * VelocityScale(),
          VSH()(r, attribute) * VelocityScale(), Eta()(r, attribute));
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

std::function<Real(Real, Int)> RadialModel::VPV() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(C()(r, attribute) / Rho()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::VSV() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(L()(r, attribute) / Rho()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::VPH() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(A()(r, attribute) / Rho()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::VSH() const {
  return [this](Real r, Real attribute) {
    return std::sqrt(N()(r, attribute) / Rho()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::Eta() const {
  return [this](Real r, Real attribute) {
    return F()(r, attribute) / (A()(r, attribute) - 2 * L()(r, attribute));
  };
}

std::function<Real(Real, Int)> RadialModel::Kappa() const {
  return [this](Real r, Real attribute) {
    constexpr auto ninth = static_cast<Real>(1) / static_cast<Real>(9);
    return ninth *
           (C()(r, attribute) +
            4 * (A()(r, attribute) - N()(r, attribute) + F()(r, attribute)));
  };
}

std::function<Real(Real, Int)> RadialModel::Mu() const {
  return [this](Real r, Real attribute) {
    constexpr auto fifteenth = static_cast<Real>(1) / static_cast<Real>(15);

    return fifteenth *
           (C()(r, attribute) + A()(r, attribute) + 6 * L()(r, attribute) +
            5 * N()(r, attribute) - 2 * F()(r, attribute));
  };
}

void RadialModel::WriteGridFunction(const mfem::GridFunction &f,
                                    const std::string &file, Real scale) const {

  // Loop over the mesh storing values.
  auto pairs = std::vector<std::pair<Real, Real>>();
  auto point = mfem::Vector();
  auto *fes = f.FESpace();
  for (auto i = 0; i < Mesh().GetNE(); i++) {
    auto *el = fes->GetFE(i);
    auto *eltrans = fes->GetElementTransformation(i);
    auto &ir = el->GetNodes();
    auto values = mfem::Vector();
    f.GetValues(*eltrans, ir, values);
    for (auto j = 0; j < ir.GetNPoints(); j++) {
      auto &ip = ir.IntPoint(j);
      eltrans->SetIntPoint(&ip);
      eltrans->Transform(ip, point);
      pairs.push_back({point[0], values[j]});
    }
  }

  // Sort the radii to insure the proper ordering.
  std::ranges::sort(pairs,
                    [](auto p1, auto p2) { return p1.first < p2.first; });

  // Write the values to the file.
  auto fout = std::ofstream(file);
  for (auto [r, v] : pairs) {
    fout << r * LengthScale() << " " << v * scale << std::endl;
  }
}

void RadialModel::WriteCoefficient(mfem::Coefficient &&f,
                                   const std::string &file, Real scale) const {

  // Loop over the mesh storing values.
  auto pairs = std::vector<std::pair<Real, Real>>();
  auto point = mfem::Vector();
  auto *fes = H1Space();
  for (auto i = 0; i < Mesh().GetNE(); i++) {
    auto *el = fes->GetFE(i);
    auto *eltrans = fes->GetElementTransformation(i);
    auto &ir = el->GetNodes();
    for (auto j = 0; j < ir.GetNPoints(); j++) {
      auto &ip = ir.IntPoint(j);
      eltrans->SetIntPoint(&ip);
      eltrans->Transform(ip, point);
      auto value = f.Eval(*eltrans, ip);
      pairs.push_back({point[0], value});
    }
  }

  // Sort the radii to insure the proper ordering.
  std::ranges::stable_sort(
      pairs, [](auto p1, auto p2) { return p1.first < p2.first; });

  // Write the values to the file.
  auto fout = std::ofstream(file);
  for (auto [r, v] : pairs) {
    fout << r * LengthScale() << " " << v * scale << std::endl;
  }
}

} // namespace LoveNumbers