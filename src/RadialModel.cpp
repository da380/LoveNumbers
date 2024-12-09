#include "LoveNumbers/RadialModel.hpp"
#include <cstddef>

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

mfem::Array<Int> RadialModel::FreeSurfaceMarker() const {
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

Int RadialModel::NumberOfBoundaries() const { return NumberOfLayers() + 1; }

Real RadialModel::JeanLength(Int degree) const {
  return 2 * std::numbers::pi_v<Real> * SurfaceRadius() /
         static_cast<Real>(degree + 1);
}

Real RadialModel::SurfaceRadius() const {
  return LayerRadii(NumberOfLayers() - 1).second;
}

std::function<Real(Real)> RadialModel::VPV(Int i) const {
  return [this, i](Real r) -> Real { return std::sqrt(C(i)(r) / Rho(i)(r)); };
}

std::function<Real(Real)> RadialModel::VPH(Int i) const {
  return [this, i](Real r) -> Real { return std::sqrt(A(i)(r) / Rho(i)(r)); };
}

std::function<Real(Real)> RadialModel::VSV(Int i) const {
  return [this, i](Real r) -> Real { return std::sqrt(L(i)(r) / Rho(i)(r)); };
}

std::function<Real(Real)> RadialModel::VSH(Int i) const {
  return [this, i](Real r) -> Real { return std::sqrt(N(i)(r) / Rho(i)(r)); };
}

std::function<Real(Real)> RadialModel::Eta(Int i) const {
  return
      [this, i](Real r) -> Real { return F(i)(r) / (A(i)(r) - 2 * L(i)(r)); };
}

std::function<Real(Real)> RadialModel::Kappa(Int i) const {
  return [this, i](Real r) -> Real {
    constexpr auto ninth = static_cast<Real>(1) / static_cast<Real>(9);
    return ninth * (C(i)(r) + 4 * (A(i)(r) - N(i)(r) + F(i)(r)));
  };
}

std::function<Real(Real)> RadialModel::Mu(Int i) const {
  return [this, i](Real r) -> Real {
    constexpr auto fifteenth = static_cast<Real>(1) / static_cast<Real>(15);
    return fifteenth *
           (C(i)(r) + A(i)(r) + 6 * L(i)(r) + 5 * N(i)(r) - 2 * F(i)(r));
  };
}

void RadialModel::BuildMesh(Real maximumElementSize) {
  auto maximumElementSizes =
      std::vector<Real>(NumberOfLayers(), maximumElementSize);
  BuildMesh(maximumElementSizes);
}

void RadialModel::BuildMesh(const std::vector<Real> &maximumElementSizes) {

  // Clear the mesh in case already set up.
  _mesh.Clear();

  // Work out number of elements in each layer.
  auto elementsInLayer = std::vector<Int>();
  for (auto i : LayerIndices()) {
    auto [r1, r2] = LayerRadii(i);
    elementsInLayer.push_back(
        std::max<int>(1, std::round((r2 - r1) / maximumElementSizes[i])));
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

  // Set up the function coefficients and pointers.
  auto rhoPointers = mfem::Array<mfem::Coefficient *>(NumberOfLayers());
  auto APointers = mfem::Array<mfem::Coefficient *>(NumberOfLayers());
  auto CPointers = mfem::Array<mfem::Coefficient *>(NumberOfLayers());
  auto FPointers = mfem::Array<mfem::Coefficient *>(NumberOfLayers());
  auto LPointers = mfem::Array<mfem::Coefficient *>(NumberOfLayers());
  auto NPointers = mfem::Array<mfem::Coefficient *>(NumberOfLayers());
  for (auto i : LayerIndices()) {
    _rhoFunctionCoefficients.push_back(mfem::FunctionCoefficient(
        [this, i](mfem::Vector x) -> Real { return Rho(i)(x[0]); }));
    _AFunctionCoefficients.push_back(mfem::FunctionCoefficient(
        [this, i](mfem::Vector x) -> Real { return A(i)(x[0]); }));
    _CFunctionCoefficients.push_back(mfem::FunctionCoefficient(
        [this, i](mfem::Vector x) -> Real { return C(i)(x[0]); }));
    _FFunctionCoefficients.push_back(mfem::FunctionCoefficient(
        [this, i](mfem::Vector x) -> Real { return F(i)(x[0]); }));
    rhoPointers[i] = &_rhoFunctionCoefficients[i];
    APointers[i] = &_AFunctionCoefficients[i];
    CPointers[i] = &_CFunctionCoefficients[i];
    FPointers[i] = &_FFunctionCoefficients[i];
    if (LayerIsSolid(i)) {
      _LFunctionCoefficients.push_back(mfem::FunctionCoefficient(
          [this, i](mfem::Vector x) -> Real { return L(i)(x[0]); }));
      _NFunctionCoefficients.push_back(mfem::FunctionCoefficient(
          [this, i](mfem::Vector x) -> Real { return N(i)(x[0]); }));
      LPointers[i] = &_LFunctionCoefficients[i];
      NPointers[i] = &_NFunctionCoefficients[i];
    } else {
      LPointers[i] = nullptr;
      NPointers[i] = nullptr;
    }
  }

  // Set the piecewise coefficients.
  auto attributes = LayerAttributes();
  _rhoCoefficient = mfem::PWCoefficient(attributes, rhoPointers);
  _ACoefficient = mfem::PWCoefficient(attributes, APointers);
  _CCoefficient = mfem::PWCoefficient(attributes, CPointers);
  _FCoefficient = mfem::PWCoefficient(attributes, FPointers);
  _LCoefficient = mfem::PWCoefficient(attributes, LPointers);
  _NCoefficient = mfem::PWCoefficient(attributes, NPointers);
}

void RadialModel::SetFiniteElementSpace(Int order) {
  _FECollection = std::make_unique<mfem::H1_FECollection>(order, 1);
  _FESpace =
      std::make_unique<mfem::FiniteElementSpace>(&_mesh, _FECollection.get());
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
    auto [r0, r1] = LayerRadii(i);
    auto n = std::max<Int>(
        2, std::round((r1 - r0) / (maximumKnotSpacing / LengthScale())));
    auto dr = (r1 - r0) / (n - 1);
    auto rho = Rho(i);
    auto vpv = VPV(i);
    auto vsv = VSV(i);
    auto qKappa = QKappa(i);
    auto qMu = QMu(i);
    auto vph = VPH(i);
    auto vsh = VSH(i);
    auto eta = Eta(i);
    for (auto j = 0; j < n; j++) {
      auto r = r0 + j * dr;
      modelFile << std::format(
          "{:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f} "
          "{:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f}\n",
          r * LengthScale(), rho(r) * DensityScale(), vpv(r) * VelocityScale(),
          vsv(r) * VelocityScale(), qKappa(r), qMu(r), vph(r) * VelocityScale(),
          vsh(r) * VelocityScale(), eta(r));
    }
  }
}

void RadialModel::WriteAsDeckModel(const std::string &fileName,
                                   Real maximumKnotSpacing) {
  auto message = std::string("Header line is ignored!");
  auto header = std::array<std::string, 3>{message, message, message};
  WriteAsDeckModel(fileName, header, maximumKnotSpacing);
}

} // namespace LoveNumbers