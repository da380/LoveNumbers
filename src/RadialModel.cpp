
#include "LoveNumbers/RadialModel.hpp"

namespace LoveNumbers {

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
          "{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} "
          "{:10.3f} {:10.3f} {:10.3f} {:10.3f}\n",
          r * LengthScale(), rho(r) * DensityScale(), vpv(r) * VelocityScale(),
          vsv(r) * VelocityScale(), qKappa(r), qMu(r), vph(r) * VelocityScale(),
          vsh(r) * VelocityScale(), eta(r));
    }
  }
}

} // namespace LoveNumbers