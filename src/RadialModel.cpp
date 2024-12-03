
#include "LoveNumbers/RadialModel.hpp"
#include <numeric>

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
  for (auto i = 1; i < NumberOfLayers(); i++) {
    marker[i] = !LayerIsSolid(i - 1) && LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModel::SolidFluidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  for (auto i = 1; i < NumberOfLayers(); i++) {
    marker[i] = LayerIsSolid(i - 1) && !LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModel::SolidSolidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  for (auto i = 1; i < NumberOfLayers(); i++) {
    marker[i] = LayerIsSolid(i - 1) && LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

mfem::Array<Int> RadialModel::FluidFluidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaries());
  for (auto i = 1; i < NumberOfLayers(); i++) {
    marker[i] = !LayerIsSolid(i - 1) && !LayerIsSolid(i) ? 1 : 0;
  }
  return marker;
}

void RadialModel::BuildMesh(Real maximumElementSize) {

  // Clear the mesh in case already set up.
  _mesh.Clear();

  // Work out number of elements in each layer.
  auto elementsInLayer = std::vector<Int>();
  for (auto i : LayerIndices()) {
    auto [r1, r2] = LayerRadii(i);
    elementsInLayer.push_back(
        std::max<int>(1, std::round((r2 - r1) / maximumElementSize)));
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

} // namespace LoveNumbers