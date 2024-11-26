
#include "LoveNumbers/RadialPlanetaryModel.h"
#include <iterator>

namespace LoveNumbers {

void RadialPlanetaryModel::BuildMesh(
    const std::vector<mfem::real_t> &maximumElementSizes) {

  // Determine the number of elements per layer.
  auto elementsPerLayer =
      std::ranges::views::zip(Layers(),
                              std::ranges::views::all(maximumElementSizes)) |
      std::ranges::views::transform([](auto triple) {
        auto [layer, drMax] = triple;
        auto [r1, r2] = layer;
        return std::max<int>(1, std::round((r2 - r1) / drMax));
      });

  // Set the numbers of elements, vertices, and boundary elements.
  auto numberOfElements =
      std::ranges::fold_left_first(elementsPerLayer, std::plus<>()).value();
  auto numberOfVertices = numberOfElements + 1;
  auto numberOfBoundaryElements = NumberOfLayers() + 1;

  // Initiallise the mesh.
  _mesh = mfem::Mesh(1, numberOfVertices, numberOfElements,
                     numberOfBoundaryElements);
  _domainAttributes.SetSize(NumberOfLayers());

  // Set the domain attributes for the layers.
  for (auto i : LayerIndices()) {
    _domainAttributes[i] = i + 1;
  }

  // Loop over the layers building up the mesh.
  auto vertex1 = 0;
  auto domainAttribute = 1;
  auto boundaryAttribute = 1;
  for (auto triple : std::ranges::views::zip(Layers(), elementsPerLayer)) {
    auto [radii, n] = triple;
    auto [r1, r2] = radii;

    // Set the vertices within the layer.
    auto vertices = std::ranges::views::iota(0, n + 1) |
                    std::ranges::views::transform([&](auto j) {
                      auto dr = (r2 - r1) / static_cast<Real>(n);
                      return r1 + j * dr;
                    }) |
                    std::ranges::views::drop(domainAttribute > 1 ? 1 : 0);

    // Add in the vertices.
    for (auto vertex : vertices) {
      _mesh.AddVertex(vertex);
    }

    // Add in the elements.
    auto vertex2 = vertex1 + n;
    for (auto j = vertex1; j < vertex2; j++) {
      _mesh.AddSegment(j, j + 1, domainAttribute);
    }

    // Label the boundary elements.
    if (domainAttribute == 1) {
      _mesh.AddBdrPoint(vertex1, boundaryAttribute++);
      _mesh.AddBdrPoint(vertex2, boundaryAttribute++);
    } else {
      _mesh.AddBdrPoint(vertex2, boundaryAttribute++);
    }

    vertex1 = vertex2;
  }

  // Finalise the mesh construction.
  _mesh.FinalizeMesh();
}

mfem::Array<Int> RadialPlanetaryModel::SolidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfLayers());
  std::ranges::transform(IsSolid(), marker.begin(),
                         [](auto isSolid) { return isSolid ? 1 : 0; });
  return marker;
}

mfem::Array<Int> RadialPlanetaryModel::FluidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfLayers());
  std::ranges::transform(IsSolid(), marker.begin(),
                         [](auto isSolid) { return isSolid ? 0 : 1; });
  return marker;
}

mfem::Array<Int> RadialPlanetaryModel::FreeSurfaceMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaryElements());
  marker = 0;
  marker[NumberOfBoundaryElements() - 1] = 1;
  return marker;
}

mfem::Array<Int> RadialPlanetaryModel::FluidSolidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaryElements());
  marker = 0;
  for (auto i = 0; i < NumberOfLayers() - 1; i++) {
    if (!_isSolid[i] && _isSolid[i + 1]) {
      marker[i + 1] = 1;
    }
  }
  return marker;
}

mfem::Array<Int> RadialPlanetaryModel::SolidFluidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaryElements());
  marker = 0;
  for (auto i = 0; i < NumberOfLayers() - 1; i++) {
    if (_isSolid[i] && !_isSolid[i + 1]) {
      marker[i + 1] = 1;
    }
  }
  return marker;
}

mfem::Array<Int> RadialPlanetaryModel::SolidSolidMarker() const {
  auto marker = mfem::Array<Int>(NumberOfBoundaryElements());
  marker = 0;
  for (auto i = 0; i < NumberOfLayers() - 1; i++) {
    if (_isSolid[i] && _isSolid[i + 1]) {
      marker[i + 1] = 1;
    }
  }
  return marker;
}

} // namespace LoveNumbers