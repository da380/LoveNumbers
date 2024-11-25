
#include "LoveNumbers/RadialPlanetaryModel.h"
#include <iterator>

namespace LoveNumbers {

void RadialPlanetaryModel::BuildMesh() {

  // Determine the number of elements per layer.
  auto elementsPerLayer =
      std::ranges::views::zip(Layers(), MaximumElementSizes()) |
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

    domainAttribute++;
    vertex1 = vertex2;
  }

  // Finalise the mesh construction.
  _mesh.FinalizeMesh();
}

} // namespace LoveNumbers