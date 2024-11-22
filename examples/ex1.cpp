
#include "mfem.hpp"

#include "mfem.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>

#include "LoveNumbers/RadialMesh.h"

// using namespace std;
using namespace mfem;

template <std::ranges::range Range1, std::ranges::range Range2>
Mesh RadialMesh(const Range1 &layerRadii, const Range2 &maximumElementSizes);

template <std::ranges::range Range1>
Mesh RadialMesh(const Range1 &layerRadii, real_t maximumElementSize);

int main(int argc, char *argv[]) {

  const auto mesh_file = std::string("radial.mesh");

  // Set the number of radial layers.
  auto numberOfLayers = 3;

  // Set the radii for the layer boundaries.
  auto layerRadii = std::vector<real_t>(numberOfLayers + 1);
  layerRadii[0] = 0;
  layerRadii[1] = 1;
  layerRadii[2] = 1.5;
  layerRadii[3] = 3.5;

  auto mesh = RadialMesh(layerRadii, 0.076);

  std::ofstream ofs(mesh_file);
  ofs.precision(8);
  mesh.Print(ofs);
  ofs.close();

  return 0;
}

template <std::ranges::range Range1>
Mesh RadialMesh(const Range1 &layerRadii, real_t maximumElementSize) {
  auto maximumElementSizes =
      std::vector<real_t>(layerRadii.size() - 1, maximumElementSize);
  return RadialMesh(layerRadii, maximumElementSizes);
}

template <std::ranges::range Range1, std::ranges::range Range2>
Mesh RadialMesh(const Range1 &layerRadii, const Range2 &maximumElementSizes) {

  using std::ranges::fold_left_first;
  using namespace std::ranges::views;

  auto numberOfLayers = layerRadii.size() - 1;
  assert(numberOfLayers == maximumElementSizes.size());

  auto layers = layerRadii | adjacent<2>;

  auto elementsPerLayer =
      zip(layers, maximumElementSizes) | transform([](auto triple) {
        auto [layer, drMax] = triple;
        auto [r1, r2] = layer;
        return std::max<int>(1, std::round((r2 - r1) / drMax));
      });

  auto numberOfElements =
      fold_left_first(elementsPerLayer, std::plus<>()).value();
  auto numberOfVertices = numberOfElements + 1;
  auto numberOfBoundaryElements = numberOfLayers + 1;

  auto mesh =
      Mesh(1, numberOfVertices, numberOfElements, numberOfBoundaryElements);

  auto vertex1 = 0;
  auto boundaryAttribute = 1;
  for (auto quadruple : enumerate(zip(layers, elementsPerLayer))) {
    auto [i, triple] = quadruple;
    auto [radii, n] = triple;
    auto [r1, r2] = radii;
    auto dr = (r2 - r1) / n;
    auto vertices = iota(0, n + 1) |
                    transform([r1, dr](auto j) { return r1 + j * dr; }) |
                    drop(i > 0 ? 1 : 0);

    for (auto vertex : vertices) {
      mesh.AddVertex(vertex);
    }

    auto vertex2 = vertex1 + n;
    for (auto j = vertex1; j < vertex2; j++) {
      mesh.AddSegment(j, j + 1, i + 1);
    }

    if (i == 0) {
      mesh.AddBdrPoint(vertex1, boundaryAttribute++);
      mesh.AddBdrPoint(vertex2, boundaryAttribute++);
    } else {
      mesh.AddBdrPoint(vertex2, boundaryAttribute++);
    }
    vertex1 = vertex2;
  }

  mesh.FinalizeMesh();

  return mesh;
}
