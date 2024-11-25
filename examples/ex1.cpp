
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

#include "LoveNumbers/RadialPlanetaryModel.h"

int main(int argc, char *argv[]) {

  using namespace mfem;
  using Real = mfem::real_t;

  const auto mesh_file = std::string("radial.mesh");

  // Set the number of radial layers.
  auto numberOfLayers = 3;

  // Set the radii for the layer boundaries.
  auto layerRadii = std::vector<real_t>(numberOfLayers + 1);
  layerRadii[0] = 0;
  layerRadii[1] = 1;
  layerRadii[2] = 1.5;
  layerRadii[3] = 5.5;

  auto maximumElementSize = 0.1;

  auto model =
      LoveNumbers::RadialPlanetaryModel(layerRadii, maximumElementSize);

  // model.PrintMesh(mesh_file);

  return 0;
}
