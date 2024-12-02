#pragma once

#include <algorithm>
#include <concepts>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

#include <Interpolation/CubicSpline>

namespace LoveNumbers {

template <std::floating_point Real> class DeckModel {
private:
  // Vectors storing the model values at the radial knots.
  std::vector<Real> _r;      // radius
  std::vector<Real> _rho;    // density
  std::vector<Real> _A;      // Love parameter, A
  std::vector<Real> _C;      // Love parameter, C
  std::vector<Real> _F;      // Love parameter, F
  std::vector<Real> _L;      // Love parameter, L
  std::vector<Real> _N;      // Love parameter, N
  std::vector<Real> _QKappa; // Bulk quality factor
  std::vector<Real> _QMu;    // Shear quality factor

  // Radii of internal and external boundaries (including the model's centre).
  std::vector<Real> _layerRadii;

  // Value for the ith layer is a pair of integers giving the knot
  // at the top and bottom of the layer.
  std::vector<std::pair<int, int>> _layerKnots;

  // Value for the ith layer is true if it is solid.
  std::vector<bool> _layerSolid;

  // Data for interpolating functions.
  using CubicSpline =
      Interpolation::CubicSpline<typename std::vector<Real>::iterator,
                                 typename std::vector<Real>::iterator>;
  std::vector<CubicSpline> _rhoSpline;
  std::vector<CubicSpline> _ASpline;
  std::vector<CubicSpline> _CSpline;
  std::vector<CubicSpline> _FSpline;
  std::vector<CubicSpline> _LSpline;
  std::vector<CubicSpline> _NSpline;
  std::vector<CubicSpline> _QKappaSpline;
  std::vector<CubicSpline> _QMuSpline;

public:
  // Default constructor.
  DeckModel() = default;

  // Standard construction is from a deck file.
  DeckModel(const std::string &fileName) {

    // Check the file exists.
    if (!std::filesystem::exists(fileName)) {
      throw std::runtime_error("deck file does not exist");
    }

    // Read in the deck file.
    ReadDeckModel(fileName);

    // Non-dimensionalise the values.
    NonDimensionlise();

    // Determine the layering.
    SetLayers();

    // Set up the cubic splines.
    SetSplines();
  }

  // Return the number of radial knots.
  auto NumberOfRadialKnots() const { return _r.size(); }

  // Return the number of layers.
  auto NumberOfLayers() const { return _layerRadii.size() - 1; }

  // Return the bounding knots for each layer.
  auto LayerKnots() const { return std::ranges::views::all(_layerKnots); }

  // Reuturn a range over the layer indices.
  auto LayerIndices() const {
    return std::ranges::views::iota(0, NumberOfLayers());
  }

  // Return the bounding radii for each layer.
  auto LayerRadii() const {
    return LayerKnots() | std::ranges::views::transform([this](auto knots) {
             return std::pair(_r[knots.first], _r[knots.second]);
           });
  }

  auto LayerRadii(int i) const { return LayerRadii()[i]; }

  // Return whether the layers are solid.
  auto LayerSolid() const { return std::ranges::views::all(_layerSolid); }

  // Return whether the layers are fluid.
  auto LayerFluid() const {
    return LayerSolid() |
           std::ranges::views::transform([](auto solid) { return !solid; });
  }

  // Return parameter function in the ith layer.
  auto Rho(int i) const {
    return [i, this](auto r) { return _rhoSpline[i](r); };
  }

  auto A(int i) const {
    return [i, this](auto r) { return _ASpline[i](r); };
  }

  auto C(int i) const {
    return [i, this](auto r) { return _CSpline[i](r); };
  }

  auto F(int i) const {
    return [i, this](auto r) { return _FSpline[i](r); };
  }

  auto L(int i) const {
    return [i, this](auto r) { return _LSpline[i](r); };
  }

  auto N(int i) const {
    return [i, this](auto r) { return _NSpline[i](r); };
  }

private:
  void ReadDeckModel(const std::string &fileName) {

    // Open the file.
    auto modelFile = std::ifstream(fileName);

    // Skip the first three lines
    for (auto i = 0; i < 3; i++) {
      auto line = std::string();
      std::getline(modelFile, line);
    }

    // Read in the model data.
    auto line = std::string();
    while (std::getline(modelFile, line)) {

      // Get the model parameters at the knot.
      Real r, rho, vpv, vsv, qKappa, qMu, vph, vsh, eta;
      if (!(std::stringstream(line) >> r >> rho >> vpv >> vsv >> qKappa >>
            qMu >> vph >> vsh >> eta)) {
        break;
      }

      // Covert velocities to modulii.
      auto A = rho * vph * vph;
      auto C = rho * vpv * vpv;
      auto L = rho * vsv * vsv;
      auto N = rho * vsh * vsh;
      auto F = eta * (A - 2 * L);

      // Store the values.
      _r.push_back(r);
      _rho.push_back(rho);
      _A.push_back(A);
      _C.push_back(C);
      _F.push_back(F);
      _L.push_back(L);
      _N.push_back(N);
    }
  }

  // Non-dimensionalise the model parameters.
  void NonDimensionlise() {}

  // Determine the layer radii.
  void SetLayers() {

    // Determine the knots at all external and internal boundaries.
    constexpr auto eps = std::numeric_limits<Real>::epsilon();
    auto knots = std::vector<int>();
    knots.push_back(0);
    for (auto [i, radii] : std::ranges::views::enumerate(
             std::ranges::views::all(_r) | std::ranges::views::adjacent<2>)) {
      const auto [r1, r2] = radii;
      const auto sum = r1 + r2;
      const auto difference = r2 - r1;
      if (difference / sum < eps) {
        knots.push_back(i);
        knots.push_back(i + 1);
      }
    }
    knots.push_back(NumberOfRadialKnots() - 1);

    // Store the indices at the bottom and top of each layer.
    for (auto i = 0; 2 * i < knots.size(); i++) {
      auto j = 2 * i;
      _layerKnots.push_back({knots[j], knots[j + 1]});
    }

    // Now determine if each layer is solid.
    for (auto [i1, i2] : LayerKnots()) {
      _layerSolid.push_back(_L[i1] > 0);
    }
  }

  // Compute and store cubic spline coefficients within each layer.
  void SetSplines() {
    for (auto [i1, i2] : LayerKnots()) {
      auto rStart = std::next(_r.begin(), i1);
      auto rFinish = std::next(_r.begin(), i2 + 1);
      _rhoSpline.push_back(
          CubicSpline(rStart, rFinish, std::next(_rho.begin(), i1)));
      _ASpline.push_back(
          CubicSpline(rStart, rFinish, std::next(_A.begin(), i1)));
      _CSpline.push_back(
          CubicSpline(rStart, rFinish, std::next(_C.begin(), i1)));
      _FSpline.push_back(
          CubicSpline(rStart, rFinish, std::next(_F.begin(), i1)));
      _LSpline.push_back(
          CubicSpline(rStart, rFinish, std::next(_L.begin(), i1)));
      _NSpline.push_back(
          CubicSpline(rStart, rFinish, std::next(_N.begin(), i1)));
    }
  }
};
} // namespace LoveNumbers