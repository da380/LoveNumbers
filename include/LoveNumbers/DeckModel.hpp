#pragma once

#include "Configure.hpp"
#include "LoveNumbers/Configure.hpp"
#include "LoveNumbers/RadialModel.hpp"
#include "RadialModel.hpp"
#include <Interpolation/CubicSpline>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <ostream>
#include <ranges>
#include <string>
#include <vector>

namespace LoveNumbers {

class DeckModel : public RadialModel {
private:
  // Store the header information.
  std::vector<std::string> _header;

  // Model values at the radial knots.
  std::vector<Real> _r;      // radius
  std::vector<Real> _rho;    // density
  std::vector<Real> _A;      // Love parameter, A
  std::vector<Real> _C;      // Love parameter, C
  std::vector<Real> _F;      // Love parameter, F
  std::vector<Real> _L;      // Love parameter, L
  std::vector<Real> _N;      // Love parameter, N
  std::vector<Real> _QKappa; // Bulk quality factor
  std::vector<Real> _QMu;    // Shear quality factor

  // Layering information.
  std::vector<Int> _boundaryIndex;
  std::vector<Real> _boundaryRadius;
  std::vector<bool> _layerSolid;

  // Cubic spline interpolating functions within each layer.
  using Spline = Interpolation::CubicSpline<std::vector<Real>::iterator,
                                            std::vector<Real>::iterator>;
  std::vector<Spline> _rhoSplines;
  std::vector<Spline> _ASplines;
  std::vector<Spline> _CSplines;
  std::vector<Spline> _FSplines;
  std::vector<Spline> _LSplines;
  std::vector<Spline> _NSplines;
  std::vector<Spline> _QKappaSplines;
  std::vector<Spline> _QMuSplines;

public:
  DeckModel() = delete;

  DeckModel(const Dimensions &dimensions, Int order,
            Real characteristicLengthScale, const std::string &fileName)
      : RadialModel(dimensions, order) {
    ReadModelFile(fileName);
    BuildMesh(characteristicLengthScale);
  }

  // Return the number of layers. Override of pure virtual function in
  // base class.
  Int NumberOfLayers() const override;

  // Return the bounding radii of the ith layer. Override of pure
  // virtual function in base class.
  std::pair<Real, Real> LayerRadii(Int i) const override;

  // Return true if ith layer is solid. Override of pure virtual function in
  // base class.
  bool LayerIsSolid(Int i) const override;

  // Return the number of knots in the model.
  Int NumberOfKnots() const;

  // Functions to return material parameters functions.
  std::function<Real(Real, Int)> Rho() const override;
  std::function<Real(Real, Int)> A() const override;
  std::function<Real(Real, Int)> C() const override;
  std::function<Real(Real, Int)> F() const override;
  std::function<Real(Real, Int)> L() const override;
  std::function<Real(Real, Int)> N() const override;
  std::function<Real(Real, Int)> QKappa() const override;
  std::function<Real(Real, Int)> QMu() const override;

private:
  // Read and process the model file. This method does not
  // set up the finite element mesh.
  void ReadModelFile(const std::string &fileName);
};

} // namespace LoveNumbers