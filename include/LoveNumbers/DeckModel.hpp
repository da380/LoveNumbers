#pragma once

#include "Configure.hpp"
#include "LoveNumbers/Configure.hpp"
#include "LoveNumbers/RadialModel.hpp"
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
  std::vector<std::pair<Int, Int>> _boundaryIndices;
  std::vector<std::pair<Real, Real>> _boundaryRadii;
  std::vector<bool> _layerSolid;

  // Cubic spline interpolating functions within each layer.
  using Spline = Interpolation::CubicSpline<std::vector<Real>::iterator,
                                            std::vector<Real>::iterator>;
  std::vector<std::unique_ptr<Spline>> _rhoSplines;
  std::vector<std::unique_ptr<Spline>> _ASplines;
  std::vector<std::unique_ptr<Spline>> _CSplines;
  std::vector<std::unique_ptr<Spline>> _FSplines;
  std::vector<std::unique_ptr<Spline>> _LSplines;
  std::vector<std::unique_ptr<Spline>> _NSplines;
  std::vector<std::unique_ptr<Spline>> _QKappaSplines;
  std::vector<std::unique_ptr<Spline>> _QMuSplines;

public:
  DeckModel() = delete;

  // Construct the deck model from file with given scale parameters.
  DeckModel(const std::string &fileName, Real lengthScale, Real massScale,
            Real timeScale)
      : RadialModel(lengthScale, massScale, timeScale) {
    ReadModelFile(fileName);
  }

  // Construct the deck model from file with given scale parameters.
  DeckModel(const std::string &fileName)
      : DeckModel(fileName, _LENGTH_SCALE, _MASS_SCALE, _TIME_SCALE) {}

  // Returns a deck model with its finite element mesh set up based
  // on a uniform maximum element size.
  static DeckModel FromMaximumElementSize(const std::string &fileName,
                                          Real maximumElementSize);

  // Returns a deck model with its finite element mesh set up with a
  // a uniform maximum element size based on the Jeans length for  the
  // given maximum degree. A default scale factor of 5 is used, but this
  // can optionally be set directly.
  static DeckModel FromMaximumDegree(const std::string &fileName,
                                     Int maximumDegree);

  // Return the number of layers. Override of pure virtual function in base
  // class.
  Int NumberOfLayers() const override;

  // Return the bounding radii of the ith layer. Override of pure
  // virtual function in base class.
  std::pair<Real, Real> LayerRadii(Int i) const override;

  // Return true if ith layer is solid. Override of pure virtual function in
  // base class.
  bool LayerIsSolid(Int i) const override;

  // Return the number of knots in the model.
  Int NumberOfKnots() const;

  // Return material parameter functions in each layer.
  std::function<Real(Real)> Rho(Int i) const override;

  std::function<Real(Real)> A(Int i) const override;

  std::function<Real(Real)> C(Int i) const override;

  std::function<Real(Real)> F(Int i) const override;

  std::function<Real(Real)> L(Int i) const override;

  std::function<Real(Real)> N(Int i) const override;

  std::function<Real(Real)> QKappa(Int i) const override;

  std::function<Real(Real)> QMu(Int i) const override;

private:
  // Read and process the model file. This method does not
  // set up the finite element mesh.
  void ReadModelFile(const std::string &fileName);
};

} // namespace LoveNumbers