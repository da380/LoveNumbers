#pragma once

#include "mfem.hpp"
#include <complex>

namespace LoveNumbers {

// Define type aliases for numeric types.
using Int = int;
using Real = mfem::real_t;
using Complex = std::complex<Real>;

// Set the default scale parameters.
const Real _LENGTH_SCALE = 6.371e6;
const Real _MASS_SCALE = 5.972e24;
const Real _TIME_SCALE = 3.6e3;

} // namespace LoveNumbers
