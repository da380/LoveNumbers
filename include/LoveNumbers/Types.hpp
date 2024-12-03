#pragma once

#include "mfem.hpp"
#include <complex>

namespace LoveNumbers {

// Define type aliases for numeric types.
using Int = int;
using Real = mfem::real_t;
using Complex = std::complex<Real>;

} // namespace LoveNumbers
