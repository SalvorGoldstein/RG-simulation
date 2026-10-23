#pragma once

#include <vector>
#include <functional>
#include <string>
#include <cmath>
#include <autodiff/forward/dual.hpp>

namespace rgsim {

using autodiff::dual;

using func   = std::function<dual(dual, dual, dual, dual)>;
using met    = std::vector<std::vector<func>>;
using christ = std::vector<std::vector<std::vector<func>>>;

// Right-hand side of the equations of motion: q = (x0..x3, u0..u3)
using rhs_fn = std::function<void(double (&)[8], double, double (&)[8], int,
                                  const christ&)>;

namespace constants {
inline constexpr double G    = 6.67430e-11;       // m^3 kg^-1 s^-2
inline constexpr double c    = 299792458.0;       // m/s
inline constexpr double eps0 = 8.8541878128e-12;  // F/m
inline constexpr double pi   = 3.14159265358979323846;  // M_PI is not standard C++
}

} // namespace rgsim
