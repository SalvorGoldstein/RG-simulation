#pragma once

#include <vector>
#include <functional>
#include <string>
#include <cmath>
#include <autodiff/forward/dual.hpp>

using autodiff::dual;

using func   = std::function<dual(dual, dual, dual, dual)>;
using met    = std::vector<std::vector<func>>;
using christ = std::vector<std::vector<std::vector<func>>>;

inline constexpr double G = 6.67430e-11;
inline constexpr double c = 299792458.0;

using rhs_fn = void (*)(double (&)[8], double, double (&)[8], int,
                        const christ&);