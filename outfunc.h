#pragma once

#include <vector>
#include <functional>
#include <string>
#include <cmath>
#include <autodiff/forward/dual.hpp>

using autodiff::dual;

using func = std::function<dual(dual, dual, dual, dual)>;
using met  = std::vector<std::vector<func>>;

inline constexpr double G = 6.67430e-11;
inline constexpr double c = 299792458.0;

using rhs_fn = void (*)(double (&)[8], double, double (&)[8], int,
                        const met&, const met&);

dual partial(const func& fn,
             dual t, dual r, dual theta, dual phi,
             int coordinate);

dual christoffel(dual t, dual r, dual theta, dual phi,
                 int i, int k, int l,
                 const met& g, const met& invg);

void euler(rhs_fn sd,
           double (&q)[8], double t, double dt, int n,
           const met& g, const met& invg);

double timelikeU(double t, double r, double theta, double phi,
                 double ur, double utheta, double uphi,
                 const met& g);

double lightlikeU(double t, double r, double theta, double phi,
                  double ur, double utheta, double uphi,
                  const met& g);

void geodesics(double r, double theta0, double phi0,
               double ur, double utheta, double uphi,
               double r_max, int steps, double dtau,
               const std::string& filename,
               const met& g, const met& invg,
               rhs_fn sd, double r_min,
               bool type = true);