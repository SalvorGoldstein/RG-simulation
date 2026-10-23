#pragma once
#include <cmath>
#include <rgsim/rgsim.hpp>

using namespace rgsim;
using namespace rgsim::constants;   // G, c, eps0, pi

inline constexpr double M       = 5.972e30;
inline constexpr double r_s     = 2.0 * G * M / (c * c);
inline constexpr double r_earth = 6378137.0;
inline constexpr double r_moon  = 3.844e8;
inline constexpr double Q       = 0.5 * M * std::sqrt(4 * pi * eps0 * G);

inline constexpr double J       = G * M * M / c;
inline constexpr double a_spin  = J / (M * c);
inline constexpr double r_Q2    = G * Q * Q / (4 * pi * eps0 * c * c * c * c);
constexpr double a = r_s;

func schwarzschild(int i, int j);
func kerr(int i, int j);
func kerr_newman(int i, int j);
func ellis_wormhole(int i, int j);