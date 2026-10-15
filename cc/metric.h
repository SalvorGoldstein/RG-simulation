#pragma once
#include "types.h"

inline constexpr double M       = 5.972e30;
inline constexpr double r_s     = 2.0 * G * M / (c * c);
inline constexpr double r_earth = 6378137.0;
inline constexpr double r_moon  = 3.844e8;


inline constexpr double J       = 0.1 * G * M * M / c;

func schwarzschild(int i, int j);
func kerr(int i, int j);