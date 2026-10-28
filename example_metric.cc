#include "example_metric.h"
#include <utility>            // swap

using namespace std;

// Uniform black string: 4D Schwarzschild times a flat fifth dimension w.
// Coordinates: x0 = t, x1 = r, x2 = theta, x3 = phi, x4 = w.
// Signature (+,-,-,-,-). A good first test: with u4 = 0 it must reproduce the 4D orbits.
func black_string(int i, int j)
{
    if (i > j) swap(i, j);

    if (i == 0 && j == 0)
        return [](const std::array<dual,N>& x) {
            dual r = x[1];
            dual res = c*c*(1.0 - r_s/r);
            return res;
        };
    if (i == 1 && j == 1)
        return [](const std::array<dual,N>& x) {
            dual r = x[1];
            dual res = -1.0/(1.0 - r_s/r);
            return res;
        };
    if (i == 2 && j == 2)
        return [](const std::array<dual,N>& x) {
            dual r = x[1];
            dual res = -r*r;
            return res;
        };
    if (i == 3 && j == 3)
        return [](const std::array<dual,N>& x) {
            dual r = x[1];
            dual th = x[2];
            dual res = -r*r*sin(th)*sin(th);
            return res;
        };
    if (i == 4 && j == 4)
        return [](const std::array<dual,N>& x) {
            dual res = -1.0;
            return res;
        };
    return constant(0.0);
}

// 5D Schwarzschild-Tangherlini black hole:
//   ds^2 = c^2 f dt^2 - dr^2/f - r^2 ( dchi^2 + sin^2(chi) ( dtheta^2 + sin^2(theta) dphi^2 ) )
//   f = 1 - (r_h5/r)^2
// Coordinates: x0 = t, x1 = r, x2 = chi, x3 = theta, x4 = phi.
// Avoid starting exactly on a pole (chi or theta = 0 or pi); chi = theta = pi/2 is fine.
func tangherlini(int i, int j)
{
    if (i > j) swap(i, j);

    if (i == 0 && j == 0)
        return [](const std::array<dual,N>& x) {
            dual r = x[1];
            dual res = c*c*(1.0 - r_h5*r_h5/(r*r));
            return res;
        };
    if (i == 1 && j == 1)
        return [](const std::array<dual,N>& x) {
            dual r = x[1];
            dual res = -1.0/(1.0 - r_h5*r_h5/(r*r));
            return res;
        };
    if (i == 2 && j == 2)
        return [](const std::array<dual,N>& x) {
            dual r = x[1];
            dual res = -r*r;
            return res;
        };
    if (i == 3 && j == 3)
        return [](const std::array<dual,N>& x) {
            dual r = x[1];
            dual chi = x[2];
            dual res = -r*r*sin(chi)*sin(chi);
            return res;
        };
    if (i == 4 && j == 4)
        return [](const std::array<dual,N>& x) {
            dual r = x[1];
            dual chi = x[2];
            dual th = x[3];
            dual res = -r*r*sin(chi)*sin(chi)*sin(th)*sin(th);
            return res;
        };
    return constant(0.0);
}