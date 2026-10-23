#include "example_metric.h"
#include <utility>            // swap

using namespace std;

func schwarzschild(int i, int j)
{
    if (i > j) swap(i, j);

    if (i == 0 && j == 0)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = c*c*(1.0 - r_s/r);
            return res;
        };
    if (i == 1 && j == 1)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -1.0/(1.0 - r_s/r);
            return res;
        };
    if (i == 2 && j == 2)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -r*r;
            return res;
        };
    if (i == 3 && j == 3)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -r*r*sin(th)*sin(th);
            return res;
        };
    return constant(0.0);
}

func kerr(int i, int j)
{
    if (i > j) swap(i, j);

    // a = J/(M c) (a length), sigma = r^2 + a^2 cos^2(th), delta = r^2 - r_s r + a^2

    if (i == 0 && j == 0)
        return [](dual t, dual r, dual th, dual ph) {
            constexpr double a = J / (M * c);
            dual sigma = r*r + a*a*cos(th)*cos(th);
            dual res = c*c*(1.0 - r_s*r/sigma);
            return res;
        };

    if (i == 0 && j == 3)
        return [](dual t, dual r, dual th, dual ph) {
            constexpr double a = J / (M * c);
            dual sigma = r*r + a*a*cos(th)*cos(th);
            dual res = c*r_s*r*a*sin(th)*sin(th)/sigma;
            return res;
        };

    if (i == 1 && j == 1)
        return [](dual t, dual r, dual th, dual ph) {
            constexpr double a = J / (M * c);
            dual sigma = r*r + a*a*cos(th)*cos(th);
            dual delta = r*r - r_s*r + a*a;
            dual res = -sigma/delta;
            return res;
        };

    if (i == 2 && j == 2)
        return [](dual t, dual r, dual th, dual ph) {
            constexpr double a = J / (M * c);
            dual res = -(r*r + a*a*cos(th)*cos(th));
            return res;
        };

    if (i == 3 && j == 3)
        return [](dual t, dual r, dual th, dual ph) {
            constexpr double a = J / (M * c);
            dual sigma = r*r + a*a*cos(th)*cos(th);
            dual res = -(r*r + a*a + r_s*r*a*a*sin(th)*sin(th)/sigma)
                       * sin(th)*sin(th);
            return res;
        };

    return constant(0.0);
}


func kerr_newman(int i, int j)
{
    if (i > j) swap(i, j);

    // sigma = r^2 + a^2 cos^2(th), delta = r^2 - r_s r + a^2 + r_Q^2
    // e = r_s r - r_Q^2

    if (i == 0 && j == 0)
        return [](dual t, dual r, dual th, dual ph) {
            dual sigma = r*r + a_spin*a_spin*cos(th)*cos(th);
            dual e = r_s*r - r_Q2;
            dual res = c*c*(1.0 - e/sigma);
            return res;
        };

    if (i == 0 && j == 3)
        return [](dual t, dual r, dual th, dual ph) {
            dual sigma = r*r + a_spin*a_spin*cos(th)*cos(th);
            dual e = r_s*r - r_Q2;
            dual res = c*e*a_spin*sin(th)*sin(th)/sigma;
            return res;
        };

    if (i == 1 && j == 1)
        return [](dual t, dual r, dual th, dual ph) {
            dual sigma = r*r + a_spin*a_spin*cos(th)*cos(th);
            dual delta = r*r - r_s*r + a_spin*a_spin + r_Q2;
            dual res = -sigma/delta;
            return res;
        };

    if (i == 2 && j == 2)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -(r*r + a_spin*a_spin*cos(th)*cos(th));
            return res;
        };

    if (i == 3 && j == 3)
        return [](dual t, dual r, dual th, dual ph) {
            dual sigma = r*r + a_spin*a_spin*cos(th)*cos(th);
            dual e = r_s*r - r_Q2;
            dual res = -(r*r + a_spin*a_spin + e*a_spin*a_spin*sin(th)*sin(th)/sigma)
                       * sin(th)*sin(th);
            return res;
        };

    return constant(0.0);
}

func ellis_wormhole(int i, int j)
{
    if (i > j) swap(i, j);

    if (i == 0 && j == 0)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = c*c;
            return res;
        };

    if (i == 1 && j == 1)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -1.0;
            return res;
        };

    if (i == 2 && j == 2)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -(r*r + a*a);
            return res;
        };

    if (i == 3 && j == 3)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -(r*r + a*a)*sin(th)*sin(th);
            return res;
        };

    return constant(0.0);
}