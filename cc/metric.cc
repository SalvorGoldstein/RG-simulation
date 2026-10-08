#include "metric.h"
#include "funcoperator.h"     // constant
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

    // sigma = r^2 + a^2 cos^2(th)     delta = r^2 - r_s r + a^2

    if (i == 0 && j == 0)
        return [](dual t, dual r, dual th, dual ph) {
            dual sigma = r*r + a_kerr*a_kerr*cos(th)*cos(th);
            dual res = c*c*(1.0 - r_s*r/sigma);
            return res;
        };
    if (i == 0 && j == 3)
        return [](dual t, dual r, dual th, dual ph) {
            dual sigma = r*r + a_kerr*a_kerr*cos(th)*cos(th);
            dual res = c*r_s*r*a_kerr*sin(th)*sin(th)/sigma;
            return res;
        };
    if (i == 1 && j == 1)
        return [](dual t, dual r, dual th, dual ph) {
            dual sigma = r*r + a_kerr*a_kerr*cos(th)*cos(th);
            dual delta = r*r - r_s*r + a_kerr*a_kerr;
            dual res = -sigma/delta;
            return res;
        };
    if (i == 2 && j == 2)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -(r*r + a_kerr*a_kerr*cos(th)*cos(th));
            return res;
        };
    if (i == 3 && j == 3)
        return [](dual t, dual r, dual th, dual ph) {
            dual sigma = r*r + a_kerr*a_kerr*cos(th)*cos(th);
            dual res = -(r*r + a_kerr*a_kerr
                         + r_s*r*a_kerr*a_kerr*sin(th)*sin(th)/sigma)
                       * sin(th)*sin(th);
            return res;
        };
    return constant(0.0);
}