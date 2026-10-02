#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <autodiff/forward/dual.hpp>
#include "outfunc.h"

using autodiff::dual;
using namespace std;

constexpr double M = 5.972e24; //for earth

double r_s = 2.0 * G * M / (c * c);
double r_earth = 6378137;
double r_moon = 3.844e8;

dual zero(dual, dual, dual, dual)
{
    return 0.0;
}

dual f(dual, dual r, dual, dual)
{
    return 1.0 - 2.0 * G * M / (c * c * r);
}

dual g00(dual t, dual r, dual theta, dual phi)
{
    return c * c * f(t, r, theta, phi);
}

dual g11(dual t, dual r, dual theta, dual phi)
{
    return -1.0 / f(t, r, theta, phi);
}

dual g22(dual, dual r, dual, dual)
{
    return -r * r;
}

dual g33(dual, dual r, dual theta, dual)
{
    return -r * r * sin(theta) * sin(theta);
}

dual ig00(dual t, dual r, dual theta, dual phi)
{
    return 1.0 / (c * c * f(t, r, theta, phi));
}

dual ig11(dual t, dual r, dual theta, dual phi)
{
    return -f(t, r, theta, phi);
}

dual ig22(dual, dual r, dual, dual)
{
    return -1.0 / (r * r);
}

dual ig33(dual, dual r, dual theta, dual)
{
    return -1.0 / (r * r * sin(theta) * sin(theta));
}

void sda(double (&q)[8], double tau, double (&qp)[8], int n, const met& g, const met& invg)
{
    double t = q[0], r = q[1], theta = q[2], phi = q[3];
    double u[4] = { q[4], q[5], q[6], q[7] };

    dual dt = t, dr = r, dtheta = theta, dphi = phi;

    for (int i = 0; i < 4; ++i)
        qp[i] = u[i];

    for (int i = 0; i < 4; ++i) {
        double sum = 0.0;
        for (int k = 0; k < 4; ++k) {
            for (int l = 0; l < 4; ++l) {
                dual gam = christoffel(dt, dr, dtheta, dphi, i, k, l, g, invg);
                sum += gam.val * u[k] * u[l];
            }
        }
        qp[4 + i] = -sum;
    }
}

int main()
{
    met g(4 , vector<func>(4,zero) );
    met invg(4 , vector<func>(4,zero) );
    g[0][0] = g00;
    g[1][1] = g11;
    g[2][2] = g22;
    g[3][3] = g33;

    invg[0][0] = ig00;
    invg[1][1] = ig11;
    invg[2][2] = ig22;
    invg[3][3] = ig33;
    fstream rich;
    rich.open("dataforpy.txt", ios::out);
    rich << r_s << endl;
    rich << r_moon << endl;
    rich << r_earth << endl;
    rich.close();
    double r0    = r_moon; // don't put anything smaller than r_s please.
    double v = 1.01*sqrt(G * M / r0); // be logical about what you need: if light, put c; for anything else think (for a circular orbit and a timelike particle it's sqrt(G * M / r0))
    double uphi  = v / (sqrt(2)*r0);
    double uteta  = v / (sqrt(2)*r0);
    double dtau  = 3600;   // adapt it to what you need (for example, near r_s, (r_s / c) / 200.0 works)
    int    steps = 20000;
    double r_max = 2.0 * r_moon;

    geodesics(r0, M_PI/2.0, 0.0, 0.0, uteta, uphi,
              r_max, steps, dtau, "Geodesics.dat", g, invg, sda, 1.01 * r_s, true); // put false for light-like particles and true for time-like particles.

    return 0;
}