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
double J =-G*M*M/c;
dual zero(dual, dual, dual, dual)
{
    return 0.0;
}

dual g00(dual t, dual r, dual theta, dual phi)
{
    return c * c - 2.0 * G * M * r / (r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c));
}

dual g03(dual t, dual r, dual theta, dual phi)
{
    return 2.0 * G * J * r * sin(theta) * sin(theta)
         / (c * c * (r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c)));
}

dual g11(dual t, dual r, dual theta, dual phi)
{
    return -(r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c))
         / (r * r - 2.0 * G * M * r / (c * c) + J * J / (M * M * c * c));
}

dual g22(dual t, dual r, dual theta, dual phi)
{
    return -(r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c));
}

dual g33(dual t, dual r, dual theta, dual phi)
{
    return -(r * r + J * J / (M * M * c * c)
           + 2.0 * G * J * J * r * sin(theta) * sin(theta)
             / (M * c * c * c * c * (r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c))))
           * sin(theta) * sin(theta);
}

dual ig00(dual t, dual r, dual theta, dual phi)
{
    return ((r * r + J * J / (M * M * c * c)) * (r * r + J * J / (M * M * c * c))
           - (r * r - 2.0 * G * M * r / (c * c) + J * J / (M * M * c * c))
             * (J * J / (M * M * c * c)) * sin(theta) * sin(theta))
         / (c * c * (r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c))
                  * (r * r - 2.0 * G * M * r / (c * c) + J * J / (M * M * c * c)));
}

dual ig03(dual t, dual r, dual theta, dual phi)
{
    return 2.0 * G * J * r
         / (c * c * c * c * (r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c))
                          * (r * r - 2.0 * G * M * r / (c * c) + J * J / (M * M * c * c)));
}

dual ig11(dual t, dual r, dual theta, dual phi)
{
    return -(r * r - 2.0 * G * M * r / (c * c) + J * J / (M * M * c * c))
         / (r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c));
}

dual ig22(dual t, dual r, dual theta, dual phi)
{
    return -1.0 / (r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c));
}

dual ig33(dual t, dual r, dual theta, dual phi)
{
    return -((r * r - 2.0 * G * M * r / (c * c) + J * J / (M * M * c * c))
           - (J * J / (M * M * c * c)) * sin(theta) * sin(theta))
         / ((r * r + J * J * cos(theta) * cos(theta) / (M * M * c * c))
          * (r * r - 2.0 * G * M * r / (c * c) + J * J / (M * M * c * c))
          * sin(theta) * sin(theta));
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
    g[0][3] = g03;
    g[3][0] = g03;

    invg[0][0] = ig00;
    invg[1][1] = ig11;
    invg[2][2] = ig22;
    invg[3][3] = ig33;
    invg[0][3] = ig03;
     invg[3][0] = ig03;
    fstream rich;
    rich.open("dataforpy.txt", ios::out);
    rich << r_s << endl;
    rich << r_moon << endl;
    rich << r_earth << endl;
    rich.close();
    double r0    = r_moon; // don't put anything smaller than r_s please.
    double v = sqrt(G * M / r0); // be logical about what you need: if light, put c; for anything else think (for a circular orbit and a timelike particle it's sqrt(G * M / r0))
    double uphi  = v / (sqrt(2)*r0);
    double uteta  = v / (sqrt(2)*r0);
    double dtau  = 3600/2;   // adapt it to what you need (for example, near r_s, (r_s / c) / 200.0 works)
    int    steps = 100000;
    double r_max = 2*r_moon;

    geodesics(r0, M_PI/2.0, 0.0, 0.0, uteta, uphi,
              r_max, steps, dtau, "Geodesics.dat", g, invg, sda, 1.01 * r_s, true); // put false for light-like particles and true for time-like particles.

    return 0;
}