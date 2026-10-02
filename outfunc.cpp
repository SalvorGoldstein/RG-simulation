#include "outfunc.h"
#include <iostream>
#include <iomanip>
#include <fstream>

using autodiff::at;
using autodiff::derivative;
using autodiff::wrt;
using namespace std;


// ---------- derivatives and Christoffel symbols ----------
dual partial(const func& fn,
             dual t, dual r, dual theta, dual phi,
             int coordinate)
{
    switch (coordinate) {
        case 0: return derivative(fn, wrt(t),     at(t, r, theta, phi));
        case 1: return derivative(fn, wrt(r),     at(t, r, theta, phi));
        case 2: return derivative(fn, wrt(theta), at(t, r, theta, phi));
        case 3: return derivative(fn, wrt(phi),   at(t, r, theta, phi));
        default: return 0.0;
    }
}

dual christoffel(dual t, dual r, dual theta, dual phi,
                 int i, int k, int l,
                 const met& g, const met& invg)
{
    dual sum = 0.0;
    for (int m = 0; m < 4; ++m) {
        const dual dg_mk_dl = partial(g[m][k], t, r, theta, phi, l);
        const dual dg_ml_dk = partial(g[m][l], t, r, theta, phi, k);
        const dual dg_kl_dm = partial(g[k][l], t, r, theta, phi, m);

        sum += invg[i][m](t, r, theta, phi)
             * (dg_mk_dl + dg_ml_dk - dg_kl_dm);
    }
    return 0.5 * sum;
}

// ---------- equations of motion ----------


void euler(rhs_fn sd,
           double (&q)[8], double t, double dt, int n,
           const met& g, const met& invg)
{
    double qp[8];
    sd(q, t, qp, n, g, invg);
    for (int i = 0; i < n; i++)
        q[i] += dt * qp[i];
}

// ---------- initial conditions ----------
double timelikeU(double t, double r, double theta, double phi,
                 double ur, double utheta, double uphi,
                 const met& g)
{
    dual dt = t, dr = r, dtheta = theta, dphi = phi;
    double G00 = g[0][0](dt, dr, dtheta, dphi).val;
    double G11 = g[1][1](dt, dr, dtheta, dphi).val;
    double G22 = g[2][2](dt, dr, dtheta, dphi).val;
    double G33 = g[3][3](dt, dr, dtheta, dphi).val;

    double rhs = c * c - (G11 * ur * ur + G22 * utheta * utheta + G33 * uphi * uphi);
    return sqrt(rhs / G00);   // NaN if invalid, caught by isfinite in geodesics()
}

double lightlikeU(double t, double r, double theta, double phi,
                  double ur, double utheta, double uphi,
                  const met& g)
{
    dual dt = t, dr = r, dtheta = theta, dphi = phi;
    double G00 = g[0][0](dt, dr, dtheta, dphi).val;
    double G11 = g[1][1](dt, dr, dtheta, dphi).val;
    double G22 = g[2][2](dt, dr, dtheta, dphi).val;
    double G33 = g[3][3](dt, dr, dtheta, dphi).val;

    double ut2 = -(G11 * ur * ur + G22 * utheta * utheta + G33 * uphi * uphi) / G00;
    if (ut2 < 0.0) {
        cerr << "No lightlike solution for k^t at this point/direction!" << endl;
        return NAN;
    }
    return sqrt(ut2);
}

// ---------- main routine ----------
void geodesics(double r, double theta0, double phi0,
               double ur, double utheta, double uphi,
               double r_max, int steps, double dtau,
               const string& filename,
               const met& g, const met& invg,
               rhs_fn sd, double r_min,
               bool type)                   
{
    const int n = 8;
    ofstream fich(filename);
    fich << scientific << setprecision(10);

    double q[8] = { 0.0, r, theta0, phi0, 0.0, ur, utheta, uphi };

    double ut = type ? timelikeU (q[0], q[1], q[2], q[3], ur, utheta, uphi, g)
                 : lightlikeU(q[0], q[1], q[2], q[3], ur, utheta, uphi, g);
    if (!isfinite(ut)) {
        cerr << "invalid initial velocity" << endl;
        return;
    }
    q[4] = ut;

    double tau = 0.0;
    for (int step = 0; step < steps; ++step) {
        if (q[1] < r_min|| q[0] < 0.0 || q[1] > r_max)
            break;
        fich << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << '\n';
        euler(sd, q, tau, dtau, n, g, invg);
        tau += dtau;
    }
}