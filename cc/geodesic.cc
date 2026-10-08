#include "geodesic.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

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

func christoffel(int i, int k, int l,
                 const met& g, const met& invg)
{
    return [i, k, l, g, invg](dual t, dual r, dual theta, dual phi)
    {
        dual sum = 0.0;
        for (int m = 0; m < 4; ++m) {
            const dual dg_mk_dl = partial(g[m][k], t, r, theta, phi, l);
            const dual dg_ml_dk = partial(g[m][l], t, r, theta, phi, k);
            const dual dg_kl_dm = partial(g[k][l], t, r, theta, phi, m);

            sum += invg[i][m](t, r, theta, phi)
                 * (dg_mk_dl + dg_ml_dk - dg_kl_dm);
        }
        dual res = 0.5 * sum;
        return res;
    };
}

// ---------- equations of motion ----------
void euler(rhs_fn sd,
           double (&q)[8], double t, double dt, int n,
           const christ& Gam)
{
    double qp[8];
    sd(q, t, qp, n, Gam);
    for (int i = 0; i < n; i++)
        q[i] += dt * qp[i];
}

// ---------- initial conditions (any symmetric metric) ----------
double timelikeU(double t, double r, double theta, double phi,
                 double ur, double utheta, double uphi,
                 const met& g)
{
    dual dt = t, dr = r, dtheta = theta, dphi = phi;

    double gm[4][4];
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            gm[a][b] = g[a][b](dt, dr, dtheta, dphi).val;

    double u[4] = { 0.0, ur, utheta, uphi };    // u[0] is the unknown

    double A = gm[0][0];
    double B = 0.0;
    double C = -c * c;                          // kappa = c^2
    for (int i = 1; i < 4; ++i) {
        B += 2.0 * gm[0][i] * u[i];
        for (int j = 1; j < 4; ++j)
            C += gm[i][j] * u[i] * u[j];
    }

    double D = B * B - 4.0 * A * C;
    return sqrt(D) / (2.0 * A) - B / (2.0 * A); // NaN if D < 0, caught by isfinite in geodesics()
}

double lightlikeU(double t, double r, double theta, double phi,
                  double ur, double utheta, double uphi,
                  const met& g)
{
    dual dt = t, dr = r, dtheta = theta, dphi = phi;

    double gm[4][4];
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            gm[a][b] = g[a][b](dt, dr, dtheta, dphi).val;

    double u[4] = { 0.0, ur, utheta, uphi };

    double A = gm[0][0];
    double B = 0.0;
    double C = 0.0;                             // kappa = 0
    for (int i = 1; i < 4; ++i) {
        B += 2.0 * gm[0][i] * u[i];
        for (int j = 1; j < 4; ++j)
            C += gm[i][j] * u[i] * u[j];
    }

    double D = B * B - 4.0 * A * C;
    if (D < 0.0) {
        cerr << "No lightlike solution for k^t at this point/direction!" << endl;
        return NAN;
    }
    return (-B + sqrt(D)) / (2.0 * A);
}

// ---------- main routine ----------
void geodesics(double r, double theta0, double phi0,
               double ur, double utheta, double uphi,
               double r_max, int steps, double dtau,
               const string& filename,
               const met& g, const christ& Gam,
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
        if (!isfinite(q[1]) || q[1] < r_min || q[0] < 0.0 || q[1] > r_max)
            break;
        fich << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << '\n';
        euler(sd, q, tau, dtau, n, Gam);
        tau += dtau;
    }
}