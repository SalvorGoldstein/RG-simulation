#pragma once
#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include "types.hpp"

namespace rgsim {

// ---------- derivatives and Christoffel symbols ----------
inline dual partial(const func& fn,
                    dual t, dual r, dual theta, dual phi,
                    int coordinate)
{
    using autodiff::at; using autodiff::derivative; using autodiff::wrt;
    switch (coordinate) {
        case 0: return derivative(fn, wrt(t),     at(t, r, theta, phi));
        case 1: return derivative(fn, wrt(r),     at(t, r, theta, phi));
        case 2: return derivative(fn, wrt(theta), at(t, r, theta, phi));
        case 3: return derivative(fn, wrt(phi),   at(t, r, theta, phi));
        default: return 0.0;
    }
}

inline func christoffel(int i, int k, int l,
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
inline void euler(const rhs_fn& sd,
                  double (&q)[8], double t, double dt, int n,
                  const christ& Gam)
{
    double qp[8];
    sd(q, t, qp, n, Gam);
    for (int i = 0; i < n; i++)
        q[i] += dt * qp[i];
}

// Classical 4th-order Runge-Kutta step. q = (x^mu, u^mu), advanced by dt.
inline void rk4(const rhs_fn& sd, double (&q)[8], double t, double dt, int n,
                const christ& Gam)
{
    double k1[8], k2[8], k3[8], k4[8], tmp[8];

    sd(q, t, k1, n, Gam);

    for (int i = 0; i < n; i++) tmp[i] = q[i] + 0.5 * dt * k1[i];
    sd(tmp, t + 0.5 * dt, k2, n, Gam);

    for (int i = 0; i < n; i++) tmp[i] = q[i] + 0.5 * dt * k2[i];
    sd(tmp, t + 0.5 * dt, k3, n, Gam);

    for (int i = 0; i < n; i++) tmp[i] = q[i] + dt * k3[i];
    sd(tmp, t + dt, k4, n, Gam);

    for (int i = 0; i < n; i++)
        q[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

// ---------- initial conditions (any symmetric metric) ----------
// Assumes signature (+,-,-,-): g_{mu nu} u^mu u^nu = kappa, kappa = c^2 or 0.
inline double timelikeU(double t, double r, double theta, double phi,
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
    double C = -constants::c * constants::c;    // kappa = c^2
    for (int i = 1; i < 4; ++i) {
        B += 2.0 * gm[0][i] * u[i];
        for (int j = 1; j < 4; ++j)
            C += gm[i][j] * u[i] * u[j];
    }

    double D = B * B - 4.0 * A * C;
    if (D < 0.0) return NAN;                    // caught by isfinite in geodesics()
    return (-B + std::sqrt(D)) / (2.0 * A);
}

inline double lightlikeU(double t, double r, double theta, double phi,
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
        std::cerr << "No lightlike solution for k^t at this point/direction!" << std::endl;
        return NAN;
    }
    return (-B + std::sqrt(D)) / (2.0 * A);
}

inline void geo(double (&q)[8], double tau, double (&qp)[8], int n, const christ& Gam)
{
    double u[4] = { q[4], q[5], q[6], q[7] };
    dual dt = q[0], dr = q[1], dtheta = q[2], dphi = q[3];

    for (int i = 0; i < 4; ++i)
        qp[i] = u[i];

    for (int i = 0; i < 4; ++i) {
        double sum = 0.0;
        for (int k = 0; k < 4; ++k)
            for (int l = 0; l < 4; ++l) {
                dual gam = Gam[i][k][l](dt, dr, dtheta, dphi);
                sum += gam.val * u[k] * u[l];
            }
        qp[4 + i] = -sum;
    }
}

inline const double infinity = std::numeric_limits<double>::infinity();
inline const double minusinfinity = std::numeric_limits<double>::lowest();

// ---------- main routine ----------
inline void geodesics(double x1, double x2, double x3,
                      double u1, double u2, double u3,
                      int steps, double dtau,
                      const std::string& filename,
                      const met& g, const christ& Gam,
                      const rhs_fn& sd = geo,
                      bool type = true,
                      double x1max = infinity, double x2max = infinity, double x3max = infinity,
                      double x1min = minusinfinity, double x2min = minusinfinity, double x3min = minusinfinity)
{
    const int n = 8;
    std::ofstream fich(filename);
    fich << std::scientific << std::setprecision(10);

    double q[8] = { 0.0, x1, x2, x3, 0.0, u1, u2, u3 };
    double ut;
    if (type == true) {
        ut = timelikeU(q[0], q[1], q[2], q[3], u1, u2, u3, g);
    }
    else {
        ut = lightlikeU(q[0], q[1], q[2], q[3], u1, u2, u3, g);
    }
    if (!std::isfinite(ut)) {
        std::cerr << "invalid initial velocity" << std::endl;
        return;
    }
    q[4] = ut;

    double tau = 0.0;
    for (int step = 0; step < steps; ++step) {
        if (!std::isfinite(q[1]) || !std::isfinite(q[2]) || !std::isfinite(q[3])
            || q[1] < x1min || q[1] > x1max
            || q[2] < x2min || q[2] > x2max
            || q[3] < x3min || q[3] > x3max)
            break;
        fich << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << '\n';
        rk4(sd, q, tau, dtau, n, Gam);
        tau += dtau;
    }
}

} // namespace rgsim