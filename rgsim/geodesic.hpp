#pragma once
#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <array>
#include "types.hpp"

namespace rgsim {

// ---------- derivatives and Christoffel symbols ----------
// d fn / d x[coordinate], obtained by seeding that coordinate's derivative part.
inline dual partial(const func& fn,
                    std::array<dual,N> x,
                    int coordinate)
{
    for (int a = 0; a < N; ++a) x[a].grad = 0.0;
    x[coordinate].grad = 1.0;
    dual y = fn(x);
    dual res = y.grad;
    return res;
}

inline func christoffel(int i, int k, int l,
                        const met& g, const met& invg)
{
    return [i, k, l, g, invg](const std::array<dual,N>& x)
    {
        dual sum = 0.0;
        for (int m = 0; m < N; ++m) {
            const dual dg_mk_dl = partial(g[m][k], x, l);
            const dual dg_ml_dk = partial(g[m][l], x, k);
            const dual dg_kl_dm = partial(g[k][l], x, m);

            sum += invg[i][m](x)
                 * (dg_mk_dl + dg_ml_dk - dg_kl_dm);
        }
        dual res = 0.5 * sum;
        return res;
    };
}

// ---------- equations of motion ----------
// Here x0 is the integration parameter (proper time / affine parameter)
// and dx0 the step, as in the original (t, dt).
inline void euler(const rhs_fn& sd,
                  double (&q)[2*N], double x0, double dx0, int n,
                  const christ& Gam)
{
    double qp[2*N];
    sd(q, x0, qp, n, Gam);
    for (int i = 0; i < n; i++)
        q[i] += dx0 * qp[i];
}

// Classical 4th-order Runge-Kutta step. q = (x^mu, u^mu), advanced by dx0.
inline void rk4(const rhs_fn& sd, double (&q)[2*N], double x0, double dx0, int n,
                const christ& Gam)
{
    double k1[2*N], k2[2*N], k3[2*N], k4[2*N], tmp[2*N];

    sd(q, x0, k1, n, Gam);

    for (int i = 0; i < n; i++) tmp[i] = q[i] + 0.5 * dx0 * k1[i];
    sd(tmp, x0 + 0.5 * dx0, k2, n, Gam);

    for (int i = 0; i < n; i++) tmp[i] = q[i] + 0.5 * dx0 * k2[i];
    sd(tmp, x0 + 0.5 * dx0, k3, n, Gam);

    for (int i = 0; i < n; i++) tmp[i] = q[i] + dx0 * k3[i];
    sd(tmp, x0 + dx0, k4, n, Gam);

    for (int i = 0; i < n; i++)
        q[i] += dx0 / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

// ---------- initial conditions (any symmetric metric) ----------
// Assumes signature (+,-,-,-,-): g_{mu nu} u^mu u^nu = kappa, kappa = c^2 or 0.
inline double timelikeU(double x0, double x1, double x2, double x3, double x4,
                        double u1, double u2, double u3, double u4,
                        const met& g)
{
    std::array<dual,N> x = { x0, x1, x2, x3, x4 };

    double gm[N][N];
    for (int a = 0; a < N; ++a)
        for (int b = 0; b < N; ++b)
            gm[a][b] = g[a][b](x).val;

    double u[N] = { 0.0, u1, u2, u3, u4 };    // u[0] is the unknown

    double A = gm[0][0];
    double B = 0.0;
    double C = -constants::c * constants::c;    // kappa = c^2
    for (int i = 1; i < N; ++i) {
        B += 2.0 * gm[0][i] * u[i];
        for (int j = 1; j < N; ++j)
            C += gm[i][j] * u[i] * u[j];
    }

    double D = B * B - 4.0 * A * C;
    if (D < 0.0) return NAN;                    // caught by isfinite in geodesics()
    return (-B + std::sqrt(D)) / (2.0 * A);
}

inline double lightlikeU(double x0, double x1, double x2, double x3, double x4,
                         double u1, double u2, double u3, double u4,
                         const met& g)
{
    std::array<dual,N> x = { x0, x1, x2, x3, x4 };

    double gm[N][N];
    for (int a = 0; a < N; ++a)
        for (int b = 0; b < N; ++b)
            gm[a][b] = g[a][b](x).val;

    double u[N] = { 0.0, u1, u2, u3, u4 };

    double A = gm[0][0];
    double B = 0.0;
    double C = 0.0;                             // kappa = 0
    for (int i = 1; i < N; ++i) {
        B += 2.0 * gm[0][i] * u[i];
        for (int j = 1; j < N; ++j)
            C += gm[i][j] * u[i] * u[j];
    }

    double D = B * B - 4.0 * A * C;
    if (D < 0.0) {
        std::cerr << "No lightlike solution for k^0 at this point/direction!" << std::endl;
        return NAN;
    }
    return (-B + std::sqrt(D)) / (2.0 * A);
}

inline void geo(double (&q)[2*N], double tau, double (&qp)[2*N], int n, const christ& Gam)
{
    double u[N];
    std::array<dual,N> x;
    for (int i = 0; i < N; ++i) {
        x[i] = q[i];
        u[i] = q[N + i];
    }

    for (int i = 0; i < N; ++i)
        qp[i] = u[i];

    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int k = 0; k < N; ++k)
            for (int l = 0; l < N; ++l) {
                dual gam = Gam[i][k][l](x);
                sum += gam.val * u[k] * u[l];
            }
        qp[N + i] = -sum;
    }
}

inline const double infinity = std::numeric_limits<double>::infinity();
inline const double minusinfinity = std::numeric_limits<double>::lowest();

// ---------- main routine ----------
inline void geodesics(double x1, double x2, double x3, double x4,
                      double u1, double u2, double u3, double u4,
                      int steps, double dtau,
                      const std::string& filename,
                      const met& g, const christ& Gam,
                      const rhs_fn& sd = geo,
                      bool type = true,
                      double x1max = infinity, double x2max = infinity, double x3max = infinity, double x4max = infinity,
                      double x1min = minusinfinity, double x2min = minusinfinity, double x3min = minusinfinity, double x4min = minusinfinity)
{
    const int n = 2 * N;
    std::ofstream fich(filename);
    fich << std::scientific << std::setprecision(10);

    double q[2*N] = { 0.0, x1, x2, x3, x4, 0.0, u1, u2, u3, u4 };
    double u0;
    if (type == true) {
        u0 = timelikeU(q[0], q[1], q[2], q[3], q[4], u1, u2, u3, u4, g);
    }
    else {
        u0 = lightlikeU(q[0], q[1], q[2], q[3], q[4], u1, u2, u3, u4, g);
    }
    if (!std::isfinite(u0)) {
        std::cerr << "invalid initial velocity" << std::endl;
        return;
    }
    q[N] = u0;

    double tau = 0.0;
    for (int step = 0; step < steps; ++step) {
        if (!std::isfinite(q[1]) || !std::isfinite(q[2]) || !std::isfinite(q[3]) || !std::isfinite(q[4])
            || q[1] < x1min || q[1] > x1max
            || q[2] < x2min || q[2] > x2max
            || q[3] < x3min || q[3] > x3max
            || q[4] < x4min || q[4] > x4max)
            break;
        fich << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " " << q[4] << '\n';
        rk4(sd, q, tau, dtau, n, Gam);
        tau += dtau;
    }
}

} // namespace rgsim