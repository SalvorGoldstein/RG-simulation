#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <autodiff/forward/dual.hpp>

using autodiff::dual;
using autodiff::at;
using autodiff::derivative;
using autodiff::wrt;
using namespace std;

constexpr double G = 6.67430e-11;
constexpr double M = 5.972e24;
constexpr double c = 299792458.0;
constexpr double r_s = 2.0 * G * M / (c * c);

using GFunc = dual (*)(dual, dual, dual, dual);

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

dual partial(GFunc function,
             dual t, dual r, dual theta, dual phi,
             int coordinate)
{
    switch(coordinate) {
        case 0:
            return derivative(function, wrt(t), at(t, r, theta, phi));
        case 1:
            return derivative(function, wrt(r), at(t, r, theta, phi));
        case 2:
            return derivative(function, wrt(theta), at(t, r, theta, phi));
        case 3:
            return derivative(function, wrt(phi), at(t, r, theta, phi));
        default:
            return 0.0;
    }
}

dual christoffel(dual t, dual r, dual theta, dual phi,
                 int i, int k, int l)
{
    GFunc g[4][4];
    GFunc invg[4][4];

    for(int a = 0; a < 4; ++a) {
        for(int b = 0; b < 4; ++b) {
            g[a][b] = zero;
            invg[a][b] = zero;
        }
    }

    g[0][0] = g00;
    g[1][1] = g11;
    g[2][2] = g22;
    g[3][3] = g33;

    invg[0][0] = ig00;
    invg[1][1] = ig11;
    invg[2][2] = ig22;
    invg[3][3] = ig33;

    dual sum = 0.0;

    for(int m = 0; m < 4; ++m) {
        const dual dg_mk_dl = partial(g[m][k], t, r, theta, phi, l);
        const dual dg_ml_dk = partial(g[m][l], t, r, theta, phi, k);
        const dual dg_kl_dm = partial(g[k][l], t, r, theta, phi, m);

        sum += invg[i][m](t, r, theta, phi)
             * (dg_mk_dl + dg_ml_dk - dg_kl_dm);
    }

    return 0.5 * sum;
}


void sda(double (&q)[8], double tau, double (&qp)[8], int n)
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
                dual g = christoffel(dt, dr, dtheta, dphi, i, k, l);
                sum += g.val * u[k] * u[l];
            }
        }
        qp[4 + i] = -sum;
    }
}


void euler(void(sd)(double(&)[8], double, double(&)[8], int),
           double (&q)[8], double t, double dt, int n)
{
    double qp[8];
    sd(q, t, qp, n);
    for (int i = 0; i < n; i++)
        q[i] = q[i] + dt * qp[i];
}
double timelikeU(double t, double r, double theta, double phi,
                          double ur, double utheta, double uphi){
                            dual dt = t, dr = r, dtheta = theta, dphi = phi;
    double G00 = g00(dt, dr, dtheta, dphi).val;
    double G11 = g11(dt, dr, dtheta, dphi).val;
    double G22 = g22(dt, dr, dtheta, dphi).val;
    double G33 = g33(dt, dr, dtheta, dphi).val;

    double rhs = -(G11*ur*ur + G22*utheta*utheta + G33*uphi*uphi);
    double ut2 = rhs / G00;
    return sqrt(ut2);
                          }

void geodesics(double r , double teta, double phi, double ur , double uteta, double uphi){
     cout << scientific << setprecision(10);
    const int n = 8;    
    fstream fich;
    fich.open("Geodesics.txt", ios::out);
    double q[8];
    q[0] = 0.0;             
    q[1] = r;        
    q[2] = teta ; 
    q[3] = phi; 
    double ut = timelikeU(q[0], q[1], q[2], q[3], ur, uteta, uphi);
    q[4] = ut;
    q[5] = ur;
    q[6] = utheta;
    q[7] = uphi;

   double  tau = 0.0;
    int steps = 1000000;
    double r_scale = q[1] / c;
    double dtau =r_scale / 1000.0;
    dual gamma[4][4][4];
    {
        dual t0 = q[0], r0 = q[1], theta0 = q[2], phi0 = q[3];
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    gamma[i][j][k] = christoffel(t0, r0, theta0, phi0, i, j, k);
    }
}
int main()
{
    
    

        
                             
  





    for (int n_step = 0; n_step < steps; ++n_step) {
    
        if(q[1] < r_s || q[0] < 0|| q[1] > 10*r_s){
            cout << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
                fich.close();
                break;
                return 0;
        }
        else{
        cout << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
        fich  << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << endl;
        euler(sda, q, tau, dtau, n);
        tau += dtau;}
    }

    fich.close();
    return 0;
}