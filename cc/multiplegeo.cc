#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <string>
#include "geodesic.h"
#include "funcoperator.h"
#include "metricoperator.h"
#include "metric.h"

using namespace std;

void sda(double (&q)[8], double tau, double (&qp)[8], int n, const christ& Gam)
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

int main()
{
    met g(4, vector<func>(4));
    met invg(4, vector<func>(4));

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            g[i][j] = schwarzschild(i, j);        // or kerr(i, j)

    invg = inverseMet(g);

    christ Gam(4, vector<vector<func>>(4, vector<func>(4)));
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 4; ++k)
            for (int l = 0; l < 4; ++l)
                Gam[i][k][l] = christoffel(i, k, l, g, invg);

    fstream rich;
    rich.open("dataforpy.txt", ios::out);
    rich << r_s << endl;
    rich << r_moon << endl;
    rich << r_earth << endl;
    rich << J << endl;
    rich.close();

    double v     = c;
    
    double uteta =0;
    double dtau  = (r_s / c) / 200.0;
    int    steps = 5000;
    double r_max = 4 * r_s;


    for(double i = 1.01; i< 1.9 ; i = i+0.01){
        double uphi  = v / (r_s*i);
    geodesics(r_s*i, M_PI/2.0, 0.0, 0.0, uteta, uphi,
              r_max, steps, dtau, "data/Geodesics"+to_string(i)+"plus.dat", g, Gam, sda, r_s, false);
    }
    for(double j = 1.01; j< 1.9 ; j = j+0.01){
    double uphi  = -v / (r_s*j);
    geodesics(r_s*j, M_PI/2.0, 0.0, 0.0, uteta, uphi,
              r_max, steps, dtau, "data/Geodesics"+to_string(j)+"minus.dat", g, Gam, sda,  r_s, false);
    }
    return 0;
}