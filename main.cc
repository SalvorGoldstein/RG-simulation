#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <rgsim/rgsim.hpp>
#include "example_metric.h"

using namespace std;

int main()
{
    int npre;
    double j0, dj;

    cout << "how many layers of geodesics do you want? ";
    cin >> npre;
    cout << "starting distance (in units of r_s, must be > 1): ";
    cin >> j0;
    cout << "distance between each layer (in units of r_s): ";
    cin >> dj;
    cout << "perfect!" << endl;

    if (npre < 1 || j0 <= 1.0 || dj < 0.0) {
        cerr << "invalid input: need layers >= 1, start > 1, spacing >= 0" << endl;
        return 1;
    }

    int n = 2 * npre;
    met g(4, vector<func>(4));
    met invg(4, vector<func>(4));

    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 4; ++k)
            g[i][k] = kerr_newman(i, k);        // or schwarzschild(i, k)

    invg = inverseMet(g);

    christ Gam(4, vector<vector<func>>(4, vector<func>(4)));
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 4; ++k)
            for (int l = 0; l < 4; ++l)
                Gam[i][k][l] = christoffel(i, k, l, g, invg);

    fstream rich;
    rich.open("python/dataforpy.dat", ios::out);
    rich << r_s << endl;
    rich << r_moon << endl;
    rich << r_earth << endl;
    rich << J << endl;
    rich << n << endl;
    rich.close();

    double v     = c;
    double dtau  = (r_s / c) / 200.0;
    int    steps = 10000;

    // largest starting distance, so r_max is always beyond every start point
    double j_last = j0 + (npre - 1) * dj;
    double r_max  = max(7.0, j_last + 1.0) * r_s;

    double j = j0;
    for (int i = 0; i < n; i++) {
        double uphi   = v / (r_s * j * sqrt(2));
        double utheta = 0;   // v / (r_s * j);
        double ur     = -v / sqrt(2);

        if (i % 2) {
            geodesics(r_s * j, pi / 2.0, 0.0, ur, utheta, uphi,
                      steps, dtau,
                      "python/data/Geodesics" + to_string(i) + ".dat",
                      g, Gam, geo, false,
                      r_max, infinity, infinity,   // x1max, x2max, x3max
                      r_s);                        // x1min
            j = j + dj;
        }
        else {
            geodesics(r_s * j, pi / 2.0, 0.0, ur, -utheta, -uphi,
                      steps, dtau,
                      "python/data/Geodesics" + to_string(i) + ".dat",
                      g, Gam, geo, false,
                      r_max, infinity, infinity,
                      r_s);
        }
    }
    return 0;
}