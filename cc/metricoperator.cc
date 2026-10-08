#include "metricoperator.h"
#include "funcoperator.h"     // constant, +, -, *, /

using namespace std;



met minorOf(const met& m, int row, int col)
{
    int n = m.size();
    met s(n - 1, vector<func>(n - 1));
    for (int i = 0, a = 0; i < n; ++i) {
        if (i == row) continue;
        for (int j = 0, b = 0; j < n; ++j) {
            if (j == col) continue;
            s[a][b++] = m[i][j];
        }
        ++a;
    }
    return s;
}

func determinant(const met& m)
{
    int n = m.size();
    if (n == 1) return m[0][0];
    func sum = constant(0.0);
    for (int j = 0; j < n; ++j) {
        func term = m[0][j] * determinant(minorOf(m, 0, j));
        if (j % 2 == 0) sum = sum + term;
        else            sum = sum - term;
    }
    return sum;
}


met inverseMet(const met& g)
{
    int n = g.size();
    func det = determinant(g);                 // built once, shared
    met inv(n, vector<func>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            func cof = determinant(minorOf(g, j, i));   // transpose built in: (j, i)
            if ((i + j) % 2 == 1) cof = constant(0.0) - cof;
            inv[i][j] = cof / det;
        }
    return inv;
}