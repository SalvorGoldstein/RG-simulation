# General Relativity Library

> **Test branch, can break.** This branch explores a test `main` (`main.cc`), a `Makefile` to build it, and a small Python program (`python/Geodesics_topview.py`) to plot the results.

With this library you can compute the Christoffel symbols for any 4x4 metric and use them to code your own metric (there is an example C++ file showing how). You can also compute the path of a particle or a light ray with the `geodesics` function. Most of the calculation relies on the derivatives provided by [autodiff](https://github.com/autodiff/autodiff).

## What you need

- `g++` (C++17) and `make`
- the `autodiff` and `eigen-3.4.0` folders in `rgsim/third_party/`
- Python with `numpy` and `matplotlib` (`pip install numpy matplotlib`)

## How to download

```bash
git clone -b <branch-name> https://github.com/SalvorGoldstein/RG-simulation.git
```

Replace `<branch-name>` with the name of this test branch.

## How to use it

Your main file must be placed in `RG-simulation`, next to the library folder `rgsim`. It must include the library header:

```cpp
#include <rgsim/rgsim.hpp>
```

To use the example metrics, also add:

```cpp
#include "example_metric.h"
```

You can then write your code in your main file. The main function is `geodesics`, and its arguments are:

```cpp
geodesics(double x1, double x2, double x3,
          double u1, double u2, double u3,
          int steps, double dtau,
          const std::string& filename,
          const met& g, const christ& Gam,
          const rhs_fn& sd = geo,
          bool type = true,
          double x1max = infinity, double x2max = infinity, double x3max = infinity,
          double x1min = minusinfinity, double x2min = minusinfinity, double x3min = minusinfinity);
```

- `type = true` computes a timelike path (a massive particle), `type = false` a lightlike path (a light ray).
- `x1max`, `x2max`, `x3max`, `x1min`, `x2min`, `x3min` stop the integration when the path leaves the given range. By default there is no limit.
- For a black hole, set `x1min` to the Schwarzschild radius, so the integration stops at the horizon.

The parameters (mass, spin and charge) of the Schwarzschild-type metrics are set in `example_metric.h`.

## How to compile and run

```bash
make            # compile, produces ./main
make execute    # compile and run ./main
make aff        # plot the results with Python
make clean      # remove build files
```

The test `main` asks for the number of geodesic layers, the starting distance and the spacing between layers (in units of r_s), then writes the paths to `python/data/`. `make aff` plots them.

Without `make`:

```bash
g++ -std=c++17 -O2 -I . -I rgsim/third_party/autodiff -I rgsim/third_party/eigen-3.4.0 main.cc example_metric.cc -o main
```

## Setting up your metric

At the beginning of your `main` function you need the following:

```cpp
met g(4, vector<func>(4));
met invg(4, vector<func>(4));

for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
        g[i][j] = metric_name(i, j);

invg = inverseMet(g);

christ Gam(4, vector<vector<func>>(4, vector<func>(4)));
for (int i = 0; i < 4; ++i)
    for (int k = 0; k < 4; ++k)
        for (int l = 0; l < 4; ++l)
            Gam[i][k][l] = christoffel(i, k, l, g, invg);
```

Replace `metric_name` with the name of your metric function.

## Coding your own metric

You can code your metric with the same technique as the Schwarzschild example:

```cpp
func schwarzschild(int i, int j)
{
    if (i > j) swap(i, j);

    if (i == 0 && j == 0)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = c*c*(1.0 - r_s/r);
            return res;
        };
    if (i == 1 && j == 1)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -1.0/(1.0 - r_s/r);
            return res;
        };
    if (i == 2 && j == 2)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -r*r;
            return res;
        };
    if (i == 3 && j == 3)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -r*r*sin(th)*sin(th);
            return res;
        };
    return constant(0.0);
}
```

## Credits

- Derivatives: [autodiff](https://github.com/autodiff/autodiff) (MIT)
- Linear algebra headers: [Eigen](https://eigen.tuxfamily.org) (MPL 2.0)
