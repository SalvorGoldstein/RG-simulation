# Geodesics-tracer

Computes the path of a particle or a light ray in a curved spacetime and plots
it in 3D. You write the metric once, and the program builds everything else
from it: the inverse metric, the Christoffel symbols and the equations of
motion. Derivatives are computed with
[autodiff](https://github.com/autodiff/autodiff).

The metric can be any symmetric 4x4 metric (diagonal or not) in any
coordinates `(x0, x1, x2, x3)`. Schwarzschild and Kerr (Boyer-Lindquist
coordinates `(t, r, theta, phi)`) are included, and you can define your own in
`main.cc`.

## What you need

- `g++` (C++17) and `make`
- the `autodiff` and `eigen-3.4.0` folders in `cc/third_party/`
- Python with `numpy` and `matplotlib` (`pip install numpy matplotlib`)

## How to download

Clone the repo with :

```
git clone https://github.com/SalvorGoldstein/RG-simulation.git
```
It should run fine.

## How to run

From the project folder (the one containing the `Makefile`):

```
make execute aff
```

This compiles the program, runs it (it writes `Geodesics.dat` and
`dataforpy.txt` in the project folder), then plots the result.

| Command | What it does |
|---|---|
| `make` | compile only |
| `make execute` | compile and run |
| `make aff` | run if needed, then plot |
| `make clean` | remove compiled files |

Always run from the project folder, because the program and the Python script
read and write their files there.

If you have Windows, well, switch to Linux. It's the best troubleshooting advice I have.

## How the code is organised

All the C++ is in `cc/`.

| File | What it does |
|---|---|
| `main.cc` | choice of metric, initial conditions, equations of motion (`sda`), call to `geodesics` |
| `types.h` | the types `func`, `met`, `christ` and the constants `G`, `c` |
| `metric.h/.cc` | the constants of the problem (`M`, `r_s`, `J`, Kerr spin...) and the ready-made metrics `schwarzschild` and `kerr` |
| `geodesic.h/.cc` | derivatives (`partial`), Christoffel symbols, Euler step, initial velocity, main loop (`geodesics`) |
| `funcoperator.h/.cc` | lets you add, subtract, multiply and divide functions (`f + g`, `f * g`, ...) and make constant functions |
| `metricoperator.h/.cc` | determinant and inverse of a matrix of functions |

The main types are:

- `func`: a function of the four coordinates `(t, r, theta, phi)`
- `met`: a 4x4 table of `func` (the metric `g`, or its inverse `invg`)
- `christ`: a 4x4x4 table of `func` (the Christoffel symbols)

## Change the simulation

**Initial conditions.** Edit the values in `main()` in `cc/main.cc`: the
starting radius `r0`, the initial velocity, the step `dtau` and the number of
`steps`. Put `true` at the end of the `geodesics(...)` call for a massive
particle and `false` for light. If the plot is empty, `dtau` is probably too
big.

**Constants.** The mass `M` and the Kerr `spin` (from 0 to 1) are in
`cc/metric.h`.

**Choose a metric.** In `main()`, the loop that fills `g` calls one metric
function. Change its name:

```cpp
g[i][j] = schwarzschild(i, j);   // or kerr(i, j), or mymetric(i, j)
```

**Define your own metric.** Write a function in `main.cc` with the same shape
as `schwarzschild` in `cc/metric.cc`. Each `if` gives one component, and every
component you don't write is 0. You only write `i <= j`, the symmetric one is
filled for you.

```cpp
func mymetric(int i, int j)
{
    if (i > j) swap(i, j);

    if (i == 1 && j == 1)
        return [](dual t, dual r, dual th, dual ph) {
            dual res = -1.0/(1.0 - r_s/r);
            return res;
        };

    return constant(0.0);
}
```

Keep the `dual res = ...; return res;` form inside every function.

Non-diagonal terms (like `g03` in Kerr) just need one more `if`, the inverse
and the initial `u^t` handle them automatically. For FLRW, write `a(t)` with
`dual` operations (`pow`, `sinh`, ...) so it can be differentiated.

## Things to keep in mind

- The stop conditions in `geodesics` (`r_min`, `r_max`, `t < 0`) assume that
  `q[1]` is a radius. In other coordinates, edit that line in `cc/geodesic.cc`.
- For Kerr, the outer horizon is at `r+ = r_s (1 + sqrt(1 - spin^2)) / 2`, so
  keep `r_min` above it (for example `1.05 * r+`). `1.01 * r_s` is for
  Schwarzschild only.
- Inside the Kerr ergoregion, `g00` can be negative and the initial `u^t` may
  need a manual check. With `spin = 0`, `kerr` must give the same result as
  `schwarzschild`.
- At a singular point (horizon, poles where `sin(theta) = 0`) the result
  becomes `inf` or `nan`. Don't start at `theta = 0` or `pi`.
- The integrator is a plain Euler method, so use a small `dtau`.
- The inverse metric is built as a chain of functions, so it is slow. This is
  fine for a few tens of thousands of steps.
  ## Credits

- Derivatives: [autodiff](https://github.com/autodiff/autodiff) (MIT)
- Linear algebra headers: [Eigen](https://eigen.tuxfamily.org) (MPL 2.0)