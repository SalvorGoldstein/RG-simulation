# Geodesics-tracer

Computes the path of a particle or a light ray around a mass (Schwarzschild
metric) and plots it in 3D. The curvature is computed with
[autodiff](https://github.com/autodiff/autodiff).

## What you need

- `g++` (C++17)
- the `autodiff` and `eigen-3.4.0` folders next to `main.cpp`
- Python with `numpy` and `matplotlib` (`pip install numpy matplotlib`)

## How to download

Clone the repo with :

```
git clone https://github.com/SalvorGoldstein/RG-simulation.git
```
It should run fine.
## How to run

In the project folder:

```
g++ -std=c++17 main.cpp -I autodiff -I eigen-3.4.0 -o main && ./main && python Schwartzschield_aff.py
```

This compiles the program, runs it (it writes `Geodesics.dat` and
`dataforpy.txt`), then plots the result.

If you have Windows, well, switch to Linux. It's the best troubleshooting advice I have.

## Change the simulation

Edit the values in `main()` in `main.cpp`: the mass `M`, the starting radius
`r0`, the initial velocity and the step `dtau`.

If the plot is empty, `dtau` is probably too big.
