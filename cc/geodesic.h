#pragma once
#include "types.h"

dual partial(const func& fn,
             dual t, dual r, dual theta, dual phi,
             int coordinate);

func christoffel(int i, int k, int l,
                 const met& g, const met& invg);

void euler(rhs_fn sd,
           double (&q)[8], double t, double dt, int n,
           const christ& Gam);

double timelikeU(double t, double r, double theta, double phi,
                 double ur, double utheta, double uphi,
                 const met& g);

double lightlikeU(double t, double r, double theta, double phi,
                  double ur, double utheta, double uphi,
                  const met& g);

void geodesics(double r, double theta0, double phi0,
               double ur, double utheta, double uphi,
               double r_max, int steps, double dtau,
               const std::string& filename,
               const met& g, const christ& Gam,
               rhs_fn sd, double r_min,
               bool type = true);