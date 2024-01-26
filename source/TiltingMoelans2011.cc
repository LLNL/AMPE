// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#include "TiltingMoelans2011.h"

// Interpolation function for 3 phases
// Adapted from multi-phases from Moelans (2011)
// note: p1, p2 order does not matter
// Derivatives are constrained derivatives (p0+p1+p2=1)
double TiltingMoelans2011::g0(const double p0, const double p1, const double p2)
{
   if (p0 < 0.) return 0.;

   return p0 * p0 / (p0 * p0 + p1 * p1 + p2 * p2);
}

// dg/dp0 - 1/3 * sum_j dg/dpj
// note: p1, p2 order does not matter
double TiltingMoelans2011::dg0dp0(const double p0, const double p1,
                                  const double p2)
{
   if (p0 < 0.) return 0.;

   const double a = p0 * p0 + p1 * p1 + p2 * p2;
   return 2. * p0 * (2. * p1 * p1 + 2. * p2 * p2 + p0 * p1 + p0 * p2) /
          (3. * a * a);
}

double TiltingMoelans2011::dg0dp1(const double p0, const double p1,
                                  const double p2)
{
   if (p0 < 0.) return 0.;

   const double a = p0 * p0 + p1 * p1 + p2 * p2;
   return 2. * p0 * (p0 * (p2 - 2. * p1) - (p1 * p1 + p2 * p2)) / (3. * a * a);
}
