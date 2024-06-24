// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#include "TiltingFolchPlapp2005.h"

// Interpolation function for 3 phases
// Folch and Plapp (2005)
// note: p1, p2 order does not matter
// Derivatives are constrained derivatives (p0+p1+p2=1)
double TiltingFolchPlapp2005::g0(const double p0, const double p1,
                                 const double p2)
{
   if (p0 < 0.) return 0.;
   if (p0 > 1.) return 1.;

   return 0.25 * p0 * p0 *
          (15. * (1. - p0) * (1. + p0 - (p2 - p1) * (p2 - p1)) +
           p0 * (9. * p0 * p0 - 5.));
}

// dg/dp0 - 1/3 * sum_j dg/dpj
// note: p1, p2 order does not matter
double TiltingFolchPlapp2005::dg0dp0(const double p0, const double p1,
                                     const double p2)
{
   if (p0 < 0.) return 0.;

   // p0=1 -> p1=p2=0
   if (p0 > 1.) return 0.;

   return 2.5 * p0 *
          ((p2 - p1) * (p2 - p1) * (3. * p0 - 2.) +
           (1. - p0) * (1. - p0) * (3. * p0 + 2.));
}

double TiltingFolchPlapp2005::dg0dp1(const double p0, const double p1,
                                     const double p2)
{
   if (p0 < 0.) return 0.;
   if (p0 > 1.) return 0.;

   if (p1 < 0.) return 0.;
   if (p1 > 1.) return 0.;

   return -0.5 * dg0dp0(p1, p0, p2) + 7.5 * (p1 * p1 * (1. - p1) * (p2 - p0));
}
