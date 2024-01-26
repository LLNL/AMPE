// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#ifndef included_TiltingMoelans2011
#define included_TiltingMoelans2011

class TiltingMoelans2011
{
 public:
   static double g0(const double p0, const double p1, const double p2);
   static double dg0dp0(const double p0, const double p1, const double p2);
   static double dg0dp1(const double p0, const double p1, const double p2);
};

#endif
