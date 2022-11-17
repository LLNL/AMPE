// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include <math.h>
#include <cassert>
#include <iostream>

// Given 3 phases fractions p0, p1, p2, compute weights associated with
// each pairs satisfying w01+w02+w12=1
void computeWeights3Pairs(const double p0, const double p1, const double p2,
                          double& w01, double& w02, double& w12)
{
   double pp[3] = {p0, p1, p2};

   // step 1: enforce max_i(pp[i])<=maxp
   const double maxp = 0.99;
   for (short i = 0; i < 3; i++) {
      if (pp[i] > maxp) {
         pp[i] = maxp;
         double delta = 1. - pp[0] - pp[1] - pp[2];
         pp[(i + 1) % 3] += 0.5 * delta;
         pp[(i + 2) % 3] += 0.5 * delta;
         break;
      }
   }
   assert(pp[0] + pp[1] > 0.);
   assert(pp[0] + pp[2] > 0.);
   assert(pp[1] + pp[2] > 0.);

   // set coordinate vectors for pairs of phases
   double v01[3] = {pp[0] / (pp[0] + pp[1]), pp[1] / (pp[0] + pp[1]), 0.};
   double v02[3] = {pp[0] / (pp[0] + pp[2]), 0., pp[2] / (pp[0] + pp[2])};
   double v12[3] = {0., pp[1] / (pp[1] + pp[2]), pp[2] / (pp[1] + pp[2])};

   // d0 is the norm of v02-v01
   double d0 = sqrt((v02[0] - v01[0]) * (v02[0] - v01[0]) +
                    (v02[1] - v01[1]) * (v02[1] - v01[1]) +
                    (v02[2] - v01[2]) * (v02[2] - v01[2]));

   // d1 is the norm of v12-v01
   double d1 = sqrt((v12[0] - v01[0]) * (v12[0] - v01[0]) +
                    (v12[1] - v01[1]) * (v12[1] - v01[1]) +
                    (v12[2] - v01[2]) * (v12[2] - v01[2]));

   // d2 is the norm of v12-v02
   double d2 = sqrt((v12[0] - v02[0]) * (v12[0] - v02[0]) +
                    (v12[1] - v02[1]) * (v12[1] - v02[1]) +
                    (v12[2] - v02[2]) * (v12[2] - v02[2]));

   assert(d0 > 1.e-32);
   assert(d1 > 1.e-32);
   assert(d2 > 1.e-32);

   // compute weights as angles
   double piinv = 1. / M_PI;
   double a01 = (d0 * d0 + d1 * d1 - d2 * d2) / (2. * d0 * d1);
   if (a01 < -1.) a01 = -1.;
   if (a01 > 1.) a01 = 1.;
   w01 = piinv * acos(a01);

   double a02 = (d0 * d0 + d2 * d2 - d1 * d1) / (2. * d0 * d2);
   if (a02 < -1.) a02 = -1.;
   if (a02 > 1.) a02 = 1.;
   w02 = piinv * acos(a02);

   double a12 = (d2 * d2 + d1 * d1 - d0 * d0) / (2. * d1 * d2);
   if (a12 < -1.) a12 = -1.;
   if (a12 > 1.) a12 = 1.;
   w12 = piinv * acos(a12);

   // if(w01!=w01 || w12!=w12 || w02!=w02)
   //{
   // std::cerr<<"a01="<<a01<<std::endl;
   // std::cerr<<"a02="<<a02<<std::endl;
   // std::cerr<<"a12="<<a12<<std::endl;
   // std::cerr<<"p0="<<p0<<", p1="<<p1<<", p2="<<p2<<std::endl;
   // std::cerr<<"d0="<<d0<<", d1="<<d1<<", d2="<<d2<<std::endl;
   // std::cerr<<"p0="<<pp[0]<<", p1="<<pp[1]<<", p2="<<pp[2]<<std::endl;
   //}
   assert(w01 == w01);
   assert(w02 == w02);
   assert(w12 == w12);
}
