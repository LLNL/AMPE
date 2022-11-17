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
#include "AMPE_internal.h"

#include <iostream>
#include <math.h>

bool checkInterval(const double w01, const double w02, const double w12)
{
   if (w01 > 1. || w01 < 0. || w01 > 1. || w01 < 0. || w01 > 1. || w01 < 0.) {
      return false;
   } else {
      return true;
   }
}

int main(int argc, char* argv[])
{
   std::cout << "Test weights computations for 3 pairs." << std::endl;

   double w01;
   double w02;
   double w12;

   {
      double p[3] = {0.2, 0.2, 0.6};
      computeWeights3Pairs(p[0], p[1], p[2], w01, w02, w12);
      std::cout << "w01 = " << w01 << ", w02 = " << w02 << ", w12 = " << w12
                << std::endl;
      if (!checkInterval(w01, w02, w12)) {
         std::cerr << "TEST: failed!!!" << std::endl;
         return 1;
      }
   }

   {
      double p[3] = {0., 0., 1.};
      computeWeights3Pairs(p[0], p[1], p[2], w01, w02, w12);
      std::cout << "w01 = " << w01 << ", w02 = " << w02 << ", w12 = " << w12
                << std::endl;
      if (!checkInterval(w01, w02, w12)) {
         std::cerr << "TEST: failed!!!" << std::endl;
         return 1;
      }
   }

   {
      double p[3] = {0.9998, 1.e-4, -1.e-4};
      computeWeights3Pairs(p[0], p[1], p[2], w01, w02, w12);
      std::cout << "w01 = " << w01 << ", w02 = " << w02 << ", w12 = " << w12
                << std::endl;
      if (!checkInterval(w01, w02, w12)) {
         std::cerr << "TEST: failed!!!" << std::endl;
         return 1;
      }
   }

   {
      double p[3] = {0.98, 0.02, 0.};
      computeWeights3Pairs(p[0], p[1], p[2], w01, w02, w12);
      std::cout << "w01 = " << w01 << ", w02 = " << w02 << ", w12 = " << w12
                << std::endl;
      if (!checkInterval(w01, w02, w12)) {
         std::cerr << "TEST: failed!!!" << std::endl;
         return 1;
      }
   }


   std::cout << "TEST successful!" << std::endl;
   return (0);
}
