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
#include "CALPHADFunctions.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Thermo4PFM;

int main(int argc, char *argv[])
{
   std::cout << "Test CALPHAD functions." << std::endl;

   const double epsilon = 1.e-8;
   const double tol = 1.e-6;
   int nfailures = 0;

   std::cout << std::setprecision(12);
   std::cerr << std::setprecision(12);

   {
      double l0 = 2.3;
      double l1 = 5.1;
      double l2 = 3.2;
      double l3 = -2.5;
      double c = 0.33;

      double f0 = CALPHADcomputeFMixBinary(l0, l1, l2, l3, c);
      double f1 = CALPHADcomputeFMixBinary(l0, l1, l2, l3, c + epsilon);

      double fd = (f1 - f0) / epsilon;
      std::cout << "Numerical derivative   = " << fd << std::endl;
      double ad = CALPHADcomputeFMix_derivBinary(l0, l1, l2, l3, c);
      std::cout << "Analytical derivative = " << ad << std::endl;

      if (fabs(fd - ad) > tol) {
         nfailures++;
         std::cerr << "TEST: CALPHADcomputeFMix_derivBinary failed!!!"
                   << std::endl;
         std::cerr << "Difference = " << fd - ad << std::endl;
      } else {
         std::cout << "TEST successful!" << std::endl;
      }
   }

   {
      double rt = 10.;
      double cA = 0.33;
      double cB = 0.21;
      double f0 = CALPHADcomputeFIdealMixTernary(rt, cA, cB);
      double f1 = CALPHADcomputeFIdealMixTernary(rt, cA + epsilon, cB);
      double f2 = CALPHADcomputeFIdealMixTernary(rt, cA, cB + epsilon);
      double deriv[2];
      CALPHADcomputeFIdealMix_derivTernary(rt, cA, cB, deriv);
      std::cout << "Analytical derivative = " << deriv[0] << "," << deriv[1]
                << std::endl;

      double fd1 = (f1 - f0) / epsilon;
      std::cout << "Numerical derivative   = " << fd1 << std::endl;
      if (fabs(fd1 - deriv[0]) > tol) {
         nfailures++;
         std::cerr << "TEST: CALPHADcomputeFIdealMix_derivTernary failed!!!"
                   << std::endl;
         std::cerr << "Difference = " << fd1 - deriv[0] << std::endl;
      } else {
         std::cout << "TEST CALPHADcomputeFIdealMix_derivTernary, c0, "
                      "successful!"
                   << std::endl;
      }

      double fd2 = (f2 - f0) / epsilon;
      std::cout << "Numerical derivative   = " << fd2 << std::endl;
      if (fabs(fd2 - deriv[1]) > tol) {
         nfailures++;
         std::cerr << "TEST: CALPHADcomputeFIdealMix_derivTernary failed!!!"
                   << std::endl;
         std::cerr << "Difference = " << fd2 - deriv[1] << std::endl;
      } else {
         std::cout << "TEST CALPHADcomputeFIdealMix_derivTernary, c1, "
                      "successful!"
                   << std::endl;
      }
   }

   {
      double lAB[4] = {1.2, 3.5, 6.1, -1.2};
      double lAC[4] = {3.2, 4.5, -3.1, 7.1};
      double lBC[4] = {0.2, 0.7, 4.1, -8.2};
      double lABC[3] = {3.7, 6.3, -1.1};

      double cA = 0.33;
      double cB = 0.21;

      double f0 = CALPHADcomputeFMixTernary(lAB, lAC, lBC, lABC, cA, cB);
      double f1 =
          CALPHADcomputeFMixTernary(lAB, lAC, lBC, lABC, cA + epsilon, cB);
      double f2 =
          CALPHADcomputeFMixTernary(lAB, lAC, lBC, lABC, cA, cB + epsilon);

      double deriv[2];
      CALPHADcomputeFMix_derivTernary(lAB, lAC, lBC, lABC, cA, cB, deriv);

      double fd1 = (f1 - f0) / epsilon;
      if (fabs(fd1 - deriv[0]) > tol) {
         nfailures++;
         std::cerr << "TEST: CALPHADcomputeFMix_derivTernary failed!!!"
                   << std::endl;
         std::cerr << "Difference = " << fd1 - deriv[0] << std::endl;
      } else {
         std::cout << "TEST CALPHADcomputeFMix_derivTernary, c0, successful!"
                   << std::endl;
      }
      double fd2 = (f2 - f0) / epsilon;
      if (fabs(fd2 - deriv[1]) > tol) {
         nfailures++;
         std::cerr << "TEST: CALPHADcomputeFMix_derivTernary failed!!!"
                   << std::endl;
         std::cerr << "Difference = " << fd2 - deriv[1] << std::endl;
      } else {
         std::cout << "TEST CALPHADcomputeFMix_derivTernary, c1, successful!"
                   << std::endl;
      }

      double deriv2[4];
      CALPHADcomputeFMix_deriv2Ternary(lAB, lAC, lBC, lABC, cA, cB, deriv2);

      double fderiv0[2];
      CALPHADcomputeFMix_derivTernary(lAB, lAC, lBC, lABC, cA + epsilon, cB,
                                      fderiv0);
      double fderiv1[2];
      CALPHADcomputeFMix_derivTernary(lAB, lAC, lBC, lABC, cA, cB + epsilon,
                                      fderiv1);

      double fd = (fderiv0[0] - deriv[0]) / epsilon;
      std::cout << "FD=" << fd << ", exact=" << deriv2[0] << std::endl;
      if (fabs(fd - deriv2[0]) > tol) {
         nfailures++;
         std::cerr << "TEST: CALPHADcomputeFMix_deriv2Ternary failed!!!"
                   << std::endl;
         std::cerr << "Difference = " << fd - deriv2[0] << std::endl;
      } else {
         std::cout << "TEST CALPHADcomputeFMix_deriv2Ternary, c0, c0, "
                      "successful!"
                   << std::endl;
      }

      fd = (fderiv1[1] - deriv[1]) / epsilon;
      std::cout << "FD=" << fd << ", exact=" << deriv2[3] << std::endl;
      if (fabs(fd - deriv2[3]) > tol) {
         nfailures++;
         std::cerr << "TEST: CALPHADcomputeFMix_deriv2Ternary failed!!!"
                   << std::endl;
         std::cerr << "Difference = " << fd - deriv2[3] << std::endl;
      } else {
         std::cout << "TEST CALPHADcomputeFMix_deriv2Ternary, c1, c1, "
                      "successful!"
                   << std::endl;
      }

      // cross derivatives
      fd = (fderiv0[1] - deriv[1]) / epsilon;
      std::cout << "FD=" << fd << ", exact=" << deriv2[1] << std::endl;
      if (fabs(fd - deriv2[1]) > tol) {
         nfailures++;
         std::cerr << "TEST: CALPHADcomputeFMix_deriv2Ternary failed!!!"
                   << std::endl;
         std::cerr << "Difference = " << fd - deriv2[0] << std::endl;
      } else {
         std::cout << "TEST CALPHADcomputeFMix_deriv2Ternary, c0, c1, "
                      "successful!"
                   << std::endl;
      }

      fd = (fderiv1[0] - deriv[0]) / epsilon;
      std::cout << "FD=" << fd << ", exact=" << deriv2[2] << std::endl;
      if (fabs(fd - deriv2[2]) > tol) {
         nfailures++;
         std::cerr << "TEST: CALPHADcomputeFMix_deriv2Ternary failed!!!"
                   << std::endl;
         std::cerr << "Difference = " << fd - deriv2[0] << std::endl;
      } else {
         std::cout << "TEST CALPHADcomputeFMix_deriv2Ternary, c1, c0, "
                      "successful!"
                   << std::endl;
      }
   }

   if (nfailures == 0) std::cout << "TEST PASSED" << std::endl;
   return nfailures;
}
