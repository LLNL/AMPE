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
#include "CALPHADTieLineConcSolverTernary.h"
#include "CALPHADConcSolverTernary.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Thermo4PFM;
typedef float jacobian_dtype;

int main(int argc, char* argv[])
{
   std::cout << "=======================================" << std::endl;
   std::cout << "Test CALPHADTieLineConcSolverTernary..." << std::endl;

   int nfailures = 0;

   double epsilon = 1.e-8;
   double tol = 1.e-6;

   double cA = 0.1;
   double cB = 0.2;
   CALPHADTieLineConcSolverTernary solver;

   // energies of 3 species, in two phase each
   double fA[2] = {2.3, 4.5};
   double fB[2] = {0.5, 2.6};
   double fC[2] = {0.9, 3.1};

   // L coefficients for 2 possible phases (L and S)
   double L_AB_L[4] = {2., 3., 4., 5.};
   double L_AC_L[4] = {8., 5., 1., 3.};
   double L_BC_L[4] = {6., 4., 9., 2.};

   double L_AB_S[4] = {2., 2., 3., 4.};
   double L_AC_S[4] = {6., 8., 5., 1.};
   double L_BC_S[4] = {2., 6., 4., 9.};

   double L_ABC_L[3] = {3.1, 4.1, 5.1};
   double L_ABC_S[3] = {2.1, 3.1, 4.1};

   double RTinv = 10.;

   solver.setup(cA, cB, RTinv, L_AB_L, L_AC_L, L_BC_L, L_AB_S, L_AC_S, L_BC_S,
                L_ABC_L, L_ABC_S, fA, fB, fC);

   double fvec1[5];
   double fvec2[5];
   jacobian_dtype* fjac[5];
   for (int i = 0; i < 5; i++)
      fjac[i] = new jacobian_dtype[5];

   double x[5] = {0.1, 0.2, 0.3, 0.4, 0.5};

   solver.RHS(x, fvec1);

   solver.Jacobian(x, fjac);

   std::cerr << std::setprecision(12);

   // loop over variables
   for (int j = 0; j < 5; j++) {
      std::cout << "----------------------------" << std::endl;
      std::cout << "Test variations of variable " << j << std::endl;
      x[j] += epsilon;
      solver.RHS(x, fvec2);

      double fd[5];
      // loop over equations
      for (int i = 0; i < 5; i++)
         fd[i] = (fvec2[i] - fvec1[i]) / epsilon;

      for (int i = 0; i < 5; i++) {
         if (fabs(fd[i] - fjac[i][j]) > tol) {
            nfailures++;
            std::cerr << "ERROR: Equation " << i << ", FD=" << fd[i]
                      << ", fjac=" << fjac[i][j] << std::endl;
         }
      }

      x[j] -= epsilon;
   }

   std::cout << "=========================================" << std::endl;
   std::cout << "Test CALPHADConcentrationSolverTernary..." << std::endl;
   CALPHADConcSolverTernary solver2;
   solver2.setup(cA, cB, 0.5, RTinv, L_AB_L, L_AC_L, L_BC_L, L_AB_S, L_AC_S,
                 L_BC_S, L_ABC_L, L_ABC_S, fA, fB, fC);

   solver2.RHS(x, fvec1);

   solver2.Jacobian(x, fjac);

   std::cerr << std::setprecision(12);

   // loop over variables
   for (int j = 0; j < 4; j++) {
      std::cout << "----------------------------" << std::endl;
      std::cout << "Test variations of variable " << j << std::endl;
      x[j] += epsilon;
      solver2.RHS(x, fvec2);

      double fd[4];
      // loop over equations
      for (int i = 0; i < 4; i++)
         fd[i] = (fvec2[i] - fvec1[i]) / epsilon;

      for (int i = 0; i < 4; i++) {
         if (fabs(fd[i] - fjac[i][j]) > tol) {
            nfailures++;
            std::cerr << "ERROR: Equation " << i << ", FD=" << fd[i]
                      << ", fjac=" << fjac[i][j] << std::endl;
         }
      }

      x[j] -= epsilon;
   }

   if (nfailures == 0) std::cout << "TEST PASSED" << std::endl;

   return nfailures;
}
