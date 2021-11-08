// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#ifndef HAVE_THERMO4PFM
#include "math_utilities.h"
#endif

#include <iostream>
#include <math.h>

#ifdef HAVE_THERMO4PFM
#include "Determinant.h"
using namespace Thermo4PFM;
#else
using namespace ampe_thermo;
#endif

int main(int argc, char* argv[])
{
   std::cout << "Test determinant computation." << std::endl;

   double* mat[4];
   double work[16] = {11, 2, 3, 4, 5, 16, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
   for (short i = 0; i < 4; ++i)
      mat[i] = &work[4 * i];

   std::cout << "Test function Determinant4..." << std::endl;
#ifdef HAVE_THERMO4PFM
   double d = evalDeterminant<4>(mat);
#else
   double d = Determinant4(mat);
#endif

   const double tol = 1.e-8;
   if (fabs(d + 400.) > tol) {
      std::cerr << "TEST: Determinant of 4x4 matrix failed!!!" << std::endl;
      return 1;
   } else {
      std::cout << "TEST successful!" << std::endl;
   }

   std::cout << "Test function DeterminantN..." << std::endl;
#ifdef HAVE_THERMO4PFM
   d = evalDeterminant<4>(mat);
#else
   d = DeterminantN(mat, 4);
#endif
   if (fabs(d + 400.) > tol) {
      std::cerr << "TEST: Determinant of 4x4 matrix failed!!!" << std::endl;
      std::cerr << "computed d=" << d << std::endl;
      std::cerr << "exact d=" << 400. << std::endl;
      return 1;
   } else {
      std::cout << "TEST successful!" << std::endl;
   }

   return (0);
}
