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
#include <iostream>
#include <math.h>

#include "Determinant.h"
using namespace Thermo4PFM;

int main(int argc, char* argv[])
{
   std::cout << "Test determinant computation." << std::endl;

   double* mat[4];
   double work[16] = {11, 2, 3, 4, 5, 16, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
   for (short i = 0; i < 4; ++i)
      mat[i] = &work[4 * i];

   std::cout << "Test function Determinant4..." << std::endl;
   double d = evalDeterminant<4>(mat);

   const double tol = 1.e-8;
   if (fabs(d + 400.) > tol) {
      std::cerr << "TEST: Determinant of 4x4 matrix failed!!!" << std::endl;
      return 1;
   } else {
      std::cout << "TEST successful!" << std::endl;
   }

   std::cout << "Test function DeterminantN..." << std::endl;
   d = evalDeterminant<4>(mat);
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
