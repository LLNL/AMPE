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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include "math_utilities.h"

#include <iostream>
#include <math.h>

#ifdef HAVE_THERMO4PFM
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
   double d = Determinant4(mat);

   const double tol = 1.e-8;
   if (fabs(d + 400.) > tol) {
      std::cerr << "TEST: Determinant of 4x4 matrix failed!!!" << std::endl;
      return 1;
   } else {
      std::cout << "TEST successful!" << std::endl;
   }

   std::cout << "Test function DeterminantN..." << std::endl;
   d = DeterminantN(mat, 4);
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
