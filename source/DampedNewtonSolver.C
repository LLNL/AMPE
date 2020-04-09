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
#include "DampedNewtonSolver.h"

#include <iostream>
#include <cmath>
#include <cassert>

#include <iomanip>

using namespace std;

//=======================================================================

DampedNewtonSolver::DampedNewtonSolver() : NewtonSolver(), d_alpha(1.){};


//=======================================================================
// note: sizes to accomodate up to ternary alloys
//
// c: solution to be updated
void DampedNewtonSolver::UpdateSolution(double* const c,
                                        const double* const fvec,
                                        double** const fjac)
{
   int nn = size();

   static double* mwork[5];
   static double mtmp[25];
   for (int ii = 0; ii < nn; ii++) {
      mwork[ii] = &mtmp[ii * nn];
   }

   const double D = Determinant(fjac);
   assert(fabs(D) > 1.e-15);

   const double D_inv = 1.0 / D;

   // cout<<setprecision(12);
   // cout << "DampedNewtonSolver::UpdateSolution(), N = "<<nn<<", D = " << D <<
   // endl;

   static double del[5];

   // use Cramer's rule to solve linear system
   for (int jj = 0; jj < nn; jj++) {

      CopyMatrix(mwork, fjac);

      // replace jth column with rhs
      for (int ii = 0; ii < nn; ii++) {
         mwork[ii][jj] = fvec[ii];
      }

      const double Dmwork = Determinant(mwork);
      // cout << "nn="<<nn<<", Dmwork="<<Dmwork <<endl;
      del[jj] = D_inv * Dmwork;

      const double maxdel = 0.25;
      if (fabs(del[jj]) > maxdel) del[jj] = del[jj] > 0 ? maxdel : -maxdel;

      // cout << "del[" << jj << "] = " << del[jj] << endl;
   }

   double w = d_alpha;
   for (int ii = 0; ii < nn; ii++) {
      c[ii] = c[ii] - w * del[ii];
   }
}
