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
#include "CALPHADEqConcSolverBinaryWithPenalty.h"
#include "CALPHADFunctions.h"

#include <iostream>
#include <cmath>
#include <cassert>
#include "SAMRAI/tbox/PIO.h"


using namespace SAMRAI;

//=======================================================================

void CALPHADEqConcentrationSolverBinaryWithPenalty::RHS(const double* const c,
                                                        double* const fvec)
{
   assert(d_penalty_parametersL.size() == 6);
   assert(d_penalty_parametersS.size() == 6);

   CALPHADEqConcentrationSolverBinary::RHS(c, fvec);

   double dfdci[2];
   dfdci[0] = CALPHADcomputeDerivPenalty(d_penalty_parametersL[0],
                                         d_penalty_parametersL[1],
                                         d_penalty_parametersL[2],
                                         d_penalty_parametersL[3],
                                         d_penalty_parametersL[4],
                                         d_penalty_parametersL[5], c[0]);
   dfdci[1] = CALPHADcomputeDerivPenalty(d_penalty_parametersS[0],
                                         d_penalty_parametersS[1],
                                         d_penalty_parametersS[2],
                                         d_penalty_parametersS[3],
                                         d_penalty_parametersS[4],
                                         d_penalty_parametersS[5], c[1]);

   fvec[0] += (CALPHADcomputePenalty(
                   d_penalty_parametersL[0], d_penalty_parametersL[1],
                   d_penalty_parametersL[2], d_penalty_parametersL[3],
                   d_penalty_parametersL[4], d_penalty_parametersL[5], c[0]) -
               CALPHADcomputePenalty(
                   d_penalty_parametersS[0], d_penalty_parametersS[1],
                   d_penalty_parametersS[2], d_penalty_parametersS[3],
                   d_penalty_parametersS[4], d_penalty_parametersS[5], c[1]) -
               (c[0] - c[1]) * dfdci[1]);

   fvec[1] += (dfdci[0] - dfdci[1]);
   // tbox::pout<<"Compute RHS for CALPHAD with Penalty...
   // fvec[0]="<<fvec[0]<<", fvec[1]="<<fvec[1]<<endl;
}

//=======================================================================

void CALPHADEqConcentrationSolverBinaryWithPenalty::Jacobian(
    const double* const c, double** const fjac)
{
   // tbox::pout<<"Compute Jacobian for CALPHAD with Penalty..."<<endl;
   CALPHADEqConcentrationSolverBinary::Jacobian(c, fjac);

   double dfdci[2];

   dfdci[0] = CALPHADcomputeDerivPenalty(d_penalty_parametersL[0],
                                         d_penalty_parametersL[1],
                                         d_penalty_parametersL[2],
                                         d_penalty_parametersL[3],
                                         d_penalty_parametersL[4],
                                         d_penalty_parametersL[5], c[0]);
   dfdci[1] = CALPHADcomputeDerivPenalty(d_penalty_parametersS[0],
                                         d_penalty_parametersS[1],
                                         d_penalty_parametersS[2],
                                         d_penalty_parametersS[3],
                                         d_penalty_parametersS[4],
                                         d_penalty_parametersS[5], c[1]);


   fjac[0][0] += (dfdci[0] - dfdci[1]);

   const double d2f1dc1 = CALPHADcompute2ndDerivPenalty(
       d_penalty_parametersS[0], d_penalty_parametersS[1],
       d_penalty_parametersS[2], d_penalty_parametersS[3],
       d_penalty_parametersS[4], d_penalty_parametersS[5], c[1]);
   const double d2f0dc0 = CALPHADcompute2ndDerivPenalty(
       d_penalty_parametersL[0], d_penalty_parametersL[1],
       d_penalty_parametersL[2], d_penalty_parametersL[3],
       d_penalty_parametersL[4], d_penalty_parametersL[5], c[0]);

   fjac[0][1] -= (c[0] - c[1]) * d2f1dc1;


   fjac[1][0] += d2f0dc0;
   fjac[1][1] -= d2f1dc1;
   // tbox::pout<<"Compute J for CALPHAD with Penalty...
   // fjac[0][0]="<<fjac[0][0]
   //                                              <<", fjac[0][1]="<<fjac[0][1]
   //                                              <<", fjac[1][0]="<<fjac[1][0]
   //                                              <<",
   //                                              fjac[1][1]="<<fjac[1][1]<<endl;
}

//=======================================================================

int CALPHADEqConcentrationSolverBinaryWithPenalty::
    ComputeConcentrationWithPenalty(
        double* const conc, const double RTinv, const double* const L0,
        const double* const L1, const double* const L2, const double* const L3,
        const double* const fA, const double* const fB,
        std::vector<std::vector<double> >& penalty_parameters)
{
   d_penalty_parametersL = penalty_parameters[0];
   d_penalty_parametersS = penalty_parameters[1];

   CALPHADEqConcentrationSolverBinary::ComputeConcentration(conc, RTinv, L0, L1,
                                                            L2, L3, fA, fB);

   return DampedNewtonSolver::ComputeSolution(conc, 2);
}
