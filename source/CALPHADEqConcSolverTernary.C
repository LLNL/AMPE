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
#include "CALPHADEqConcSolverTernary.h"
#include "CALPHADFunctions.h"

#include "SAMRAI/tbox/IEEE.h"

#include <iostream>
#include <cmath>
#include <cassert>


//=======================================================================

CALPHADEqConcentrationSolverTernary::CALPHADEqConcentrationSolverTernary()
{
   double def_val = SAMRAI::tbox::IEEE::getSignalingNaN();

   d_fA[0] = def_val;
   d_fA[1] = def_val;
   d_fB[0] = def_val;
   d_fB[1] = def_val;
   d_fC[0] = def_val;
   d_fC[1] = def_val;
}


void CALPHADEqConcentrationSolverTernary::RHS(const double* const c,
                                              double* const fvec)
{
   assert(d_fA[0] == d_fA[0]);
   assert(d_fC[1] == d_fC[1]);

   const double* const cL = &c[0];  // composition of Species A and B in phase L
   const double* const cS = &c[2];  // composition of Species A and B in phase S
   // tbox::pout<<"Compute RHS for CALPHAD..."<<endl;

   double derivIdealMixL[2];
   CALPHADcomputeFIdealMix_derivTernary(d_RT, cL[0], cL[1], derivIdealMixL);

   double derivFMixL[2];
   CALPHADcomputeFMix_derivTernary(d_L_AB_L, d_L_AC_L, d_L_BC_L, d_L_ABC_L,
                                   cL[0], cL[1], derivFMixL);

   double dfLdciL[2];
   // 1st species
   dfLdciL[0] = d_fA[0] - d_fC[0] + derivFMixL[0] + derivIdealMixL[0];

   // 2nd species
   dfLdciL[1] = d_fB[0] - d_fC[0] + derivFMixL[1] + derivIdealMixL[1];

   double derivIdealMixS[2];
   CALPHADcomputeFIdealMix_derivTernary(d_RT, cS[0], cS[1], derivIdealMixS);

   double derivFMixS[2];
   CALPHADcomputeFMix_derivTernary(d_L_AB_S, d_L_AC_S, d_L_BC_S, d_L_ABC_S,
                                   cS[0], cS[1], derivFMixS);

   double dfSdciS[2];
   // 1st species
   dfSdciS[0] = d_fA[1] - d_fC[1] + derivFMixS[0] + derivIdealMixS[0];

   // 2nd species
   dfSdciS[1] = d_fB[1] - d_fC[1] + derivFMixS[1] + derivIdealMixS[1];

   // equation fL-fS-(cL-cS)*dfL/dcL=0
   const double fL = cL[0] * d_fA[0] + cL[1] * d_fB[0] +
                     (1.0 - cL[0] - cL[1]) * d_fC[0] +
                     CALPHADcomputeFIdealMixTernary(d_RT, cL[0], cL[1]) +
                     CALPHADcomputeFMixTernary(d_L_AB_L, d_L_AC_L, d_L_BC_L,
                                               d_L_ABC_L, cL[0], cL[1]);
   const double fS = cS[0] * d_fA[1] + cS[1] * d_fB[1] +
                     (1.0 - cS[0] - cS[1]) * d_fC[1] +
                     CALPHADcomputeFIdealMixTernary(d_RT, cS[0], cS[1]) +
                     CALPHADcomputeFMixTernary(d_L_AB_S, d_L_AC_S, d_L_BC_S,
                                               d_L_ABC_S, cS[0], cS[1]);
   // std::cout<<"fL="<<fL<<", fS="<<fS<<endl;
   fvec[0] =
       fL - fS - (cL[0] - cS[0]) * dfLdciL[0] - (cL[1] - cS[1]) * dfLdciL[1];

   // equation: slope 0 in tangent plane for std::vector orthogonal to cL-CS
   fvec[1] = dfLdciL[0] * cL[1] - dfLdciL[0] * cS[1] - dfLdciL[1] * cL[0] +
             dfLdciL[1] * cS[0];

   fvec[2] = dfLdciL[0] - dfSdciS[0];
   fvec[3] = dfLdciL[1] - dfSdciS[1];

   // std::cout<<"fvec="<<fvec[0]<<","<<fvec[1]<<","<<fvec[2]<<","<<fvec[3]<<endl;
}

//=======================================================================

void CALPHADEqConcentrationSolverTernary::Jacobian(const double* const c,
                                                   double** const fjac)
{
   // tbox::pout<<"Compute Jacobian for CALPHAD..."<<endl;
   const double* const cL = &c[0];
   const double* const cS = &c[2];
   // tbox::pout<<"Compute RHS for CALPHAD..."<<endl;

   double derivIdealMixL[2];
   CALPHADcomputeFIdealMix_derivTernary(d_RT, cL[0], cL[1], derivIdealMixL);

   double derivFMixL[2];
   CALPHADcomputeFMix_derivTernary(d_L_AB_L, d_L_AC_L, d_L_BC_L, d_L_ABC_L,
                                   cL[0], cL[1], derivFMixL);

   double deriv2IdealMixL[4];
   CALPHADcomputeFIdealMix_deriv2Ternary(d_RT, cL[0], cL[1], deriv2IdealMixL);

   double deriv2FMixL[4];
   CALPHADcomputeFMix_deriv2Ternary(d_L_AB_L, d_L_AC_L, d_L_BC_L, d_L_ABC_L,
                                    cL[0], cL[1], deriv2FMixL);

   double dfLdciL[2];
   // 1st species
   dfLdciL[0] = d_fA[0] - d_fC[0] + derivFMixL[0] + derivIdealMixL[0];

   // 2nd species
   dfLdciL[1] = d_fB[0] - d_fC[0] + derivFMixL[1] + derivIdealMixL[1];

   double d2fLdciL2[3];  // include only one cross term (other one equal by
                         // symmetry)
   d2fLdciL2[0] = deriv2FMixL[0] + deriv2IdealMixL[0];
   d2fLdciL2[1] = deriv2FMixL[1] + deriv2IdealMixL[1];
   d2fLdciL2[2] = deriv2FMixL[3] + deriv2IdealMixL[3];

   double derivIdealMixS[2];
   CALPHADcomputeFIdealMix_derivTernary(d_RT, cS[0], cS[1], derivIdealMixS);

   double derivFMixS[2];
   CALPHADcomputeFMix_derivTernary(d_L_AB_S, d_L_AC_S, d_L_BC_S, d_L_ABC_S,
                                   cS[0], cS[1], derivFMixS);

   double deriv2IdealMixS[4];
   CALPHADcomputeFIdealMix_deriv2Ternary(d_RT, cS[0], cS[1], deriv2IdealMixS);

   double deriv2FMixS[4];
   CALPHADcomputeFMix_deriv2Ternary(d_L_AB_S, d_L_AC_S, d_L_BC_S, d_L_ABC_S,
                                    cS[0], cS[1], deriv2FMixS);

   double d2fSdciS2[3];
   d2fSdciS2[0] = deriv2FMixS[0] + deriv2IdealMixS[0];
   d2fSdciS2[1] = deriv2FMixS[1] + deriv2IdealMixS[1];
   d2fSdciS2[2] = deriv2FMixS[3] + deriv2IdealMixS[3];

   // f[i][j]=df[i]/dc[j]
   fjac[0][0] = d_fA[0] - d_fC[0] + derivIdealMixL[0] + derivFMixL[0] -
                dfLdciL[0] - (cL[0] - cS[0]) * d2fLdciL2[0] -
                (cL[1] - cS[1]) * d2fLdciL2[1];

   fjac[0][1] = d_fB[0] - d_fC[0] + derivIdealMixL[1] + derivFMixL[1] -
                dfLdciL[1] - (cL[0] - cS[0]) * d2fLdciL2[1] -
                (cL[1] - cS[1]) * d2fLdciL2[2];

   fjac[0][2] =
       -d_fA[1] + d_fC[1] - derivIdealMixS[0] - derivFMixS[0] + dfLdciL[0];

   fjac[0][3] =
       -d_fB[1] + d_fC[1] - derivIdealMixS[1] - derivFMixS[1] + dfLdciL[1];

   fjac[1][0] =
       d2fLdciL2[0] * (cL[1] - cS[1]) - dfLdciL[1] + d2fLdciL2[1] * cS[0];

   fjac[1][1] =
       dfLdciL[0] - d2fLdciL2[1] * cS[1] + d2fLdciL2[1] * (cL[0] - cS[0]);

   fjac[1][2] = dfLdciL[1];
   fjac[1][3] = -dfLdciL[0];


   fjac[2][0] = d2fLdciL2[0];
   fjac[2][1] = d2fLdciL2[1];
   fjac[2][2] = -d2fSdciS2[0];
   fjac[2][3] = -d2fSdciS2[1];

   fjac[3][0] = d2fLdciL2[1];
   fjac[3][1] = d2fLdciL2[2];
   fjac[3][2] = -d2fSdciS2[1];
   fjac[3][3] = -d2fSdciS2[2];

   // std::cout<<"Jacobian:"<<endl;
   // std::cout<<"("<<fjac[0][0]<<","<<fjac[0][1]<<","<<fjac[0][2]<<","<<fjac[0][3]<<")"<<endl;
   // std::cout<<"("<<fjac[1][0]<<","<<fjac[1][1]<<","<<fjac[1][2]<<","<<fjac[1][3]<<")"<<endl;
   // std::cout<<"("<<fjac[2][0]<<","<<fjac[2][1]<<","<<fjac[2][2]<<","<<fjac[2][3]<<")"<<endl;
   // std::cout<<"("<<fjac[3][0]<<","<<fjac[3][1]<<","<<fjac[3][2]<<","<<fjac[3][3]<<")"<<endl;
}

//=======================================================================

int CALPHADEqConcentrationSolverTernary::ComputeConcentration(
    double* const conc, const double RTinv, const double* const L_AB_L,
    const double* const L_AC_L, const double* const L_BC_L,
    const double* const L_AB_S, const double* const L_AC_S,
    const double* const L_BC_S, const double* const L_ABC_L,
    const double* const L_ABC_S, const double* const fA, const double* const fB,
    const double* const fC)
{
   d_RTinv = RTinv;
   d_RT = 1. / RTinv;

   for (int ii = 0; ii < 4; ii++) {
      d_L_AB_L[ii] = L_AB_L[ii];
      d_L_AC_L[ii] = L_AC_L[ii];
      d_L_BC_L[ii] = L_BC_L[ii];
   }
   for (int ii = 0; ii < 4; ii++) {
      d_L_AB_S[ii] = L_AB_S[ii];
      d_L_AC_S[ii] = L_AC_S[ii];
      d_L_BC_S[ii] = L_BC_S[ii];
   }
   for (int ii = 0; ii < 3; ii++) {
      d_L_ABC_L[ii] = L_ABC_L[ii];
      d_L_ABC_S[ii] = L_ABC_S[ii];
   }

   // loop over phases (L and S)
   for (int ii = 0; ii < 2; ii++) {
      d_fA[ii] = fA[ii];
      d_fB[ii] = fB[ii];
      d_fC[ii] = fC[ii];
   }

   return DampedNewtonSolver::ComputeSolution(conc, 4);
}
