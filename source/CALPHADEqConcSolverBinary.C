// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
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
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
#include "CALPHADEqConcSolverBinary.h"
#include "CALPHADFunctions.h"

#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

//=======================================================================

void CALPHADEqConcentrationSolverBinary::RHS(
   const double* const c, //composition of species A in various phases 
   double* const fvec )
{
   //tbox::pout<<"Compute RHS for CALPHAD..."<<endl;
   double dfdci[2];
   
   //loop over phases
   for ( int ii = 0; ii < 2; ii++ ) {

      dfdci[ii] = d_fA[ii] - d_fB[ii]
       + CALPHADcomputeFMix_derivBinary( d_L0[ii], d_L1[ii], d_L2[ii], d_L3[ii], c[ii] )
       + CALPHADcomputeFIdealMix_derivBinary( d_RT, c[ii] );

   }
      
   fvec[0] =
        c[0] * d_fA[0] + ( 1.0 - c[0] ) * d_fB[0]
      + CALPHADcomputeFIdealMixBinary( d_RT, c[0] )
      + CALPHADcomputeFMixBinary( d_L0[0], d_L1[0], d_L2[0], d_L3[0], c[0] )
      - c[1] * d_fA[1] - ( 1.0 - c[1] ) * d_fB[1]
      - CALPHADcomputeFIdealMixBinary( d_RT, c[1] )
      - CALPHADcomputeFMixBinary( d_L0[1], d_L1[1], d_L2[1], d_L3[1], c[1] )
      - (c[0]-c[1])* ( dfdci[1] );

   // dfL/dcL - dfS/dcS
   fvec[1] = dfdci[0] - dfdci[1];
}

//=======================================================================

void CALPHADEqConcentrationSolverBinary::Jacobian(
   const double* const c,
   double** const fjac )
{
   //tbox::pout<<"Compute Jacobian for CALPHAD..."<<endl;
   double dfdci[2];

   //loop over phases
   for ( int ii = 0; ii < 2; ii++ ) {

      dfdci[ii] = d_fA[ii] - d_fB[ii]
       + CALPHADcomputeFMix_derivBinary( d_L0[ii], d_L1[ii], d_L2[ii], d_L3[ii], c[ii] )
       + CALPHADcomputeFIdealMix_derivBinary( d_RT, c[ii] );

   }

   // f[i][j]=df[i]/dc[j]
   fjac[0][0] = dfdci[0] - dfdci[1];
   fjac[0][1] = - (c[0]-c[1])*
                     ( CALPHADcomputeFIdealMix_deriv2Binary( d_RT, c[1] )
                     + CALPHADcomputeFMix_deriv2Binary( d_L0[1], d_L1[1], d_L2[1], d_L3[1], c[1] ) );

   fjac[1][0] = CALPHADcomputeFMix_deriv2Binary( d_L0[0], d_L1[0], d_L2[0], d_L3[0], c[0] )
              + CALPHADcomputeFIdealMix_deriv2Binary( d_RT, c[0] );
         
   fjac[1][1] = - CALPHADcomputeFMix_deriv2Binary( d_L0[1], d_L1[1], d_L2[1], d_L3[1], c[1] )
                - CALPHADcomputeFIdealMix_deriv2Binary( d_RT, c[1] );
}

//=======================================================================

int CALPHADEqConcentrationSolverBinary::ComputeConcentration(
   double* const conc,
   const double RTinv,
   const double* const L0,
   const double* const L1,
   const double* const L2,
   const double* const L3,
   const double* const fA,
   const double* const fB )
{
   d_RTinv = RTinv;
   d_RT    = 1./RTinv;

   for ( int ii = 0; ii < 2; ii++ ) {
      d_L0[ii] = L0[ii];
      d_L1[ii] = L1[ii];
      d_L2[ii] = L2[ii];
      d_L3[ii] = L3[ii];
      d_fA[ii] = fA[ii];
      d_fB[ii] = fB[ii];
   }

   return DampedNewtonSolver::ComputeSolution( conc, 2 );
}
