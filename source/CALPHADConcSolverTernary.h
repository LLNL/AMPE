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
#ifndef included_CALPHADConcSolverTernary
#define included_CALPHADConcSolverTernary

#include "DampedNewtonSolver.h"

class CALPHADConcentrationSolverTernary :
   public DampedNewtonSolver
{
public :

   CALPHADConcentrationSolverTernary();
      
   virtual ~CALPHADConcentrationSolverTernary() {};
      
   void setup(const double c0,
      const double c1,
      const double hphi,
      const double RTinv,
      const double* const L_AB_L,
      const double* const L_AC_L,
      const double* const L_BC_L,
      const double* const L_AB_S,
      const double* const L_AC_S,
      const double* const L_BC_S,
      const double* const L_ABC_L,
      const double* const L_ABC_S,
      const double* const fA,
      const double* const fB,
      const double* const fC );

   int ComputeConcentration(
      double* const conc,
      const double c0,
      const double c1,
      const double hphi,
      const double RTinv,
      const double* const L_AB_L,
      const double* const L_AC_L,
      const double* const L_BC_L,
      const double* const L_AB_S,
      const double* const L_AC_S,
      const double* const L_BC_S,
      const double* const L_ABC_L,
      const double* const L_ABC_S,
      const double* const fA,
      const double* const fB,
      const double* const fC );

   void RHS(
      const double* const c,
      double* const fvec );

   void Jacobian(
      const double* const c,
      double** const fjac );

private:

   //energies of 3 species, in two phase each
   double d_fA[2];
   double d_fB[2];
   double d_fC[2];

   // L coefficients for phase L
   double d_L_AB_L[4];
   double d_L_AC_L[4];
   double d_L_BC_L[4];
   double d_L_ABC_L[3];
   
   // L coefficients for phase S
   double d_L_AB_S[4];
   double d_L_AC_S[4];
   double d_L_BC_S[4];
   double d_L_ABC_S[3];
 
   double d_c0[2];
   double d_hphi;

   double d_RT;
};

#endif
