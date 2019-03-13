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
#include "KKSdiluteBinaryConcentrationSolver.h"
#include "xlogx.h"

#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

//=======================================================================

KKSdiluteBinaryConcentrationSolver::KKSdiluteBinaryConcentrationSolver()
{
   d_N = 2;
}

//=======================================================================

// solve for c=(c_L, c_A)
void KKSdiluteBinaryConcentrationSolver::RHS(
   const double* const c,
   double* const fvec )
{
   fvec[0] =
      -d_c0 + ( 1.0 - d_hphi ) * c[0] + 
      d_hphi * c[1];
   fvec[1] = xlogx_deriv(c[0])-xlogx_deriv(1.-c[0])
            -xlogx_deriv(c[1])+xlogx_deriv(1.-c[1])
            - (d_fA-d_fB);
   cout<< "d_fA="<<d_fA<<", d_fB="<<d_fB<<endl;
}

//=======================================================================

void KKSdiluteBinaryConcentrationSolver::Jacobian(
   const double* const c,
   double** const fjac )
{
   fjac[0][0] = ( 1.0 - d_hphi );
   fjac[0][1] = d_hphi;

   fjac[1][0] =  xlogx_deriv2(c[0])+xlogx_deriv2(1.-c[0]);
   fjac[1][1] = -xlogx_deriv2(c[1])-xlogx_deriv2(1.-c[1]);
}

/*
 ********************************************************************
 * conc: initial guess and final solution (concentration in each phase)
 * c0: local composition
 ********************************************************************
 */
int KKSdiluteBinaryConcentrationSolver::ComputeConcentration(
   double* const conc,
   const double c0,
   const double hphi,
   const double RTinv,
   const double fA,
   const double fB )
{
   (void) RTinv;

   //std::cout<<"KKSdiluteBinaryConcentrationSolver::ComputeConcentration()"<<endl;
   d_c0 = c0;
   d_hphi = hphi;
   d_fA = fA;
   d_fB = fB;

   int ret= NewtonSolver::ComputeSolution( conc, d_N );
   return ret;
}
