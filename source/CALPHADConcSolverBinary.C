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
#include "CALPHADConcSolverBinary.h"
#include "CALPHADFunctions.h"

#include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

static const double s_smallc = 1.0e-8;
static const double s_inv_smallc = 1./s_smallc;

//=======================================================================

CALPHADConcentrationSolverBinary::CALPHADConcentrationSolverBinary(
   const bool with_third_phase )
{
   d_with_third_phase = with_third_phase;

   if ( d_with_third_phase ) {
      d_N = 3;
   }
   else {
      d_N = 2;
   }
}

//=======================================================================

void CALPHADConcentrationSolverBinary::computeXi(const double* const c, double xi[3])const
{
   //std::cout<<"CALPHADConcentrationSolverBinary::computeXi()"<<endl;
   for ( int ii = 0; ii < d_N; ii++ ) {

      double omega = CALPHADcomputeFMix_derivBinary(d_L0[ii],d_L1[ii],d_L2[ii],d_L3[ii],c[ii]);

      double eps = d_fA[ii] - d_fB[ii];

      xi[ii] = d_RTinv * (eps + omega);

      //cout << "d_L2["<<ii<<"] = " << d_L2[ii] << endl;
   }
}
//=======================================================================

// solve for c=(c_L, c_A, c_B)
void CALPHADConcentrationSolverBinary::RHS(
   const double* const c,
   double* const fvec )
{
   double xi[3]={0.,0.,0.};

   computeXi(c,xi);

   fvec[0] =
      -d_c0 + ( 1.0 - d_hphi ) * c[0] + 
      d_hphi * c[1];
#if 1
   fvec[1] = xlogx_deriv(c[0])-xlogx_deriv(1.-c[0])
            -xlogx_deriv(c[1])+xlogx_deriv(1.-c[1])
            +(xi[0] - xi[1] );

   if ( d_N > 2 ) {
      fvec[2] = xlogx_deriv(c[0])-xlogx_deriv(1.-c[0])
               -xlogx_deriv(c[2])+xlogx_deriv(1.-c[2])
               +(xi[0] - xi[2] );
   }
#else      
   if ( d_N > 2 ) {
      fvec[0] += d_hphi * d_heta * ( - c[1] + c[2] );
   }

   const double c0 = c[0] > s_smallc ? c[0] : s_smallc*exp(c[0]*s_inv_smallc-1.);
   const double c1 = c[1] > s_smallc ? c[1] : s_smallc*exp(c[1]*s_inv_smallc-1.);
   const double onemc0 = (1.-c[0]) > s_smallc ? 1.-c[0] : s_smallc*exp((1.-c[0])*s_inv_smallc-1.);
   const double onemc1 = (1.-c[1]) > s_smallc ? 1.-c[1] : s_smallc*exp((1.-c[1])*s_inv_smallc-1.);
   
   fvec[1] =
      ( onemc1 * c0 ) * exp( xi[0] - xi[1] ) - 
      ( c1 * onemc0 );

   if ( d_N > 2 ) {
      const double c2 = c[2] > s_smallc ? c[2] : s_smallc*exp(c[2]*s_inv_smallc-1.);
      const double onemc2 = (1.-c[2]) > s_smallc ? 1.-c[2] : s_smallc*exp((1.-c[2])*s_inv_smallc-1.);
      fvec[2] =
         ( onemc2 * c0 ) * exp( xi[0] - xi[2] ) -
         ( c2 * onemc0 );
      //cout << "xi[0] = " << xi[0] << endl;
      //cout << "xi[1] = " << xi[1] << endl;
      //cout << "xi[2] = " << xi[2] << endl;
      //cout << "c2 = " << c2 << endl;
      //cout << "fvec[2] = " << fvec[2] << endl;
   }
#endif
}

//=======================================================================

void CALPHADConcentrationSolverBinary::computeDxiDc(const double* const c, double xi[3], double dxidc[3])const
{
   //std::cout<<"CALPHADConcentrationSolverBinary::computeDxiDc()"<<endl;
   computeXi(c,xi);
   
   for ( int ii = 0; ii < d_N; ii++ ) {
      dxidc[ii] = d_RTinv * CALPHADcomputeFMix_deriv2Binary(d_L0[ii],d_L1[ii],d_L2[ii],d_L3[ii],c[ii]);
   }
}

//=======================================================================

void CALPHADConcentrationSolverBinary::Jacobian(
   const double* const c,
   double** const fjac )
{
   double xi[3];
   double dxidc[3];
   double expxi[3];
   
   computeDxiDc(c,xi,dxidc);

   fjac[0][0] = ( 1.0 - d_hphi );
   fjac[0][1] = d_hphi;
   if ( d_N > 2 ) {
      fjac[0][1] -= d_hphi * d_heta;
      fjac[0][2] = d_hphi * d_heta;
   }

#if 1
   fjac[1][0] =  dxidc[0]+xlogx_deriv2(c[0])+xlogx_deriv2(1.-c[0]);

   fjac[1][1] = -dxidc[1]-xlogx_deriv2(c[1])-xlogx_deriv2(1.-c[1]);
   
#else
   for ( int ii = 1; ii < d_N; ii++ ) {  // Not using [0]
      expxi[ii] = exp( xi[0] - xi[ii] );
   }

   double c0 = c[0];
   double expc0xs = 1.;
   if( c[0] < s_smallc )
   {
      expc0xs = exp(c[0]*s_inv_smallc-1.); // dc0/dc[0]
      c0      = s_smallc*expc0xs;
   }
   double c1 = c[1];
   double expc1xs = 1.;
   if( c[1] < s_smallc )
   {
      expc1xs = exp(c[1]*s_inv_smallc-1.); // dc1/dc[1]
      c1      = s_smallc*expc1xs;
   }
   const double onemc0 = (1.-c[0]) > s_smallc ? 1.-c[0] : s_smallc*exp((1.-c[0])*s_inv_smallc-1.);
   const double onemc1 = (1.-c[1]) > s_smallc ? 1.-c[1] : s_smallc*exp((1.-c[1])*s_inv_smallc-1.);

   fjac[1][0] = c1*expc0xs + 
      onemc1 * expxi[1] * ( c0 * dxidc[0] + 1.0*expc0xs );
         
   fjac[1][1] = -onemc0*expc1xs +
      c0 * expxi[1] * ( -onemc1 * dxidc[1] - 1.0*expc1xs );

   if ( d_N > 2 ) {
      double c2 = c[2];
      double expc2xs = 1.;
      if( c[2] < s_smallc )
      {
         expc2xs = exp(c[2]*s_inv_smallc-1.);
         c2      = s_smallc*expc2xs;
      }
      const double onemc2 = (1.-c[2]) > s_smallc ? 1.-c[2] : s_smallc*exp((1.-c[2])*s_inv_smallc-1.);
            
      fjac[1][2] = 0.0;

      fjac[2][0] = c2*expc0xs +
         onemc2 * expxi[2] * ( c0 * dxidc[0] + 1.0*expc0xs );

      fjac[2][1] = 0.0;

      fjac[2][2] = -onemc0*expc2xs +
         c0 * expxi[2] * ( -onemc2 * dxidc[2] - 1.0*expc2xs );
   }
#endif
}

//=======================================================================

int CALPHADConcentrationSolverBinary::ComputeConcentration(
   double* const conc,
   const double c0,
   const double hphi,
   const double heta,
   const double RTinv,
   const double* const L0,
   const double* const L1,
   const double* const L2,
   const double* const L3,
   const double* const fA,
   const double* const fB )
{
   //std::cout<<"CALPHADConcentrationSolverBinary::ComputeConcentration()"<<endl;
   d_c0 = c0;
   d_hphi = hphi;
   d_heta = heta;
   d_RTinv = RTinv;

   for ( int ii = 0; ii < d_N; ii++ ) {
      d_L0[ii] = L0[ii];
      d_L1[ii] = L1[ii];
      d_L2[ii] = L2[ii];
      d_L3[ii] = L3[ii];
      d_fA[ii] = fA[ii];
      d_fB[ii] = fB[ii];
   }

   int ret= NewtonSolver::ComputeSolution( conc, d_N );
   return ret;
}
