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
#include "CALPHADFunctions.h"

#include <cmath>
#include <cassert>
using namespace std;

static const double s_smallx = 1.0e-8;
static const double s_inv_smallx = 1./s_smallx;
static const double s_log_smallx = log( s_smallx );
static const double s_smallx_log_smallx = s_smallx * s_log_smallx;

//-----------------------------------------------------------------------
// C2 extension of x(log(x) function for x<=smallx
//
double xlogx( const double x )
{
   double r;
   
   if ( x > s_smallx ) {
      r = x * log( x );
   }
   else {
      r = s_smallx_log_smallx +
         ( x - s_smallx ) * s_log_smallx +
         0.5 * ( x * x * s_inv_smallx - s_smallx );
   }

   return r;
}

double xlogx_deriv( const double x )
{
   double r;
   
   if ( x > s_smallx ) {
      r = log( x ) + 1.0;
   }
   else {
      r = s_log_smallx + x * s_inv_smallx;
   }

   return r;
}

double xlogx_deriv2( const double x )
{
   double r;
   
   if ( x > s_smallx ) {
      r = 1./x;
   }
   else {
      r = s_inv_smallx;
   }

   return r;
}

double CALPHADcomputeFMixBinary(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const double conc )
{
   double two_c_minus_one = 2.0 * conc - 1.0;

   double fmix =
      conc * ( 1.0 - conc ) *
      ( l0 +
        l1 * two_c_minus_one +
        l2 * two_c_minus_one * two_c_minus_one +
        l3 * two_c_minus_one * two_c_minus_one * two_c_minus_one );

   return fmix;
}   

double CALPHADcomputeFMix_derivBinary(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const double conc )
{
   double two_c_minus_one = 2.0 * conc - 1.0;

   double fmix_deriv =
      conc * ( 1.0 - conc ) * (
         2.0 * l1 +
         4.0 * l2 * two_c_minus_one +
         6.0 * l3 * two_c_minus_one * two_c_minus_one )
      -
      two_c_minus_one * (
         l0 +
         l1 * two_c_minus_one +
         l2 * two_c_minus_one * two_c_minus_one +
         l3 * two_c_minus_one * two_c_minus_one * two_c_minus_one );

   return fmix_deriv;
}

// double CALPHADcomputeFMix_deriv(
//    const double l0,
//    const double l1,
//    const double l2,
//    const double conc )
// {
//    double fmix_deriv =
//       -l0*(2.*conc-1.)
//       -l1*(6.*conc*conc-6.*conc+1.)
//       -l2*(2.*conc-1.)*(8.*conc*conc-8.*conc+1.);

//    return fmix_deriv;
// }

double CALPHADcomputeFMix_deriv2Binary(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const double conc )
{
   const double conc2=conc*conc;
   double fmix_deriv2 =
      -2.*( l0
           +l1*(-3.+6.*conc)
           +l2*(5.-24.*conc+24.*conc2)
           +l3*(-7.+54.*conc-120.*conc2+80.*conc2*conc) );

   return fmix_deriv2;
}

double CALPHADcomputeFIdealMixBinary(
   const double rt,
   const double conc )
{
   double fmix =
      rt * ( xlogx( conc ) + xlogx( 1.0 - conc ) );

   return fmix;
}

double CALPHADcomputeFIdealMix_derivBinary(
   const double rt,
   const double conc )
{
   double fmix_deriv =
      rt * ( xlogx_deriv( conc ) - xlogx_deriv( 1.0 - conc ) );

   return fmix_deriv;
}

double CALPHADcomputeFIdealMix_deriv2Binary(
   const double rt,
   const double conc )
{
   double fmix_deriv2 =
      rt * ( xlogx_deriv2( conc ) + xlogx_deriv2( 1.0 - conc ) );

   return fmix_deriv2;
}

void CALPHADcomputeFIdealMix_deriv2Ternary(
   const double rt,
   const double cA,
   const double cB,
   double* deriv )
{
   deriv[0] =
      rt * ( xlogx_deriv2( cA ) 
           + xlogx_deriv2( 1.0 - cA - cB ) );
   deriv[1] =
      rt * ( xlogx_deriv2( cB ) 
           + xlogx_deriv2( 1.0 - cA - cB ) );
}

void CALPHADcomputeFIdealMix_deriv2(
   const double rt,
   const vector<double>& conc,
   vector<double>& d2fdc2 )
{
   const short nc=(int)conc.size();
   double cN = 1.;
   for(short ic=0;ic<nc;ic++) cN -= conc[ic];
   
   //diagonal terms
   for(short ic=0;ic<nc;ic++)
   {
      d2fdc2[ic+nc*ic]= rt * ( xlogx_deriv2( conc[ic] ) + xlogx_deriv2( cN ) );
   }
   
   const double val = rt * xlogx_deriv2( cN );
   for(short ic=0;ic<nc;ic++)
   for(short jc=0;jc<ic;jc++){
      d2fdc2[ic+nc*jc]= val;
      d2fdc2[jc+nc*ic]= val;
   }
}

double CALPHADcomputeGMix_mixDeriv2(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const double c0,
   const double c1 )
{
   const double dc=(c0-c1);
   const double dc2=dc*dc;
   
   return l0+2.*l1*dc+l2*(3.*dc2-2.*c0*c1)+l3*(4.*dc2*dc-6.*c0*c1*dc);
}

// multicomponents
double CALPHADcomputeGMix_deriv2(
   const double l1,
   const double l2,
   const double l3,
   const vector<double>& conc,
   const int ic )
{
   assert( l3<=1.e-15 );
   
   const int n=conc.size();
   
   double gmix_deriv2 = 0.;
   for(short i=0;i<ic;i++)
      gmix_deriv2 -= 2.*conc[i]*
         (l1+2.*l2*(conc[i]-conc[ic]));
   for(short i=ic+1;i<n;i++)
      gmix_deriv2 += 2.*conc[i]*
         (l1+2.*l2*(conc[ic]-conc[i]));
   for(short i=0;i<ic;i++)
      gmix_deriv2 += conc[i]*conc[ic]*l2*2.;
   for(short i=ic+1;i<n;i++)
      gmix_deriv2 += conc[i]*conc[ic]*l2*2.;

   return gmix_deriv2;
}

double CALPHADcomputeGMix_mixDeriv2(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const vector<double>& conc,
   const int ic0,
   const int ic1 )
{
   if( ic0==ic1 )
      return CALPHADcomputeGMix_deriv2(l1, l2, l3, conc, ic0 );
   else
      return CALPHADcomputeGMix_mixDeriv2(l0, l1, l2, l3, conc[ic0], conc[ic1] );
}

double CALPHADcomputeFMix_mixDeriv2(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const vector<double>& concf,
   const int ic0,
   const int ic1 )
{
   const int icN=(int)concf.size();

   vector<double> concg(concf);
   double cN=1.;
   for(short i=0;i<icN;i++)cN-=concf[i];
   concg.push_back(cN);

   double df2=
      CALPHADcomputeGMix_mixDeriv2(l0,l1,l2,l3,concg,ic0,ic1)
     -CALPHADcomputeGMix_mixDeriv2(l0,l1,l2,l3,concg,ic0,icN)
     -CALPHADcomputeGMix_mixDeriv2(l0,l1,l2,l3,concg,ic1,icN)
     +CALPHADcomputeGMix_deriv2(l1,l2,l3,concg,icN);

   return df2;
}

double CALPHADcomputePenalty(
   const double alpha1, const double p12, const double p13,
   const double alpha2, const double p22, const double p23,
   const double conc)
{
   if( conc<alpha1 )
   {
      const double dd=(alpha1-conc);
      return (p12+p13*dd)*dd*dd;
   }
   if( conc>=alpha2 )
   {
      const double dd=(conc-alpha2);
      return (p22+p23*dd)*dd*dd;
   }
   return 0.;
}

double CALPHADcomputeDerivPenalty(
   const double alpha1, const double p12, const double p13,
   const double alpha2, const double p22, const double p23,
   const double conc)
{
   if( conc<alpha1 )
   {
      const double dd=(alpha1-conc);
      return (-2.*p12-3.*p13*dd)*dd;
   }
   if( conc>=alpha2 )
   {
      const double dd=(conc-alpha2);
      return (2.*p22+3.*p23*dd)*dd;
   }
   return 0.;
}

double CALPHADcompute2ndDerivPenalty(
   const double alpha1, const double p12, const double p13,
   const double alpha2, const double p22, const double p23,
   const double conc)
{
   if( conc<alpha1 )
   {
      const double dd=(alpha1-conc);
      return (2.*p12+6.*p13*dd);
   }
   if( conc>=alpha2 )
   {
      const double dd=(conc-alpha2);
      return (2.*p22+6.*p23*dd);
   }
   return 0.;
}

double CALPHADcomputeFMixTernary(
   const double* lAB,
   const double* lAC,
   const double* lBC,
   const double cA,
   const double cB )
{
   double cC=1.-cA-cB;
   
   double fmix =
      cA*cB*
      ( lAB[0] +
        lAB[1] * (cA-cB) +
        lAB[2] * (cA-cB)*(cA-cB) +
        lAB[3] * (cA-cB)*(cA-cB)*(cA-cB) );

   fmix +=
      cA*cC*
      ( lAC[0] +
        lAC[1] * (cA-cC) +
        lAC[2] * (cA-cC)*(cA-cC) +
        lAC[3] * (cA-cC)*(cA-cC)*(cA-cC) );

   fmix +=
      cB*cC*
      ( lBC[0] +
        lBC[1] * (cB-cC) +
        lBC[2] * (cB-cC)*(cB-cC) +
        lBC[3] * (cB-cC)*(cB-cC)*(cB-cC) );

   return fmix;
}   

double CALPHADcomputeFIdealMixTernary(
   const double rt,
   const double conc0,
   const double conc1 )
{
   double fmix =
      rt * ( xlogx( conc0 ) 
           + xlogx( conc1 ) 
           + xlogx( 1.0 - conc0-conc1 ) );

   return fmix;
}

void CALPHADcomputeFMix_derivTernary(
   const double* lAB,
   const double* lAC,
   const double* lBC,
   const double cA,
   const double cB,
   double* deriv )
{
   double cC=1.-cA-cB;
   
   deriv[0] =
      cB*
      ( lAB[0] +
        lAB[1] * (cA-cB) +
        lAB[2] * (cA-cB)*(cA-cB) +
        lAB[3] * (cA-cB)*(cA-cB)*(cA-cB) );

   deriv[0] +=
      cA*cB*
      ( lAB[1] +
        lAB[2] * 2.*(cA-cB) +
        lAB[3] * 3.*(cA-cB)*(cA-cB) );

   deriv[0] +=
      cC*
      ( lAC[0] +
        lAC[1] * (cA-cC) +
        lAC[2] * (cA-cC)*(cA-cC) +
        lAC[3] * (cA-cC)*(cA-cC)*(cA-cC) );

   deriv[0] +=
      cA*cC*
      ( lAC[1] +
        lAC[2] * 2. *(cA-cC) +
        lAC[3] * 3. *(cA-cC)*(cA-cC) );


   deriv[1] =
      cA*
      ( lAB[0] +
        lAB[1] * (cA-cB) +
        lAB[2] * (cA-cB)*(cA-cB) +
        lAB[3] * (cA-cB)*(cA-cB)*(cA-cB) );

   deriv[1] +=
      -1.*cA*cB*
      ( lAB[1] +
        lAB[2] * 2. *(cA-cB) +
        lAB[3] * 3. *(cA-cB)*(cA-cB) );

   deriv[1] +=
      cC*
      ( lBC[0] +
        lBC[1] * (cB-cC) +
        lBC[2] * (cB-cC)*(cB-cC) +
        lBC[3] * (cB-cC)*(cB-cC)*(cB-cC) );

   deriv[1] +=
      cB*cC*
      ( lBC[1]  +
        lBC[2] * 2. *(cB-cC) +
        lBC[3] * 3. *(cB-cC)*(cB-cC) );
}

void CALPHADcomputeFMix_deriv2Ternary(
   const double* lAB,
   const double* lAC,
   const double* lBC,
   const double cA,
   const double cB,
   double* deriv )
{
   assert( deriv!=NULL );
   
   double cC=1.-cA-cB;
   
   // /dcA*dcA
   deriv[0] =
      cB*
      ( lAB[1] +
        lAB[2] * 2.*(cA-cB) +
        lAB[3] * 3.*(cA-cB)*(cA-cB) );

   deriv[0] +=
      cB*
      ( lAB[1] +
        lAB[2] * 2.*(cA-cB) +
        lAB[3] * 3.*(cA-cB)*(cA-cB) );

   deriv[0] +=
      cA*cB*
      ( lAB[2] * 2. +
        lAB[3] * 6.*(cA-cB) );

   deriv[0] +=
      cC*
      ( lAC[1]  +
        lAC[2] * 2.*(cA-cC) +
        lAC[3] * 3.*(cA-cC)*(cA-cC) );

   deriv[0] +=
      cC*
      ( lAC[1] +
        lAC[2] * 2. *(cA-cC) +
        lAC[3] * 3. *(cA-cC)*(cA-cC) );

   deriv[0] +=
      cA*cC*
      ( lAC[2] * 2.  +
        lAC[3] * 6. *(cA-cC) );

   // /dcA*dcB
   deriv[1] =
      ( lAB[0] +
        lAB[1] * (cA-cB) +
        lAB[2] * (cA-cB)*(cA-cB) +
        lAB[3] * (cA-cB)*(cA-cB)*(cA-cB) );

   deriv[1] +=
      cB*
      ( -lAB[1]  +
        -lAB[2] * 2.*(cA-cB) +
        -lAB[3] * 3.*(cA-cB)*(cA-cB) );

   deriv[1] +=
      cA*
      ( lAB[1] +
        lAB[2] * 2.*(cA-cB) +
        lAB[3] * 3.*(cA-cB)*(cA-cB) );

   deriv[1] +=
      cA*cB*
      ( -lAB[2] * 2. +
        -lAB[3] * 6.*(cA-cB) );

   // /dcB*dcA
   deriv[2] =
      ( lAB[0] +
        lAB[1] * (cA-cB) +
        lAB[2] * (cA-cB)*(cA-cB) +
        lAB[3] * (cA-cB)*(cA-cB)*(cA-cB) );

   deriv[2] +=
      cA*
      ( lAB[1]  +
        lAB[2] * 2.*(cA-cB) +
        lAB[3] * 3.*(cA-cB)*(cA-cB) );

   deriv[2] +=
      -1.*cB*
      ( lAB[1] +
        lAB[2] * 2. *(cA-cB) +
        lAB[3] * 3. *(cA-cB)*(cA-cB) );

   deriv[2] +=
      -1.*cA*cB*
      ( lAB[2] * 2.  +
        lAB[3] * 6. *(cA-cB) );

   deriv[2] +=
      cC*
      ( lBC[0] +
        lBC[1] * (cB-cC) +
        lBC[2] * (cB-cC)*(cB-cC) +
        lBC[3] * (cB-cC)*(cB-cC)*(cB-cC) );
   
   // /dcB*dcB

   deriv[3] =
      cA*
      ( -lAB[1] +
        -lAB[2] * 2.*(cA-cB) +
        -lAB[3] * 3.*(cA-cB)*(cA-cB)*(cA-cB) );

   deriv[3] +=
      -1.*cA*
      ( lAB[1] +
        lAB[2] * 2. *(cA-cB) +
        lAB[3] * 3. *(cA-cB)*(cA-cB) );

   deriv[3] +=
      -1.*cA*cB*
      ( -lAB[2] * 2. +
        -lAB[3] * 6. *(cA-cB) );

   deriv[3] +=
      cC*
      ( lBC[1]  +
        lBC[2] * 2.*(cB-cC) +
        lBC[3] * 3.*(cB-cC)*(cB-cC) );

   deriv[3] +=
      cC*
      ( lBC[1]  +
        lBC[2] * 2. *(cB-cC) +
        lBC[3] * 3. *(cB-cC)*(cB-cC) );

   deriv[3] +=
      cB*cC*
      ( lBC[2] * 2. +
        lBC[3] * 6. *(cB-cC) );
}

void CALPHADcomputeFIdealMix_derivTernary(
   const double rt,
   const double cA,
   const double cB,
   double* deriv )
{
   deriv[0] =
      rt * ( xlogx_deriv( cA ) - xlogx_deriv( 1.0 - cA-cB ) );

   deriv[1] =
      rt * ( xlogx_deriv( cB ) - xlogx_deriv( 1.0 - cA-cB ) );

}
