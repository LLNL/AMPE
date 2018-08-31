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
#ifndef included_CALPHADFunctions
#define included_CALPHADFunctions

#include <vector>

double xlogx( const double x );
double xlogx_deriv( const double x );
double xlogx_deriv2( const double x );

double CALPHADcomputeFMixBinary(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const double conc );
double CALPHADcomputeFMix_derivBinary(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const double conc );
double CALPHADcomputeFMix_deriv2Binary(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const double conc );
double CALPHADcomputeFIdealMixBinary(
   const double rt,
   const double conc );
double CALPHADcomputeFIdealMix_derivBinary(
   const double rt,
   const double conc );
double CALPHADcomputeFIdealMix_deriv2Binary(
   const double rt,
   const double conc );
void CALPHADcomputeFIdealMix_deriv2(
   const double rt,
   const std::vector<double>& conc,
   std::vector<double>& d2fdc2 );
double CALPHADcomputeGMix_deriv2(
   const double l1,
   const double l2,
   const double l3,
   const std::vector<double>& conc,
   const int ic );
double CALPHADcomputeGMix_mixDeriv2(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const std::vector<double>& conc,
   const int ic0,
   const int ic1 );
double CALPHADcomputeFMix_mixDeriv2(
   const double l0,
   const double l1,
   const double l2,
   const double l3,
   const std::vector<double>& concf,
   const int ic0,
   const int ic1 );

double CALPHADcomputePenalty(
   const double alpha1, const double p12, const double p13,
   const double alpha2, const double p22, const double p23,
   const double conc);
double CALPHADcomputeDerivPenalty(
   const double alpha1, const double p12, const double p13,
   const double alpha2, const double p22, const double p23,
   const double conc);
double CALPHADcompute2ndDerivPenalty(
   const double alpha1, const double p12, const double p13,
   const double alpha2, const double p22, const double p23,
   const double conc);

double CALPHADcomputeFMixTernary(
   const double* lAB,
   const double* lAC,
   const double* lBC,
   const double* lABC,
   const double cA,
   const double cB );
double CALPHADcomputeFIdealMixTernary(
   const double rt,
   const double conc0,
   const double conc1 );
void CALPHADcomputeFIdealMix_derivTernary(
   const double rt,
   const double cA,
   const double cB,
   double* deriv );
void CALPHADcomputeFIdealMix_deriv2Ternary(
   const double rt,
   const double cA,
   const double cB,
   double* deriv );

void CALPHADcomputeFMix_derivTernary(
   const double* lAB,
   const double* lAC,
   const double* lBC,
   const double* lABC,
   const double cA,
   const double cB,
   double* deriv );
void CALPHADcomputeFMix_deriv2Ternary(
   const double* lAB,
   const double* lAC,
   const double* lBC,
   const double* lABC,
   const double cA,
   const double cB,
   double* deriv );
#endif
