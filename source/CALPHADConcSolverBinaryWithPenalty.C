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
#include "CALPHADConcSolverBinaryWithPenalty.h"
#include "CALPHADFunctions.h"
#include <cassert>
#include <iostream>

void CALPHADConcentrationSolverBinaryWithPenalty::computeXi(const double* const c, 
                                                      double xi[3])const
{
   //std::cout<<"CALPHADConcentrationSolverBinaryWithPenalty::computeXi()"<<std::endl;
   CALPHADConcentrationSolverBinary::computeXi(c,xi);
   
   for ( short ii = 0; ii < d_N; ii++ ) {

      xi[ii] += d_RTinv * computeDerivPenalty(ii, c[ii]);

   }
}

void CALPHADConcentrationSolverBinaryWithPenalty::computeDxiDc(const double* const c, 
                                                         double dxidc[3])const
{
   //std::cout<<"CALPHADConcentrationSolverBinaryWithPenalty::computeDxiDc()"<<std::endl;
   CALPHADConcentrationSolverBinary::computeDxiDc(c,dxidc);

   for ( short ii = 0; ii < d_N; ii++ ) {
      dxidc[ii] += d_RTinv * compute2ndDerivPenalty(ii, c[ii]);
   }
   //std::cout<<"CALPHADConcentrationSolverBinaryWithPenalty::computeDxiDc() done..."<<std::endl;
}

double CALPHADConcentrationSolverBinaryWithPenalty::computeDerivPenalty(const short phase_index, 
                                                                  const double conc)const
{
   assert( d_penalty_parameters.size()>0 );
   assert( d_penalty_parameters.size()>phase_index );
   const std::vector<double>& penalty_parameters(d_penalty_parameters[phase_index]);
   assert( penalty_parameters.size()==6 );
   
   return CALPHADcomputeDerivPenalty(penalty_parameters[0], penalty_parameters[1],  penalty_parameters[2],
                                     penalty_parameters[3], penalty_parameters[4],  penalty_parameters[5],
                                     conc);
}

double CALPHADConcentrationSolverBinaryWithPenalty::compute2ndDerivPenalty(const short phase_index, 
                                                                     const double conc)const
{
   assert( d_penalty_parameters.size()>0 );
   assert( d_penalty_parameters.size()>phase_index );
   const std::vector<double>& penalty_parameters(d_penalty_parameters[phase_index]);
   assert( penalty_parameters.size()==6 );
   
   return CALPHADcompute2ndDerivPenalty(penalty_parameters[0], penalty_parameters[1],  penalty_parameters[2],
                                        penalty_parameters[3], penalty_parameters[4],  penalty_parameters[5],
                                        conc);
}
