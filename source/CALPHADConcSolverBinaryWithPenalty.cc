// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
//
#include "CALPHADConcSolverBinaryWithPenalty.h"
#include "CALPHADFunctions.h"
#include <cassert>
#include <iostream>

void CALPHADConcentrationSolverBinaryWithPenalty::computeXi(
    const double* const c, double xi[3]) const
{
   // std::cout<<"CALPHADConcentrationSolverBinaryWithPenalty::computeXi()"<<std::endl;
   CALPHADConcentrationSolverBinary::computeXi(c, xi);

   for (short ii = 0; ii < d_N; ii++) {

      xi[ii] += d_RTinv * computeDerivPenalty(ii, c[ii]);
   }
}

void CALPHADConcentrationSolverBinaryWithPenalty::computeDxiDc(
    const double* const c, double dxidc[3]) const
{
   // std::cout<<"CALPHADConcentrationSolverBinaryWithPenalty::computeDxiDc()"<<std::endl;
   CALPHADConcentrationSolverBinary::computeDxiDc(c, dxidc);

   for (short ii = 0; ii < d_N; ii++) {
      dxidc[ii] += d_RTinv * compute2ndDerivPenalty(ii, c[ii]);
   }
   // std::cout<<"CALPHADConcentrationSolverBinaryWithPenalty::computeDxiDc()
   // done..."<<std::endl;
}

double CALPHADConcentrationSolverBinaryWithPenalty::computeDerivPenalty(
    const short phase_index, const double conc) const
{
   assert(d_penalty_parameters.size() > 0);
   assert(d_penalty_parameters.size() > phase_index);
   const std::vector<double>& penalty_parameters(
       d_penalty_parameters[phase_index]);
   assert(penalty_parameters.size() == 6);

   return CALPHADcomputeDerivPenalty(penalty_parameters[0],
                                     penalty_parameters[1],
                                     penalty_parameters[2],
                                     penalty_parameters[3],
                                     penalty_parameters[4],
                                     penalty_parameters[5], conc);
}

double CALPHADConcentrationSolverBinaryWithPenalty::compute2ndDerivPenalty(
    const short phase_index, const double conc) const
{
   assert(d_penalty_parameters.size() > 0);
   assert(d_penalty_parameters.size() > phase_index);
   const std::vector<double>& penalty_parameters(
       d_penalty_parameters[phase_index]);
   assert(penalty_parameters.size() == 6);

   return CALPHADcompute2ndDerivPenalty(penalty_parameters[0],
                                        penalty_parameters[1],
                                        penalty_parameters[2],
                                        penalty_parameters[3],
                                        penalty_parameters[4],
                                        penalty_parameters[5], conc);
}
