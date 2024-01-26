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
#ifndef included_CALPHADConcSolverBinaryWithPenalty
#define included_CALPHADConcSolverBinaryWithPenalty

#include "CALPHADConcSolverBinary.h"

#include <vector>

using namespace ampe_thermo;

class CALPHADConcentrationSolverBinaryWithPenalty
    : public CALPHADConcentrationSolverBinary
{
 public:
   CALPHADConcentrationSolverBinaryWithPenalty(
       const bool with_third_phase,
       const std::vector<std::vector<double> >& penalty_parameters)
       : CALPHADConcentrationSolverBinary(with_third_phase),
         d_penalty_parameters(penalty_parameters){};

   virtual ~CALPHADConcentrationSolverBinaryWithPenalty(){};

   virtual void computeXi(const double* const c, double xi[3]) const;

   virtual void computeDxiDc(const double* const c, double dxidc[3]) const;

 private:
   std::vector<std::vector<double> > d_penalty_parameters;

   double computeDerivPenalty(const short phase_index, const double conc) const;
   double compute2ndDerivPenalty(const short phase_index,
                                 const double conc) const;
};

#endif
