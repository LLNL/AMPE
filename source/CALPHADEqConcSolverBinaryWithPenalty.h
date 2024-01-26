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
#ifndef included_CALPHADEqConcSolverBinaryWithPenalty
#define included_CALPHADEqConcSolverBinaryWithPenalty

#include "CALPHADEqConcSolverBinary.h"
#include "Phases.h"

#include <vector>

using namespace ampe_thermo;

class CALPHADEqConcentrationSolverBinaryWithPenalty
    : public CALPHADEqConcentrationSolverBinary
{
 public:
   CALPHADEqConcentrationSolverBinaryWithPenalty(){};

   virtual ~CALPHADEqConcentrationSolverBinaryWithPenalty(){};

   int ComputeConcentrationWithPenalty(
       double* const conc, const double RTinv, const double* const L0,
       const double* const L1, const double* const L2, const double* const L3,
       const double* const fA, const double* const fB,
       std::vector<std::vector<double> >& penalty_parameters);

 private:
   virtual void RHS(const double* const x, double* const fvec);

   virtual void Jacobian(const double* const x, double** const fjac);

   std::vector<double> d_penalty_parametersL;
   std::vector<double> d_penalty_parametersS;
};

#endif
