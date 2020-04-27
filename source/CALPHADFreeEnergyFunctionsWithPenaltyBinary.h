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
#ifndef included_CALPHADFreeEnergyFunctionsWithPenalty
#define included_CALPHADFreeEnergyFunctionsWithPenalty

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFunctions.h"
#include <vector>

class CALPHADFreeEnergyFunctionsWithPenaltyBinary
    : public CALPHADFreeEnergyFunctionsBinary
{
 public:
   CALPHADFreeEnergyFunctionsWithPenaltyBinary(
       std::shared_ptr<SAMRAI::tbox::Database> input_db,
       std::shared_ptr<SAMRAI::tbox::Database> newton_db,
       const EnergyInterpolationType energy_interp_func_type,
       const ConcInterpolationType conc_interp_func_type,
       const bool with_third_phase);

   ~CALPHADFreeEnergyFunctionsWithPenaltyBinary(){};

   virtual double computeFreeEnergy(const double temperature,
                                    const double* const conc,
                                    const PhaseIndex pi, const bool gp = false);
   virtual void computeDerivFreeEnergy(const double temperature,
                                       const double* const conc,
                                       const PhaseIndex pi, double* deriv);
   virtual void computeSecondDerivativeFreeEnergy(
       const double temp, const std::vector<double>& conc, const PhaseIndex pi,
       std::vector<double>& d2fdc2);
   virtual bool computeCeqT(const double temperature, const PhaseIndex pi0,
                            const PhaseIndex pi1, double* ceq);

 private:
   std::vector<std::vector<double> > d_penalty_parameters;

   void readParameters(std::shared_ptr<SAMRAI::tbox::Database> calphad_db);
   void setupSolver(std::shared_ptr<tbox::Database> newton_db);

   double computePenalty(const PhaseIndex index, const double conc)
   {
      assert(d_penalty_parameters.size() > 0);
      assert(d_penalty_parameters.size() > index);
      const std::vector<double>& penalty_parameters(
          d_penalty_parameters[static_cast<int>(index)]);
      assert(penalty_parameters.size() == 6);

      double dd =
          CALPHADcomputePenalty(penalty_parameters[0], penalty_parameters[1],
                                penalty_parameters[2], penalty_parameters[3],
                                penalty_parameters[4], penalty_parameters[5],
                                conc);

      return dd;
   }

   double computeDerivPenalty(const PhaseIndex index, const double conc)
   {
      assert(d_penalty_parameters.size() > 0);
      assert(d_penalty_parameters.size() > index);
      const std::vector<double>& penalty_parameters(
          d_penalty_parameters[static_cast<int>(index)]);
      assert(penalty_parameters.size() == 6);

      double dd = CALPHADcomputeDerivPenalty(penalty_parameters[0],
                                             penalty_parameters[1],
                                             penalty_parameters[2],
                                             penalty_parameters[3],
                                             penalty_parameters[4],
                                             penalty_parameters[5], conc);

      return dd;
   }

   double compute2ndDerivPenalty(const PhaseIndex index, const double conc)
   {
      assert(d_penalty_parameters.size() > 0);
      assert(d_penalty_parameters.size() > index);
      const std::vector<double>& penalty_parameters(
          d_penalty_parameters[static_cast<int>(index)]);
      assert(penalty_parameters.size() == 6);

      double dd = CALPHADcompute2ndDerivPenalty(penalty_parameters[0],
                                                penalty_parameters[1],
                                                penalty_parameters[2],
                                                penalty_parameters[3],
                                                penalty_parameters[4],
                                                penalty_parameters[5], conc);
      return dd;
   }
};

#endif
