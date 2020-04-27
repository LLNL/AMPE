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
#ifndef included_CompositionRHSStrategy
#define included_CompositionRHSStrategy

#include "FreeEnergyStrategy.h"
#include "FuncFort.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/Patch.h"

#include <string>

using namespace SAMRAI;

class CompositionRHSStrategy
{
 public:
   CompositionRHSStrategy(const std::string& avg_func_type);

   virtual ~CompositionRHSStrategy(){};

   virtual void computeFluxOnPatch(hier::Patch& patch, const int flux_id) = 0;
   virtual void addFluxFromGradTonPatch(hier::Patch& patch,
                                        const int temperature_id,
                                        const int flux_id);
   virtual void printDiagnostics(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy)
   {
      (void)hierarchy;
   };
   void setZeroFluxAtBoundaryOnPatch(hier::Patch& patch, const int flux_id);
   virtual void addFluxFromAntitrappingonPatch(hier::Patch& patch,
                                               const int phase_scratch_id,
                                               const int dphidt_id,
                                               const double alpha,
                                               const int flux_id);

 protected:
   double average(const double a, const double b) const
   {
      return FORT_AVERAGE_FUNC(a, b, d_avg_func_type.c_str());
   }

   const char* average() const { return d_avg_func_type.c_str(); }

 private:
   /*!
    * function to use to take averages between two phase values
    * at cell centers and define quantity at cell boundary
    */
   std::string d_avg_func_type;
};

#endif
