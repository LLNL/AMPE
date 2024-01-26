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
#ifndef included_CompositionRHSStrategy
#define included_CompositionRHSStrategy

#include "FreeEnergyStrategy.h"
#include "FuncAvgFort.h"

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
      return AVERAGE_FUNC(a, b, d_avg_func_type.c_str());
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
