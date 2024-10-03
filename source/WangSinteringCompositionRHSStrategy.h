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
#ifndef included_WangSinteringCompositionRHSStrategy
#define included_WangSinteringCompositionRHSStrategy

#include "CompositionRHSStrategy.h"
#include "InterpolationType.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"

#include <string>

/*!
 * D_conc is a phase weighted average, and is temperature dependent based
 * on D*exp(-Q/RT) for each phase
 */
class WangSinteringCompositionRHSStrategy : public CompositionRHSStrategy
{
 public:
   WangSinteringCompositionRHSStrategy(const int conc_id, const int phase_id,
                                       const double mobility,
                                       const double parameter_a,
                                       const double parameter_b,
                                       const double kappa,
                                       const std::string& avg_func_type);

   ~WangSinteringCompositionRHSStrategy(){};

   void computeFluxOnPatch(hier::Patch& patch, const int flux_id) override;
   void setDiffusionCoeff(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                          const double time);

 private:
   int d_conc_id;
   int d_phase_id;

   double d_mobility;
   double d_A;
   double d_B;
   double d_kappa;

   // Timers
   std::shared_ptr<tbox::Timer> t_set_diffcoeff_timer;
};

#endif
