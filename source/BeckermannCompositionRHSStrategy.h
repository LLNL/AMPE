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
#ifndef included_BeckermannCompositionRHSStrategy
#define included_BeckermannCompositionRHSStrategy

#include "CompositionRHSStrategy.h"
#include "InterpolationType.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"

#include <string>

using namespace Thermo4PFM;

class QuatModel;

class BeckermannCompositionRHSStrategy : public CompositionRHSStrategy
{
 public:
   BeckermannCompositionRHSStrategy(
       QuatModel* quat_model, const int conc_scratch_id,
       const int phase_scratch_id, const int partition_coeff_scratch_id,
       const int diffusion0_id, const int phase_coupling_diffusion_id,
       const double D_liquid, const double D_solid_A,
       const ConcInterpolationType phase_interp_func_type,
       const std::string& avg_func_type);
   ~BeckermannCompositionRHSStrategy(){};

   void computeFluxOnPatch(hier::Patch& patch, const int flux_id);

   void setDiffusionCoeff(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                          const double time);

 private:
   QuatModel* d_quat_model;

   int d_conc_scratch_id;
   int d_phase_scratch_id;

   int d_partition_coeff_scratch_id;

   int d_diffusion0_id;
   int d_conc_phase_coupling_diffusion_id;

   ConcInterpolationType d_phase_interp_func_type;
   std::string d_avg_func_type;

   double d_D_liquid;
   double d_D_solid_A;

   // Timers
   std::shared_ptr<tbox::Timer> t_set_diffcoeff_timer;

   void setDiffusionCoeffForConcentration(
       const std::shared_ptr<hier::PatchHierarchy>, const int concentration_id,
       const int phase_id, const int conc_tilde_diffusion_id,
       const int conc_phase_coupling_diffusion_id);

   void setDiffusionCoeffForPhaseOnPatch(
       std::shared_ptr<pdat::SideData<double> > sd_phi_diff_coeff,
       std::shared_ptr<pdat::SideData<double> > sd_d0_coeff,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_c,
       std::shared_ptr<pdat::CellData<double> > cd_k, const hier::Box& pbox);
   void computeDiffusionOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::SideData<double> > sd_diffusion0,
       std::shared_ptr<pdat::SideData<double> > sd_diffusion,
       const hier::Box& pbox);
};

#endif
