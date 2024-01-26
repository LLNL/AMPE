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
#ifndef included_SpinodalCompositionRHSStrategy
#define included_SpinodalCompositionRHSStrategy

#include "CompositionRHSStrategyWithMobilities.h"

#include <vector>
#include <string>

class SpinodalCompositionRHSStrategy
    : public CompositionRHSStrategyWithMobilities
{
 public:
   SpinodalCompositionRHSStrategy(
       std::shared_ptr<tbox::Database> input_db, const int conc_scratch_id,
       const int phase_scratch_id, const int eta_scratch_id,
       const unsigned int ncompositions, const int conc_a_scratch_id,
       const int conc_b_scratch_id, const int temperature_scratch_id,
       const int diffusion_id, const double kappa, const int Mq_id,
       const std::vector<double>& Q_heat_transport,
       const std::string& phase_interp_func_type,
       const std::string& avg_func_type,
       FreeEnergyStrategy* free_energy_strategy);

   ~SpinodalCompositionRHSStrategy(){};

   void computeFluxOnPatch(hier::Patch& patch, const int flux_id);

   void setDiffusionCoeff(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                          const double time);

 private:
   int d_conc_scratch_id;
   int d_eta_scratch_id;
   int d_diffusion_id;
   int d_temperature_scratch_id;
   int d_phase_scratch_id;
   int d_conc_a_scratch_id;
   int d_conc_b_scratch_id;

   unsigned int d_ncompositions;

   double d_kappa;

   void setDiffusionForConc(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void setDiffusionCoeffForConcOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_c,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff,
       const hier::Box& pbox);
};

#endif
