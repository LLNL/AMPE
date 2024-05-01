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
#ifndef included_MobilityCompositionDiffusionStrategy
#define included_MobilityCompositionDiffusionStrategy

#include "FuncAvgFort.h"
#include "CompositionDiffusionStrategy.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include <string>
using namespace SAMRAI;

class CompositionStrategyMobilities;
class FreeEnergyStrategy;

class MobilityCompositionDiffusionStrategy : public CompositionDiffusionStrategy
{
 public:
   MobilityCompositionDiffusionStrategy(
       const unsigned short ncompositions, const int conc_l_scratch_id,
       const int conc_a_scratch_id, const int pfm_diffusion_l_id,
       const int pfm_diffusion_a_id, const int diffusion_coeff_l_id,
       const int diffusion_coeff_a_id, const std::string& avg_func_type,
       DiffusionInterpolationType diff_interp_type,
       CompositionStrategyMobilities* mobilities_strategy,
       std::shared_ptr<FreeEnergyStrategy> free_energy_strategy);

   ~MobilityCompositionDiffusionStrategy(){};

   /*
    * compute actual diffusion in each phase by weighting diffusion coefficients
    * in each phase with phase variable
    */
   void setDiffusion(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                     const int temperature_id, const int phase_id);

 private:
   double average(const double a, const double b) const
   {
      return AVERAGE_FUNC(a, b, d_avg_func_type.c_str());
   }

   /*
    * Compute diffusion coefficient in each phase
    */
   void setDiffCoeffInEachPhase(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id);

   /*
    * compute PFM diffusion in each phase based on diffusion coefficients
    * in each phase
    */
   void setPFMDiffOnPatch(std::shared_ptr<pdat::CellData<double> > cd_phi,
                          std::shared_ptr<pdat::SideData<double> > sd_d_coeff_l,
                          std::shared_ptr<pdat::SideData<double> > sd_d_coeff_a,
                          std::shared_ptr<pdat::SideData<double> > sd_d_l,
                          std::shared_ptr<pdat::SideData<double> > sd_d_a,
                          const hier::Box& pbox);

   void setDiffCoeffInEachPhaseOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff_l,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff_a,
       const hier::Box& pbox);

   void computeLocalDiffusionMatrixL(const double temperature,
                                     const std::vector<double>& c);
   void computeLocalDiffusionMatrixA(const double temperature,
                                     const std::vector<double>& c);

   unsigned short d_ncompositions;

   int d_conc_l_scratch_id;
   int d_conc_a_scratch_id;

   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_pfm_diffusion_l_id;
   int d_pfm_diffusion_a_id;

   /*!
    * holds data for diffusion coefficients in each phase
    */
   int d_diffusion_coeff_l_id;
   int d_diffusion_coeff_a_id;

   CompositionStrategyMobilities* d_mobilities_strategy;

   // free energy needed to compute diffusion in each phase
   std::shared_ptr<FreeEnergyStrategy> d_free_energy_strategy;

   /*!
    * function to use to take averages between two phase values
    * at cell centers and define quantity at cell boundary
    */
   std::string d_avg_func_type;

   /*!
    * small work arrays
    */
   std::vector<double> d_d2f;
   std::vector<double> d_mobmat;
   std::vector<double> d_local_dmat;
};

#endif
