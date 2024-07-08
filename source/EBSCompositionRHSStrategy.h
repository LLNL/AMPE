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
// Ref: Eiken, Boettger, Steinbach, PRE 73, 066122 (2006)
#ifndef included_EBSCompositionRHSStrategy
#define included_EBSCompositionRHSStrategy

#include "CompositionRHSStrategy.h"
#include "CALPHADMobility.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CompositionDiffusionStrategy.h"

#include <vector>
#include <string>

class CompositionDiffusionStrategy;

class EBSCompositionRHSStrategy : public CompositionRHSStrategy
{
 public:
   EBSCompositionRHSStrategy(
       const int phase_scratch_id, const unsigned short ncompositions,
       const int conc_l_scratch_id, const int conc_a_scratch_id,
       const int conc_b_scratch_id, const int temperature_scratch_id,
       const int diffusion_l_id, const int diffusion_a_id,
       const int diffusion_b_id, const std::vector<int> diffusion_precond_id,
       const std::string& avg_func_type,
       std::shared_ptr<FreeEnergyStrategy> free_energy_strategy,
       std::shared_ptr<CompositionDiffusionStrategy>
           diffusion_for_conc_in_phase,
       const std::string flux_type);

   ~EBSCompositionRHSStrategy(){};

   void addFluxFromAntitrappingonPatch(hier::Patch& patch,
                                       const int phase_scratch_id,
                                       const int dphidt_id, const double alpha,
                                       const int flux_id);

   void computeFluxOnPatch(hier::Patch& patch, const int flux_id);

   /*
    * Take sum of diffusion coefficients in each phase
    * to get diffusion used in preconditioner
    * (same diffusion coefficient as used in KKS equations)
    * Assumes coefficients in each pahse includes a phase fraction weight
    */
   void setDiffusionCoeffForPreconditioner(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

 protected:
   virtual void setDiffusionCoeff(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time);

 private:
   unsigned short d_ncompositions;

   int d_phase_scratch_id;

   int d_conc_l_scratch_id;
   int d_conc_a_scratch_id;
   int d_conc_b_scratch_id;

   int d_temperature_scratch_id;

   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_diffusion_l_id;
   int d_diffusion_a_id;
   int d_diffusion_b_id;

   bool d_with_phaseB;

   bool d_with_diffusion_for_preconditioner;
   std::vector<int> d_diffusion_precond_id;

   std::shared_ptr<CompositionDiffusionStrategy> d_diffusion_for_conc_in_phase;

   // free energy needed to compute diffusion in each phase
   std::shared_ptr<FreeEnergyStrategy> d_free_energy_strategy;

   void setDiffusionCoeffForPreconditionerOnPatch(
       std::shared_ptr<pdat::SideData<double> > sd_d_l,
       std::shared_ptr<pdat::SideData<double> > sd_d_a,
       std::shared_ptr<pdat::SideData<double> > sd_d_b,
       const int depth_in_Dmatrix,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff,
       const hier::Box& pbox, const int depth);

#if (NDIM == 3)
   void (*d_add_flux)(const int& ifirst0, const int& ilast0, const int& ifirst1,
                      const int& ilast1, const int& ifirst2, const int& ilast2,
                      const double* dx, const double* conc, const int& ngconc,
                      const int& ncomp, const double* diffconc0,
                      const double* diffconc1, const double* diffconc2,
                      const int& ngdiffconc, const double* flux0,
                      const double* flux1, const double* flux2,
                      const int& ngflux, const int* physbc);
#else
   void (*d_add_flux)(const int& ifirst0, const int& ilast0, const int& ifirst1,
                      const int& ilast1, const double* dx, const double* conc,
                      const int& ngconc, const int& ncomp,
                      const double* diffconc0, const double* diffconc1,
                      const int& ngdiffconc, const double* flux0,
                      const double* flux1, const int& ngflux,
                      const int* physbc);
#endif
};

#endif
