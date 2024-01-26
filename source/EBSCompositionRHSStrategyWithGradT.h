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
#ifndef included_EBSCompositionRHSStrategyWithGradT
#define included_EBSCompositionRHSStrategyWithGradT

#include "EBSCompositionRHSStrategy.h"

class CompositionDiffusionStrategy;
class CompositionStrategyMobilities;

class EBSCompositionRHSStrategyWithGradT : public EBSCompositionRHSStrategy
{
 public:
   EBSCompositionRHSStrategyWithGradT(
       const int phase_scratch_id, const int eta_scratch_id,
       const unsigned short ncompositions, const int conc_l_scratch_id,
       const int conc_a_scratch_id, const int conc_b_scratch_id,
       const int temperature_scratch_id, const int diffusion_l_id,
       const int diffusion_a_id, const int diffusion_b_id, const int Mq_id,
       const std::vector<double>& Q_heat_transport,
       const std::vector<int> diffusion_precond_id,
       const std::string& avg_func_type,
       std::shared_ptr<FreeEnergyStrategy> free_energy_strategy,
       CompositionStrategyMobilities* mobilities_strategy,
       std::shared_ptr<CompositionDiffusionStrategy>
           diffusion_for_conc_in_phase);

   ~EBSCompositionRHSStrategyWithGradT(){};

   virtual void setDiffusionCoeff(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time);

   // void printDiagnostics(const std::shared_ptr<hier::PatchHierarchy >
   // hierarchy);

 private:
   unsigned short d_ncompositions;

   std::vector<double> d_Q_heat_transport;

   int d_phase_scratch_id;
   int d_eta_scratch_id;

   int d_conc_l_scratch_id;
   int d_conc_a_scratch_id;
   int d_conc_b_scratch_id;

   int d_Mq_id;

   int d_temperature_scratch_id;

   bool d_with_three_phases;

   bool d_with_diffusion_for_preconditioner;
   std::vector<int> d_diffusion_precond_id;

   CompositionStrategyMobilities* d_mobilities_strategy;

   void addFluxFromGradTonPatch(hier::Patch& patch, const int temperature_id,
                                const int flux_id);

   void setDiffusionCoeffForTOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_c,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff,
       const hier::Box& pbox);
   void setDiffusionCoeffForTOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::SideData<double> > mq, const hier::Box& pbox);
   void setDiffusionCoeffForT(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int concentration_l_id, const int concentration_a_id,
       const int concentration_b_id, const int temperature_id);
};

#endif
