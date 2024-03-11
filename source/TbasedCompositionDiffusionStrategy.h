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
#ifndef TbasedCompositionDiffusionStrategy_H
#define TbasedCompositionDiffusionStrategy_H

#include "CompositionDiffusionStrategy.h"

#include <cstring>

class TbasedCompositionDiffusionStrategy : public CompositionDiffusionStrategy
{
 public:
   TbasedCompositionDiffusionStrategy(
       const short norderp, const short norderpA, const bool with3phases,
       const int pfm_diffusion_l_id, const int pfm_diffusion_a_id,
       const int pfm_diffusion_b_id, const double D_liquid,
       const double Q0_liquid, const double D_solid_A, const double Q0_solid_A,
       const double D_solid_B, const double Q0_solid_B,
       DiffusionInterpolationType interp_func_type,
       const std::string& avg_func_type);


   /*
    * compute actual diffusion in each phase by weighting diffusion coefficients
    * in each phase with phase variable
    */
   virtual void setDiffusion(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int phase_id, const int eta_id);

 private:
   const short d_norderp;
   const short d_norderpA;

   // distinguish 3 phases implementation (Folch-Plapp)
   const bool d_with3phases;

   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_pfm_diffusion_l_id;
   int d_pfm_diffusion_a_id;
   int d_pfm_diffusion_b_id;

   double d_D_liquid;
   double d_Q0_liquid;

   double d_D_solid_A;
   double d_Q0_solid_A;

   double d_D_solid_B;
   double d_Q0_solid_B;

   std::string d_avg_func_type;

   bool d_with_phaseB;
};

#endif
