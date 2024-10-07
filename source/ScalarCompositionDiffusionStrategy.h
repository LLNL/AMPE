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
#ifndef ScalarCompositionDiffusionStrategy_H
#define ScalarCompositionDiffusionStrategy_H

#include "CompositionDiffusionStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include <cstring>

/*!
 * Class to set composition diffusion coefficient to a scalar field that depends
 * on local phase and temperature
 */
class ScalarCompositionDiffusionStrategy : public CompositionDiffusionStrategy
{
 public:
   ScalarCompositionDiffusionStrategy(
       const short norderp, const short norderpA, const short norderpB,
       const bool with3phases, const int pfm_diffusion_id,
       const double D_liquid, const double D_solidA, const double D_solidB,
       const double D0_LA, const double D0_LB, const double D0_AA,
       const double D0_AB, const double D0_BB,
       DiffusionInterpolationType interp_func_type,
       const std::string& avg_func_type);

   /*
    * compute actual diffusion in each phase by weighting diffusion coefficients
    * in each phase with phase variable
    */
   virtual void setDiffusion(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int phase_id) override;

 private:
   const short d_norderp;
   const short d_norderpA;
   const short d_norderpB;

   // distinguish 3 phases implementation (Folch-Plapp)
   const bool d_with3phases;

   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_pfm_diffusion_id;

   double d_D0_liquid;
   double d_D0_solidA;
   double d_D0_solidB;

   /*!
    * additional interfacial diffusion
    */
   double d_d0_LA;
   double d_d0_LB;
   double d_d0_AA;
   double d_d0_AB;
   double d_d0_BB;

   std::string d_avg_func_type;

   bool d_with_phaseB;

   void setDiffusion(std::shared_ptr<hier::Patch> patch,
                     std::shared_ptr<pdat::CellData<double>> phi,
                     std::shared_ptr<pdat::SideData<double>> pfm_diffusion);
   void setDiffusionInterfaces(
       std::shared_ptr<hier::Patch> patch,
       std::shared_ptr<pdat::CellData<double>> phi,
       std::shared_ptr<pdat::SideData<double>> pfm_diffusion);
};

#endif
