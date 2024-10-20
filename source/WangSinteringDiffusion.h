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
#ifndef WangSinteringDiffusion_H
#define WangSinteringDiffusion_H

#include "CompositionDiffusionStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include <cstring>

/*!
 * Class to set composition diffusion coefficient to a scalar field that depends
 * on local phase and temperature
 */
class WangSinteringDiffusion : public CompositionDiffusionStrategy
{
 public:
   WangSinteringDiffusion(const int conc_id, const int diffusion_id,
                          const double D_liquid, const double D_solidA,
                          const double D0_LA, const double D0_AA,
                          DiffusionInterpolationType interp_func_type,
                          const std::string& avg_func_type);

   /*
    * compute actual diffusion in each phase by weighting diffusion coefficients
    * in each phase with phase variable
    */
   virtual void setDiffusion(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int phase_id) override;

 protected:
   const int d_conc_id;

   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_diffusion_id;

   double d_D0_liquid;
   double d_D0_solidA;

   /*!
    * additional interfacial diffusion
    */
   double d_d0_LA;
   double d_d0_AA;

   std::string d_avg_func_type;

   virtual void setDiffusion(std::shared_ptr<hier::Patch> patch,
                             std::shared_ptr<pdat::CellData<double>> phi,
                             std::shared_ptr<pdat::CellData<double>> conc,
                             std::shared_ptr<pdat::SideData<double>> diffusion);
   virtual void setDiffusionInterfaces(
       std::shared_ptr<hier::Patch> patch,
       std::shared_ptr<pdat::CellData<double>> phi,
       std::shared_ptr<pdat::CellData<double>> conc,
       std::shared_ptr<pdat::SideData<double>> diffusion);
};

#endif
