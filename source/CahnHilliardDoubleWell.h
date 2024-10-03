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
#ifndef included_CahnHilliardDoubleWell
#define included_CahnHilliardDoubleWell

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
class CahnHilliardDoubleWell : public CompositionRHSStrategy
{
 public:
   CahnHilliardDoubleWell(const int conc_id, const int temperature_id,
                          const int diffusion_id, const double mobility,
                          const double ca, const double cb, const double kappa,
                          const double well_scale,
                          const std::string& avg_func_type);

   ~CahnHilliardDoubleWell(){};

   void computeFluxOnPatch(hier::Patch& patch, const int flux_id) override;
   void setDiffusionCoeff(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                          const double time);

 private:
   const int d_conc_id;
   const int d_temperature_id;

   /*!
    * pre-computed diffusion coeffient
    */
   const int d_diffusion_id;

   /*!
    * coefficients of D = D0*exp(-Q0/RT)
    */
   const double d_d0;
   const double d_q0;

   const double d_mobility;
   const double d_ca;
   const double d_cb;
   const double d_kappa;
   const double d_well_scale;

   void addDiffusionCoeff(std::shared_ptr<pdat::CellData<double> > temperature,
                          std::shared_ptr<pdat::SideData<double> > diffusion,
                          const hier::Box& pbox);


   // Timers
   std::shared_ptr<tbox::Timer> t_set_diffcoeff_timer;
};

#endif
