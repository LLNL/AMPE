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
#include "BiasDoubleWellUTRCFreeEnergyStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "QuatFort.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include <cassert>

using namespace SAMRAI;


BiasDoubleWellUTRCFreeEnergyStrategy::BiasDoubleWellUTRCFreeEnergyStrategy(
    const double alpha, const double gamma,
    MeltingTemperatureStrategy* meltingTstrat)
    : d_alpha(alpha), d_gamma(gamma), d_meltingTstrat(meltingTstrat)
{
   tbox::plog << "BiasDoubleWellFreeFreeEnergyStrategy:" << std::endl;
   tbox::plog << "alpha =" << d_alpha << std::endl;
   tbox::plog << "gamma =" << d_gamma << std::endl;
}

//=======================================================================

void BiasDoubleWellUTRCFreeEnergyStrategy::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int conc_id, const int f_l_id, const int f_a_id,
    const int f_b_id, const int rhs_id)
{
   assert(phase_id >= 0);
   assert(rhs_id >= 0);
   assert(temperature_id >= 0);
   assert(d_meltingTstrat != NULL);

   (void)time;     // unused
   (void)conc_id;  // unused
   (void)f_l_id;   // unused
   (void)f_a_id;   // unused
   (void)f_b_id;

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);

   std::shared_ptr<pdat::CellData<double> > temp(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(temp);

   std::shared_ptr<pdat::CellData<double> > rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(rhs_id)));
   assert(rhs);

   assert(rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   // evaluate melting temperature field
   //(which may depend on composition in linear model for instance)
   d_meltingTstrat->evaluate(patch);

   std::shared_ptr<pdat::CellData<double> > eq_temp(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_meltingTstrat->equilibrium_temperature_id())));
   assert(eq_temp);
#ifdef DEBUG_CHECK_ASSERTIONS
   SAMRAI::math::PatchCellDataNormOpsReal<double> ops;
   double l2rhs = ops.L2Norm(eq_temp, patch.getBox());
   assert(l2rhs == l2rhs);
#endif


   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   COMPUTERHSBIASWELL(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                      ifirst(2), ilast(2),
#endif
                      phase->getPointer(), phase->getGhostCellWidth()[0],
                      temp->getPointer(), temp->getGhostCellWidth()[0], d_alpha,
                      d_gamma, eq_temp->getPointer(), 0, rhs->getPointer(), 0);
}
