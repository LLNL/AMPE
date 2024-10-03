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
#include "WangSinteringCompositionRHSStrategy.h"
#include "QuatFort.h"
#include "ConcFort.h"
#include "FuncFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <cassert>

const double gas_constant_R_JpKpmol = GASCONSTANT_R_JPKPMOL;

WangSinteringCompositionRHSStrategy::WangSinteringCompositionRHSStrategy(
    const int conc_id, const int phase_id, const double mobility,
    const double parameter_a, const double parameter_b, const double kappa,
    const std::string& avg_func_type)
    : CompositionRHSStrategy(avg_func_type),
      d_conc_id(conc_id),
      d_phase_id(phase_id),
      d_mobility(mobility),
      d_A(parameter_a),
      d_B(parameter_b),
      d_kappa(kappa)
{
   assert(conc_id >= 0);
   assert(mobility > 0.);
   assert(d_kappa > 0.);
   assert(> 0.);

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_set_diffcoeff_timer = tman->getTimer(
       "AMPE::WangSinteringCompositionRHSStrategy::setDiffusionCoeff()");
}

//-----------------------------------------------------------------------

void WangSinteringCompositionRHSStrategy::setDiffusionCoeff(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   // mobility is a scalar value -> nothing to do
   (void)time;

   t_set_diffcoeff_timer->start();

   // tbox::pout<<"QuatIntegrator::setDiffCoeffForConcentration"<<endl;
   assert(hierarchy);

   t_set_diffcoeff_timer->stop();
}

//-----------------------------------------------------------------------

void WangSinteringCompositionRHSStrategy::computeFluxOnPatch(hier::Patch& patch,
                                                             const int flux_id)
{
   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > conc(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_id)));
   assert(conc);
   assert(conc->getGhostCellWidth()[0] > 1);

   std::shared_ptr<pdat::CellData<double> > phi(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_phase_id)));
   assert(phase);
   assert(phase->getGhostCellWidth()[0] > 1);

   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);
   assert(flux->getDepth() == conc->getDepth());
   assert(flux->getGhostCellWidth()[0] > 0);

   // M*grad delta F/delta c
   flux->fillAll(0.);
   std::vector<double> tmp1(pbox.size());
   std::vector<double> tmp2(pbox.size());

   ADD_WANG_SINTERING_FLUX(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                           ifirst(2), ilast(2),
#endif
                           dx, conc->getPointer(), conc->getGhostCellWidth()[0],
                           phi->getPointer(), phi->getGhostCellWidth()[0],
                           phi->getDepth(), d_mobility, d_A, d_B, d_kappa,
                           flux->getPointer(0), flux->getPointer(1),
#if (NDIM == 3)
                           flux->getPointer(2),
#endif
                           flux->getGhostCellWidth()[0], tmp1.data(),
                           tmp2.data());
}
