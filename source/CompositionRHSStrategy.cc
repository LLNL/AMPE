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
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include "ConcFort.h"
#include "QuatFort.h"
#include "CompositionRHSStrategy.h"

#include <cassert>


CompositionRHSStrategy::CompositionRHSStrategy(const std::string& avg_func_type)
    : d_avg_func_type(avg_func_type){};

void CompositionRHSStrategy::addFluxFromGradTonPatch(hier::Patch& patch,
                                                     const int temperature_id,
                                                     const int flux_id)
{
   (void)patch;
   (void)temperature_id;
   (void)flux_id;

   TBOX_ERROR(
       "CompositionRHSStrategy::addFluxFromGradTonPatch() not implemented..."
       << std::endl);
}

void CompositionRHSStrategy::addFluxFromAntitrappingonPatch(
    hier::Patch& patch, const int phase_scratch_id, const int dphidt_id,
    const double alpha, const int flux_id)
{
   (void)patch;
   (void)phase_scratch_id;
   (void)dphidt_id;
   (void)alpha;
   (void)flux_id;

   TBOX_ERROR(
       "CompositionRHSStrategy::addFluxFromAntitrappingonPatch() not "
       "implemented..."
       << std::endl);
}

//-----------------------------------------------------------------------

void CompositionRHSStrategy::setZeroFluxAtBoundaryOnPatch(hier::Patch& patch,
                                                          const int flux_id)
{
   std::shared_ptr<geom::CartesianPatchGeometry> pg(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));

   // Get face boundary boxes.
   const std::vector<hier::BoundaryBox> bdry = pg->getCodimensionBoundaries(1);

   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);

   const hier::Box& gbox = flux->getGhostBox();
   const hier::Box& pbox = patch.getBox();


   const hier::Index plo(pbox.lower());
   const hier::Index pup(pbox.upper());

   const hier::Index glo(gbox.lower());
   for (int i = 0; i < NDIM; i++)
      assert(glo(i) == plo(i));

   for (size_t i = 0; i < bdry.size(); i++) {

      const int locind = bdry[i].getLocationIndex();
      const int dir = locind >> 1;
      const int side = locind & 1;
      if (pg->getTouchesRegularBoundary(dir, side)) {
         // std::cout<<"CompositionRHSStrategy::setZeroFluxAtBoundaryOnPatch"
         //          <<" for side "<<side<<" in direction "<<dir<<endl;

         hier::Index gup(gbox.upper());
         gup(dir)++;

         hier::Index lo(plo);
         hier::Index up(pup);

         // set direction perpendicular to boundary
         if (side == 0) {  // lower side
            up(dir) = plo(dir);
            lo(dir) = plo(dir);
         } else {  // upper side
            lo(dir) = pup(dir) + 1;
            up(dir) = pup(dir) + 1;
         }

#if (NDIM == 3)
         SETTOZERO(lo(0), lo(1), lo(2), up(0), up(1), up(2), glo(0), glo(1),
                   glo(2), gup(0), gup(1), gup(2), flux->getPointer(dir));
#endif
#if (NDIM == 2)
         SETTOZERO(lo(0), lo(1), up(0), up(1), glo(0), glo(1), gup(0), gup(1),
                   flux->getPointer(dir));
#endif
      }
   }  // loop over bdry boxes
}
