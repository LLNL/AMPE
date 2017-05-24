#include "PhaseFluxStrategyIsotropic.h"
#include "QuatFort.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

void PhaseFluxStrategyIsotropic::computeFluxes(const boost::shared_ptr<hier::PatchLevel> level,
                   const int phase_id,
                   const int quat_id,
                   const int flux_id)
{
   //this strategy is independent of grain orientation
   (void) quat_id;
   
   //  Compute phase "flux" on patches in level.
   for ( hier::PatchLevel::Iterator ip(level->begin()); ip != level->end(); ++ip ) {
      boost::shared_ptr< hier::Patch > patch = *ip;

      const boost::shared_ptr<geom::CartesianPatchGeometry > patch_geom (
         patch->getPatchGeometry(), boost::detail::dynamic_cast_tag());
      const double* dx  = patch_geom->getDx();

      boost::shared_ptr<pdat::CellData<double> >
         phase ( patch->getPatchData(phase_id), boost::detail::dynamic_cast_tag());

      assert( phase->getGhostCellWidth()[0]>0 );
      
      boost::shared_ptr<pdat::SideData<double> >
         phase_flux ( patch->getPatchData( flux_id ), boost::detail::dynamic_cast_tag());

      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast  = pbox.upper();

      FORT_COMPUTE_FLUX_ISOTROPIC(
         ifirst(0),ilast(0),
         ifirst(1),ilast(1),
#if (NDIM == 3)
         ifirst(2),ilast(2),
#endif
         dx,
         d_epsilon_phase,
         phase->getPointer(),
         phase->getGhostCellWidth()[0],
         phase_flux->getPointer(0),
         phase_flux->getPointer(1),
#if (NDIM == 3)
         phase_flux->getPointer(2),
#endif
         phase_flux->getGhostCellWidth()[0]
         );
   }
}
