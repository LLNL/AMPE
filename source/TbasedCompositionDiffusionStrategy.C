#include "TbasedCompositionDiffusionStrategy.h"
#include "toolsSAMRAI.h"

#include "ConcFort.h"
#include "PhysicalConstants.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

TbasedCompositionDiffusionStrategy::TbasedCompositionDiffusionStrategy(
   const int diffusion_l_id,
   const int diffusion_a_id,
   const double D_liquid, const double Q0_liquid,
   const double D_solid_A, const double Q0_solid_A,
   const string& phase_interp_func_type,
   const string& avg_func_type):
      d_diffusion_l_id(diffusion_l_id),
      d_diffusion_a_id(diffusion_a_id),
      d_D_liquid(D_liquid),
      d_D_solid_A(D_solid_A),
      d_phase_interp_func_type(phase_interp_func_type),
      d_avg_func_type(avg_func_type)
{
}

void TbasedCompositionDiffusionStrategy::setDiffCoeff(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int phase_id,
   const int eta_id)
{
   //tbox::pout<<"TbasedCompositionDiffusionStrategy::setDiffCoeff()"<<endl;
   assert( temperature_id >= 0 );
   assert( phase_id >= 0 );
   assert( d_diffusion_l_id >= 0 );
   assert( d_diffusion_a_id >= 0 );

   const int maxl = hierarchy->getNumberOfLevels();

   for ( int amr_level = 0; amr_level < maxl; amr_level++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( amr_level );

      for ( hier::PatchLevel::Iterator p(level->begin()); 
            p!=level->end(); ++p ) {
         boost::shared_ptr<hier::Patch > patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast  = pbox.upper();

         boost::shared_ptr< pdat::CellData<double> > phi (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( phase_id) ) );

         boost::shared_ptr< pdat::CellData<double> > temperature (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( temperature_id) ) );

         boost::shared_ptr< pdat::SideData<double> > diffusionL (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData( d_diffusion_l_id ) ) );

         boost::shared_ptr< pdat::SideData<double> > diffusionA (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(
               patch->getPatchData( d_diffusion_a_id ) ) );

         assert( diffusionA->getDepth()==diffusionL->getDepth() );
         assert( diffusionA->getDepth()==1 || // binary
                 diffusionA->getDepth()==4 ); // ternary

         //compute depth 0 of diffusion variables
         FORT_CONCENTRATIONDIFFUSION_OF_T(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            phi->getPointer(), phi->getGhostCellWidth()[0],
            diffusionL->getPointer(0),
            diffusionL->getPointer(1),
#if (NDIM == 3)
            diffusionL->getPointer(2),
#endif
            diffusionA->getPointer(0),
            diffusionA->getPointer(1),
#if (NDIM == 3)
            diffusionA->getPointer(2),
#endif
            0, //assuming no ghosts for diffusion data
            temperature->getPointer(),
            temperature->getGhostCellWidth()[0],
            d_D_liquid, d_Q0_liquid,
            d_D_solid_A, d_Q0_solid_A,
            gas_constant_R_JpKpmol,
            d_phase_interp_func_type.c_str(),
            d_avg_func_type.c_str());

         //fill other diagonal value with same value for ternaries for now
         if( diffusionL->getDepth()>1 ){
            copyDepthSideData(diffusionL,3,diffusionL,0);
            copyDepthSideData(diffusionA,3,diffusionA,0);
         }
      }
   }
}


