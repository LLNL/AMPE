#include "DeltaTemperatureFreeEnergyStrategy.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"

#include <cassert>

using namespace SAMRAI;
using namespace std;

DeltaTemperatureFreeEnergyStrategy::DeltaTemperatureFreeEnergyStrategy(
   const double Tm, const double latent_heat):
   d_Tm(Tm),
   d_L(latent_heat)
{
   tbox::plog<<"DeltaTemperatureFreeEnergyStrategy..."<<endl;
   tbox::plog<<"Tm="<<d_Tm<<endl;
   tbox::plog<<"L="<<d_L<<endl;

   assert( d_L==d_L );
   assert( d_Tm==d_Tm );
}

//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::addComponentRhsPhi(
   hier::Patch& patch,
   const int temperature_id,
   const int phase_id,
   const int eta_id,
   const int conc_id,
   const int f_l_id,
   const int f_a_id,
   const int f_b_id,
   const int rhs_id )
{
   assert( phase_id >= 0 );
   assert( rhs_id >= 0 );
   assert( temperature_id >= 0 );
   assert( d_Tm>0. );
   assert( d_L>0. );

   (void) conc_id;  // unused
   (void) f_l_id; // unused
   (void) f_a_id; // unused

   //tbox::pout<<"DeltaTemperatureFreeEnergyStrategy::addComponentRhsPhi()..."<<endl;

   boost::shared_ptr< pdat::CellData<double> > phase (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData(phase_id) ) );
   assert( phase );

   boost::shared_ptr< pdat::CellData<double> > temp (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData(temperature_id) ) );
   assert( temp );

   boost::shared_ptr< pdat::CellData<double> > rhs (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( rhs_id) ) );
   assert( rhs );

   assert( rhs->getGhostCellWidth() == hier::IntVector(tbox::Dimension(NDIM),0) );

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast  = pbox.upper();

   FORT_COMP_RHS_DELTA_TEMPERATURE(
      ifirst(0),ilast(0),
      ifirst(1),ilast(1),
#if (NDIM == 3)
      ifirst(2),ilast(2),
#endif
      phase->getPointer(), phase->getGhostCellWidth()[0],
      temp->getPointer(), temp->getGhostCellWidth()[0],
      d_Tm, d_L,
      rhs->getPointer(), 0 );
}

