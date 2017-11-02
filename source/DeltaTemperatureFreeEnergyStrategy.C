#include "DeltaTemperatureFreeEnergyStrategy.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"

#include <cassert>

using namespace SAMRAI;
using namespace std;

#define FORT_COMP_DPHIDTEMPERATURE_DELTA_TEMPERATURE computedphidtemperaturedeltatemperature_

extern "C" {

   void FORT_COMP_DPHIDTEMPERATURE_DELTA_TEMPERATURE(
      const int& ifirst0, const int& ilast0,
      const int& ifirst1, const int& ilast1,
#if (NDIM == 3)
      const int& ifirst2, const int& ilast2,
#endif
      const double*, const int&,
      const double*, const int&,
      const double&,const double&,const double&,
      double* rhs, const int&
      );
}

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
//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::applydPhidTBlock(const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int phase_id,
   const int rhs_id,
   const double phase_mobility)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for ( int ln = 0; ln <= maxln; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );

      for ( hier::PatchLevel::Iterator p(level->begin()); p!=level->end(); ++p ) {
         boost::shared_ptr<hier::Patch > patch = *p;

         const hier::Box& box = patch->getBox();
         const hier::Index& ifirst = box.lower();
         const hier::Index& ilast = box.upper();

         boost::shared_ptr< pdat::CellData<double> > temp(
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( temperature_id ) ) );
         boost::shared_ptr< pdat::CellData<double> > phase(
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( phase_id ) ) );
         boost::shared_ptr< pdat::CellData<double> > rhs(
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( rhs_id ) ) );

         FORT_COMP_DPHIDTEMPERATURE_DELTA_TEMPERATURE(
            ifirst(0),ilast(0),
            ifirst(1),ilast(1),
#if (NDIM == 3)
            ifirst(2),ilast(2),
#endif
            phase->getPointer(), phase->getGhostCellWidth()[0],
            temp->getPointer(), temp->getGhostCellWidth()[0],
            d_Tm, d_L, phase_mobility,
            rhs->getPointer(), rhs->getGhostCellWidth()[0]);
      }
   }
}

