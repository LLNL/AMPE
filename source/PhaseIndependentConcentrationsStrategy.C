#include "PhaseIndependentConcentrationsStrategy.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"

void PhaseIndependentConcentrationsStrategy::computePhaseConcentrationsOnPatch(
   boost::shared_ptr< pdat::CellData<double> > cd_temperature,
   boost::shared_ptr< pdat::CellData<double> > cd_phi,
   boost::shared_ptr< pdat::CellData<double> > cd_eta,
   boost::shared_ptr< pdat::CellData<double> > cd_concentration,
   boost::shared_ptr< pdat::CellData<double> > cd_c_l,
   boost::shared_ptr< pdat::CellData<double> > cd_c_a,
   boost::shared_ptr< pdat::CellData<double> > cd_c_b,
   boost::shared_ptr<hier::Patch > patch )
{
   (void)cd_temperature;
   (void)cd_phi;
   (void)cd_eta;
   (void)cd_c_b;

   assert( cd_concentration );
   assert( cd_c_l );
   assert( cd_c_a );

   const hier::Box& pbox = patch->getBox();

   math::PatchCellDataOpsReal<double> cellops;

   cellops.copyData(cd_c_l,cd_concentration,pbox);
   cellops.copyData(cd_c_a,cd_concentration,pbox); 
}

