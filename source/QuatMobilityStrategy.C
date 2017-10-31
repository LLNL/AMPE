#include "QuatMobilityStrategy.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

void QuatMobilityStrategy::computePhaseTemperatureMobility(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& phase_mobility_id, int& cp_id,
   int& mobility_id)
{
   assert( phase_mobility_id>=0 );
   assert( cp_id>=0 );
   assert( mobility_id>=0 );

   math::HierarchyCellDataOpsReal<double> ops(hierarchy);

   ops.divide(mobility_id,phase_mobility_id,cp_id,false);
}

