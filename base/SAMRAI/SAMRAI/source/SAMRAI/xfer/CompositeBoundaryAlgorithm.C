/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:  Manages data on stencils at coarse-fine boundaries
 *
 ************************************************************************/
#include "SAMRAI/xfer/CompositeBoundaryAlgorithm.h"

#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/HierarchyNeighbors.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/xfer/PatchInteriorVariableFillPattern.h"
#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */

CompositeBoundaryAlgorithm::CompositeBoundaryAlgorithm(
   boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int stencil_width)
: d_hierarchy(hierarchy),
  d_stencil_width(stencil_width),
  d_schedule_created(false)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(stencil_width >= 1);
}

/*
 *************************************************************************
 *
 * Destructor explicitly deallocates PatchData
 *
 *************************************************************************
 */

CompositeBoundaryAlgorithm::~CompositeBoundaryAlgorithm()
{
}

/*
 *************************************************************************
 *
 * Creates a stencil at the coarse-fine boundary of the given level
 *
 *************************************************************************
 */

boost::shared_ptr<CompositeBoundarySchedule>
CompositeBoundaryAlgorithm::createSchedule(int level_num)
{
   d_schedule_created = true;

   return boost::make_shared<CompositeBoundarySchedule>(
             d_hierarchy,
             level_num,
             d_stencil_width,
             d_data_ids);

}

void CompositeBoundaryAlgorithm::addDataId(int data_id)
{
   TBOX_ASSERT(data_id >= 0);

   if (d_schedule_created) {
      TBOX_ERROR("CompositeBoundaryAlgorithm::addDataId error:  No data ids may be added after creating a CompositeBoundarySchedule." << std::endl);
   }

   d_data_ids.insert(data_id);
}




}
}
