/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch level objects
 *
 ************************************************************************/

#ifndef included_hier_PatchLevelFactory_C
#define included_hier_PatchLevelFactory_C

#include "SAMRAI/hier/PatchLevelFactory.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace hier {

PatchLevelFactory::PatchLevelFactory()
{
}

PatchLevelFactory::~PatchLevelFactory()
{
}

boost::shared_ptr<PatchLevel>
PatchLevelFactory::allocate(
   const BoxLevel& mapped_box_level,
   const boost::shared_ptr<BaseGridGeometry>& grid_geometry,
   const boost::shared_ptr<PatchDescriptor>& descriptor,
   const boost::shared_ptr<PatchFactory>& factory) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(mapped_box_level, *grid_geometry);
   boost::shared_ptr<PatchLevel> pl(
      boost::make_shared<PatchLevel>(
         mapped_box_level,
         grid_geometry,
         descriptor,
         factory));
   return pl;
}

boost::shared_ptr<PatchLevel>
PatchLevelFactory::allocate(
   const boost::shared_ptr<tbox::Database>& database,
   const boost::shared_ptr<BaseGridGeometry>& grid_geometry,
   const boost::shared_ptr<PatchDescriptor>& descriptor,
   const ComponentSelector& component_selector,
   const boost::shared_ptr<PatchFactory>& factory,
   const bool defer_boundary_box_creation) const
{
   boost::shared_ptr<PatchLevel> pl(
      boost::make_shared<PatchLevel>(
         database,
         grid_geometry,
         descriptor,
         factory,
         component_selector,
         defer_boundary_box_creation));
   return pl;
}

}
}
#endif
