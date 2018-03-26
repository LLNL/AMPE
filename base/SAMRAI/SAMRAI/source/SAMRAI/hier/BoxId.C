/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Identifier for a Box.
 *
 ************************************************************************/

#ifndef included_hier_BoxId_C
#define included_hier_BoxId_C

#include "SAMRAI/hier/BoxId.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

/*
 ************************************************************************
 * Constructors
 ************************************************************************
 */
BoxId::BoxId():
   d_global_id(),
   d_periodic_id()
{
}

BoxId::BoxId(
   const GlobalId& id,
   const PeriodicId& periodic_id):
   d_global_id(id),
   d_periodic_id(periodic_id)
{
   TBOX_ASSERT(periodic_id.isValid());
}

BoxId::BoxId(
   const LocalId& local_id,
   const int owner,
   const PeriodicId& periodic_id):
   d_global_id(local_id, owner),
   d_periodic_id(periodic_id)
{
   TBOX_ASSERT(periodic_id.isValid());
}

BoxId::BoxId(
   const BoxId& r):
   d_global_id(r.d_global_id),
   d_periodic_id(r.d_periodic_id)
{
   TBOX_ASSERT(r.d_periodic_id.isValid());
}

/*
 ************************************************************************
 * Destructor
 ************************************************************************
 */
BoxId::~BoxId()
{
}

/*
 ******************************************************************************
 * Stream-insert operator.
 ******************************************************************************
 */
std::ostream&
operator << (
   std::ostream& co,
   const BoxId& r)
{
   co << r.d_global_id.getOwnerRank()
   << '#' << r.d_global_id.getLocalId()
   << '/' << r.d_periodic_id;
   return co;
}

}
}
#endif
