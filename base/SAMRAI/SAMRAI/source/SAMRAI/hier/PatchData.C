/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for patch data objects
 *
 ************************************************************************/

#ifndef included_hier_PatchData_C
#define included_hier_PatchData_C

#include "SAMRAI/hier/PatchData.h"

namespace SAMRAI {
namespace hier {

const int PatchData::HIER_PATCH_DATA_VERSION = 2;

PatchData::PatchData(
   const Box& domain,
   const IntVector& ghosts):
   d_box(domain),
   d_ghost_box(domain),
   d_ghosts(ghosts),
   d_timestamp(0.0)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(domain, ghosts);

   d_ghost_box.grow(ghosts);
}

PatchData::~PatchData()
{
}

/*
 *************************************************************************
 *
 * Checks that clas and restart file version number are same.  If so,
 * reads in data members common to all patch data and then invoke
 * getSpecializedFromDatabase() to read in data particular to the
 * specific derived class.
 *
 *************************************************************************
 */

void
PatchData::getFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("HIER_PATCH_DATA_VERSION");
   if (ver != HIER_PATCH_DATA_VERSION) {
      TBOX_ERROR("PatchData::getFromDatabase() error...\n"
         << "  Restart file version different than class version" << std::endl);
   }

   d_box = database->getDatabaseBox("d_box");
   d_ghost_box = database->getDatabaseBox("d_ghost_box");
   database->getIntegerArray("d_ghosts", &d_ghosts[0], d_ghosts.getDim().getValue());
   d_timestamp = database->getDouble("d_timestamp");

   getSpecializedFromDatabase(database);
}

/*
 *************************************************************************
 *
 * Write out data members common to all patch data and then invoke
 * putSpecializedToDatabase() to write out data particular to the
 * specific derived class.
 *
 *************************************************************************
 */

void
PatchData::putUnregisteredToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{
   TBOX_ASSERT(database);

   database->putInteger("HIER_PATCH_DATA_VERSION", HIER_PATCH_DATA_VERSION);
   database->putDatabaseBox("d_box", d_box);
   database->putDatabaseBox("d_ghost_box", d_ghost_box);
   database->putDouble("d_timestamp", d_timestamp);
   database->putIntegerArray("d_ghosts", &d_ghosts[0], d_ghosts.getDim().getValue());

   putSpecializedToDatabase(database);
}

}
}
#endif
