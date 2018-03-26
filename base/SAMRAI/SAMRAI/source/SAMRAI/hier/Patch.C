/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Patch container class for patch data objects
 *
 ************************************************************************/

#ifndef included_hier_Patch_C
#define included_hier_Patch_C

#include "SAMRAI/hier/Patch.h"

#include <typeinfo>
#include <string>

namespace SAMRAI {
namespace hier {

const int Patch::HIER_PATCH_VERSION = 2;

/*
 *************************************************************************
 *
 * Allocate a patch container but do not instantiate any components.
 *
 *************************************************************************
 */

Patch::Patch(
   const Box& mapped_box,
   const boost::shared_ptr<PatchDescriptor>& descriptor):
   d_mapped_box(mapped_box),
   d_descriptor(descriptor),
   d_patch_data(d_descriptor->getMaxNumberRegisteredComponents()),
   d_patch_level_number(-1),
   d_patch_in_hierarchy(false)
{
   TBOX_ASSERT(mapped_box.getLocalId() >= 0);
}

/*
 *************************************************************************
 *
 * The virtual destructor does nothing; all memory deallocation is
 * managed automatically by the pointer and array classes.
 *
 *************************************************************************
 */

Patch::~Patch()
{
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space required to allocate the
 * specified component(s).  This information can then be used by a
 * fixed-size memory allocator.
 *
 *************************************************************************
 */

size_t
Patch::getSizeOfPatchData(
   const ComponentSelector& components) const
{
   size_t size = 0;
   const int max_set_component = components.getMaxIndex();

   for (int i = 0; i < max_set_component && components.isSet(i); i++) {
      size += d_descriptor->getPatchDataFactory(i)->getSizeOfMemory(
            d_mapped_box);
   }

   return size;
}

/*
 *************************************************************************
 *
 * Allocate the specified patch data object(s) on the patch.
 *
 *************************************************************************
 */

void
Patch::allocatePatchData(
   const int id,
   const double time)
{
   const int ncomponents = d_descriptor->getMaxNumberRegisteredComponents();

   TBOX_ASSERT((id >= 0) && (id < ncomponents));

   if (ncomponents > d_patch_data.getSize()) {
      d_patch_data.resizeArray(ncomponents);
   }

   if (!checkAllocated(id)) {
      d_patch_data[id] =
         d_descriptor->getPatchDataFactory(id)->allocate(*this);
   }
   d_patch_data[id]->setTime(time);
}

void
Patch::allocatePatchData(
   const ComponentSelector& components,
   const double time)
{
   const int ncomponents = d_descriptor->getMaxNumberRegisteredComponents();
   if (ncomponents > d_patch_data.getSize()) {
      d_patch_data.resizeArray(ncomponents);
   }

   for (int i = 0; i < ncomponents; i++) {
      if (components.isSet(i)) {
         if (!checkAllocated(i)) {
            d_patch_data[i] =
               d_descriptor->getPatchDataFactory(i)->allocate(*this);
         }
         d_patch_data[i]->setTime(time);
      }
   }
}

/*
 *************************************************************************
 *
 * Deallocate (or set to null) the specified component(s).
 *
 *************************************************************************
 */

void
Patch::deallocatePatchData(
   const ComponentSelector& components)
{
   const int ncomponents = d_patch_data.getSize();
   for (int i = 0; i < ncomponents; i++) {
      if (components.isSet(i)) {
         d_patch_data[i].reset();
      }
   }
}

/*
 *************************************************************************
 *
 * Set the time stamp for the specified components in the patch.
 *
 *************************************************************************
 */

void
Patch::setTime(
   const double timestamp,
   const ComponentSelector& components)
{
   const int ncomponents = d_patch_data.getSize();
   for (int i = 0; i < ncomponents; i++) {
      if (components.isSet(i) && d_patch_data[i]) {
         d_patch_data[i]->setTime(timestamp);
      }
   }
}

void
Patch::setTime(
   const double timestamp)
{
   const int ncomponents = d_patch_data.getSize();
   for (int i = 0; i < ncomponents; i++) {
      if (d_patch_data[i]) {
         d_patch_data[i]->setTime(timestamp);
      }
   }
}

/*
 *************************************************************************
 *
 * Checks that class and restart file version numbers are equal.  If so,
 * reads in data from database and have each patch_data item read
 * itself in from the database
 *
 *************************************************************************
 */

void
Patch::getFromDatabase(
   const boost::shared_ptr<tbox::Database>& database,
   const ComponentSelector& component_selector)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("HIER_PATCH_VERSION");
   if (ver != HIER_PATCH_VERSION) {
      TBOX_ERROR("Patch::getFromDatabase() error...\n"
         << "   Restart file version different than class version" << std::endl);
   }

   Box box(database->getDatabaseBox("d_box"));
   const LocalId patch_local_id(database->getInteger("d_patch_local_id"));
   int patch_owner = database->getInteger("d_patch_owner");
   int block_id = database->getInteger("d_block_id");
   box.setBlockId(BlockId(block_id));
   d_mapped_box.initialize(box,
      patch_local_id,
      patch_owner);

   d_patch_level_number = database->getInteger("d_patch_level_number");
   d_patch_in_hierarchy = database->getBool("d_patch_in_hierarchy");

   d_patch_data.resizeArray(d_descriptor->getMaxNumberRegisteredComponents());

   int namelist_count = database->getInteger("patch_data_namelist_count");
   tbox::Array<std::string> patch_data_namelist;
   if (namelist_count) {
      patch_data_namelist = database->getStringArray("patch_data_namelist");
   }

   ComponentSelector local_selector(component_selector);

   for (int i = 0; i < patch_data_namelist.getSize(); i++) {
      std::string patch_data_name;
      int patch_data_index;

      patch_data_name = patch_data_namelist[i];

      if (!database->isDatabase(patch_data_name)) {
         TBOX_ERROR("Patch::getFromDatabase() error...\n"
            << "   patch data" << patch_data_name
            << " not found in database" << std::endl);
      }
      boost::shared_ptr<tbox::Database> patch_data_database(
         database->getDatabase(patch_data_name));

      patch_data_index = d_descriptor->mapNameToIndex(patch_data_name);

      if ((patch_data_index >= 0) &&
          (local_selector.isSet(patch_data_index))) {
         boost::shared_ptr<PatchDataFactory> patch_data_factory(
            d_descriptor->getPatchDataFactory(patch_data_index));
         d_patch_data[patch_data_index] = patch_data_factory->allocate(*this);
         d_patch_data[patch_data_index]->getFromDatabase(patch_data_database);

         local_selector.clrFlag(patch_data_index);
      }
   }

   if (local_selector.any()) {
      TBOX_WARNING("Patch::getFromDatabase() warning...\n"
         << "   Some requested patch data components not "
         << "found in database" << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Write out the class version number to database.  Then,
 * writes out data to database and have each patch_data item write
 * itself out to the database.  The following data
 * members are written out: d_mapped_box, d_patch_number,
 * d_patch_level_number,
 * d_patch_in_hierarchy, d_patch_data[].
 * The database key for all data members is identical to the
 * name of the data member except for the d_patch_data.  These have
 * keys of the form "variable##context" which is the form that they
 * are stored by the patch descriptor.  In addition a list of the
 * patch_data names ("patch_data_namelist") and the number of patch data
 * items saved ("namelist_count") are also written to the database.
 * The patchdata_write_table determines which patchdata are written to
 * the database.
 *
 *************************************************************************
 */
void
Patch::putUnregisteredToDatabase(
   const boost::shared_ptr<tbox::Database>& database,
   const ComponentSelector& patchdata_write_table) const
{
   TBOX_ASSERT(database);

   int i;

   database->putInteger("HIER_PATCH_VERSION", HIER_PATCH_VERSION);
   database->putDatabaseBox("d_box", d_mapped_box);
   database->putInteger("d_patch_local_id",
      d_mapped_box.getLocalId().getValue());
   database->putInteger("d_patch_owner",
      d_mapped_box.getOwnerRank());
   database->putInteger("d_block_id",
      d_mapped_box.getBlockId().getBlockValue());
   database->putInteger("d_patch_level_number", d_patch_level_number);
   database->putBool("d_patch_in_hierarchy", d_patch_in_hierarchy);

   int namelist_count = 0;
   for (i = 0; i < d_patch_data.getSize(); i++) {
      if (patchdata_write_table.isSet(i) && checkAllocated(i)) {
         namelist_count++;
      }
   }

   std::string patch_data_name;
   tbox::Array<std::string> patch_data_namelist(namelist_count);
   namelist_count = 0;
   for (i = 0; i < d_patch_data.getSize(); i++) {
      if (patchdata_write_table.isSet(i) && checkAllocated(i)) {
         patch_data_namelist[namelist_count++] =
            patch_data_name = d_descriptor->mapIndexToName(i);
         boost::shared_ptr<tbox::Database> patch_data_database(
            database->putDatabase(patch_data_name));
         (d_patch_data[i])->putUnregisteredToDatabase(patch_data_database);
      }
   }

   database->putInteger("patch_data_namelist_count", namelist_count);
   if (namelist_count > 0) {
      database->putStringArray("patch_data_namelist", patch_data_namelist);
   }
}

/*
 *************************************************************************
 *
 * Print information about the patch.
 *
 *************************************************************************
 */

int
Patch::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int depth) const
{
   NULL_USE(depth);

   const tbox::Dimension& dim(d_mapped_box.getDim());

   os << border
      << d_mapped_box
      << "\tdims: " << d_mapped_box.numberCells(0)
   ;
   for (int i = 1; i < dim.getValue(); ++i) {
      os << " X " << d_mapped_box.numberCells(i);
   }
   os << "\tsize: " << d_mapped_box.size()
      << "\n";
   return 0;
}

std::ostream&
operator << (
   std::ostream& s,
   const Patch& patch)
{
   s << "Patch::mapped_box = "
   << patch.d_mapped_box << std::endl << std::flush;
   s << "Patch::patch_level_number = " << patch.d_patch_level_number
   << std::endl << std::flush;
   s << "Patch::patch_in_hierarchy = " << patch.d_patch_in_hierarchy
   << std::endl << std::flush;
   s << "Patch::number_components = " << patch.d_patch_data.getSize()
   << std::endl << std::flush;
   const int ncomponents = patch.d_patch_data.getSize();
   for (int i = 0; i < ncomponents; i++) {
      s << "Component(" << i << ")=";
      if (!patch.d_patch_data[i]) {
         s << "NULL\n";
      } else {
         s << typeid(*patch.d_patch_data[i]).name()
         << " [GCW=" << patch.d_patch_data[i]->getGhostCellWidth() << "]\n";
      }
   }
   return s;
}

}
}
#endif
