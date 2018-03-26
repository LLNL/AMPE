/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Factory class for patch data objects that live on a patch
 *
 ************************************************************************/

#ifndef included_hier_PatchDescriptor_C
#define included_hier_PatchDescriptor_C

#include "SAMRAI/hier/PatchDescriptor.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include <typeinfo>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

const int PatchDescriptor::INDEX_UNDEFINED = -1;

/*
 *************************************************************************
 *
 * The constructor sets the max number of registered components to zero
 * and allocates the factory and name arrays to the fixed length set
 * by the SAMRAIManager utility.  The free list of indices
 * is initialized to the full set of potentially used indices.
 *
 * The destructor clears the free index list and implicitly
 * deallocates the arrays of name strings and factory pointers.
 *
 *************************************************************************
 */

PatchDescriptor::PatchDescriptor():
   d_min_gcw(tbox::Dimension::MAXIMUM_DIMENSION_VALUE)
{
   const int max_num_patch_data_components_allowed =
      tbox::SAMRAIManager::getMaxNumberPatchDataEntries();
   d_max_number_registered_components = 0;
   d_names.resizeArray(max_num_patch_data_components_allowed);
   d_factories.resizeArray(max_num_patch_data_components_allowed);
   for (int i = 0; i < max_num_patch_data_components_allowed; i++) {
      d_free_indices.push_back(i);
   }
   for (unsigned short d = 0; d < d_min_gcw.size(); ++d) {
      d_min_gcw[d] = IntVector::getZero(tbox::Dimension(static_cast<unsigned short>(d + 1)));
   }
}

PatchDescriptor::~PatchDescriptor()
{
   d_free_indices.clear();
}

/*
 *************************************************************************
 *
 * Add the new factory to the list of patch data factories and assign
 * it an integer index identifier.  Use a free list item if possible.
 *
 *************************************************************************
 */

int
PatchDescriptor::definePatchDataComponent(
   const std::string& name,
   const boost::shared_ptr<PatchDataFactory>& factory)
{
   TBOX_ASSERT(!name.empty());
   TBOX_ASSERT(factory);

   int ret_index = INDEX_UNDEFINED;
   if (d_free_indices.empty()) {
      TBOX_ERROR(
         "PatchDescriptor::definePatchDataComponent error...\n"
         << "No available patch data component indices left.\n"
         << "Application must be restarted and size must be increased.\n"
         << "See tbox::SAMRAIManager utility for more information."
         << std::endl);
   } else {
      ret_index = d_free_indices.front();
      d_free_indices.pop_front();
      if (d_max_number_registered_components < ret_index + 1) {
         d_max_number_registered_components = ret_index + 1;
      }
      d_factories[ret_index] = factory;
      d_names[ret_index] = name;
   }
   return ret_index;
}

/*
 *************************************************************************
 *
 * Remove the specified patch data factory index and place the index on
 * the list of free indices.
 *
 *************************************************************************
 */

void
PatchDescriptor::removePatchDataComponent(
   const int id)
{
   if ((id >= 0) && (id < d_max_number_registered_components)) {
      if (!d_names[id].empty()) {
         d_names[id] = std::string();
      }
      if (d_factories[id]) {
         d_factories[id].reset();
         d_free_indices.push_front(id);
      }
   }
}

/*
 *************************************************************************
 *
 * Look up the factory by name; if no matching factory exists, then a
 * pointer to null is returned.  The first matching factory is returned.
 *
 *************************************************************************
 */

boost::shared_ptr<PatchDataFactory>
PatchDescriptor::getPatchDataFactory(
   const std::string& name) const
{
   boost::shared_ptr<PatchDataFactory> factory;
   const int id = mapNameToIndex(name);
   if (id >= 0) {
      factory = d_factories[id];
   }
   return factory;
}

/*
 *************************************************************************
 *
 * Search the factory list for a match and return the associated
 * factory.  If no match exists, return a negative identifier.
 *
 *************************************************************************
 */

int
PatchDescriptor::mapNameToIndex(
   const std::string& name) const
{
   int ret_index = INDEX_UNDEFINED;
   int id = 0;
   while ((ret_index == INDEX_UNDEFINED) &&
          (id < d_max_number_registered_components)) {
      if (name == d_names[id]) {
         ret_index = id;
      }
      id++;
   }
   return ret_index;
}

/*
 *************************************************************************
 *
 * Print index, name, and factory data for the patch descriptor.
 *
 *************************************************************************
 */

void
PatchDescriptor::printClassData(
   std::ostream& stream) const
{
   stream << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          << std::endl;
   stream << "Printing PatchDescriptor state ..." << std::endl;
   stream << "this = " << (PatchDescriptor *)this << std::endl;
   stream << "d_max_number_registered_components = "
          << d_max_number_registered_components << std::endl;
   stream << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          << std::endl;
   for (int i = 0; i < d_max_number_registered_components; i++) {
      stream << "Patch Data Index=" << i << std::endl;
      if (d_factories[i]) {
         stream << "   Patch Data Factory Name = "
                << d_names[i] << std::endl;
         stream << "   Patch Data Factory = "
                << typeid(*d_factories[i]).name() << std::endl;
      } else {
         stream << "   Patch Data Factory = NULL" << std::endl;
      }
   }
   stream << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          << std::endl;
}

/*
 *************************************************************************
 * Return the maximum ghost cell width across all factories and the
 * user-specified minimum value.
 *************************************************************************
 */

IntVector
PatchDescriptor::getMaxGhostWidth(
   const tbox::Dimension& dim) const
{
   IntVector max_gcw(d_min_gcw[dim.getValue() - 1]);
   for (int i = 0; i < d_max_number_registered_components; i++) {
      if (d_factories[i] && (d_factories[i]->getDim() == dim)) {
         max_gcw.max(d_factories[i]->getGhostCellWidth());
      }
   }
   return max_gcw;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
