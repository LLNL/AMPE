/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   tbox
 *
 ************************************************************************/

#ifndef included_hier_ProcessorMapping
#define included_hier_ProcessorMapping

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

/**
 * Class ProcessorMapping represents the processor assignments of
 * patches to processors.  It makes sure that all processor assignments
 * are in the range from 0 through NODES-1 and answers whether a particular
 * assignment is local to the processor.
 */

class ProcessorMapping
{
public:
   /**
    * Create a default processor mapping array with 0 elements.  Before
    * the mapping can be used, its size should be set using the function
    * setMappingSize() and each element of the mapping should be set
    * by setProcessorAssignment().
    */
   ProcessorMapping();

   /**
    * Create a processor mapping array with enough space for n elements.
    * All elements of the mapping are initialized to processor zero, but
    * they should be set by setProcessorAssignment() later.
    */
   explicit ProcessorMapping(
      const int n);

   /**
    * Create a new processor mapping and copy the processor assignments
    * from the argument.
    */
   ProcessorMapping(
      const ProcessorMapping& mapping);

   /**
    * Create a new processor mapping and get processor assignments
    * from the the tbox::Array<int> argument.
    */
   explicit ProcessorMapping(
      const tbox::Array<int>& mapping);

   /**
    * The destructor simply releases the storage for the mapping.
    */
   ~ProcessorMapping();

   /**
    * Resize the mapping so that it has n elements.  Before it can be
    * used, each element should be set using setProcessorAssignment().
    */
   void
   setMappingSize(
      const int n);

   /**
    * Sets the number of mapped_boxes to n.
    * IMPORTANT NOTE: This method should only be used for
    * testing purposes.  Under normal circumstances, the number of
    * mapped_boxes is set by a call to tbox::SAMRAI_MPI::getNodes() and should NOT
    * be changed.
    */
   void
   setNumberNodes(
      const int n)
   {
      d_nodes = n;
   }

   /**
    * Return the processor assignment for the specified patch index.
    */
   int
   getProcessorAssignment(
      const int i) const
   {
      TBOX_ASSERT((i >= 0) && (i < d_mapping.getSize()));
      return d_mapping[i];
   }

   /**
    * Set the processor assignment (second argument) for the specified
    * patch index (first argument).
    */
   void
   setProcessorAssignment(
      const int i,
      const int p)
   {
      TBOX_ASSERT((i >= 0) && (i < d_mapping.getSize()));
      TBOX_ASSERT((p >= 0) && (p < d_nodes));
      d_mapping[i] = p % d_nodes;
   }

   /**
    * Return an tbox::Array<int> of the processor mappings.
    */
   tbox::Array<int>
   getProcessorMapping() const
   {
      return d_mapping;
   }

   /**
    * Sets the processor mappings from an tbox::Array<int>.  Remaps the
    * processors so that patches are not accidentally mapped to
    * non-existent mapped_boxes.
    */
   void
   setProcessorMapping(
      const tbox::Array<int>& mapping);

   /**
    * Return the number of local indices (that is, those indices mapped to
    * the local processor).
    */
   int
   getNumberOfLocalIndices() const
   {
      computeLocalIndices();
      return d_local_id_count;
   }

   /**
    * Return an array containing the local indices (that is,
    * those indices mapped to the local processor).
    */
   const tbox::Array<int>&
   getLocalIndices() const
   {
      computeLocalIndices();
      return d_local_indices;
   }

   /**
    * Return the total number of indices in the mapping array.
    */
   int
   getSizeOfMappingArray() const
   {
      return d_mapping.getSize();
   }

   /**
    * Check whether the specified index is a local index (that is, mapped
    * to the local processor).
    */
   bool
   isMappingLocal(
      const int i) const
   {
      TBOX_ASSERT((i >= 0) && (i < d_mapping.getSize()));
      return d_mapping[i] == d_my_rank;
   }

private:
   /**
    * Fills in the array d_local_indices, and sets d_local_id_count.
    */
   void
   computeLocalIndices() const;

   void
   operator = (
      const ProcessorMapping&);                 // not implemented

   int d_my_rank;
   int d_nodes;
   tbox::Array<int> d_mapping;
   mutable int d_local_id_count;
   mutable tbox::Array<int> d_local_indices;
};

}
}

#endif
