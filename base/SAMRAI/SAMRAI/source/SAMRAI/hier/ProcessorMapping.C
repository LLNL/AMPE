/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   tbox
 *
 ************************************************************************/

#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

namespace SAMRAI {
namespace hier {

ProcessorMapping::ProcessorMapping():
   d_my_rank(-1),
   d_nodes(-1),
   d_mapping(0),
   d_local_id_count(-1)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_my_rank = mpi.getRank();
   d_nodes = mpi.getSize();
}

ProcessorMapping::ProcessorMapping(
   const int n):
   d_my_rank(-1),
   d_nodes(-1),
   d_mapping(n),
   d_local_id_count(-1)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_my_rank = mpi.getRank();
   d_nodes = mpi.getSize();
   for (int i = 0; i < n; i++) {
      d_mapping[i] = 0;
   }
}

ProcessorMapping::ProcessorMapping(
   const ProcessorMapping& mapping):
   d_my_rank(-1),
   d_nodes(-1),
   d_mapping(mapping.d_mapping.getSize()),
   d_local_id_count(-1)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_my_rank = mpi.getRank();
   d_nodes = mpi.getSize();
   const int n = d_mapping.getSize();
   for (int i = 0; i < n; i++) {
      d_mapping[i] = mapping.d_mapping[i];
   }
}

ProcessorMapping::ProcessorMapping(
   const tbox::Array<int>& mapping):
   d_my_rank(-1),
   d_nodes(-1),
   d_local_id_count(-1)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_my_rank = mpi.getRank();
   d_nodes = mpi.getSize();
   setProcessorMapping(mapping);
}

ProcessorMapping::~ProcessorMapping()
{
}

void
ProcessorMapping::setMappingSize(
   const int n)
{
   d_mapping.resizeArray(n);

   for (int i = 0; i < n; i++) {
      d_mapping[i] = 0;
   }
   d_local_id_count = -1;
}

void
ProcessorMapping::setProcessorMapping(
   const tbox::Array<int>& mapping)
{
   d_mapping.resizeArray(mapping.getSize());

   for (int i = 0; i < d_mapping.getSize(); i++) {
      //  (mapping[i] % d_nodes) keeps patches from being assigned
      //  non-existent processors.
      setProcessorAssignment(i, mapping[i] % d_nodes);
   }
   d_local_id_count = -1;
}

void
ProcessorMapping::computeLocalIndices() const
{
   if (d_local_id_count != -1) {
      return;
   }

   /*
    * first, count the number of local indices,
    * so we can set the array size.
    */
   const int n = d_mapping.getSize();
   d_local_id_count = 0;

   for (int i = 0; i < n; i++) {
      if (d_mapping[i] == d_my_rank) {
         ++d_local_id_count;
      }
   }

   /*
    * second, resize the array and fill in the data
    */
   d_local_indices.resizeArray(d_local_id_count);
   int idx = 0;
   for (int i = 0; i < n; i++) {
      if (d_mapping[i] == d_my_rank) {
         d_local_indices[idx++] = i;
      }
   }
}

}
}
