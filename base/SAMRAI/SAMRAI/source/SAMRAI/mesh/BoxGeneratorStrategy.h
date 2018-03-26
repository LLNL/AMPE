/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for box generation routines.
 *
 ************************************************************************/

#ifndef included_mesh_BoxGeneratorStrategy
#define included_mesh_BoxGeneratorStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchLevel.h"

namespace SAMRAI {
namespace mesh {

/**
 * Class BoxGeneratorStrategy is an abstract base class that defines
 * a Strategy pattern interface for operations to build boxes that cover a
 * collection of tagged cells on a single AMR patch hierarchy level.
 *
 * @see hier::PatchLevel
 */

class BoxGeneratorStrategy
{
public:
   /**
    * Default constructor.
    */
   BoxGeneratorStrategy();

   /**
    * Virtual destructor.
    */
   virtual ~BoxGeneratorStrategy();

   /*!
    * @brief Cluster tags using the DLBG interfaces.
    */
   virtual void
   findBoxesContainingTags(
      hier::BoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const boost::shared_ptr<hier::PatchLevel>& tag_level,
      const int tag_data_index,
      const int tag_val,
      const hier::Box& bound_box,
      const hier::IntVector& min_box,
      const double efficiency_tol,
      const double combine_tol,
      const hier::IntVector& max_gcw,
      const hier::BlockId& block_id,
      const hier::LocalId& first_local_id) const = 0;

private:
   // The following are not implemented:
   BoxGeneratorStrategy(
      const BoxGeneratorStrategy&);
   void
   operator = (
      const BoxGeneratorStrategy&);

};

}
}
#endif
