/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch level objects
 *
 ************************************************************************/

#ifndef included_hier_PatchLevelFactory
#define included_hier_PatchLevelFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchFactory.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Database.h"

#include <boost/shared_ptr.hpp>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Factory used to create new patch levels.
 *
 * New types of patch level objects can be introduced into SAMRAI by deriving
 * from PatchLevelFactory and re-defining allocate.
 *
 * @see hier::PatchLevel
 */
class PatchLevelFactory
{
public:
   /*!
    * @brief Construct a patch level factory object.
    */
   PatchLevelFactory();

   /*!
    * @brief Virtual destructor for patch level factory objects.
    */
   virtual ~PatchLevelFactory();

   /*!
    * @brief Allocate a patch level with the specified boxes and processor mappings.
    *
    * Redefine this function to change the method for creating patch levels.
    *
    * @return A boost::shared_ptr to the newly created PatchLevel.
    *
    * @param[in]  mapped_box_level
    * @param[in]  grid_geometry
    * @param[in]  descriptor
    * @param[in]  factory @b Default: a boost::shared_ptr to the standard
    *             PatchFactory
    */
   virtual boost::shared_ptr<PatchLevel>
   allocate(
      const BoxLevel& mapped_box_level,
      const boost::shared_ptr<BaseGridGeometry>& grid_geometry,
      const boost::shared_ptr<PatchDescriptor>& descriptor,
      const boost::shared_ptr<PatchFactory>& factory =
         boost::shared_ptr<PatchFactory>()) const;

   /*!
    * @brief Allocate a patch level using the data from the database to
    * initialize it.
    *
    * The component_selector argument is used to specify which patch data
    * components to allocate and read in from the database.
    * @note
    * If desired, pass a ComponentSelector with all bits set to false to
    * indicate that no patch data components are read/allocated.
    *
    * Redefine this function to change the method for creating
    * patch levels from a database.
    *
    * @return A boost::shared_ptr to the newly created PatchLevel.
    *
    * @param[in]  database
    * @param[in]  grid_geometry
    * @param[in]  descriptor
    * @param[in]  component_selector
    * @param[in]  factory @b Default: a boost::shared_ptr to the standard
    *             PatchFactory
    * @param[in]  defer_boundary_box_creation @b Default: false
    */
   virtual boost::shared_ptr<PatchLevel>
   allocate(
      const boost::shared_ptr<tbox::Database>& database,
      const boost::shared_ptr<BaseGridGeometry>& grid_geometry,
      const boost::shared_ptr<PatchDescriptor>& descriptor,
      const ComponentSelector& component_selector,
      const boost::shared_ptr<PatchFactory>& factory =
         boost::shared_ptr<PatchFactory>(),
      const bool defer_boundary_box_creation = false) const;

private:
   /*
    * Copy constructor and assignment are not implemented.
    */
   PatchLevelFactory(
      const PatchLevelFactory&);
   void
   operator = (
      const PatchLevelFactory&);

};

}
}

#endif
