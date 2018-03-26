/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Common Box operations for Box containers.
 *
 ************************************************************************/
#ifndef included_hier_BoxContainerUtils
#define included_hier_BoxContainerUtils

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Connector.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Utilities for performing simple common tasks on a container
 * of Boxes.
 *
 * TODO: Arguments should be re-ordered to the SAMRAI standard, output
 * before input.
 *
 * TODO: There are some very similar methods in this class.  In many cases,
 * one version supports input and output being the same object, but the other
 * does not.  For uniformity, all these methods should support input and
 * output containers being the same object.  It's simple to implement.
 */
class BoxContainerUtils
{

public:
   //@{

   //! @name I/O operations for containers that lack built-in versions.

   /*!
    * @brief Print a vector of Boxes to an output stream.
    *
    * @param[in] mapped_boxes
    *
    * @param[in] output_stream
    *
    * @param[in] border
    *
    * @param[in] detail_depth
    */
   static void
   recursivePrintBoxVector(
      const std::vector<Box>& mapped_boxes,
      std::ostream& output_stream = tbox::plog,
      const std::string& border = std::string(),
      int detail_depth = 0);

   //@}

private:
   // Disabled constructor.  No need for objects of this class.
   BoxContainerUtils();

};

}
}

#endif  // included_hier_BoxContainerUtils
