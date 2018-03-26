/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   TreeLoadBalancer test.
 *
 ************************************************************************/
#ifndef included_DerivedVisOwnerData
#define included_DerivedVisOwnerData

#include <string>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"

#include "SAMRAI/tbox/Database.h"

#include <boost/shared_ptr.hpp>

/*
 * SAMRAI classes
 */
#include "SAMRAI/appu/VisDerivedDataStrategy.h"

using namespace SAMRAI;

/*!
 * @brief Write owner rank using VisDerivedDataStrategy.
 */
class DerivedVisOwnerData:
   public appu::VisDerivedDataStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   DerivedVisOwnerData();

   ~DerivedVisOwnerData();

   //@{ @name SAMRAI::appu::VisDerivedDataStrategy virtuals

   virtual bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_id) const;

   //@}

private:
};

#endif  // included_DerivedVisOwnerData
