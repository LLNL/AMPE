/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   VisItDerivedData class to write patch owner value.
 *
 ************************************************************************/
#ifndef included_VisItDerivedData
#define included_VisItDerivedData

/*
 * SAMRAI classes
 */
#include "SAMRAI/appu/VisDerivedDataStrategy.h"

using namespace SAMRAI;

/*!
 * @brief Class to implemement patch owners as VisIt derived data.
 */
class VisItDerivedData:
   public appu::VisDerivedDataStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   VisItDerivedData();

   ~VisItDerivedData();

   //@{ @name SAMRAI::appu::VisDerivedDataStrategy virtuals

   virtual bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_id) const;

   //@}

};

#endif  // included_VisItDerivedData
