/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceVariable_C
#define included_pdat_OuterfaceVariable_C

#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OuterfaceDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for face variable objects
 *
 *************************************************************************
 */

template<class TYPE>
OuterfaceVariable<TYPE>::OuterfaceVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  boost::make_shared<OuterfaceDataFactory<TYPE> >(dim, depth))
{
}

template<class TYPE>
OuterfaceVariable<TYPE>::~OuterfaceVariable()
{
}

template<class TYPE>
int OuterfaceVariable<TYPE>::getDepth() const
{
   boost::shared_ptr<OuterfaceDataFactory<TYPE> > factory(
      BOOST_CAST<OuterfaceDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
