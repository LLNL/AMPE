/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Variable class for defining outernode centered variables
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeVariable_C
#define included_pdat_OuternodeVariable_C

#include "SAMRAI/pdat/OuternodeVariable.h"
#include "SAMRAI/pdat/OuternodeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for side variable objects
 *
 *************************************************************************
 */

template<class TYPE>
OuternodeVariable<TYPE>::OuternodeVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  boost::make_shared<OuternodeDataFactory<TYPE> >(dim, depth))
{
}

template<class TYPE>
OuternodeVariable<TYPE>::~OuternodeVariable()
{
}

template<class TYPE>
int OuternodeVariable<TYPE>::getDepth() const
{
   boost::shared_ptr<OuternodeDataFactory<TYPE> > factory(
      BOOST_CAST<OuternodeDataFactory<TYPE>, hier::PatchDataFactory>(
         getPatchDataFactory()));
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

}
}
#endif
