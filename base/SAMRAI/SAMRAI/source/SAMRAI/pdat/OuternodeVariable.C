/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Variable<DIM> class for defining outernode centered variables
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeVariable_C
#define included_pdat_OuternodeVariable_C

#include "SAMRAI/pdat/OuternodeVariable.h"
#include "SAMRAI/pdat/OuternodeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/make_shared.hpp>

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
      getPatchDataFactory());
   TBOX_ASSERT(factory);
   return factory->getDepth();
}

/*
 *************************************************************************
 *
 * These are private and should not be used.  They are defined here
 * because some template instantiation methods fail if some member
 * functions are left undefined.
 *
 *************************************************************************
 */

template<class TYPE>
OuternodeVariable<TYPE>::OuternodeVariable(
   const OuternodeVariable<TYPE>& foo):
   hier::Variable(NULL,
                  boost::shared_ptr<hier::PatchDataFactory>())
{
   NULL_USE(foo);
}

template<class TYPE>
void OuternodeVariable<TYPE>::operator = (
   const OuternodeVariable<TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
