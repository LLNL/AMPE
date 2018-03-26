/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Variable class for defining outeredge centered variables
 *
 ************************************************************************/

#ifndef included_pdat_OuteredgeVariable_C
#define included_pdat_OuteredgeVariable_C

#include "SAMRAI/pdat/OuteredgeVariable.h"
#include "SAMRAI/pdat/OuteredgeDataFactory.h"
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
OuteredgeVariable<TYPE>::OuteredgeVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth):
   hier::Variable(name,
                  boost::make_shared<OuteredgeDataFactory<TYPE> >(dim, depth))
{
}

template<class TYPE>
OuteredgeVariable<TYPE>::~OuteredgeVariable()
{
}

template<class TYPE>
int OuteredgeVariable<TYPE>::getDepth() const
{
   boost::shared_ptr<OuteredgeDataFactory<TYPE> > factory(
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
OuteredgeVariable<TYPE>::OuteredgeVariable(
   const OuteredgeVariable<TYPE>& foo):
   hier::Variable(NULL,
                  boost::shared_ptr<hier::PatchDataFactory>())
{
   NULL_USE(foo);
}

template<class TYPE>
void OuteredgeVariable<TYPE>::operator = (
   const OuteredgeVariable<TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
