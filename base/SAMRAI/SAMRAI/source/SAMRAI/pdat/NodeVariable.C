/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeVariable_C
#define included_pdat_NodeVariable_C

#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * Constructor and destructor for node variable objects
 *
 *************************************************************************
 */

template<class TYPE>
NodeVariable<TYPE>::NodeVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth,
   bool fine_boundary_represents_var):
   hier::Variable(name,
                  boost::make_shared<NodeDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var)),

   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
}

template<class TYPE>
NodeVariable<TYPE>::~NodeVariable()
{
}

template<class TYPE>
int NodeVariable<TYPE>::getDepth() const
{
   boost::shared_ptr<NodeDataFactory<TYPE> > factory(getPatchDataFactory());
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
NodeVariable<TYPE>::NodeVariable(
   const NodeVariable<TYPE>& foo):
   hier::Variable(NULL,
                  boost::shared_ptr<hier::PatchDataFactory>())
{
   NULL_USE(foo);
}

template<class TYPE>
void NodeVariable<TYPE>::operator = (
   const NodeVariable<TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
