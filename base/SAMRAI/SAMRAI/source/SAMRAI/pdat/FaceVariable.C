/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_FaceVariable_C
#define included_pdat_FaceVariable_C

#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/FaceDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/make_shared.hpp>

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
FaceVariable<TYPE>::FaceVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth,
   const bool fine_boundary_represents_var):
   hier::Variable(name,
                  boost::make_shared<FaceDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var)),
   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
}

template<class TYPE>
FaceVariable<TYPE>::~FaceVariable()
{
}

template<class TYPE>
int FaceVariable<TYPE>::getDepth() const
{
   boost::shared_ptr<FaceDataFactory<TYPE> > factory(getPatchDataFactory());
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
FaceVariable<TYPE>::FaceVariable(
   const FaceVariable<TYPE>& foo):
   hier::Variable(NULL,
                  boost::shared_ptr<hier::PatchDataFactory>())
{
   NULL_USE(foo);
}

template<class TYPE>
void FaceVariable<TYPE>::operator = (
   const FaceVariable<TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
