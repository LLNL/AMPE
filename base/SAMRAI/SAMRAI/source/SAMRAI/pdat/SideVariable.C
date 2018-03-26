/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_SideVariable_C
#define included_pdat_SideVariable_C

#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/SideDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int SideVariable<TYPE>::ALL_DIRECTIONS = -1;

/*
 *************************************************************************
 *
 * Constructor and destructor for side variable objects
 *
 *************************************************************************
 */

template<class TYPE>
SideVariable<TYPE>::SideVariable(
   const tbox::Dimension& dim,
   const std::string& name,
   int depth,
   bool fine_boundary_represents_var,
   int direction):
   hier::Variable(name,
                  boost::make_shared<SideDataFactory<TYPE> >(
                     depth,
                     // default zero ghost cells
                     hier::IntVector::getZero(dim),
                     fine_boundary_represents_var)),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(dim)
{
   TBOX_ASSERT((direction >= ALL_DIRECTIONS) && (direction < getDim().getValue()));

   d_directions = hier::IntVector::getOne(getDim());
   if ((direction != ALL_DIRECTIONS)) {
      // SGS this loop seems stupid, why not just set directions(direction) = 1?
      for (int id = 0; id < getDim().getValue(); id++) {
         d_directions(id) = ((direction == id) ? 1 : 0);
      }
      const hier::IntVector& zero_vector(hier::IntVector::getZero(getDim()));
      setPatchDataFactory(
         boost::make_shared<SideDataFactory<TYPE> >(
            depth,
            zero_vector,
            fine_boundary_represents_var,
            d_directions));
   }
}

template<class TYPE>
SideVariable<TYPE>::~SideVariable()
{
}

template<class TYPE>
const hier::IntVector& SideVariable<TYPE>::getDirectionVector() const
{
   return d_directions;
}

template<class TYPE>
int SideVariable<TYPE>::getDepth() const
{
   boost::shared_ptr<SideDataFactory<TYPE> > factory(getPatchDataFactory());
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
SideVariable<TYPE>::SideVariable(
   const SideVariable<TYPE>& foo):
   hier::Variable(NULL, boost::shared_ptr<hier::PatchDataFactory>()),
   d_directions(hier::IntVector(foo.getDim()))
{
   NULL_USE(foo);
}

template<class TYPE>
void SideVariable<TYPE>::operator = (
   const SideVariable<TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
