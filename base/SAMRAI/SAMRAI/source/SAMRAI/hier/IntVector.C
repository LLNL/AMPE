/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A n-dimensional integer vector
 *
 ************************************************************************/

#ifndef included_hier_IntVector_C
#define included_hier_IntVector_C

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

namespace SAMRAI {
namespace hier {

IntVector * IntVector::s_zeros[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
IntVector * IntVector::s_ones[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

tbox::StartupShutdownManager::Handler
IntVector::s_initialize_finalize_handler(
   IntVector::initializeCallback,
   0,
   0,
   IntVector::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

IntVector::IntVector():
   d_dim(tbox::Dimension::getInvalidDimension())
{
#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = 0; i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; i++) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   // an explicit setting Invalid is allowed.
   TBOX_DIM_ASSERT((!d_dim.isValid()) ||
      (d_dim >= tbox::Dimension(1) && d_dim <= tbox::Dimension::getMaxDimension()));

#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = 0; i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; i++) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif

}

IntVector::IntVector(
   const tbox::Dimension& dim,
   const int value):
   d_dim(dim)
{
   // an explicit setting Invalid is allowed.
   TBOX_DIM_ASSERT((!d_dim.isValid()) ||
      (d_dim >= tbox::Dimension(1) && d_dim <= tbox::Dimension::getMaxDimension()));

   if (d_dim.isValid()) {
      for (int i = 0; i < d_dim.getValue(); i++)
         d_vector[i] = value;

#ifdef DEBUG_INITIALIZE_UNDEFINED
      for (int i = d_dim.getValue(); i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE;
           i++) {
         d_vector[i] = tbox::MathUtilities<int>::getMin();
      }
#endif
   } else {
      for (int i = 0; i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; i++) {
         d_vector[i] = tbox::MathUtilities<int>::getMin();
      }
   }
}

IntVector::IntVector(
   const tbox::Array<int>& a):
   d_dim(static_cast<unsigned short>(a.getSize()))
{
   TBOX_DIM_ASSERT(a.getSize() > 1 &&
      a.getSize() <= tbox::Dimension::MAXIMUM_DIMENSION_VALUE);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] = a[i];

#ifdef DEBUG_INITIALIZE_UNDEFINED
   for (int i = d_dim.getValue(); i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; i++) {
      d_vector[i] = tbox::MathUtilities<int>::getMin();
   }
#endif
}

IntVector::IntVector(
   const IntVector& rhs):
   d_dim(rhs.getDim())
{
   /*
    * STL needs to be able to copy invalid values.
    */
   if (rhs.getDim().isValid()) {
      TBOX_DIM_ASSERT_CHECK_DIM(rhs.getDim());

      for (int i = 0; i < d_dim.getValue(); i++)
         d_vector[i] = rhs.d_vector[i];
   }
}

IntVector::IntVector(
   const tbox::Dimension& dim,
   const int array[]):
   d_dim(dim)
{
   TBOX_DIM_ASSERT_CHECK_DIM(dim);

   for (int i = 0; i < d_dim.getValue(); i++)
      d_vector[i] = array[i];
}

IntVector::~IntVector()
{
}

std::istream&
operator >> (
   std::istream& s,
   IntVector& rhs)
{
   while (s.get() != '(') ;

   for (int i = 0; i < rhs.getDim().getValue(); i++) {
      s >> rhs(i);
      if (i < rhs.getDim().getValue() - 1)
         while (s.get() != ',') ;
   }

   while (s.get() != ')') ;

   return s;
}

std::ostream& operator << (
   std::ostream& s,
   const IntVector& rhs)
{
   s << '(';

   for (int i = 0; i < rhs.getDim().getValue(); i++) {
      s << rhs(i);
      if (i < rhs.getDim().getValue() - 1)
         s << ",";
   }
   s << ')';

   return s;
}

void
IntVector::putUnregisteredToDatabase(
   tbox::Database& database,
   const std::string& name) const
{
   database.putIntegerArray(name, d_vector, d_dim.getValue());
}

void
IntVector::getFromDatabase(
   tbox::Database& database,
   const std::string& name)
{
   int d = database.getArraySize(name);
   d_dim = tbox::Dimension(static_cast<unsigned short>(d));
   database.getIntegerArray(name, d_vector, d_dim.getValue());
}

/*
 *************************************************************************
 * Sort the sizes of the given IntVector from smallest to largest value.
 *************************************************************************
 */
void
IntVector::sortIntVector(
   const IntVector& values)
{
   const IntVector num_cells = values;

   for (int d = 0; d < d_dim.getValue(); d++) {
      d_vector[d] = d;
   }
   for (int d0 = 0; d0 < d_dim.getValue() - 1; d0++) {
      for (int d1 = d0 + 1; d1 < d_dim.getValue(); d1++) {
         if (values(d_vector[d0]) > values(d_vector[d1])) {
            int tmp_d = d_vector[d0];
            d_vector[d0] = d_vector[d1];
            d_vector[d1] = tmp_d;
         }
      }
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int d = 0; d < d_dim.getValue() - 1; d++) {
      TBOX_ASSERT(values(d_vector[d]) <= values(d_vector[d + 1]));
   }
#endif
}

void
IntVector::initializeCallback()
{
   for (unsigned short d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      s_zeros[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 0);
   }

   for (unsigned short d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      s_ones[d] = new IntVector(tbox::Dimension(static_cast<unsigned short>(d + 1)), 1);
   }
}

void
IntVector::finalizeCallback()
{
   for (int d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      delete s_zeros[d];
      delete s_ones[d];
   }
}

}
}

#endif
