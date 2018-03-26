/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Utilities class to access common POSIX constants and math ops
 *
 ************************************************************************/

#ifndef included_tbox_MathUtilities_C
#define included_tbox_MathUtilities_C

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Routines to initialize arrays to signaling NaNs.
 *
 *************************************************************************
 */

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToSignalingNaN(
   Array<TYPE>& array)
{
   for (int i = 0; i < array.getSize(); i++) {
      array[i] = getSignalingNaN();
   }
}

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToSignalingNaN(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; i++) {
      array[i] = getSignalingNaN();
   }
}

/*
 *************************************************************************
 *
 * Routines to initialize arrays to max value for type.
 *
 *************************************************************************
 */

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToMax(
   Array<TYPE>& array)
{
   for (int i = 0; i < array.getSize(); i++) {
      array[i] = getMax();
   }
}

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToMax(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; i++) {
      array[i] = getMax();
   }
}

/*
 *************************************************************************
 *
 * Routines to initialize arrays to min value for type.
 *
 *************************************************************************
 */

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToMin(
   Array<TYPE>& array)
{
   for (int i = 0; i < array.getSize(); i++) {
      array[i] = getMin();
   }
}

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToMin(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; i++) {
      array[i] = getMin();
   }
}

/*
 *************************************************************************
 *
 * Routines to initialize arrays to epsilon value for type.
 *
 *************************************************************************
 */

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToEpsilon(
   Array<TYPE>& array)
{
   for (int i = 0; i < array.getSize(); i++) {
      array[i] = getEpsilon();
   }
}

template<class TYPE>
void
MathUtilities<TYPE>::setArrayToEpsilon(
   TYPE* array,
   int n)
{
   for (int i = 0; i < n; i++) {
      array[i] = getEpsilon();
   }
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getZero()
{
   return s_zero;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getOne()
{
   return s_one;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getSignalingNaN()
{
   return s_signaling_nan;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getMax()
{
   return s_max;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getMin()
{
   return s_min;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::getEpsilon()
{
   return s_epsilon;
}

template<class TYPE>
bool
MathUtilities<TYPE>::isNaN(
   const TYPE& value)
{
   NULL_USE(value);
   return false;
}

template<class TYPE>
bool
MathUtilities<TYPE>::equalEps(
   const TYPE& a,
   const TYPE& b)
{
   return a == b;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::Min(
   TYPE a,
   TYPE b)
{
   return a < b ? a : b;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::Max(
   TYPE a,
   TYPE b)
{
   return a > b ? a : b;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::Abs(
   TYPE value)
{
   return value;
}

template<class TYPE>
TYPE
MathUtilities<TYPE>::round(
   TYPE x)
{
   return x;
}

}
}

#endif
