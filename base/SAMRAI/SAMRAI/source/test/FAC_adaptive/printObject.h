/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Misc printing functions in FAC solver test.
 *
 ************************************************************************/
#ifndef included_printObject_h
#define included_printObject_h

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/Box.h"

#include <boost/shared_ptr.hpp>
#include <string>

using namespace std;
using namespace SAMRAI;

/*!
 * @brief Print a box
 */
int
printObject(
   std::ostream& os,
   const hier::Box& box,
   const std::string& border = std::string(),
   unsigned short depth = 0);

/*!
 * @brief Print a patch data object
 */
int
printObject(
   std::ostream& os,
   const hier::PatchData& pdat,
   const std::string& border = std::string(),
   unsigned short depth = 0);

/*!
 * @brief Print an array
 */
int
printObject(
   std::ostream& os,
   const tbox::Dimension& dim,
   const double* a_ptr,
   const int* a_lower,
   const int* a_upper,
   const int* lower = NULL,
   const int* upper = NULL);

/*!
 * @brief Print an array data object
 */
template<class T>
int
printObject(
   std::ostream& os,
   const tbox::Dimension& dim,
   const pdat::ArrayData<T>& adat,
   const int depth = 0,
   const std::string& border = std::string());

#endif  // included_printObject_h
