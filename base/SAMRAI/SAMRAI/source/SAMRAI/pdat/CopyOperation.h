/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Copy operation on single array data elements templated on data type
 *
 ************************************************************************/

#ifndef included_pdat_CopyOperation
#define included_pdat_CopyOperation

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class CopyOperation<TYPE> encapsulates a copy
 * operation into an object.
 */

template<class TYPE>
class CopyOperation
{
public:
   /*!
    * The default constructor does nothing interesting.
    */
   CopyOperation();

   /*!
    * The destructor does nothing interesting.
    */
   ~CopyOperation();

   /*!
    * The operator copies the source value to the destination.
    */
   void
   operator () (
      TYPE& vdst,
      const TYPE& vsrc) const;

private:
   CopyOperation(
      const CopyOperation&);              // not implemented
   CopyOperation&
   operator = (
      const CopyOperation&);              // not implemented
};

}
}

#include "SAMRAI/pdat/CopyOperation.C"

#endif
