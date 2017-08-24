/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Sum operation on single array data elements templated on data type
 *
 ************************************************************************/

#ifndef included_pdat_SumOperation
#define included_pdat_SumOperation

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class SumOperation<TYPE> encapsulates a summation
 * operation into an object.
 */

template<class TYPE>
class SumOperation
{
public:
   /*!
    * The default constructor does nothing interesting.
    */
   SumOperation();

   /*!
    * The destructor does nothing interesting.
    */
   ~SumOperation();

   /*!
    * The operator adds the source value to the destination.
    */
   void
   operator () (
      TYPE& vdst,
      const TYPE& vsrc) const;

private:
   SumOperation(
      const SumOperation&);             // not implemented
   SumOperation&
   operator = (
      const SumOperation&);              // not implemented
};

}
}

#include "SAMRAI/pdat/SumOperation.C"

#endif
