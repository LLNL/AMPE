/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/

#include "ModifiedBratuProblem.h"

#include <boost/shared_ptr.hpp>

#if !defined(HAVE_PETSC) || !defined(HAVE_SUNDIALS) || !defined(HAVE_HYPRE)

/*
 *************************************************************************
 * If the library is not compiled with PETSC -and- KINSOL, print an error.
 * If we're running autotests, skip the error
 *************************************************************************
 */
#if (TESTING != 1)
#error \
   "This example requires SAMRAI be compiled with KINSOL, PETSC, and HYPRE."
#endif

#else

template class boost::shared_ptr<ModifiedBratuProblem>;

#endif
