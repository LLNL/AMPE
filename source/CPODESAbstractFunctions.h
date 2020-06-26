// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#ifndef included_CPODESAbstractFunctions
#define included_CPODESAbstractFunctions

#include "SAMRAI/SAMRAI_config.h"
#include "SundialsAbstractVector.h"

using namespace SAMRAI;

/**
 * Class CPODESAbstractFunctions is an abstract base class that defines
 * an interface for the user-supplied RHSFunction and preconditioner
 * routines to be used with CPODES and CPSpgmr via the C++ wrapper
 * class CPODESSolver.  To use CPODES with the C++ wrapper one must
 * derive a subclass of this base class and pass it into the CPODESSolver
 * constructor.  The pure virtual member functions in this interface are
 * used by CPODES and CPSpgmr during the ODE integration process.  The
 * complete argument lists in the function signatures defined by CPODES
 * for the user-supplied routines have been preserved for the most part.
 * In a few cases, some arguments do not appear in the function signatures
 * below since they are superfluous via this interface.
 *
 * @see solv::CPODESSolver
 * @see SundialsAbstractVector
 */

class CPODESAbstractFunctions
{
 public:
   /**
    * The constructor and destructor for CPODESAbstractFunctions
    * is empty.
    */
   CPODESAbstractFunctions();
   virtual ~CPODESAbstractFunctions();

   /**
    * User-supplied right-hand side function evaluation.
    *
    * The function arguments are:
    *


    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b y_dot   (OUTPUT){current value of the derivative of y}
    *


    *
    * IMPORTANT: This function must not modify the vector y.
    */
   virtual int evaluateRHSFunction(double t, SundialsAbstractVector* y,
                                   SundialsAbstractVector* y_dot,
                                   int fd_flag) = 0;

   /**
    * User-supplied function to project the current solution and estimated
    * error onto the constraint manifold.
    *
    * The function arguments are:
    *


    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b corr     (OUTPUT){correction such that y+corr is on the constraint
    manifold}
    * - \b epsProj  (INPUT){WRMS norm tolerance for a nonlinear solver
    iteration}
    * - \b err      (INPUT/OUTPUT){projection of err constraint manifold}
    *


    *
    * IMPORTANT: This function must not modify the vector y.
    */
   virtual int applyProjection(double t, SundialsAbstractVector* y,
                               SundialsAbstractVector* corr, double epsProj,
                               SundialsAbstractVector* err) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CPSpgmrPrecondSet(double t, SundialsAbstractVector* y,
                                 SundialsAbstractVector* fy, int jok,
                                 int* jcurPtr, double gamma,
                                 SundialsAbstractVector* vtemp1,
                                 SundialsAbstractVector* vtemp2,
                                 SundialsAbstractVector* vtemp3) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CPSpgmrPrecondSolve(double t, SundialsAbstractVector* y,
                                   SundialsAbstractVector* fy,
                                   SundialsAbstractVector* r,
                                   SundialsAbstractVector* z, double gamma,
                                   double delta, int lr,
                                   SundialsAbstractVector* vtemp) = 0;
};

#endif
