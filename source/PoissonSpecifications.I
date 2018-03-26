// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
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
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
/*
 * Description: Specifications for the scalar Poisson equation
 * adapted from SAMRAI library class
 */

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/IEEE.h"
#include <cassert>

using namespace SAMRAI;
#define SAMRAI_INLINE_KEYWORD inline

/*
*******************************************************************
* Default constructor                                             *
*******************************************************************
*/

SAMRAI_INLINE_KEYWORD
PoissonSpecifications::PoissonSpecifications(
   const std::string &object_name
   ) : d_object_name(object_name),
       d_D_id(-1),
       d_D_constant(1.0),
       d_C_zero(true),
       d_C_id(-1),
       d_C_constant(0.0),
       d_M_id(-1),
       d_M_constant(1.0) {
   return;
}


/*
*******************************************************************
* Copy constructor                                                *
*******************************************************************
*/

SAMRAI_INLINE_KEYWORD
PoissonSpecifications::PoissonSpecifications(
   const std::string &object_name,
   const PoissonSpecifications &r
   ) : d_object_name(object_name),
       d_D_id(r.d_D_id),
       d_D_constant(r.d_D_constant),
       d_C_zero(r.d_C_zero),
       d_C_id(r.d_C_id),
       d_C_constant(r.d_C_constant),
       d_M_id(r.d_M_id),
       d_M_constant(r.d_M_constant) {
   return;
}

/*
*******************************************************************
* Destructor (does nothing).                                      *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
PoissonSpecifications::~PoissonSpecifications() {
   return;
}


/*
*******************************************************************
* Assignment operator                                             *
*******************************************************************
*/

SAMRAI_INLINE_KEYWORD
const PoissonSpecifications& PoissonSpecifications::operator=(
   const PoissonSpecifications &r ) {

   d_D_id = r.d_D_id;
   d_D_constant = r.d_D_constant;

   d_M_id = r.d_M_id;
   d_M_constant = r.d_M_constant;

   d_C_zero = r.d_C_zero;
   d_C_id = r.d_C_id;
   d_C_constant = r.d_C_constant;

   return *this;
}

/*
*******************************************************************
* Set the patch data index for variable D.                        *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setDPatchDataId( int id ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if( id < 0 ) {
      TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
   }
#endif
   d_D_id = id;
   d_D_constant = tbox::IEEE::getSignalingNaN();
   return;
}

/*
*******************************************************************
* Set the constant value variable D.                              *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setDConstant( double constant ) {
   d_D_id = -1;
   d_D_constant = constant;
   return;
}

/*
*******************************************************************
* Whether D is variable.                                          *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::dIsVariable() const {
   return d_D_id != -1;
}

/*
*******************************************************************
* Whether D is constant.                                          *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::dIsConstant() const {
   return d_D_id == -1;
}

/*
*******************************************************************
* Get the patch data index for variable D.                        *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
int PoissonSpecifications::getDPatchDataId() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_D_id == -1 ) {
      TBOX_ERROR(d_object_name << ": D not prepresented by a patch data.\n");
   }
#endif
   return d_D_id;
}

/*
*******************************************************************
* Get the constant D value.                                       *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
double PoissonSpecifications::getDConstant() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_D_id != -1 ) {
      TBOX_ERROR(d_object_name << ": D not prepresented by a constant.\n");
   }
#endif
   return d_D_constant;
}

/*
*******************************************************************
* Set the patch data index for variable M.                        *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setMPatchDataId( int id ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if( id < 0 ) {
      TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
   }
#endif
   d_M_id = id;
   d_M_constant = tbox::IEEE::getSignalingNaN();
   return;
}

/*
*******************************************************************
* Set the constant value variable M.                              *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setMConstant( double constant ) {
   d_M_id = -1;
   d_M_constant = constant;
   return;
}

/*
*******************************************************************
* Whether M is variable.                                          *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::mIsVariable() const {
   return d_M_id != -1;
}

/*
*******************************************************************
* Whether M is constant.                                          *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::mIsConstant() const {
   return d_M_id == -1;
}

/*
*******************************************************************
* Get the patch data index for variable D.                        *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
int PoissonSpecifications::getMPatchDataId() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_M_id == -1 ) {
      TBOX_ERROR(d_object_name << ": M not prepresented by a patch data.\n");
   }
#endif
   return d_M_id;
}

/*
*******************************************************************
* Get the constant M value.                                       *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
double PoissonSpecifications::getMConstant() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_M_id != -1 ) {
      TBOX_ERROR(d_object_name << ": M not prepresented by a constant.\n");
   }
#endif
   return d_M_constant;
}

/*
*******************************************************************
* Set the constant value variable C.                              *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setCPatchDataId( int id ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if( id < 0 ) {
      TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
   }
#endif
   d_C_zero = false;
   d_C_id = id;
   d_C_constant = tbox::IEEE::getSignalingNaN();
   return;
}

/*
*******************************************************************
* Set the patch data index for variable C.                        *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setCConstant( double constant ) {
   assert( constant>0. || constant<0. );
   d_C_zero = false;
   d_C_id = -1;
   d_C_constant = constant;
   return;
}

/*
*******************************************************************
* Set the value of C to zero.                                     *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setCZero() {
   d_C_zero = true;
   d_C_id = -1;
   d_C_constant = 0.0;
   return;
}

/*
*******************************************************************
* Whether C is variable.                                          *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::cIsVariable() const {
   return d_C_id != -1;
}

/*
*******************************************************************
* Whether C is zero.                                              *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::cIsZero() const {
   return d_C_zero;
}

/*
*******************************************************************
* Whether C is constant.                                          *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::cIsConstant() const {
   return (d_C_id == -1);
}

/*
*******************************************************************
* Get the patch data index for variable C.                        *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
int PoissonSpecifications::getCPatchDataId() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_C_id == -1 ) {
      TBOX_ERROR(d_object_name << ": C not prepresented by a an index.\n");
   }
#endif
   return d_C_id;
}

/*
*******************************************************************
* Get the constant C value.                                       *
*******************************************************************
*/
SAMRAI_INLINE_KEYWORD
double PoissonSpecifications::getCConstant() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_C_id != -1 ) {
     TBOX_ERROR(d_object_name << ": C is not prepresented by a constant.\n");
   }
#endif
   return d_C_constant;
}
