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
/*
 * Description: Specifications for the scalar Poisson equation
 * adapted from SAMRAI library class
 */

#ifndef included_PoissonSpecifications
#define included_PoissonSpecifications

#ifndef included_SAMRAI_config
#include "SAMRAI/SAMRAI_config.h"
#endif

#include "SAMRAI/tbox/Utilities.h"

#include <string>

/*!
 * @brief Light class holding specifications for cell-centered
 * implementation of the scalar Poisson equation.
 *
 * The scalar Poisson equation is
 * @f$ M \nabla ( D \nabla u ) + C u = f @f$,
 * where C is a scalar field, D is the diffusion coefficient.
 * and u and f are scalar quantities.
 *
 * This class describes the things you can set: M, C, D.
 *
 * Note that the storage and alignment of u, f, M, C and D depend
 * on the implementation of the solver.  For example, if the
 * solver is cell centered, u, f and M, C are cell-centered while
 * D is side-centered.
 */

class PoissonSpecifications
{
public:

   /*!
    * @brief Constructor.
    *
    * Sets the specifications to their default state:
    * - C is zero
    * - D is uniformly 1
    *
    * @param object_name Name of object.
    */
   PoissonSpecifications( const std::string &object_name );

   /*!
    * @brief Copy constructor.
    */
   PoissonSpecifications(
      const std::string &object_name,
      const PoissonSpecifications &r );

   /*!
    * @brief Destructor (does nothing).
    */
   virtual ~PoissonSpecifications();

   /*!
    * @brief Assignment operator
    *
    * Assign everything except name.
    */
   PoissonSpecifications&
   operator = (
      const PoissonSpecifications& r)
   {
      d_D_id = r.d_D_id;
      d_D_constant = r.d_D_constant;
      d_C_zero = r.d_C_zero;
      d_C_id = r.d_C_id;
      d_C_constant = r.d_C_constant;
      d_M_id = r.d_M_id;
      d_M_constant = r.d_M_constant;
      return *this;
   }

   /*!
    * @brief Print out class data.
    */
   virtual void printClassData( std::ostream &stream ) const;

   //@{
   //! @name Functions for setting and getting D

   /*!
    * @brief Set the patch data index for variable D.
    *
    * In addition, disregard any previous value
    * specified by setDConstant().
    */
   void
   setDPatchDataId(
      int id)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (id < 0) {
         TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
      }
#endif
      d_D_id = id;
      d_D_constant = 0.0;
   }

   /*!
    * @brief Set the constant value variable D.
    *
    * In addition, disregard any previous patch data index
    * specified by setDPatchDataId().
    */
   void
   setDConstant(
      double constant)
   {
      d_D_id = -1;
      d_D_constant = constant;
   }

   /*!
    * @brief Whether D is variable (described by a patch data id).
    *
    * @return True if D is variable, described by the patch data
    *         id given in setCPatchDataId().
    */
   bool dIsVariable() const;

   /*!
    * @brief Whether D is constant.
    *
    * @return True if D is constant, as specified by setCConstant().
    */
   bool dIsConstant() const;

   /*!
    * @brief Get D's patch data id
    *
    * Error if D is not represented by a patch data id.
    *
    * @return D's id
    */
   int getDPatchDataId() const;

   /*!
    * @brief Get D constant value
    *
    * Error if D is not represented by a constant.
    *
    * @return D's constant value
    */
   double getDConstant() const;

   //@{
   //! @name Functions for setting and getting M

   /*!
    * @brief Set the patch data index for variable M.
    *
    * In addition, disregard any previous value
    * specified by setMConstant().
    */
   void setMPatchDataId( int id );

   /*!
    * @brief Set the constant value variable M.
    *
    * In addition, disregard any previous patch data index
    * specified by setMPatchDataId().
    */
   void setMConstant( double constant );

   /*!
    * @brief Whether M is variable (described by a patch data id).
    *
    * @return True if M is variable, described by the patch data
    *         id given in setMPatchDataId().
    */
   bool mIsVariable() const;

   /*!
    * @brief Whether M is constant.
    *
    * @return True if M is constant, as specified by setMConstant().
    */
   bool mIsConstant() const;

   /*!
    * @brief Get M's patch data id
    *
    * Error if M is not represented by a patch data id.
    *
    * @return M's id
    */
   int getMPatchDataId() const;

   /*!
    * @brief Get M constant value
    *
    * Error if M is not represented by a constant.
    *
    * @return M's constant value
    */
   double getMConstant() const;

   //@}

   //@{
   //! @name Functions for setting and getting C

   /*!
    * @brief Set the patch data index for C.
    *
    * In addition, disregard any previous values
    * specified by setCConstant() or setCZero().
    */
   void
   setCPatchDataId(
      int id)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (id < 0) {
         TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
      }
#endif
      d_C_zero = false;
      d_C_id = id;
      d_C_constant = 0.0;
   }

   /*!
    * @brief Set C to a constant.
    *
    * In addition, disregard any previous value
    * specified by setCPatchDataId() or setCZero().
    *
    * If you want to set C to zero, use setCZero() instead.
    * This allows solvers to take advantage of fact C is absent.
    */
   void setCConstant( double constant );

   /*!
    * @brief Set the value of C to zero.
    *
    * In addition, disregard any previous patch data index
    * specified by setCPatchDataId() and any previous constant
    * specified by setCConstant().
    */
   void setCZero();

   /*!
    * @brief Whether C is variable (described by a patch data id).
    *
    * @return True if C is variable, described by the patch data
    *         id given in setCPatchDataId().
    */
   bool cIsVariable() const;

   /*!
    * @brief Whether C is zero.
    *
    * As it pertains to what this function returns,
    * C is zero @em only by calling setCZero().
    * Calling setCConstant() does @em not make C zero,
    * even if you pass in the value of zero.
    *
    * @return True if C is exactly zero, as set by setCZero().
    */
   bool cIsZero() const;

   /*!
    * @brief Whether C is constant.
    *
    * As it pertains to what this function returns,
    * C is constant @em only by calling setCConstant().
    * Calling setCZero() does @em not make C a constant.
    *
    * @return True if C is constant, as specified by setCConstant().
    */
   bool cIsConstant() const;

   /*!
    * @brief Get C's patch data id
    *
    * Error if C is not represented by a patch data id.
    *
    * @return C's patch data id
    */
   int getCPatchDataId() const;

   /*!
    * @brief Get C as a constant value.
    *
    * Error if C is not represented by a constant.
    *
    * @return C's constant value
    */
   double
   getCConstant() const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_C_id != -1 || d_C_zero) {
         TBOX_ERROR(d_object_name << ": C is not prepresented by a constant.\n");
      }
#endif
      return d_C_constant;
   }

   //@}

   /**
    * @brief Get the name of this object.
    *
    * @return The name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:

   /*!
    * @brief Object name.
    */
   std::string d_object_name;

   int d_D_id;
   double d_D_constant;

   bool d_C_zero;
   int d_C_id;
   double d_C_constant;

   int d_M_id;
   double d_M_constant;

};

#include "PoissonSpecifications.I"

#endif
