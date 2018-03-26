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
#ifndef included_PhaseConcentrationRefinePatchStrategy
#define included_PhaseConcentrationRefinePatchStrategy

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/solv/CartesianRobinBcHelper.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <boost/make_shared.hpp>

#include <string>
using namespace std;
using namespace SAMRAI;

/**
 * Class PhaseConcentrationRefinePatchStrategy is a light weight class implementing the
 * function of RefinePatchStrategy for phase concentrations
 */

class PhaseConcentrationRefinePatchStrategy : 
   public xfer::RefinePatchStrategy
{
public: 

   PhaseConcentrationRefinePatchStrategy(
      const string& object_name,
      boost::shared_ptr< tbox::Database > input_db,
      const int conc_l_id,
      const int conc_a_id,
      const int conc_b_id );

   /**
    * Destructor for PhaseConcentrationRefinePatchStrategy class does nothing.
    */
   ~PhaseConcentrationRefinePatchStrategy(){};

   /*!
    * @brief Set solution ghost cell values along physical boundaries.
    *
    * Function is overloaded from RefinePatchStrategy.
    */

   void setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double time,
      const hier::IntVector& ghost_width_to_fill );

   //@{
   /*!
    * @name Empty functions for applying user-defined data refine operations
    */

   /*
    * These are overloaded from RefinePatchStrategy.
    * There are no such user-defined operations here.
    */

   void preprocessRefine(hier::Patch& fine,
                         const hier::Patch& coarse,
                         const hier::Box& fine_box,
                         const hier::IntVector& ratio)
   {
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
   }

   void postprocessRefine(hier::Patch& fine,
                          const hier::Patch& coarse,
                          const hier::Box& fine_box,
                          const hier::IntVector& ratio)
   {
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
   }

   hier::IntVector getRefineOpStencilWidth() const 
   {
       return( hier::IntVector(tbox::Dimension(NDIM),1) );
   }

   //@}

   /**
    * Write class data to given output stream.
    */
   void printClassData( ostream& os ) const;

private:

   string d_object_name;

   int d_conc_l_id;
   solv::CartesianRobinBcHelper* d_conc_l_refine_strategy;
   solv::LocationIndexRobinBcCoefs* d_conc_l_bc_coefs;

   int d_conc_a_id;
   solv::CartesianRobinBcHelper* d_conc_a_refine_strategy;
   solv::LocationIndexRobinBcCoefs* d_conc_a_bc_coefs;

   int d_conc_b_id;
   solv::CartesianRobinBcHelper* d_conc_b_refine_strategy;
   solv::LocationIndexRobinBcCoefs* d_conc_b_bc_coefs;
};

#endif
