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
#include "PartitionCoeffRefinePatchStrategy.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include <boost/make_shared.hpp>

#include <cassert>
using namespace std;

PartitionCoeffRefinePatchStrategy::PartitionCoeffRefinePatchStrategy(
   const string& object_name,
   boost::shared_ptr< tbox::Database > input_bc_db,
   const int partition_coeff_id )
   : xfer::RefinePatchStrategy(),
     d_object_name( object_name ),
     d_partition_coeff_id( partition_coeff_id )
{
   assert(!object_name.empty());
   assert( partition_coeff_id>=0 );

   d_partition_coeff_refine_strategy =
      new solv::CartesianRobinBcHelper(tbox::Dimension(NDIM), "partition_coeffBcHelper" );
   d_partition_coeff_refine_strategy->setTargetDataId( d_partition_coeff_id );
   boost::shared_ptr<tbox::Database> partition_coeff_bc_db =
      input_bc_db->getDatabase( "PartitionCoeff" );
   d_partition_coeff_bc_coefs = new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
      "PartitionCoeffBcCoefs", partition_coeff_bc_db );
   d_partition_coeff_refine_strategy->setCoefImplementation( d_partition_coeff_bc_coefs );
}

//=======================================================================

void PartitionCoeffRefinePatchStrategy::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill )
{
   assert( d_partition_coeff_id>=0 );
   
   d_partition_coeff_refine_strategy->setPhysicalBoundaryConditions(
      patch,
      fill_time,
      ghost_width_to_fill );
   
}

//=======================================================================
//
// Print all class data members to given output stream.

void PartitionCoeffRefinePatchStrategy::printClassData(ostream& os) const
{
   os << "\nPartitionCoeffRefinePatchStrategy::printClassData..." << endl;
   os << "PartitionCoeffRefinePatchStrategy: this = " << (PartitionCoeffRefinePatchStrategy*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_partition_coeff_id =   " << d_partition_coeff_id << endl;
}
