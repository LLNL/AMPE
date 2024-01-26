// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
//
#include "PartitionCoeffRefinePatchStrategy.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include <cassert>


PartitionCoeffRefinePatchStrategy::PartitionCoeffRefinePatchStrategy(
    const std::string& object_name, std::shared_ptr<tbox::Database> input_bc_db,
    const int partition_coeff_id)
    : xfer::RefinePatchStrategy(),
      d_object_name(object_name),
      d_partition_coeff_id(partition_coeff_id)
{
   assert(!object_name.empty());
   assert(partition_coeff_id >= 0);

   d_partition_coeff_refine_strategy =
       new solv::CartesianRobinBcHelper(tbox::Dimension(NDIM),
                                        "partition_coeffBcHelper");
   d_partition_coeff_refine_strategy->setTargetDataId(d_partition_coeff_id);
   std::shared_ptr<tbox::Database> partition_coeff_bc_db =
       input_bc_db->getDatabase("PartitionCoeff");
   d_partition_coeff_bc_coefs =
       new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                           "PartitionCoeffBcCoefs",
                                           partition_coeff_bc_db);
   d_partition_coeff_refine_strategy->setCoefImplementation(
       d_partition_coeff_bc_coefs);
}

//=======================================================================

void PartitionCoeffRefinePatchStrategy::setPhysicalBoundaryConditions(
    hier::Patch& patch, const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
   assert(d_partition_coeff_id >= 0);

   d_partition_coeff_refine_strategy->setPhysicalBoundaryConditions(
       patch, fill_time, ghost_width_to_fill);
}

//=======================================================================
//
// Print all class data members to given output stream.

void PartitionCoeffRefinePatchStrategy::printClassData(std::ostream& os) const
{
   os << "\nPartitionCoeffRefinePatchStrategy::printClassData..." << std::endl;
   os << "PartitionCoeffRefinePatchStrategy: this = "
      << (PartitionCoeffRefinePatchStrategy*)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_partition_coeff_id =   " << d_partition_coeff_id << std::endl;
}
