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
#ifndef included_QuatRefinePatchStrategy
#define included_QuatRefinePatchStrategy

#include "CartesianRobinBcHelperWithDepth.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/solv/CartesianRobinBcHelper.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <string>

using namespace SAMRAI;

/**
 * Class QuatRefinePatchStrategy is a light weight class implementing the
 * function of RefinePatchStrategy for quaternions
 */

class QuatRefinePatchStrategy : public xfer::RefinePatchStrategy
{
 public:
   QuatRefinePatchStrategy(const std::string& object_name,
                           std::shared_ptr<tbox::Database> input_db,
                           const int phase_id, const int eta_id,
                           const int quat_id, const int conc_id,
                           const int temperature_id, const int nghosts,
                           const double rescaled_temperature_coeff = -1.);

   /**
    * Destructor for QuatRefinePatchStrategy class does nothing.
    */
   ~QuatRefinePatchStrategy();

   /*************************************************************************
    *
    * Methods inherited from RefinePatchStrategy.
    *
    ************************************************************************/

   /*!
    * @brief Set solution ghost cell values along physical boundaries.
    *
    * Function is overloaded from RefinePatchStrategy.
    */

   void setPhysicalBoundaryConditions(
       hier::Patch& patch, const double time,
       const hier::IntVector& ghost_width_to_fill);

   //@{
   /*!
    * @name Empty functions for applying user-defined data refine operations
    */

   /*
    * These are overloaded from RefinePatchStrategy.
    * There are no such user-defined operations here.
    */

   void preprocessRefine(hier::Patch& fine, const hier::Patch& coarse,
                         const hier::Box& fine_box,
                         const hier::IntVector& ratio)
   {
      (void)fine;
      (void)coarse;
      (void)fine_box;
      (void)ratio;
   }

   void postprocessRefine(hier::Patch& fine, const hier::Patch& coarse,
                          const hier::Box& fine_box,
                          const hier::IntVector& ratio)
   {
      (void)fine;
      (void)coarse;
      (void)fine_box;
      (void)ratio;
   }

   /**
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.
    */
   hier::IntVector getRefineOpStencilWidth(const tbox::Dimension& dim) const
   {
      return (hier::IntVector(dim, d_nghosts));
   }

   //@}

   /**
    * Write class data to given output stream.
    */
   void printClassData(std::ostream& os) const;

 private:
   std::string d_object_name;

   int d_phase_id;
   CartesianRobinBcHelperWithDepth* d_phase_refine_strategy;
   solv::LocationIndexRobinBcCoefs* d_phase_bc_coefs;

   int d_eta_id;
   solv::CartesianRobinBcHelper* d_eta_refine_strategy;
   solv::LocationIndexRobinBcCoefs* d_eta_bc_coefs;

   int d_quat_id;
   CartesianRobinBcHelperWithDepth* d_quat_refine_strategy;
   solv::LocationIndexRobinBcCoefs* d_quat_bc_coefs;

   int d_conc_id;
   CartesianRobinBcHelperWithDepth* d_conc_refine_strategy;
   solv::LocationIndexRobinBcCoefs* d_conc_bc_coefs;

   int d_temperature_id;
   solv::CartesianRobinBcHelper* d_temp_refine_strategy;
   solv::RobinBcCoefStrategy* d_temp_bc_coefs;

   const int d_nghosts;
};

#endif
