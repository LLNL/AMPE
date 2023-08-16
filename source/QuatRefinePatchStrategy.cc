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
#include "QuatRefinePatchStrategy.h"
#include "TimeLocationIndexRobinBcCoefs.h"

#include <cassert>


QuatRefinePatchStrategy::QuatRefinePatchStrategy(
    const std::string& object_name, std::shared_ptr<tbox::Database> input_bc_db,
    const int phase_id, const int eta_id, const int quat_id, const int conc_id,
    const int temperature_id, const int nghosts,
    const double rescaled_temperature_coeff)
    : xfer::RefinePatchStrategy(),
      d_object_name(object_name),
      d_phase_id(phase_id),
      d_eta_id(eta_id),
      d_quat_id(quat_id),
      d_conc_id(conc_id),
      d_temperature_id(temperature_id),
      d_nghosts(nghosts)
{
   assert(!object_name.empty());

   if (d_phase_id >= 0) {
      d_phase_refine_strategy =
          new CartesianRobinBcHelperWithDepth(tbox::Dimension(NDIM),
                                              "PhaseBcHelper");
      d_phase_refine_strategy->setTargetDataId(d_phase_id);
      std::shared_ptr<tbox::Database> phase_bc_db =
          input_bc_db->getDatabase("Phase");
      d_phase_bc_coefs =
          new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                              "QRPSPhaseBcCoefs", phase_bc_db);
      d_phase_refine_strategy->setCoefImplementation(d_phase_bc_coefs);
   }

   if (d_eta_id >= 0) {
      d_eta_refine_strategy =
          new solv::CartesianRobinBcHelper(tbox::Dimension(NDIM),
                                           "EtaBcHelpe"
                                           "r");
      d_eta_refine_strategy->setTargetDataId(d_eta_id);
      std::shared_ptr<tbox::Database> eta_bc_db =
          input_bc_db->getDatabase("Eta");
      d_eta_bc_coefs =
          new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                              "QRPSEtaBcCoefs", eta_bc_db);
      d_eta_refine_strategy->setCoefImplementation(d_eta_bc_coefs);
   }

   if (d_quat_id >= 0) {
      d_quat_refine_strategy =
          new CartesianRobinBcHelperWithDepth(tbox::Dimension(NDIM),
                                              "QuatBcHelper");
      d_quat_refine_strategy->setTargetDataId(d_quat_id);
      std::shared_ptr<tbox::Database> quat_bc_db =
          input_bc_db->getDatabase("Quat");
      d_quat_bc_coefs =
          new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                              "QRPSQuatBcCoefs", quat_bc_db);
      d_quat_refine_strategy->setCoefImplementation(d_quat_bc_coefs);
   }

   if (d_conc_id >= 0) {
      d_conc_refine_strategy =
          new CartesianRobinBcHelperWithDepth(tbox::Dimension(NDIM),
                                              "ConcBcHelper");
      d_conc_refine_strategy->setTargetDataId(d_conc_id);
      std::shared_ptr<tbox::Database> conc_bc_db =
          input_bc_db->getDatabase("Conc");
      d_conc_bc_coefs =
          new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                              "QRPSConcBcCoefs", conc_bc_db);
      d_conc_refine_strategy->setCoefImplementation(d_conc_bc_coefs);
   }

   if (d_temperature_id >= 0 && input_bc_db->keyExists("Temperature")) {
      d_temp_refine_strategy =
          new solv::CartesianRobinBcHelper(tbox::Dimension(NDIM),
                                           "TemperatureBcHelper");
      d_temp_refine_strategy->setTargetDataId(d_temperature_id);
      std::shared_ptr<tbox::Database> temp_bc_db =
          input_bc_db->getDatabase("Temperature");
      bool flag = false;
      std::string name("boundary_0");
      if (temp_bc_db->isString(name)) {
         std::vector<std::string> specs = temp_bc_db->getStringVector(name);
         if (specs[0] == "file") flag = true;
      }
      if (flag) {
         tbox::plog << "Use TimeLocationIndexRobinBcCoefs for temperature"
                    << std::endl;
         d_temp_bc_coefs =
             new TimeLocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                               "QRPSTemperatureBcCoefs",
                                               temp_bc_db);
         if (rescaled_temperature_coeff > 0.) {
            tbox::plog << "Rescale Temperature boundary conditions by factor "
                       << rescaled_temperature_coeff << std::endl;
            TimeLocationIndexRobinBcCoefs* bc_coefs =
                dynamic_cast<TimeLocationIndexRobinBcCoefs*>(d_temp_bc_coefs);
            bc_coefs->rescaleGcoefficients(rescaled_temperature_coeff);
         }
      } else {
         tbox::plog << "Use solv::LocationIndexRobinBcCoefs for temperature"
                    << std::endl;
         d_temp_bc_coefs =
             new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                                 "QRPSTemperatureBcCoefs",
                                                 temp_bc_db);
         if (rescaled_temperature_coeff > 0.) {
            tbox::plog << "Rescale Temperature boundary conditions by factor "
                       << rescaled_temperature_coeff << std::endl;
            double a, b, g;
            solv::LocationIndexRobinBcCoefs* bc_coefs =
                dynamic_cast<solv::LocationIndexRobinBcCoefs*>(d_temp_bc_coefs);
            for (int n = 0; n < 2 * NDIM; n++) {
               bc_coefs->getCoefficients(n, a, b, g);
               tbox::plog << "old values: " << a << "," << b << "," << g
                          << std::endl;
               g *= rescaled_temperature_coeff;
               tbox::plog << "new values: " << a << "," << b << "," << g
                          << std::endl;
               bc_coefs->setRawCoefficients(n, a, b, g);
            }
         }
      }
      d_temp_refine_strategy->setCoefImplementation(d_temp_bc_coefs);
   } else {
      d_temp_refine_strategy = NULL;
   }
}

//=======================================================================

QuatRefinePatchStrategy::~QuatRefinePatchStrategy() {}

//=======================================================================

void QuatRefinePatchStrategy::setPhysicalBoundaryConditions(
    hier::Patch& patch, const double fill_time,
    const hier::IntVector& ghost_width_to_fill)
{
   if (d_phase_id >= 0) {
      d_phase_refine_strategy->setPhysicalBoundaryConditions(
          patch, fill_time, ghost_width_to_fill);
   }

   if (d_eta_id >= 0) {
      d_eta_refine_strategy->setPhysicalBoundaryConditions(patch, fill_time,
                                                           ghost_width_to_fill);
   }

   if (d_quat_id >= 0) {
      d_quat_refine_strategy->setPhysicalBoundaryConditions(
          patch, fill_time, ghost_width_to_fill);
   }

   if (d_conc_id >= 0) {
      d_conc_refine_strategy->setPhysicalBoundaryConditions(
          patch, fill_time, ghost_width_to_fill);
   }

   if (d_temperature_id >= 0 && d_temp_refine_strategy != NULL) {
      d_temp_refine_strategy->setPhysicalBoundaryConditions(
          patch, fill_time, ghost_width_to_fill);
   }
}

//=======================================================================
//
// Print all class data members to given output stream.

void QuatRefinePatchStrategy::printClassData(std::ostream& os) const
{
   os << "\nQuatRefinePatchStrategy::printClassData..." << std::endl;
   os << "QuatRefinePatchStrategy: this = " << (QuatRefinePatchStrategy*)this
      << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_phase_id =   " << d_phase_id << std::endl;
   os << "d_eta_id =   " << d_eta_id << std::endl;
   os << "d_quat_id =   " << d_quat_id << std::endl;
   os << "d_conc_id =   " << d_conc_id << std::endl;
   os << "d_temperature_id = " << d_temperature_id << std::endl;
}
