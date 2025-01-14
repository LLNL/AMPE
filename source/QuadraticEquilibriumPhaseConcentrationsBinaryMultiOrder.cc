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
#include "QuadraticEquilibriumPhaseConcentrationsBinaryMultiOrder.h"
#include "QuadraticFreeEnergyFunctionsBinary.h"
#include "FuncFort.h"


QuadraticEquilibriumPhaseConcentrationsBinaryMultiOrder::
    QuadraticEquilibriumPhaseConcentrationsBinaryMultiOrder(
        const int conc_l_id, const int conc_a_id,
        const QuatModelParameters& model_parameters,
        std::shared_ptr<tbox::Database> conc_db)
    : EquilibriumPhaseConcentrationsBinaryMultiOrder(conc_l_id, conc_a_id,
                                                     conc_db)
{
   std::shared_ptr<tbox::Database> quad_db = conc_db->getDatabase("Quadratic");

   double Tref = quad_db->getDouble("T_ref");

   double A_liquid = quad_db->getDouble("A_liquid");
   double Ceq_liquid = quad_db->getDouble("Ceq_liquid");
   double m_liquid = quad_db->getDouble("m_liquid");

   double A_solid_A = quad_db->getDouble("A_solid");
   double Ceq_solid_A = quad_db->getDouble("Ceq_solid");
   double m_solid = quad_db->getDouble("m_solid");

   d_fenergy.reset(new Thermo4PFM::QuadraticFreeEnergyFunctionsBinary(
       Tref, A_liquid, Ceq_liquid, m_liquid, A_solid_A, Ceq_solid_A, m_solid,
       model_parameters.energy_interp_func_type(),
       model_parameters.conc_interp_func_type()));
}
