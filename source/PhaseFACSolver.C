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
#include "PhaseFACSolver.h"
#include "PhaseFACOps.h"

using namespace std;

PhaseFACSolver::PhaseFACSolver (
   const std::string &object_name,
   boost::shared_ptr<PhaseFACOps> fac_ops,
   boost::shared_ptr<tbox::Database> database )
   :
   EllipticFACSolver( object_name, fac_ops, database )
{
   t_set_op_coef = tbox::TimerManager::getManager()->getTimer("PhaseFACSolver::setOperatorCoefficients");
   
   return;
}

void PhaseFACSolver::setOperatorCoefficients(
   const int phase_id,
   const int eta_id,
   const int phase_mobility_id,
   const double epsilon_phase,
   const double gamma,
   const string phase_interp_func_type,
   const double phase_well_scale,
   const string phase_well_func_type,
   const double eta_well_scale,
   const string eta_well_func_type )
{
   t_set_op_coef->start();

   boost::shared_ptr<PhaseFACOps> phase_fac_ops (
      boost::dynamic_pointer_cast< PhaseFACOps, EllipticFACOps >( d_fac_ops ) );

   phase_fac_ops->setOperatorCoefficients(
      phase_id,
      eta_id,
      phase_mobility_id,
      epsilon_phase, 
      gamma,
      phase_interp_func_type,
      phase_well_scale,
      phase_well_func_type,
      eta_well_scale,
      eta_well_func_type );

   finalizeCoefficients();

   t_set_op_coef->stop();
}
