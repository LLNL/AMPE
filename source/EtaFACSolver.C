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
#include "EtaFACSolver.h"
#include "EtaFACOps.h"

using namespace std;

EtaFACSolver::EtaFACSolver (
   const std::string &object_name,
   boost::shared_ptr<EtaFACOps> fac_ops,
   boost::shared_ptr<tbox::Database> database )
   :
   EllipticFACSolver( object_name, fac_ops, database )
{
   return;
}

void EtaFACSolver::setOperatorCoefficients(
   const int phase_id,
   const int eta_id,
   const int eta_mobility_id,
   const double epsilon_eta,
   const double gamma,
   const string phase_interp_func_type,
   const double eta_well_scale,
   const string eta_well_func_type )
{
   boost::shared_ptr<EtaFACOps> eta_fac_ops (
      boost::dynamic_pointer_cast< EtaFACOps, EllipticFACOps >( d_fac_ops ) );

   eta_fac_ops->setOperatorCoefficients(
      phase_id,
      eta_id,
      eta_mobility_id,
      epsilon_eta, 
      gamma,
      phase_interp_func_type,
      eta_well_scale,
      eta_well_func_type );

   finalizeCoefficients();
}
