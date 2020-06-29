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
#ifndef included_KimMobilityStrategy
#define included_KimMobilityStrategy

#include "CALPHADFreeEnergyFunctions.h"
#include "SimpleQuatMobilityStrategy.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Patch.h"

class QuatModel;

#include <string>

using namespace SAMRAI;

/*
 * Based on S.G. Kim, Acta Mat. 55 (2007), p. 4391-4399
 */
class KimMobilityStrategy : public SimpleQuatMobilityStrategy
{
 public:
   KimMobilityStrategy(QuatModel* quat_model, const int conc_l_id,
                       const int conc_s_id, const int temp_id,
                       const EnergyInterpolationType energy_interp_func_type,
                       const ConcInterpolationType conc_interp_func_type,
                       std::shared_ptr<tbox::Database> conc_db,
                       const unsigned ncompositions);

   void computePhaseMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual double evaluateMobility(const double temp,
                                   const std::vector<double>& phaseconc) = 0;

 protected:
   const int d_conc_l_id;
   const int d_conc_s_id;
   const int d_temp_id;

   const unsigned d_ncompositions;

   FreeEnergyFunctions* d_fenergy;

 private:
   void update(std::shared_ptr<pdat::CellData<double> > cd_te,
               std::shared_ptr<pdat::CellData<double> > cd_cl,
               std::shared_ptr<pdat::CellData<double> > cd_cs,
               std::shared_ptr<pdat::CellData<double> > cd_mobility,
               std::shared_ptr<hier::Patch> patch);

   /*
    * Timers for performance measurement.
    */
   std::shared_ptr<tbox::Timer> t_compute;
};

#endif
