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
#ifndef included_KimMobilityStrategy
#define included_KimMobilityStrategy

#include "SimpleQuatMobilityStrategy.h"
#include "InterpolationType.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Patch.h"

class QuatModel;

#include <string>

using namespace SAMRAI;

/*
 * Based on S.G. Kim, Acta Mat. 55 (2007), p. 4391-4399
 */
template <class FreeEnergyType>
class KimMobilityStrategy : public SimpleQuatMobilityStrategy
{
 public:
   KimMobilityStrategy(
       QuatModel* quat_model, const int conc_l_id, const int conc_a_id,
       const int conc_b_id, const int temp_id,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions);

   void computePhaseMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual double evaluateMobility(const double temp,
                                   const std::vector<double>& phaseconc,
                                   const std::vector<double>& phi) = 0;

 protected:
   const int d_conc_l_id;
   const int d_conc_a_id;
   const int d_conc_b_id;

   const int d_temp_id;

   const unsigned d_ncompositions;

   FreeEnergyType* d_fenergy;

 private:
   void update(std::shared_ptr<pdat::CellData<double> > cd_te,
               std::shared_ptr<pdat::CellData<double> > cd_phi,
               std::shared_ptr<pdat::CellData<double> > cd_cl,
               std::shared_ptr<pdat::CellData<double> > cd_ca,
               std::shared_ptr<pdat::CellData<double> > cd_cb,
               std::shared_ptr<pdat::CellData<double> > cd_mobility,
               std::shared_ptr<hier::Patch> patch);

   /*
    * Timers for performance measurement.
    */
   std::shared_ptr<tbox::Timer> t_compute;
};

#endif
