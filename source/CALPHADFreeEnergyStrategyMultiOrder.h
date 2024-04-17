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
#ifndef included_CALPHADFreeEnergyStrategyMultiOrder
#define included_CALPHADFreeEnergyStrategyMultiOrder

#include "CALPHADFreeEnergyStrategyBinary.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
class MolarVolumeStrategy;

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

template <class FreeEnergyFunctionType>
class CALPHADFreeEnergyStrategyMultiOrder
    : public CALPHADFreeEnergyStrategyBinary<FreeEnergyFunctionType>
{
 public:
   CALPHADFreeEnergyStrategyMultiOrder(
       boost::property_tree::ptree calphad_db,
       std::shared_ptr<tbox::Database> newton_db,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       const short _norderp_A, MolarVolumeStrategy* mvstrategy,
       const int conc_l_id, const int conc_a_id, const int conc_b_id);

   ~CALPHADFreeEnergyStrategyMultiOrder(){};

   void addDrivingForce(std::shared_ptr<pdat::CellData<double> > cd_rhs,
                        std::shared_ptr<pdat::CellData<double> > cd_temperature,
                        std::shared_ptr<pdat::CellData<double> > cd_phi,
                        std::shared_ptr<pdat::CellData<double> > cd_eta,
                        std::shared_ptr<pdat::CellData<double> > cd_f_l,
                        std::shared_ptr<pdat::CellData<double> > cd_f_a,
                        std::shared_ptr<pdat::CellData<double> > cd_f_b,
                        std::shared_ptr<pdat::CellData<double> > cd_c_l,
                        std::shared_ptr<pdat::CellData<double> > cd_c_a,
                        std::shared_ptr<pdat::CellData<double> > cd_c_b,
                        const hier::Box& pbox) override;

 private:
   // number of order parameters associated with phase A
   const short d_norderp_A;
};

#endif
