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
#ifndef included_CALPHADFreeEnergyStrategyBinaryFolchPlapp
#define included_CALPHADFreeEnergyStrategyBinaryFolchPlapp

#include "ConcFreeEnergyStrategy.h"
#include "InterpolationType.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
class MolarVolumeStrategy;

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

template <class FreeEnergyFunctionType>
class CALPHADFreeEnergyStrategyBinaryFolchPlapp : public ConcFreeEnergyStrategy
{
 public:
   CALPHADFreeEnergyStrategyBinaryFolchPlapp(
       boost::property_tree::ptree calphad_db,
       std::shared_ptr<tbox::Database> newton_db,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       MolarVolumeStrategy* mvstrategy, const int conc_l_id,
       const int conc_a_id, const int conc_b_id);

   virtual ~CALPHADFreeEnergyStrategyBinaryFolchPlapp(){};

   void setup(boost::property_tree::ptree calphad_db,
              std::shared_ptr<tbox::Database> newton_db);

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int conc_id, const int f_l_id, const int f_a_id,
                        const int f_b_id, const int rhs_id);

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true)
   {
      defaultComputeSecondDerivativeEnergyPhaseL(temperature, c, d2fdc2,
                                                 use_internal_units);
      // if( d2fdc2[0]<0. )
      //   tbox::pout<<"CALPHADFreeEnergyStrategy, WARNING: fcc<0. in phase L
      //   for c="<<c[0]<<"!!!"<<std::endl;
   }
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true)
   {
      defaultComputeSecondDerivativeEnergyPhaseA(temperature, c, d2fdc2,
                                                 use_internal_units);
      // if( d2fdc2[0]<0. )
      //   tbox::pout<<"CALPHADFreeEnergyStrategy, WARNING: fcc<0. in phase A
      //   for c="<<c[0]<<"!!!"<<std::endl;
   }
   void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true)
   {
      defaultComputeSecondDerivativeEnergyPhaseB(temperature, c, d2fdc2,
                                                 use_internal_units);
      // if( d2fdc2[0]<0. )
      //   tbox::pout<<"CALPHADFreeEnergyStrategy, WARNING: fcc<0. in phase B
      //   for c="<<c[0]<<"!!!"<<std::endl;
   }

   void computeSecondDerivativeEnergyPhase(const char phase, const double temp,
                                           const std::vector<double>& c,
                                           std::vector<double>& d2fdc2,
                                           const bool use_internal_units)
   {
      switch (phase) {
         case 'l':
            computeSecondDerivativeEnergyPhaseL(temp, c, d2fdc2,
                                                use_internal_units);
            break;

         case 'a':
            computeSecondDerivativeEnergyPhaseA(temp, c, d2fdc2,
                                                use_internal_units);
            break;

         case 'b':
            computeSecondDerivativeEnergyPhaseB(temp, c, d2fdc2,
                                                use_internal_units);
            break;

         default:
            tbox::pout << "undefined phase=" << phase << "!!!" << std::endl;
            tbox::SAMRAI_MPI::abort();
      }
   }

   void preRunDiagnostics(const double temperature)
   {
      d_calphad_fenergy->preRunDiagnostics(temperature);
   }

 protected:
   Thermo4PFM::EnergyInterpolationType d_energy_interp_func_type;
   Thermo4PFM::ConcInterpolationType d_conc_interp_func_type;

 private:
   std::shared_ptr<FreeEnergyFunctionType> d_calphad_fenergy;

   MolarVolumeStrategy* d_mv_strategy;

   void defaultComputeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units);
   void defaultComputeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units);
   void defaultComputeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units);

   double computeMuA(const double t, const double c);

   double computeMuL(const double t, const double c);

   double computeMuB(const double t, const double c);

   void addDrivingForceOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_f_b,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b, const hier::Box& pbox);

   // implements pure virtual functions of base class ConcFreeEnergyStrategy
   void computeFreeEnergy(hier::Patch& patch, const int temperature_id,
                          const int f_id, const int c_i_id,
                          const Thermo4PFM::PhaseIndex pi, const bool gp);

   void computeDerivFreeEnergy(hier::Patch& patch, const int temperature_id,
                               const int f_id, const int c_i_id,
                               const Thermo4PFM::PhaseIndex pi);


   void computeFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       const Thermo4PFM::PhaseIndex pi, const bool gp);

   void computeDerivFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       const Thermo4PFM::PhaseIndex pi);
};

#endif
