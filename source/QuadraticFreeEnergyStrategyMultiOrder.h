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
#ifndef included_QuadraticFreeEnergyStrategyMultiOrder
#define included_QuadraticFreeEnergyStrategyMultiOrder

#include "QuadraticFreeEnergyFunctionsBinary.h"
#include "ConcFreeEnergyStrategy.h"
#include "InterpolationType.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
class MolarVolumeStrategy;

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

class QuadraticFreeEnergyStrategyMultiOrder : public FreeEnergyStrategy
{
 public:
   QuadraticFreeEnergyStrategyMultiOrder(
       std::shared_ptr<tbox::Database> input_db,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const double vml, const double vma, const int conc_l_id,
       const int conc_a_id);

   ~QuadraticFreeEnergyStrategyMultiOrder(){};

   // implement pure virtual functions of FreeEnergyStrategy
   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp) override;

   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp) override;

   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp) override
   {
      tbox::plog << "QuadraticFreeEnergyStrategyMultiOrder::"
                    "computeFreeEnergySolidB..."
                 << std::endl;
   };

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int eta_id, const int conc_id, const int f_l_id,
                        const int f_a_id, const int f_b_id,
                        const int rhs_id) override;

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override
   {
      (void)temperature;

      defaultComputeSecondDerivativeEnergyPhaseL(c, d2fdc2, use_internal_units);
      // if( d2fdc2[0]<0. )
      //   tbox::pout<<"CALPHADFreeEnergyStrategy, WARNING: fcc<0. in phase L
      //   for c="<<c[0]<<"!!!"<<std::endl;
   }
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override
   {
      (void)temperature;

      defaultComputeSecondDerivativeEnergyPhaseA(c, d2fdc2, use_internal_units);
      // if( d2fdc2[0]<0. )
      //   tbox::pout<<"CALPHADFreeEnergyStrategy, WARNING: fcc<0. in phase A
      //   for c="<<c[0]<<"!!!"<<std::endl;
   }

   void computeSecondDerivativeEnergyPhase(const char phase,
                                           const std::vector<double>& c,
                                           std::vector<double>& d2fdc2,
                                           const bool use_internal_units)
   {
      switch (phase) {
         case 'l':
            computeSecondDerivativeEnergyPhaseL(-1., c, d2fdc2,
                                                use_internal_units);
            break;

         case 'a':
            computeSecondDerivativeEnergyPhaseA(-1., c, d2fdc2,
                                                use_internal_units);
            break;

         default:
            tbox::pout << "undefined phase=" << phase << "!!!" << std::endl;
            tbox::SAMRAI_MPI::abort();
      }
   }

   void preRunDiagnostics(const double temperature) override{};


 private:
   std::shared_ptr<Thermo4PFM::QuadraticFreeEnergyFunctionsBinary>
       d_quadratic_fenergy;

   double d_energy_conv_factor_L;  // molar volume
   double d_energy_conv_factor_A;  // molar volume

   int d_conc_l_id;
   int d_conc_a_id;

   void defaultComputeSecondDerivativeEnergyPhaseL(
       const std::vector<double>& c, std::vector<double>& d2fdc2,
       const bool use_internal_units);
   void defaultComputeSecondDerivativeEnergyPhaseA(
       const std::vector<double>& c, std::vector<double>& d2fdc2,
       const bool use_internal_units);
   void defaultComputeSecondDerivativeEnergyPhaseB(
       const std::vector<double>& c, std::vector<double>& d2fdc2,
       const bool use_internal_units);

   double computeMu(const double t, const double c);

   void addDrivingForceOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a, const hier::Box& pbox);

   void computeFreeEnergy(hier::Patch& patch, const int temperature_id,
                          const int f_id, const int c_i_id,
                          Thermo4PFM::PhaseIndex pi,
                          const double energy_factor);

   void computeDerivFreeEnergy(hier::Patch& patch, const int temperature_id,
                               const int f_id, const int c_i_id,
                               Thermo4PFM::PhaseIndex pi,
                               const double energy_factor);

   void computeFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       Thermo4PFM::PhaseIndex pi, const double energy_factor);

   void computeDerivFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       Thermo4PFM::PhaseIndex pi, const double energy_factor);
};

#endif
