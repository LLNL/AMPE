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
#ifndef included_CALPHADFreeEnergyStrategyTernary
#define included_CALPHADFreeEnergyStrategyTernary

#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "ConcFreeEnergyStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
class MolarVolumeStrategy;

#include <string>
#include <vector>

class CALPHADFreeEnergyStrategyTernary : public ConcFreeEnergyStrategy
{
 public:
   CALPHADFreeEnergyStrategyTernary(
       std::shared_ptr<tbox::Database> input_db,
       std::shared_ptr<tbox::Database> newton_db,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       MolarVolumeStrategy* mvstrategy, const int conc_l_id,
       const int conc_a_id);

   virtual ~CALPHADFreeEnergyStrategyTernary() { delete d_calphad_fenergy; };

   virtual void setup(std::shared_ptr<tbox::Database> calphad_db,
                      std::shared_ptr<tbox::Database> newton_db);

   virtual void addDrivingForce(const double time, hier::Patch& patch,
                                const int temperature_id, const int phase_id,
                                const int conc_id, const int f_l_id,
                                const int f_a_id, const int f_b_id,
                                const int rhs_id);

   virtual void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true)
   {
      defaultComputeSecondDerivativeEnergyPhaseL(temperature, c, d2fdc2,
                                                 use_internal_units);
      // if( d2fdc2[0]<0. )
      //   tbox::pout<<"CALPHADFreeEnergyStrategy, WARNING: fcc<0. in phase L
      //   for c="<<c[0]<<"!!!"<<std::endl;
   }
   virtual void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true)
   {
      defaultComputeSecondDerivativeEnergyPhaseA(temperature, c, d2fdc2,
                                                 use_internal_units);
      // if( d2fdc2[0]<0. )
      //   tbox::pout<<"CALPHADFreeEnergyStrategy, WARNING: fcc<0. in phase A
      //   for c="<<c[0]<<"!!!"<<std::endl;
   }
   virtual void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true)
   {
      (void)temperature;
      (void)c;
      (void)d2fdc2;
      (void)use_internal_units;

      tbox::pout << "Function not implemented" << std::endl;
      tbox::SAMRAI_MPI::abort();
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

   bool computeCeqT(const double temperature, const Thermo4PFM::PhaseIndex pi0,
                    const Thermo4PFM::PhaseIndex pi1, double* ceq)
   {
      return d_calphad_fenergy->computeCeqT(temperature, &ceq[0], 50, true);
   }

   void energyVsPhiAndC(const double temperature, const double* const ceq,
                        const bool found_ceq, const double phi_well_scale,
                        const std::string& phi_well_type, const int npts_phi,
                        const int npts_c)
   {
      d_calphad_fenergy->energyVsPhiAndC(temperature, ceq, found_ceq,
                                         phi_well_scale, npts_phi, npts_c);
   }

 protected:
   void defaultComputeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units);
   void defaultComputeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units);

 protected:
   MolarVolumeStrategy* d_mv_strategy;

   Thermo4PFM::CALPHADFreeEnergyFunctionsTernary* d_calphad_fenergy;

   Thermo4PFM::EnergyInterpolationType d_energy_interp_func_type;
   Thermo4PFM::ConcInterpolationType d_conc_interp_func_type;

   void computeMuA(const double t, const double c0, const double c1,
                   double* mu);

   void computeMuL(const double t, const double c0, const double c1,
                   double* mu);

 private:
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
