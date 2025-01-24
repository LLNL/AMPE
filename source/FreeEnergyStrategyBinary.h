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
#ifndef included_FreeEnergyStrategyBinary
#define included_FreeEnergyStrategyBinary

#include "ConcFreeEnergyStrategy.h"
#include "FuncFort.h"
#include "InterpolationType.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"

#include <boost/property_tree/ptree.hpp>

#include <string>
#include <vector>

/*!
 * Implements free energy and driving force components that do not depend
 * on specific free energy form
 */
class FreeEnergyStrategyBinary : public ConcFreeEnergyStrategy
{
 public:
   FreeEnergyStrategyBinary(
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       const int conc_l_id, const int conc_a_id, const int conc_b_id,
       const bool with_third_phase);

   virtual ~FreeEnergyStrategyBinary(){};

   virtual void addDrivingForce(const double time, hier::Patch& patch,
                                const int temperature_id, const int phase_id,
                                const int conc_id, const int f_l_id,
                                const int f_a_id, const int f_b_id,
                                const int rhs_id);

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

         default:
            tbox::pout << "undefined phase=" << phase << "!!!" << std::endl;
            tbox::SAMRAI_MPI::abort();
      }
   }

   virtual void addDrivingForce(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_f_b,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b, const hier::Box& pbox);

   virtual double computeMuA(const double t, const double c) = 0;
   virtual double computeMuL(const double t, const double c) = 0;
   // virtual double computeMuB(const double t, const double c) = 0;

   void computeDrivingForce(const double time, hier::Patch& patch,
                            const int temperature_id, const int phase_id,
                            const int conc_id, const int f_l_id,
                            const int f_a_id, const int f_b_id,
                            const int rhs_id);

 protected:
   Thermo4PFM::EnergyInterpolationType d_energy_interp_func_type;
   Thermo4PFM::ConcInterpolationType d_conc_interp_func_type;

   int d_conc_l_id;
   int d_conc_a_id;
   int d_conc_b_id;

   bool d_with_third_phase;

   virtual double computeFreeEnergy(const double t, double* c,
                                    const Thermo4PFM::PhaseIndex pi,
                                    const bool gp) = 0;
   virtual double computeDerivFreeEnergy(const double t, double* c,
                                         const Thermo4PFM::PhaseIndex pi) = 0;

   virtual void computeSecondDerivativeEnergyPhaseL(
       const double temp, const std::vector<double>& c_l,
       std::vector<double>& d2fdc2, const bool use_internal_units) = 0;
   virtual void computeSecondDerivativeEnergyPhaseA(
       const double temp, const std::vector<double>& c_l,
       std::vector<double>& d2fdc2, const bool use_internal_units) = 0;

 private:
   double hprime(const double phi)
   {
      const char interp =
          Thermo4PFM::energyInterpChar(d_energy_interp_func_type);
      return DERIV_INTERP_FUNC(phi, &interp);
   }

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
