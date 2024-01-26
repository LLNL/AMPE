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
#ifndef included_KKSdiluteBinary
#define included_KKSdiluteBinary

#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "FreeEnergyStrategy.h"
#include "InterpolationType.h"
#include "Phases.h"
#include "FuncFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
class MolarVolumeStrategy;

#include <string>
#include <vector>

class KKSdiluteBinary : public FreeEnergyStrategy
{
 public:
   KKSdiluteBinary(std::shared_ptr<tbox::Database> conc_db,
                   const EnergyInterpolationType energy_interp_func_type,
                   const ConcInterpolationType conc_interp_func_type,
                   MolarVolumeStrategy* mvstrategy, const int conc_l_id,
                   const int conc_a_id);

   ~KKSdiluteBinary(){};

   virtual void setup(std::shared_ptr<tbox::Database> calphad_db);

   void computeDerivFreeEnergyLiquid(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fl_id);

   void computeDerivFreeEnergySolidA(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fs_id);

   void computeDerivFreeEnergySolidB(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fs_id);

   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp);

   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp);

   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp);

   virtual void addDrivingForce(const double time, hier::Patch& patch,
                                const int temperature_id, const int phase_id,
                                const int eta_id, const int conc_id,
                                const int f_l_id, const int f_a_id,
                                const int f_b_id, const int rhs_id);

   void computeDrivingForce(const double time, hier::Patch& patch,
                            const int temperature_id, const int phase_id,
                            const int eta_id, const int conc_id,
                            const int f_l_id, const int f_a_id,
                            const int f_b_id, const int rhs_id);

   void addDrivingForceEta(const double time, hier::Patch& patch,
                           const int temperature_id, const int phase_id,
                           const int eta_id, const int conc_id,
                           const int f_l_id, const int f_a_id, const int f_b_id,
                           const int rhs_id);

   virtual void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   virtual void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
#ifndef HAVE_THERMO4PFM
   virtual void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
#endif
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

#ifndef HAVE_THERMO4PFM
         case 'b':
            computeSecondDerivativeEnergyPhaseB(temp, c, d2fdc2,
                                                use_internal_units);
            break;
#endif

         default:
            tbox::pout << "undefined phase=" << phase << "!!!" << std::endl;
            tbox::SAMRAI_MPI::abort();
      }
   }

   void preRunDiagnostics(const double temperature)
   {
      d_kksdilute_fenergy->preRunDiagnostics(temperature);
   }

 private:
   EnergyInterpolationType d_energy_interp_func_type;
   ConcInterpolationType d_conc_interp_func_type;

   int d_conc_l_id;
   int d_conc_a_id;

   double computeMuA(const double t, const double c);

   double computeMuL(const double t, const double c);

   void defaultComputeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units);

   void defaultComputeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units);

   void defaultComputeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units);

   MolarVolumeStrategy* d_mv_strategy;

   KKSFreeEnergyFunctionDiluteBinary* d_kksdilute_fenergy;

   double hprime(const double phi)
   {
      const char interp = energyInterpChar(d_energy_interp_func_type);
      return DERIV_INTERP_FUNC(phi, &interp);
   }

   void addDrivingForceOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a, const hier::Box& pbox);

   void computeFreeEnergy(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                          const int temperature_id, const int f_id,
                          const PhaseIndex pi, const bool gp);

   void computeDerivFreeEnergy(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_id, const PhaseIndex pi);

   void computeFreeEnergy(hier::Patch& patch, const int temperature_id,
                          const int f_id, const PhaseIndex pi, const bool gp);

   void computeDerivFreeEnergy(hier::Patch& patch, const int temperature_id,
                               const int f_id, const PhaseIndex pi);

   void computeFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i, const PhaseIndex pi,
       const bool gp);

   void computeDerivFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i, const PhaseIndex pi);

   void addDrivingForceEtaOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a, const hier::Box& pbox);
};

#endif
