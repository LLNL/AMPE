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
#ifndef included_HBSMFreeEnergyStrategy
#define included_HBSMFreeEnergyStrategy

#include "FreeEnergyStrategy.h"
#include "InterpolationType.h"

#include <string>

class HBSMFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   HBSMFreeEnergyStrategy(std::shared_ptr<tbox::Database> input_db,
                          const EnergyInterpolationType energy_interp_func_type,
                          const double vml, const double vma, const double vmb,
                          const double D_liquid, const double D_solid_A,
                          const double D_solid_B, const double Q0_liquid,
                          const double Q0_solid_A, const double Q0_solid_B,
                          const int conc_l_id, const int conc_a_id,
                          const int conc_b_id, const bool with_third_phase);

   ~HBSMFreeEnergyStrategy(){};

   //   double computeValFreeEnergyLiquid(
   //      const double temperature, const double conc,
   //      const bool gp = false );

   //   double computeValFreeEnergySolidA(
   //      const double temperature, const double conc,
   //      const bool gp = false );

   //   double computeValFreeEnergySolidB(
   //      const double temperature, const double conc,
   //      const bool gp = false );

   void computeFreeEnergyLiquid(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fl_id, const bool gp = false);

   void computeDerivFreeEnergyLiquid(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fl_id);

   void computeFreeEnergySolidA(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fs_id, const bool gp = false);

   void computeDerivFreeEnergySolidA(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fs_id);

   void computeFreeEnergySolidB(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fs_id, const bool gp = false);

   void computeDerivFreeEnergySolidB(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fs_id);

   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp = false);

   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp = false);

   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp = false);

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int eta_id, const int conc_id, const int f_l_id,
                        const int f_a_id, const int f_b_id, const int rhs_id);

   void addDrivingForceEta(const double time, hier::Patch& patch,
                           const int temperature_id, const int phase_id,
                           const int eta_id, const int conc_id,
                           const int f_l_id, const int f_a_id, const int f_b_id,
                           const int rhs_id);

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);

   double computeLiquidConcentration(const double hphi, const double heta,
                                     const double c) const;

   double computeSolidAConcentration(const double hphi, const double heta,
                                     const double c) const;

   double computeSolidBConcentration(const double hphi, const double heta,
                                     const double cB) const;

 private:
   double computeFreeEnergyPrivate(const double temperature, const double conc,
                                   const double A, const double Ceq,
                                   const double energy_factor,
                                   const bool gp = false) const;

   void computeFreeEnergyPrivate(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const double A, const double Ceq,
       const int f_id, const int c_i_id, const double energy_factor);

   void computeFreeEnergyPrivate(hier::Patch& patch, const int temperature_id,
                                 const double A, const double Ceq,
                                 const int f_id, const int c_i_id,
                                 const double energy_factor);

   void computeDerivFreeEnergyPrivate(hier::Patch& patch,
                                      const int temperature_id, const double A,
                                      const double Ceq, const int f_id,
                                      const int c_i_id,
                                      const double energy_factor);

   double computeLiquidConcentration(const double hphi, const double heta,
                                     const double c, const double Al,
                                     const double Aa, const double Ab,
                                     const double Ceql, const double CeqA,
                                     const double CeqB) const;

   double computeSolidAConcentration(const double hphi, const double heta,
                                     const double c, const double Al,
                                     const double Aa, const double Ab,
                                     const double Ceql, const double CeqA,
                                     const double CeqB) const;

   double computeSolidBConcentration(const double hphi, const double heta,
                                     const double c, const double Al,
                                     const double Aa, const double Ab,
                                     const double Ceql, const double CeqA,
                                     const double CeqB) const;

   void computeFreeEnergyPrivatePatch(
       const hier::Box& pbox,
       std::shared_ptr<pdat::CellData<double> > cd_temp, const double A,
       const double Ceq,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       const double energy_factor);

   void computeDerivFreeEnergyPrivatePatch(
       const hier::Box& pbox,
       std::shared_ptr<pdat::CellData<double> > cd_temp, const double A,
       const double Ceq,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       const double energy_factor);

   void addDrivingForceOnPatchPrivate(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_f_b,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       const hier::Box& pbox);

   void addDrivingForceEtaOnPatchPrivate(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_f_b,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       const hier::Box& pbox);

   double computeMu(const double t, const double c);

   double computeDerivFreeEnergyPrivate(const double temperature,
                                        const double conc, const double A,
                                        const double Ceq,
                                        const double energy_factor) const;

   EnergyInterpolationType d_energy_interp_func_type;

   double d_vm_L;  // molar volume
   double d_vm_A;  // molar volume
   double d_vm_B;  // molar volume
   // double d_vm; // molar volume
   // double d_jpmol2pjpmumcube;

   double d_energy_conv_factor_L;  // molar volume
   double d_energy_conv_factor_A;  // molar volume
   double d_energy_conv_factor_B;  // molar volume

   double d_D_liquid;
   double d_D_solid_A;
   double d_D_solid_B;
   double d_Q0_liquid;
   double d_Q0_solid_A;
   double d_Q0_solid_B;

   double d_A_liquid;
   double d_A_solid_A;
   double d_A_solid_B;
   double d_Ceq_liquid;
   double d_Ceq_solid_A;
   double d_Ceq_solid_B;

   bool d_with_third_phase;

   int d_conc_l_id;
   int d_conc_a_id;
   int d_conc_b_id;
};

#endif
