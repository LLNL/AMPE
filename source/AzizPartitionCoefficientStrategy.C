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

#include "AzizPartitionCoefficientStrategy.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "Phases.h"

template <>
AzizPartitionCoefficientStrategy<CALPHADFreeEnergyFunctionsBinary>::
    AzizPartitionCoefficientStrategy(
        const int velocity_id, const int temperature_id,
        const int partition_coeff_id, const double vd, const double keq,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db)
    : PartitionCoefficientStrategy(velocity_id, temperature_id,
                                   partition_coeff_id),
      d_inv_vd(1. / vd),
      d_keq(keq)
{
   assert(d_inv_vd == d_inv_vd);

   std::shared_ptr<tbox::Database> conc_calphad_db =
       conc_db->getDatabase("Calphad");
   std::string calphad_filename = conc_calphad_db->getString("filename");
   std::shared_ptr<tbox::MemoryDatabase> calphad_db(
       new tbox::MemoryDatabase("calphad_db"));
   tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                    calphad_db);

   std::shared_ptr<tbox::Database> newton_db;
   if (conc_db->isDatabase("NewtonSolver"))
      newton_db = conc_db->getDatabase("NewtonSolver");
#ifdef HAVE_THERMO4PFM
   pt::ptree calphad_pt;
   copyDatabase(calphad_db, calphad_pt);
   pt::ptree newton_pt;
   copyDatabase(newton_db, newton_pt);
#endif

   d_fenergy = std::unique_ptr<CALPHADFreeEnergyFunctionsBinary>(
       new CALPHADFreeEnergyFunctionsBinary(
#ifdef HAVE_THERMO4PFM
           calphad_pt, newton_pt,
#else
           calphad_db, newton_db,
#endif
           energy_interp_func_type, conc_interp_func_type
#ifndef HAVE_THERMO4PFM
           ,
           false  // no 3rd phase
#endif
           ));
}

template <>
AzizPartitionCoefficientStrategy<CALPHADFreeEnergyFunctionsTernary>::
    AzizPartitionCoefficientStrategy(
        const int velocity_id, const int temperature_id,
        const int partition_coeff_id, const double vd, const double keq,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db)
    : PartitionCoefficientStrategy(velocity_id, temperature_id,
                                   partition_coeff_id),
      d_inv_vd(1. / vd),
      d_keq(keq)
{
   assert(d_inv_vd == d_inv_vd);

   std::shared_ptr<tbox::Database> conc_calphad_db =
       conc_db->getDatabase("Calphad");
   std::string calphad_filename = conc_calphad_db->getString("filename");
   std::shared_ptr<tbox::MemoryDatabase> calphad_db(
       new tbox::MemoryDatabase("calphad_db"));
   tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                    calphad_db);

   std::shared_ptr<tbox::Database> newton_db;
   if (conc_db->isDatabase("NewtonSolver"))
      newton_db = conc_db->getDatabase("NewtonSolver");
#ifdef HAVE_THERMO4PFM
   pt::ptree calphad_pt;
   copyDatabase(calphad_db, calphad_pt);
   pt::ptree newton_pt;
   copyDatabase(newton_db, newton_pt);
#endif

   d_fenergy = std::unique_ptr<CALPHADFreeEnergyFunctionsTernary>(
       new CALPHADFreeEnergyFunctionsTernary(
#ifdef HAVE_THERMO4PFM
           calphad_pt, newton_pt,
#else
           calphad_db, newton_db,
#endif
           energy_interp_func_type, conc_interp_func_type));
}

template <>
AzizPartitionCoefficientStrategy<KKSFreeEnergyFunctionDiluteBinary>::
    AzizPartitionCoefficientStrategy(
        const int velocity_id, const int temperature_id,
        const int partition_coeff_id, const double vd, const double keq,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db)
    : PartitionCoefficientStrategy(velocity_id, temperature_id,
                                   partition_coeff_id),
      d_inv_vd(1. / vd),
      d_keq(keq)
{
   assert(d_inv_vd == d_inv_vd);

#ifdef HAVE_THERMO4PFM
   pt::ptree conc_pt;
   copyDatabase(conc_db, conc_pt);
   d_fenergy = std::unique_ptr<KKSFreeEnergyFunctionDiluteBinary>(
       new KKSFreeEnergyFunctionDiluteBinary(conc_pt, energy_interp_func_type,
                                             conc_interp_func_type));
#else
   d_fenergy = std::unique_ptr<KKSFreeEnergyFunctionDiluteBinary>(
       new KKSFreeEnergyFunctionDiluteBinary(conc_db, energy_interp_func_type,
                                             conc_interp_func_type));
#endif
}

template <class FreeEnergyType>
void AzizPartitionCoefficientStrategy<FreeEnergyType>::evaluate(
    hier::Patch& patch, std::shared_ptr<pdat::CellData<double> > cd_velocity,
    std::shared_ptr<pdat::CellData<double> > cd_temperature,
    std::shared_ptr<pdat::CellData<double> > cd_partition_coeff)
{
   const hier::Box& patch_box = patch.getBox();

   pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
   for (pdat::CellIterator i(pdat::CellGeometry::begin(patch_box)); i != iend;
        ++i) {
      const pdat::CellIndex ccell = *i;

      double temperature = (*cd_temperature)(ccell);
      double vel = fabs((*cd_velocity)(ccell));
      assert(vel == vel);

      const double ke = computeKeq(temperature);
      const double scaledv = vel * d_inv_vd;
      if (ke < 1.)
         (*cd_partition_coeff)(ccell) = (ke + scaledv) / (1. + scaledv);
      else
         (*cd_partition_coeff)(ccell) = (1. + scaledv) / (1. / ke + scaledv);

      assert((*cd_partition_coeff)(ccell) > 0.);
   }
}

template <class FreeEnergyType>
double AzizPartitionCoefficientStrategy<FreeEnergyType>::computeKeq(
    const double temperature)
{
   if (d_keq > 0.) return d_keq;

   double ceq[2] = {0.5, 0.5};
   bool flag = d_fenergy->computeCeqT(temperature,
#ifndef HAVE_THERMO4PFM
                                      PhaseIndex::phaseL, PhaseIndex::phaseA,
#endif
                                      &ceq[0], 20);

   assert(flag);

   return ceq[1] / ceq[0];
}

template class AzizPartitionCoefficientStrategy<
    CALPHADFreeEnergyFunctionsBinary>;
