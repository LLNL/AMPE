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
//

#include "AzizPartitionCoefficientStrategy.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "Phases.h"

#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/InputManager.h"

#include "Database2JSON.h"
namespace pt = boost::property_tree;

template <>
AzizPartitionCoefficientStrategy<Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>::
    AzizPartitionCoefficientStrategy(
        const int velocity_id, const int temperature_id,
        const int partition_coeff_id, const double vd, const double keq,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
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
   pt::ptree calphad_pt;
   copyDatabase(calphad_db, calphad_pt);
   pt::ptree newton_pt;
   copyDatabase(newton_db, newton_pt);

   d_fenergy = std::unique_ptr<Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>(
       new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary(calphad_pt, newton_pt,
                                                        energy_interp_func_type,
                                                        conc_interp_func_type));
}

template <>
AzizPartitionCoefficientStrategy<
    Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>::
    AzizPartitionCoefficientStrategy(
        const int velocity_id, const int temperature_id,
        const int partition_coeff_id, const double vd, const double keq,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
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
   pt::ptree calphad_pt;
   copyDatabase(calphad_db, calphad_pt);
   pt::ptree newton_pt;
   copyDatabase(newton_db, newton_pt);

   d_fenergy = std::unique_ptr<Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>(
       new Thermo4PFM::CALPHADFreeEnergyFunctionsTernary(
           calphad_pt, newton_pt, energy_interp_func_type,
           conc_interp_func_type));
}

template <>
AzizPartitionCoefficientStrategy<
    Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary>::
    AzizPartitionCoefficientStrategy(
        const int velocity_id, const int temperature_id,
        const int partition_coeff_id, const double vd, const double keq,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db)
    : PartitionCoefficientStrategy(velocity_id, temperature_id,
                                   partition_coeff_id),
      d_inv_vd(1. / vd),
      d_keq(keq)
{
   assert(d_inv_vd == d_inv_vd);

   pt::ptree conc_pt;
   copyDatabase(conc_db, conc_pt);
   d_fenergy = std::unique_ptr<Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary>(
       new Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary(
           conc_pt, energy_interp_func_type, conc_interp_func_type));
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
   bool flag = d_fenergy->computeCeqT(temperature, &ceq[0], 20);

   assert(flag);

   if (flag)
      return ceq[1] / ceq[0];
   else
      return -1.;
}

template class AzizPartitionCoefficientStrategy<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>;
