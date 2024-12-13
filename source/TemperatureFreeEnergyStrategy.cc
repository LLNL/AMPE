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
#include "TemperatureFreeEnergyStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "ConcFort.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include <cassert>

using namespace SAMRAI;


TemperatureFreeEnergyStrategy::TemperatureFreeEnergyStrategy(
    const Thermo4PFM::EnergyInterpolationType phase_interp_func_type,
    const double fa, const double vma, const double latent_heat,
    const double meltingT)
{
   assert(meltingT > 0.);
   assert(meltingT < 100000.);
   assert(latent_heat > 0.);

   d_phase_interp_func_type = phase_interp_func_type;
   d_latent_heat = latent_heat;
   d_meltingT = meltingT;
   d_invMeltingT = 1. / d_meltingT;
   d_f_a = fa;

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;
   d_f_a *= 1.e-6 / vma;

   tbox::plog << "TemperatureFreeEnergyStrategy:" << std::endl;
   tbox::plog << "Molar volume A = " << vma << std::endl;
   tbox::plog << "Latent heat    = " << d_latent_heat << std::endl;
   tbox::plog << "meltingT       = " << d_meltingT << std::endl;
   tbox::plog << "f solid A      = " << d_f_a << std::endl;
}

double TemperatureFreeEnergyStrategy::computeValFreeEnergyLiquid(
    const double temperature, const double conc, const bool gp)
{
   (void)conc;  // unused
   (void)gp;
   return d_f_a + d_latent_heat * (d_meltingT - temperature) * d_invMeltingT;
}

//=======================================================================

void TemperatureFreeEnergyStrategy::computeFreeEnergyLiquid(
    hier::Patch& patch, const int temperature_id, const int fl_id,
    const bool gp)
{
   assert(fl_id >= 0);
   (void)temperature_id;  // unused
   (void)gp;

   const hier::Box& pbox = patch.getBox();

   std::shared_ptr<pdat::CellData<double> > fl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fl_id)));
   assert(fl);
   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(temperature);

   pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
   for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend; ++i) {
      pdat::CellIndex cell = *i;
      (*fl)(cell) = computeValFreeEnergyLiquid((*temperature)(cell), 0.);
   }
}

//=======================================================================

void TemperatureFreeEnergyStrategy::computeFreeEnergySolidA(
    hier::Patch& patch, const int temperature_id, const int fs_id,
    const bool gp)
{
   assert(fs_id >= 0);
   (void)temperature_id;  // unused
   (void)gp;

   std::shared_ptr<pdat::CellData<double> > fs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fs_id)));
   assert(fs);

   fs->fillAll(d_f_a);
}

//=======================================================================

void TemperatureFreeEnergyStrategy::computeFreeEnergySolidB(
    hier::Patch& patch, const int temperature_id, const int fs_id,
    const bool gp)
{
   (void)temperature_id;  // unused
   (void)gp;
}

//=======================================================================

void TemperatureFreeEnergyStrategy::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int conc_id, const int f_l_id, const int f_a_id,
    const int f_b_id, const int rhs_id)
{
   assert(phase_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(rhs_id >= 0);

   (void)time;
   (void)conc_id;         // unused
   (void)temperature_id;  // unused
   (void)f_b_id;

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);

   std::shared_ptr<pdat::CellData<double> > fl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_l_id)));
   assert(fl);

   std::shared_ptr<pdat::CellData<double> > fa(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_a_id)));
   assert(fa);

   std::shared_ptr<pdat::CellData<double> > rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(rhs_id)));
   assert(rhs);

   assert(rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   const char interpf = energyInterpChar(d_phase_interp_func_type);

   PHASERHS_FENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                    ifirst(2), ilast(2),
#endif
                    fl->getPointer(), fa->getPointer(), phase->getPointer(),
                    phase->getGhostCellWidth()[0], rhs->getPointer(), 0,
                    &interpf);
}

//=======================================================================

void TemperatureFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseL(
    const double temp, const std::vector<double>& c_l,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;
   (void)c_l;
   (void)d2fdc2;
   (void)use_internal_units;

   return;
}

//=======================================================================

void TemperatureFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseA(
    const double temp, const std::vector<double>& c_a,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;
   (void)c_a;
   (void)d2fdc2;
   (void)use_internal_units;

   return;
}

//=======================================================================

void TemperatureFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseB(
    const double temp, const std::vector<double>& c_b,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;
   (void)c_b;
   (void)d2fdc2;
   (void)use_internal_units;

   return;
}
