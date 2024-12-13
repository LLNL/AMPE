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
#include "DeltaTemperatureFreeEnergyStrategy.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"

#include <cassert>

using namespace SAMRAI;


extern "C" {

void COMPUTEDPHIDTEMPERATUREDELTATEMPERATURE(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1,
#if (NDIM == 3)
    const int& ifirst2, const int& ilast2,
#endif
    const double*, const int&, const double&, const double&, const double&,
    double* rhs, const int&, const char* const);
}

DeltaTemperatureFreeEnergyStrategy::DeltaTemperatureFreeEnergyStrategy(
    const double Tm, const double latent_heat,
    const Thermo4PFM::EnergyInterpolationType phase_interp_func_type)
    : d_Tm(Tm),
      d_L(latent_heat),
      d_phase_interp_func_type(phase_interp_func_type)
{
   tbox::plog << "DeltaTemperatureFreeEnergyStrategy..." << std::endl;
   tbox::plog << "Tm=" << d_Tm << std::endl;
   tbox::plog << "L=" << d_L << std::endl;

   assert(d_L == d_L);
   assert(d_Tm == d_Tm);
}

//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseL(
    const double temp, const std::vector<double>& c_l,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;
   (void)c_l;
   (void)use_internal_units;

   d2fdc2.assign(d2fdc2.size(), 0.);
}

//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseA(
    const double temp, const std::vector<double>& c_a,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;
   (void)c_a;
   (void)use_internal_units;

   d2fdc2.assign(d2fdc2.size(), 0.);
}

void DeltaTemperatureFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseB(
    const double temp, const std::vector<double>& c_b,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;
   (void)c_b;
   (void)use_internal_units;

   d2fdc2.assign(d2fdc2.size(), 0.);
}

void DeltaTemperatureFreeEnergyStrategy::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int conc_id, const int f_l_id, const int f_a_id,
    const int f_b_id, const int rhs_id)
{
   assert(phase_id >= 0);
   assert(rhs_id >= 0);
   assert(temperature_id >= 0);
   assert(d_Tm > 0.);
   assert(d_L > 0.);

   (void)time;
   (void)conc_id;  // unused
   (void)f_l_id;   // unused
   (void)f_a_id;   // unused
   (void)f_b_id;   // unused

   // tbox::pout<<"DeltaTemperatureFreeEnergyStrategy::addDrivingForce()..."
   //           <<endl;

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);

   std::shared_ptr<pdat::CellData<double> > temp(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(temp);

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

   COMPUTERHSDELTATEMPERATURE(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                              ifirst(2), ilast(2),
#endif
                              phase->getPointer(),
                              phase->getGhostCellWidth()[0], temp->getPointer(),
                              temp->getGhostCellWidth()[0], d_Tm, d_L,
                              rhs->getPointer(), rhs->getGhostCellWidth()[0],
                              &interpf);
}
//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::applydPhidTBlock(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id, const int rhs_id,
    const double phase_mobility)
{
   const char interpf = energyInterpChar(d_phase_interp_func_type);

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& box = patch->getBox();
         const hier::Index& ifirst = box.lower();
         const hier::Index& ilast = box.upper();

         std::shared_ptr<pdat::CellData<double> > temp(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(temperature_id)));
         std::shared_ptr<pdat::CellData<double> > phase(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));
         std::shared_ptr<pdat::CellData<double> > rhs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(rhs_id)));

         COMPUTEDPHIDTEMPERATUREDELTATEMPERATURE(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             phase->getPointer(), phase->getGhostCellWidth()[0], d_Tm, d_L,
             phase_mobility, rhs->getPointer(), rhs->getGhostCellWidth()[0],
             &interpf);
      }
   }
}

//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::computeFreeEnergyLiquid(
    hier::Patch& patch, const int temperature_id, const int fl_id,
    const bool gp)
{
   (void)gp;

   const hier::Box& box = patch.getBox();
   const hier::Index& ifirst = box.lower();
   const hier::Index& ilast = box.upper();

   std::shared_ptr<pdat::CellData<double> > temp(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   std::shared_ptr<pdat::CellData<double> > fl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fl_id)));

   double factor = 0.5 * d_L;

   TEMPERATURE_ENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                      ifirst(2), ilast(2),
#endif
                      temp->getPointer(), temp->getGhostCellWidth()[0],
                      fl->getPointer(), d_Tm, factor);
}

//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::computeFreeEnergySolidA(
    hier::Patch& patch, const int temperature_id, const int fa_id,
    const bool gp)
{
   (void)gp;

   const hier::Box& box = patch.getBox();
   const hier::Index& ifirst = box.lower();
   const hier::Index& ilast = box.upper();

   std::shared_ptr<pdat::CellData<double> > temp(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));

   std::shared_ptr<pdat::CellData<double> > fa(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fa_id)));

   double factor = -0.5 * d_L;

   TEMPERATURE_ENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                      ifirst(2), ilast(2),
#endif
                      temp->getPointer(), temp->getGhostCellWidth()[0],
                      fa->getPointer(), d_Tm, factor);
}

void DeltaTemperatureFreeEnergyStrategy::computeFreeEnergySolidB(
    hier::Patch& patch, const int temperature_id, const int fb_id,
    const bool gp)
{
   computeFreeEnergySolidA(patch, temperature_id, fb_id, gp);
}
