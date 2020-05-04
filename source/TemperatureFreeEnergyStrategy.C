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
#include "TemperatureFreeEnergyStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "ConcFort.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include <cassert>

using namespace SAMRAI;


TemperatureFreeEnergyStrategy::TemperatureFreeEnergyStrategy(
    const EnergyInterpolationType phase_interp_func_type,
    const EnergyInterpolationType eta_interp_func_type, const double fa,
    const double fb, const double vma, const double vmb,
    const double latent_heat, const double meltingT,
    const bool with_third_phase)
{
   assert(meltingT > 0.);
   assert(meltingT < 100000.);
   assert(latent_heat > 0.);

   d_phase_interp_func_type = phase_interp_func_type;
   d_eta_interp_func_type = eta_interp_func_type;
   d_latent_heat = latent_heat;
   d_meltingT = meltingT;
   d_invMeltingT = 1. / d_meltingT;
   d_f_a = fa;
   d_f_b = fb;

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;
   d_f_a *= 1.e-6 / vma;
   d_f_b *= 1.e-6 / vmb;

   d_with_third_phase = with_third_phase;

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
    const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fl_id, const bool gp)
{
   assert(fl_id >= 0);
   (void)temperature_id;  // unused
   (void)gp;

   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
      boost::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ip++) {
         boost::shared_ptr<hier::Patch> patch = *ip;
         const hier::Box& pbox = patch->getBox();

         boost::shared_ptr<pdat::CellData<double> > fl(
             BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(fl_id)));
         assert(fl);
         boost::shared_ptr<pdat::CellData<double> > temperature(
             BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(temperature_id)));
         assert(temperature);

         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend;
              ++i) {
            pdat::CellIndex cell = *i;
            (*fl)(cell) = computeValFreeEnergyLiquid((*temperature)(cell), 0.);
         }
      }
   }
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

   boost::shared_ptr<pdat::CellData<double> > fl(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fl_id)));
   assert(fl);
   boost::shared_ptr<pdat::CellData<double> > temperature(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
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
    const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fs_id, const bool gp)
{
   assert(fs_id >= 0);
   (void)temperature_id;  // unused
   (void)gp;

   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
      boost::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ip++) {
         boost::shared_ptr<hier::Patch> patch = *ip;

         boost::shared_ptr<pdat::CellData<double> > fs(
             BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(fs_id)));
         assert(fs);

         fs->fillAll(d_f_a);
      }
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

   boost::shared_ptr<pdat::CellData<double> > fs(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fs_id)));
   assert(fs);

   fs->fillAll(d_f_a);
}

//=======================================================================

void TemperatureFreeEnergyStrategy::computeFreeEnergySolidB(
    const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fs_id, const bool gp)
{
   assert(fs_id >= 0);
   (void)temperature_id;  // unused
   (void)gp;

   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
      boost::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ++ip) {
         boost::shared_ptr<hier::Patch> patch = *ip;

         boost::shared_ptr<pdat::CellData<double> > fs(
             BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(fs_id)));
         assert(fs);

         fs->fillAll(d_f_b);
      }
   }
}

//=======================================================================

void TemperatureFreeEnergyStrategy::computeFreeEnergySolidB(
    hier::Patch& patch, const int temperature_id, const int fs_id,
    const bool gp)
{
   assert(fs_id >= 0);
   (void)temperature_id;  // unused
   (void)gp;

   boost::shared_ptr<pdat::CellData<double> > fs(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fs_id)));
   assert(fs);

   fs->fillAll(d_f_b);
}

//=======================================================================

void TemperatureFreeEnergyStrategy::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int eta_id, const int conc_id, const int f_l_id,
    const int f_a_id, const int f_b_id, const int rhs_id)
{
   assert(phase_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(rhs_id >= 0);
   if (d_with_third_phase) {
      assert(eta_id >= 0);
   }
   (void)time;
   (void)conc_id;         // unused
   (void)temperature_id;  // unused

   boost::shared_ptr<pdat::CellData<double> > phase(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);

   boost::shared_ptr<pdat::CellData<double> > fl(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_l_id)));
   assert(fl);

   boost::shared_ptr<pdat::CellData<double> > fa(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_a_id)));
   assert(fa);

   boost::shared_ptr<pdat::CellData<double> > rhs(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(rhs_id)));
   assert(rhs);

   assert(rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   int three_phase = 0;
   double* ptr_fb = nullptr;
   double* ptr_eta = nullptr;
   if (d_with_third_phase) {
      three_phase = 1;
      boost::shared_ptr<pdat::CellData<double> > eta(
          BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(eta_id)));
      ptr_eta = eta->getPointer();

      boost::shared_ptr<pdat::CellData<double> > fb(
          BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(f_b_id)));
      ptr_fb = fb->getPointer();
   }

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   const char interpf = energyInterpChar(d_phase_interp_func_type);
   const char interpe = energyInterpChar(d_eta_interp_func_type);

   FORT_PHASERHS_FENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                         ifirst(2), ilast(2),
#endif
                         fl->getPointer(), fa->getPointer(), ptr_fb,
                         phase->getPointer(), phase->getGhostCellWidth()[0],
                         ptr_eta, phase->getGhostCellWidth()[0],
                         rhs->getPointer(), 0, &interpf, &interpe, three_phase);
}

//=======================================================================

void TemperatureFreeEnergyStrategy::addDrivingForceEta(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int eta_id, const int conc_id, const int f_l_id,
    const int f_a_id, const int f_b_id, const int rhs_id)
{
   assert(phase_id >= 0);
   assert(eta_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(f_b_id >= 0);
   assert(rhs_id >= 0);
   (void)time;
   (void)conc_id;         // unused
   (void)temperature_id;  // unused

   boost::shared_ptr<pdat::CellData<double> > eta(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(eta_id)));
   assert(eta);

   boost::shared_ptr<pdat::CellData<double> > phase(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);

   boost::shared_ptr<pdat::CellData<double> > fl(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_l_id)));
   assert(fl);

   boost::shared_ptr<pdat::CellData<double> > fa(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_a_id)));
   assert(fa);

   boost::shared_ptr<pdat::CellData<double> > fb(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_b_id)));
   assert(fb);

   boost::shared_ptr<pdat::CellData<double> > rhs(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(rhs_id)));
   assert(rhs);

   assert(rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   const char interpf = energyInterpChar(d_phase_interp_func_type);
   const char interpe = energyInterpChar(d_eta_interp_func_type);

   FORT_ETARHS_FENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                       ifirst(2), ilast(2),
#endif
                       fl->getPointer(), fa->getPointer(), fb->getPointer(),
                       phase->getPointer(), phase->getGhostCellWidth()[0],
                       eta->getPointer(), eta->getGhostCellWidth()[0],
                       rhs->getPointer(), 0, &interpf, &interpe);
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
