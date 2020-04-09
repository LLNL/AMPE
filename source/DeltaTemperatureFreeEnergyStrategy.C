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
#include "DeltaTemperatureFreeEnergyStrategy.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"

#include <cassert>

using namespace SAMRAI;
using namespace std;

#define FORT_COMP_DPHIDTEMPERATURE_DELTA_TEMPERATURE \
   computedphidtemperaturedeltatemperature_

extern "C" {

void FORT_COMP_DPHIDTEMPERATURE_DELTA_TEMPERATURE(
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
    const EnergyInterpolationType phase_interp_func_type)
    : d_Tm(Tm),
      d_L(latent_heat),
      d_phase_interp_func_type(phase_interp_func_type)
{
   tbox::plog << "DeltaTemperatureFreeEnergyStrategy..." << endl;
   tbox::plog << "Tm=" << d_Tm << endl;
   tbox::plog << "L=" << d_L << endl;

   assert(d_L == d_L);
   assert(d_Tm == d_Tm);
}

//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseL(
    const double temp, const vector<double>& c_l, vector<double>& d2fdc2,
    const bool use_internal_units)
{
   (void)temp;
   (void)c_l;
   (void)use_internal_units;

   d2fdc2.assign(d2fdc2.size(), 0.);
}

//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseA(
    const double temp, const vector<double>& c_a, vector<double>& d2fdc2,
    const bool use_internal_units)
{
   (void)temp;
   (void)c_a;
   (void)use_internal_units;

   d2fdc2.assign(d2fdc2.size(), 0.);
}

void DeltaTemperatureFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseB(
    const double temp, const vector<double>& c_b, vector<double>& d2fdc2,
    const bool use_internal_units)
{
   (void)temp;
   (void)c_b;
   (void)use_internal_units;

   d2fdc2.assign(d2fdc2.size(), 0.);
}

void DeltaTemperatureFreeEnergyStrategy::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int eta_id, const int conc_id, const int f_l_id,
    const int f_a_id, const int f_b_id, const int rhs_id)
{
   assert(phase_id >= 0);
   assert(rhs_id >= 0);
   assert(temperature_id >= 0);
   assert(d_Tm > 0.);
   assert(d_L > 0.);

   (void)time;
   (void)eta_id;
   (void)conc_id;  // unused
   (void)f_l_id;   // unused
   (void)f_a_id;   // unused
   (void)f_b_id;   // unused

   // tbox::pout<<"DeltaTemperatureFreeEnergyStrategy::addDrivingForce()..."
   //           <<endl;

   boost::shared_ptr<pdat::CellData<double> > phase(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);

   boost::shared_ptr<pdat::CellData<double> > temp(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(temp);

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

   FORT_COMP_RHS_DELTA_TEMPERATURE(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                                   ifirst(2), ilast(2),
#endif
                                   phase->getPointer(),
                                   phase->getGhostCellWidth()[0],
                                   temp->getPointer(),
                                   temp->getGhostCellWidth()[0], d_Tm, d_L,
                                   rhs->getPointer(),
                                   rhs->getGhostCellWidth()[0], &interpf);
}
//=======================================================================

void DeltaTemperatureFreeEnergyStrategy::applydPhidTBlock(
    const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id, const int rhs_id,
    const double phase_mobility)
{
   const char interpf = energyInterpChar(d_phase_interp_func_type);

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         boost::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& box = patch->getBox();
         const hier::Index& ifirst = box.lower();
         const hier::Index& ilast = box.upper();

         boost::shared_ptr<pdat::CellData<double> > temp(
             BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(temperature_id)));
         boost::shared_ptr<pdat::CellData<double> > phase(
             BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));
         boost::shared_ptr<pdat::CellData<double> > rhs(
             BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(rhs_id)));

         FORT_COMP_DPHIDTEMPERATURE_DELTA_TEMPERATURE(
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

   boost::shared_ptr<pdat::CellData<double> > temp(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   boost::shared_ptr<pdat::CellData<double> > fl(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fl_id)));

   double factor = 0.5 * d_L;

   FORT_TEMPERATURE_ENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
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

   boost::shared_ptr<pdat::CellData<double> > temp(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));

   boost::shared_ptr<pdat::CellData<double> > fa(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fa_id)));

   double factor = -0.5 * d_L;

   FORT_TEMPERATURE_ENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
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
