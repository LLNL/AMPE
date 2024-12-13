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
#include "FreeEnergyStrategyThreePhase.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

#include <cassert>

//=======================================================================

FreeEnergyStrategyThreePhase::FreeEnergyStrategyThreePhase(
    std::shared_ptr<tbox::Database> input_db, const double vml,
    const double vma, const double vmb, const int conc_l_id,
    const int conc_a_id, const int conc_b_id)
    : d_conc_l_id(conc_l_id), d_conc_a_id(conc_a_id), d_conc_b_id(conc_b_id)
{
   tbox::plog << "FreeEnergyStrategyThreePhase..." << std::endl;

   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
   assert(d_conc_b_id >= 0);
   assert(vml > 0.);
   assert(vma > 0.);
   assert(vmb > 0.);

   d_energy_conv_factor_L = 1.e-6 / vml;
   d_energy_conv_factor_A = 1.e-6 / vma;
   d_energy_conv_factor_B = 1.e-6 / vmb;

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1
   // tbox::plog << "FreeEnergyStrategyThreePhase:" <<
   // std::endl; tbox::plog << "Molar volume L =" << vml << std::endl;
   // tbox::plog << "Molar volume A =" << vma << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;
}

//=======================================================================

void FreeEnergyStrategyThreePhase ::computeFreeEnergyLiquid(
    hier::Patch& patch, const int temperature_id, const int fl_id,
    const bool gp)
{
   assert(fl_id >= 0);
   assert(temperature_id >= 0.);
   assert(d_conc_l_id >= 0);

   computeFreeEnergy(patch, temperature_id, fl_id, d_conc_l_id,
                     Thermo4PFM::PhaseIndex::phaseL, d_energy_conv_factor_L);
}

//=======================================================================

void FreeEnergyStrategyThreePhase ::computeFreeEnergySolidA(
    hier::Patch& patch, const int temperature_id, const int fa_id,
    const bool gp)
{
   assert(fa_id >= 0);
   assert(temperature_id >= 0.);
   assert(d_conc_a_id >= 0);

   computeFreeEnergy(patch, temperature_id, fa_id, d_conc_a_id,
                     Thermo4PFM::PhaseIndex::phaseA, d_energy_conv_factor_A);
}

//=======================================================================

void FreeEnergyStrategyThreePhase ::computeFreeEnergySolidB(
    hier::Patch& patch, const int temperature_id, const int fb_id,
    const bool gp)
{
   assert(fb_id >= 0);
   assert(temperature_id >= 0.);
   assert(d_conc_b_id >= 0);

   computeFreeEnergy(patch, temperature_id, fb_id, d_conc_b_id,
                     Thermo4PFM::PhaseIndex::phaseB, d_energy_conv_factor_B);
}

//=======================================================================

void FreeEnergyStrategyThreePhase ::computeFreeEnergy(
    hier::Patch& patch, const int temperature_id, const int f_id,
    const int conc_i_id, Thermo4PFM::PhaseIndex pi, const double energy_factor)
{
   assert(temperature_id >= 0);
   assert(f_id >= 0);
   assert(conc_i_id >= 0);

   const hier::Box& pbox = patch.getBox();

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));

   std::shared_ptr<pdat::CellData<double> > f(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_id)));

   std::shared_ptr<pdat::CellData<double> > c_i(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(conc_i_id)));

   computeFreeEnergy(pbox, temperature, f, c_i, pi, energy_factor);
}

//=======================================================================

void FreeEnergyStrategyThreePhase::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int conc_id, const int f_l_id, const int f_a_id,
    const int f_b_id, const int rhs_id)
{
   (void)time;

   assert(conc_id >= 0);
   assert(phase_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(f_b_id >= 0);
   assert(rhs_id >= 0);
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
   assert(d_conc_b_id >= 0);
   assert(temperature_id >= 0);

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);
   assert(phase->getDepth() > 1);

   std::shared_ptr<pdat::CellData<double> > t(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(t);

   std::shared_ptr<pdat::CellData<double> > fl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_l_id)));
   assert(fl);

   std::shared_ptr<pdat::CellData<double> > fa(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_a_id)));
   assert(fa);

   std::shared_ptr<pdat::CellData<double> > fb(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_b_id)));
   assert(fb);

   std::shared_ptr<pdat::CellData<double> > c_l(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_l_id)));
   assert(c_l);

   std::shared_ptr<pdat::CellData<double> > c_a(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_a_id)));
   assert(c_a);

   std::shared_ptr<pdat::CellData<double> > c_b(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_b_id)));
   assert(c_b);

   std::shared_ptr<pdat::CellData<double> > rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(rhs_id)));

   assert(rhs);
   assert(rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   const hier::Box& pbox(patch.getBox());

   addDrivingForceOnPatch(rhs, t, phase, fl, fa, fb, c_l, c_a, c_b, pbox);
}
