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
#include "PhaseFreeEnergyStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "ConcFort.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include <cassert>

using namespace SAMRAI;


PhaseFreeEnergyStrategy::PhaseFreeEnergyStrategy(
    const EnergyInterpolationType phase_interp_func_type,
    const EnergyInterpolationType eta_interp_func_type, const double fl,
    const double fa, const double fb, const double vml, const double vma,
    const double vmb, const bool with_third_phase)
{
   assert(vml == vml);
   assert(vma == vma);
   assert(vml > 0.);
   assert(vma > 0.);

   d_phase_interp_func_type = phase_interp_func_type;
   d_eta_interp_func_type = eta_interp_func_type;
   d_f_l = fl;
   d_f_a = fa;
   d_f_b = fb;

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;
   d_f_l *= 1.e-6 / vml;
   d_f_a *= 1.e-6 / vma;
   d_f_b *= 1.e-6 / vmb;

   d_with_third_phase = with_third_phase;

   tbox::plog << "FreeEnergyStrategy:" << std::endl;
   tbox::plog << "Molar volume L =" << vml << std::endl;
   tbox::plog << "Molar volume A =" << vma << std::endl;
   tbox::plog << "f liquid  =" << d_f_l << std::endl;
   tbox::plog << "f solid A =" << d_f_a << std::endl;
}

//=======================================================================

void PhaseFreeEnergyStrategy::computeFreeEnergyLiquid(hier::Patch& patch,
                                                      const int temperature_id,
                                                      const int fl_id,
                                                      const bool gp)
{
   assert(fl_id >= 0);
   (void)temperature_id;  // unused

   std::shared_ptr<pdat::CellData<double> > fl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fl_id)));
   assert(fl);

   fl->fillAll(d_f_l);
}

//=======================================================================

void PhaseFreeEnergyStrategy::computeFreeEnergySolidA(hier::Patch& patch,
                                                      const int temperature_id,
                                                      const int fs_id,
                                                      const bool gp)
{
   assert(fs_id >= 0);
   (void)temperature_id;  // unused

   std::shared_ptr<pdat::CellData<double> > fs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fs_id)));
   assert(fs);

   fs->fillAll(d_f_a);
}

//=======================================================================

void PhaseFreeEnergyStrategy::computeFreeEnergySolidB(hier::Patch& patch,
                                                      const int temperature_id,
                                                      const int fs_id,
                                                      const bool gp)
{
   assert(fs_id >= 0);
   (void)temperature_id;  // unused

   std::shared_ptr<pdat::CellData<double> > fs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(fs_id)));
   assert(fs);

   fs->fillAll(d_f_b);
}

//=======================================================================

void PhaseFreeEnergyStrategy::addDrivingForce(
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
   (void)conc_id;         // unused
   (void)temperature_id;  // unused

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

   int three_phase = 0;
   double* ptr_fb = nullptr;
   double* ptr_eta = nullptr;
   if (d_with_third_phase) {
      three_phase = 1;
      std::shared_ptr<pdat::CellData<double> > eta(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(eta_id)));
      ptr_eta = eta->getPointer();

      std::shared_ptr<pdat::CellData<double> > fb(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(f_b_id)));
      ptr_fb = fb->getPointer();
   }

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   const char interpf = energyInterpChar(d_phase_interp_func_type);
   const char interpe = energyInterpChar(d_eta_interp_func_type);

   PHASERHS_FENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                    ifirst(2), ilast(2),
#endif
                    fl->getPointer(), fa->getPointer(), ptr_fb,
                    phase->getPointer(), phase->getGhostCellWidth()[0], ptr_eta,
                    phase->getGhostCellWidth()[0], rhs->getPointer(), 0,
                    &interpf, &interpe, three_phase);
}

//=======================================================================

void PhaseFreeEnergyStrategy::addDrivingForceEta(
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
   (void)conc_id;         // unused
   (void)temperature_id;  // unused

   std::shared_ptr<pdat::CellData<double> > eta(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(eta_id)));
   assert(eta);

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

   std::shared_ptr<pdat::CellData<double> > fb(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_b_id)));
   assert(fb);

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
   const char interpe = energyInterpChar(d_eta_interp_func_type);

   ETARHS_FENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                  ifirst(2), ilast(2),
#endif
                  fl->getPointer(), fa->getPointer(), fb->getPointer(),
                  phase->getPointer(), phase->getGhostCellWidth()[0],
                  eta->getPointer(), eta->getGhostCellWidth()[0],
                  rhs->getPointer(), 0, &interpf, &interpe);
}

//=======================================================================

void PhaseFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseL(
    const double temp, const std::vector<double>& c_l,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   d2fdc2.assign(d2fdc2.size(), 0.);
}

//=======================================================================

void PhaseFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseA(
    const double temp, const std::vector<double>& c_a,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   d2fdc2.assign(d2fdc2.size(), 0.);
}

//=======================================================================

void PhaseFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseB(
    const double temp, const std::vector<double>& c_b,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   d2fdc2.assign(d2fdc2.size(), 0.);
}
