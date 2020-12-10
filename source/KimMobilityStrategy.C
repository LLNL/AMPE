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
#include "KimMobilityStrategy.h"

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "QuatModel.h"
#ifdef HAVE_THERMO4PFM
#include "Database2JSON.h"
namespace pt = boost::property_tree;
#endif

#include "SAMRAI/tbox/InputManager.h"

KimMobilityStrategy::KimMobilityStrategy(
    QuatModel* quat_model, const int conc_l_id, const int conc_s_id,
    const int temp_id, const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type,
    std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions)
    : SimpleQuatMobilityStrategy(quat_model),
      d_conc_l_id(conc_l_id),
      d_conc_s_id(conc_s_id),
      d_temp_id(temp_id),
      d_ncompositions(ncompositions)
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_s_id >= 0);
   assert(d_temp_id >= 0);
   assert(d_ncompositions > 0);

   t_compute = tbox::TimerManager::getManager()->getTimer(
       "AMPE::KimMobilityStrategy::compute");

   std::string conc_model = conc_db->getStringWithDefault("model", "undefined");

   if (conc_model[0] == 'c') {

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

      if (ncompositions == 1) {
         d_fenergy = new CALPHADFreeEnergyFunctionsBinary(
#ifdef HAVE_THERMO4PFM
             calphad_pt, newton_pt,
#else
             calphad_db, newton_db,
#endif
             energy_interp_func_type, conc_interp_func_type,
             false);  // no 3rd phase
      } else {
         d_fenergy = new CALPHADFreeEnergyFunctionsTernary(
#ifdef HAVE_THERMO4PFM
             calphad_pt, newton_pt,
#else
             calphad_db, newton_db,
#endif
             energy_interp_func_type, conc_interp_func_type);
      }
   } else if (conc_model[0] == 'd') {
#ifdef HAVE_THERMO4PFM
      pt::ptree conc_pt;
      copyDatabase(conc_db, conc_pt);
      d_fenergy = new KKSFreeEnergyFunctionDiluteBinary(conc_pt,
                                                        energy_interp_func_type,
                                                        conc_interp_func_type);
#else
      d_fenergy = new KKSFreeEnergyFunctionDiluteBinary(conc_db,
                                                        energy_interp_func_type,
                                                        conc_interp_func_type);
#endif
   } else {
      TBOX_ERROR("Error: unknown concentration model in KimMobilityStrategy");
   }
}

void KimMobilityStrategy::computePhaseMobility(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& mobility_id, const double time, const CACHE_TYPE cache)
{
   (void)phase_id;
   (void)time;
   (void)cache;

   t_compute->start();

   const int maxl = hierarchy->getNumberOfLevels();

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {
         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<pdat::CellData<double> > temperature(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temp_id)));
         assert(temperature);

         std::shared_ptr<pdat::CellData<double> > concl(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_l_id)));
         assert(concl);

         std::shared_ptr<pdat::CellData<double> > concs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_s_id)));
         assert(concs);

         std::shared_ptr<pdat::CellData<double> > mobility(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(mobility_id)));
         assert(mobility);

         update(temperature, concl, concs, mobility, patch);
      }
   }

   t_compute->stop();
}

void KimMobilityStrategy::update(
    std::shared_ptr<pdat::CellData<double> > cd_te,
    std::shared_ptr<pdat::CellData<double> > cd_cl,
    std::shared_ptr<pdat::CellData<double> > cd_cs,
    std::shared_ptr<pdat::CellData<double> > cd_mob,
    std::shared_ptr<hier::Patch> patch)
{
   assert(cd_te);
   assert(cd_cl);
   assert(cd_cs);
   assert(cd_mob);
   assert(cd_mob->getGhostCellWidth()[0] <= cd_te->getGhostCellWidth()[0]);
   assert(cd_mob->getGhostCellWidth()[0] <= cd_cl->getGhostCellWidth()[0]);

   const hier::Box& temp_gbox = cd_te->getGhostBox();
   int min_te[3] = {temp_gbox.lower(0), temp_gbox.lower(1), 0};
   int inc_j_te = temp_gbox.numberCells(0);
#if (NDIM == 3)
   min_te[2] = temp_gbox.lower(2);
   int inc_k_te = inc_j_te * temp_gbox.numberCells(1);
#else
   int inc_k_te = 0;
#endif

   // assume cs and cl have the same number of ghost values
   const hier::Box& ci_gbox = cd_cl->getGhostBox();
   int min_ci[3] = {ci_gbox.lower(0), ci_gbox.lower(1), 0};
   int inc_j_ci = ci_gbox.numberCells(0);
#if (NDIM == 3)
   min_ci[2] = ci_gbox.lower(2);
   int inc_k_ci = inc_j_ci * ci_gbox.numberCells(1);
#else
   int inc_k_ci = 0;
#endif

   const hier::Box& mo_gbox = cd_mob->getGhostBox();
   int min_mo[3] = {mo_gbox.lower(0), mo_gbox.lower(1), 0};
#if (NDIM == 3)
   min_mo[2] = mo_gbox.lower(2);
#endif

   int max_mo[3] = {mo_gbox.upper(0), mo_gbox.upper(1), 0};
#if (NDIM == 3)
   max_mo[2] = mo_gbox.upper(2);
#endif

   int idx_te = (min_mo[0] - min_te[0]) + (min_mo[1] - min_te[1]) * inc_j_te +
                (min_mo[2] - min_te[2]) * inc_k_te;
   int idx_ci = (min_mo[0] - min_ci[0]) + (min_mo[1] - min_ci[1]) * inc_j_ci +
                (min_mo[2] - min_ci[2]) * inc_k_ci;
   int idx_mo = 0;

   std::vector<double> phaseconc(2 * d_ncompositions);  // 2 for two phases
   double* cl = &phaseconc[0];
   double* cs = &phaseconc[d_ncompositions];

   for (int kk = min_mo[2]; kk <= max_mo[2]; kk++) {
      for (int jj = min_mo[1]; jj <= max_mo[1]; jj++) {
         for (int ii = min_mo[0]; ii <= max_mo[0]; ii++) {

            const double temp = cd_te->getPointer()[idx_te];
            for (unsigned ic = 0; ic < d_ncompositions; ic++)
               cl[ic] = cd_cl->getPointer(ic)[idx_ci];
            for (unsigned ic = 0; ic < d_ncompositions; ic++)
               cs[ic] = cd_cs->getPointer(ic)[idx_ci];

            cd_mob->getPointer()[idx_mo] = evaluateMobility(temp, phaseconc);

            idx_te++;
            idx_ci++;
            idx_mo++;
         }
         idx_te += 2 * (cd_te->getGhostCellWidth()[0] -
                        cd_mob->getGhostCellWidth()[0]);
         idx_ci += 2 * (cd_cl->getGhostCellWidth()[0] -
                        cd_mob->getGhostCellWidth()[0]);
      }
      idx_te +=
          2 * inc_j_te *
          (cd_te->getGhostCellWidth()[1] - cd_mob->getGhostCellWidth()[1]);
      idx_ci +=
          2 * inc_j_ci *
          (cd_cl->getGhostCellWidth()[1] - cd_mob->getGhostCellWidth()[1]);
   }
}
