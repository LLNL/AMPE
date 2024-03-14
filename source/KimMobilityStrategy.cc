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
#include "KimMobilityStrategy.h"

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "QuatModel.h"
#include "Database2JSON.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"

#include "SAMRAI/tbox/InputManager.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

template <class FreeEnergyType>
KimMobilityStrategy<FreeEnergyType>::KimMobilityStrategy(
    QuatModel* quat_model, const int conc_l_id, const int conc_a_id,
    const int conc_b_id, const int temp_id,
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
    std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions)
    : SimpleQuatMobilityStrategy(quat_model),
      d_conc_l_id(conc_l_id),
      d_conc_a_id(conc_a_id),
      d_conc_b_id(conc_b_id),
      d_temp_id(temp_id),
      d_ncompositions(ncompositions)
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
   assert(d_temp_id >= 0);
   assert(d_ncompositions > 0);

   t_compute = tbox::TimerManager::getManager()->getTimer(
       "AMPE::KimMobilityStrategy::compute");

   std::shared_ptr<tbox::Database> conc_calphad_db =
       conc_db->getDatabase("Calphad");
   std::string calphad_filename = conc_calphad_db->getString("filename");

   std::shared_ptr<tbox::MemoryDatabase> calphad_db;
   boost::property_tree::ptree calphad_pt;

   if (calphad_filename.compare(calphad_filename.size() - 4, 4, "json") == 0) {
      boost::property_tree::read_json(calphad_filename, calphad_pt);
   } else {
      calphad_db.reset(new tbox::MemoryDatabase("calphad_db"));
      tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                       calphad_db);
      copyDatabase(calphad_db, calphad_pt);
   }

   std::shared_ptr<tbox::Database> newton_db;
   if (conc_db->isDatabase("NewtonSolver"))
      newton_db = conc_db->getDatabase("NewtonSolver");
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);

   d_fenergy =
       new FreeEnergyType(calphad_pt, newton_pt, energy_interp_func_type,
                          conc_interp_func_type);
}

template <>
KimMobilityStrategy<Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>::
    KimMobilityStrategy(
        QuatModel* quat_model, const int conc_l_id, const int conc_a_id,
        const int conc_b_id, const int temp_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions)
    : SimpleQuatMobilityStrategy(quat_model),
      d_conc_l_id(conc_l_id),
      d_conc_a_id(conc_a_id),
      d_conc_b_id(conc_b_id),
      d_temp_id(temp_id),
      d_ncompositions(ncompositions)
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
   assert(d_temp_id >= 0);
   assert(d_ncompositions > 0);

   t_compute = tbox::TimerManager::getManager()->getTimer(
       "AMPE::KimMobilityStrategy::compute");

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
   if (newton_db) copyDatabase(newton_db, newton_pt);

   d_fenergy = new Thermo4PFM::CALPHADFreeEnergyFunctionsTernary(
       calphad_pt, newton_pt, energy_interp_func_type, conc_interp_func_type);
}

template <>
KimMobilityStrategy<Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary>::
    KimMobilityStrategy(
        QuatModel* quat_model, const int conc_l_id, const int conc_a_id,
        const int conc_b_id, const int temp_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions)
    : SimpleQuatMobilityStrategy(quat_model),
      d_conc_l_id(conc_l_id),
      d_conc_a_id(conc_a_id),
      d_conc_b_id(conc_b_id),
      d_temp_id(temp_id),
      d_ncompositions(ncompositions)
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
   assert(d_temp_id >= 0);
   assert(d_ncompositions > 0);

   t_compute = tbox::TimerManager::getManager()->getTimer(
       "AMPE::KimMobilityStrategy::compute");

   pt::ptree conc_pt;
   copyDatabase(conc_db, conc_pt);
   d_fenergy = new Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary(
       conc_pt, energy_interp_func_type, conc_interp_func_type);
}

template <class FreeEnergyType>
void KimMobilityStrategy<FreeEnergyType>::computePhaseMobility(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& mobility_id, const double time, const CACHE_TYPE cache)
{
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

         std::shared_ptr<pdat::CellData<double> > phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));

         std::shared_ptr<pdat::CellData<double> > concl(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_l_id)));
         assert(concl);

         std::shared_ptr<pdat::CellData<double> > conca(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_a_id)));
         assert(conca);

         std::shared_ptr<pdat::CellData<double> > concb;
         if (d_conc_b_id > -1) {
            concb =
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_conc_b_id));
            assert(concb);
         }

         std::shared_ptr<pdat::CellData<double> > mobility(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(mobility_id)));
         assert(mobility);

         update(temperature, phi, concl, conca, concb, mobility, patch);
      }
   }

   t_compute->stop();
}

template <class FreeEnergyType>
void KimMobilityStrategy<FreeEnergyType>::update(
    std::shared_ptr<pdat::CellData<double> > cd_te,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_cl,
    std::shared_ptr<pdat::CellData<double> > cd_ca,
    std::shared_ptr<pdat::CellData<double> > cd_cb,
    std::shared_ptr<pdat::CellData<double> > cd_mob,
    std::shared_ptr<hier::Patch> patch)
{
   assert(cd_te);
   assert(cd_cl);
   assert(cd_ca);
   assert(cd_mob);
   assert(cd_mob->getGhostCellWidth()[0] <= cd_te->getGhostCellWidth()[0]);
   assert(cd_mob->getGhostCellWidth()[0] <= cd_cl->getGhostCellWidth()[0]);
   if (d_conc_b_id > -1)
      assert(cd_mob->getGhostCellWidth()[0] <= cd_phi->getGhostCellWidth()[0]);

   const hier::Box& temp_gbox = cd_te->getGhostBox();
   int min_te[3] = {temp_gbox.lower(0), temp_gbox.lower(1), 0};
   int inc_j_te = temp_gbox.numberCells(0);
#if (NDIM == 3)
   min_te[2] = temp_gbox.lower(2);
   int inc_k_te = inc_j_te * temp_gbox.numberCells(1);
#else
   int inc_k_te = 0;
#endif

   // assume ca and cl have the same number of ghost values
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

   std::vector<double> phaseconc(3 *
                                 d_ncompositions);  // 3 for up to three phases
   double* cl = &phaseconc[0];
   double* ca = &phaseconc[d_ncompositions];
   double* cb = &phaseconc[2 * d_ncompositions];
   std::vector<double> phi(3);

   for (int kk = min_mo[2]; kk <= max_mo[2]; kk++) {
      for (int jj = min_mo[1]; jj <= max_mo[1]; jj++) {
         for (int ii = min_mo[0]; ii <= max_mo[0]; ii++) {

            const double temp = cd_te->getPointer()[idx_te];
            for (unsigned ic = 0; ic < d_ncompositions; ic++)
               cl[ic] = cd_cl->getPointer(ic)[idx_ci];
            for (unsigned ic = 0; ic < d_ncompositions; ic++)
               ca[ic] = cd_ca->getPointer(ic)[idx_ci];
            if (d_conc_b_id > -1) {
               for (unsigned ic = 0; ic < d_ncompositions; ic++)
                  cb[ic] = cd_cb->getPointer(ic)[idx_ci];
               for (unsigned i = 0; i < 3; i++)
                  phi[i] = cd_phi->getPointer(i)[idx_ci];
            }
            cd_mob->getPointer()[idx_mo] =
                evaluateMobility(temp, phaseconc, phi);

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

template class KimMobilityStrategy<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>;
template class KimMobilityStrategy<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl>;
template class KimMobilityStrategy<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl>;
