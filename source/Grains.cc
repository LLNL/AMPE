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
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include "Grains.h"
#include "MinIntCoarsen.h"
#include "tools.h"

#include <cassert>
#include <set>


Grains::Grains(const int qlen, const bool visit_output,
               std::shared_ptr<tbox::Database> input_db)
    : d_qlen(qlen), d_visit_output(visit_output)
{
   d_grain_number_id = -1;
   d_grain_number_scr_id = -1;
   d_number_of_grains = -1;
   d_grain_diag_isActive = false;
   d_grain_extend_isActive = false;
   d_grain_phase_threshold = 0.85;
   d_grain_size_minimum = 0.0;

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_findAndNumberGrains_timer =
       tman->getTimer("AMPE::Grains::findAndNumberGrains()");
   t_extendGrainOrientation_timer =
       tman->getTimer("AMPE::Grains::extendGrainOrientation()");

   if (input_db->isDatabase("GrainDiagnostics")) {
      std::shared_ptr<tbox::Database> g_diag_db =
          input_db->getDatabase("GrainDiagnostics");
      d_grain_diag_isActive = true;

      if (g_diag_db->keyExists("phase_threshold")) {
         d_grain_phase_threshold = g_diag_db->getDouble("phase_threshold");
      }

      if (g_diag_db->keyExists("minimum_size")) {
         d_grain_size_minimum = g_diag_db->getDouble("minimum_size");
      }
   }
   if (input_db->isDatabase("GrainExtension")) {
      d_grain_extend_isActive = true;
   }
}

//=======================================================================

void Grains::initialize(std::shared_ptr<tbox::Database> input_db,
                        const bool all_periodic)
{
   (void)input_db;
   (void)all_periodic;
}

//=======================================================================

void Grains::registerVariables(void)
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::VariableContext> current =
       variable_db->getContext("CURRENT");
   std::shared_ptr<hier::VariableContext> scratch =
       variable_db->getContext("SCRATCH");

   if (d_grain_diag_isActive) {
      d_grain_number_var.reset(
          new pdat::CellVariable<int>(tbox::Dimension(NDIM), "grain_number"));
      assert(d_grain_number_var);
      d_grain_number_id = variable_db->registerVariableAndContext(
          d_grain_number_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 1));
      d_grain_number_scr_id = variable_db->registerVariableAndContext(
          d_grain_number_var, scratch,
          hier::IntVector(tbox::Dimension(NDIM), 1));
      d_local_data.setFlag(d_grain_number_id);
      d_local_data.setFlag(d_grain_number_scr_id);

      if (d_visit_output) {
         d_grain_volume_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            "grain_volume"));
         assert(d_grain_volume_var);
         d_grain_volume_id = variable_db->registerVariableAndContext(
             d_grain_volume_var, current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         d_local_data.setFlag(d_grain_volume_id);
      }
   }

   if (d_grain_extend_isActive) {
      d_grain_extend_var.reset(
          new pdat::CellVariable<int>(tbox::Dimension(NDIM), "grain_extend"));
      assert(d_grain_extend_var);
      d_grain_extend_id = variable_db->registerVariableAndContext(
          d_grain_extend_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 1));
      d_grain_extend_scr_id = variable_db->registerVariableAndContext(
          d_grain_extend_var, scratch,
          hier::IntVector(tbox::Dimension(NDIM), 1));

      d_grain_quat_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM), "grain_quat",
                                         d_qlen));
      assert(d_grain_quat_var);
      d_grain_quat_id = variable_db->registerVariableAndContext(
          d_grain_quat_var, current, hier::IntVector(tbox::Dimension(NDIM), 1));
      d_grain_quat_scr_id = variable_db->registerVariableAndContext(
          d_grain_quat_var, scratch, hier::IntVector(tbox::Dimension(NDIM), 1));
   }
}

//=======================================================================

void Grains::initializeRefineCoarsenAlgorithms(
    std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
    std::shared_ptr<hier::CoarsenOperator> quat_coarsen_op)
{
   tbox::plog << "Grains::initializeRefineCoarsenAlgorithms()" << std::endl;

   if (d_grain_diag_isActive) {
      assert(d_grain_number_id >= 0);
      assert(grid_geom);

      d_grain_number_refine_op =
          grid_geom->lookupRefineOperator(d_grain_number_var,
                                          "CONSTANT_"
                                          "REFINE");

      d_grain_number_refine_alg.reset(new xfer::RefineAlgorithm());

      d_grain_number_refine_alg->registerRefine(
          d_grain_number_id,      // destination
          d_grain_number_id,      // source
          d_grain_number_scr_id,  // temporary
          d_grain_number_refine_op);

      d_grain_number_coarsen_alg.reset(
          new xfer::CoarsenAlgorithm(tbox::Dimension(NDIM)));

      d_grain_number_coarsen_op.reset(new MinIntCoarsen());

      d_grain_number_coarsen_alg->registerCoarsen(d_grain_number_id,
                                                  d_grain_number_id,
                                                  d_grain_number_coarsen_op);
   }

   if (d_grain_extend_isActive) {
      assert(d_grain_quat_scr_id >= 0);
      assert(d_grain_extend_scr_id >= 0);
      assert(d_grain_extend_id >= 0);

      d_grain_quat_refine_op =
          grid_geom->lookupRefineOperator(d_grain_quat_var, "CONSTANT_REFINE");

      d_grain_extend_refine_op =
          grid_geom->lookupRefineOperator(d_grain_extend_var,
                                          "CONSTANT_"
                                          "REFINE");

      d_grain_quat_refine_alg.reset(new xfer::RefineAlgorithm());

      d_grain_quat_refine_alg->registerRefine(
          d_grain_quat_scr_id,  // destination
          d_grain_quat_id,      // source
          d_grain_quat_scr_id,  // temporary
          d_grain_quat_refine_op);

      d_grain_extend_refine_alg.reset(new xfer::RefineAlgorithm());

      d_grain_extend_refine_alg->registerRefine(
          d_grain_extend_scr_id,  // destination
          d_grain_extend_id,      // source
          d_grain_extend_scr_id,  // temporary
          d_grain_extend_refine_op);

      d_grain_quat_coarsen_op = quat_coarsen_op;
      d_grain_quat_coarsen_alg.reset(
          new xfer::CoarsenAlgorithm(tbox::Dimension(NDIM)));
      d_grain_quat_coarsen_alg->registerCoarsen(d_grain_quat_id,
                                                d_grain_quat_id,
                                                d_grain_quat_coarsen_op);

      d_grain_extend_coarsen_op.reset(new MinIntCoarsen());
      d_grain_extend_coarsen_alg.reset(
          new xfer::CoarsenAlgorithm(tbox::Dimension(NDIM)));
      d_grain_extend_coarsen_alg->registerCoarsen(d_grain_extend_id,
                                                  d_grain_extend_id,
                                                  d_grain_extend_coarsen_op);
   }
}

//=======================================================================

void Grains::initializeLevelData(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number, const double time, const bool can_be_refined,
    const bool initial_time, const std::shared_ptr<hier::PatchLevel>& old_level,
    const bool allocate_data)
{
   (void)time;
   (void)can_be_refined;
   (void)initial_time;

   tbox::pout << "Grains::initializeLevelData()" << std::endl;
   std::shared_ptr<hier::PatchLevel> level =
       hierarchy->getPatchLevel(level_number);

   if (level_number > 0) {
      if (old_level) old_level->deallocatePatchData(d_local_data);
   }

   level->allocatePatchData(d_local_data);

   // initialize scratch data to -1 so that values at physical boundaries
   // are at -1 after call to fillData
   if (d_grain_diag_isActive)
      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {

         std::shared_ptr<hier::Patch> patch = *p;

         assert(d_grain_number_scr_id >= 0);
         std::shared_ptr<pdat::CellData<int> > gs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                 patch->getPatchData(d_grain_number_scr_id)));
         assert(gs);

         // initialize data with -1, including ghost values
         gs->fillAll(-1, gs->getGhostBox());
      }
}

//=======================================================================

void Grains::registerWithVisit(
    std::shared_ptr<appu::VisItDataWriter> visit_data_writer)
{
   visit_data_writer->registerPlotQuantity("grain_number", "SCALAR",
                                           d_grain_number_id, 0);

   visit_data_writer->registerPlotQuantity("grain_volume", "SCALAR",
                                           d_grain_volume_id, 0);
}

//-----------------------------------------------------------------------

// This method relies on the weight_id to be correctly set with the
// cell volumes.

void Grains::findAndNumberGrains(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int phase_id,
    const int weight_id, const double time)
{
   assert(d_grain_diag_isActive);
   assert(d_grain_phase_threshold > 0.);

   tbox::plog << "Grains::findAndNumberGrains() with threshold "
              << d_grain_phase_threshold << std::endl;

   t_findAndNumberGrains_timer->start();

   assert(phase_id >= 0);
   assert(weight_id >= 0);
   assert(d_grain_number_id >= 0);

   const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

   int maxln = hierarchy->getFinestLevelNumber();

   int max_iteration_count = 0;

   //------------------------------------------------------------
   // Give each cell a unique number across entire hierarchy.

   int nlevels = maxln + 1;
   int* level_stride = new int[nlevels * NDIM];
   int* level_offset = new int[nlevels];
   level_offset[0] = 0;

   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      const hier::BoxContainer& bl(
          level->getPhysicalDomain(hier::BlockId::zero()));
      bl.print();

      hier::Box bbox = bl.getBoundingBox();
      int bbox_size = bl.size();

      if (ln < maxln) {
         level_offset[ln + 1] = level_offset[ln] + bbox_size;
      }

      int nn = NDIM * ln;
      level_stride[nn] = 1;
      for (tbox::Dimension::dir_t dd = 0; dd < NDIM; dd++) {
         level_stride[nn + dd + 1] =
             level_stride[nn + dd] * bbox.numberCells(dd);
      }

      for (tbox::Dimension::dir_t dd = 0; dd < NDIM; dd++) {
         int width = bbox.numberCells(dd);
         if (width > max_iteration_count) {
            max_iteration_count =
                2 * width;  // need factor 2 for nonperiodic BC
         }
      }
   }

   // safety factor to take into account non-convex grains
   max_iteration_count *= 2;

   int count_cells = 0;
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {

         std::shared_ptr<hier::Patch> patch = *p;
         const hier::Box& pbox = patch->getBox();

         std::shared_ptr<pdat::CellData<double> > phase(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));

         std::shared_ptr<pdat::CellData<int> > g(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                 patch->getPatchData(d_grain_number_id)));
         assert(g);

         // initialize data with -1, including ghost values
         g->fillAll(-1, g->getGhostBox());

         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend;
              ++i) {
            pdat::CellIndex cell(*i);

            if ((*phase)(cell) >= d_grain_phase_threshold) {
               int nn = level_offset[ln];
               for (tbox::Dimension::dir_t dd = 0; dd < NDIM; ++dd) {
                  nn += level_stride[NDIM * ln + dd] * cell(dd);
               }
               (*g)(cell) = nn;
               assert((*g)(cell) >= 0);
               count_cells++;
            }
         }
      }
   }

   delete[] level_stride;
   delete[] level_offset;

   int count_cells_sum;
   mpi.Allreduce(&count_cells, &count_cells_sum, 1, MPI_INT, MPI_SUM);
   if (count_cells_sum == 0)
      tbox::pout << "WARNING: not cell above threshold of "
                 << d_grain_phase_threshold << std::endl;

   for (int ln = maxln - 1; ln >= 0; ln--) {
      d_grain_number_coarsen_sched[ln]->coarsenData();
   }

   //------------------------------------------------------------
   // Bucket fill grains with lowest grain number found.
   assert(d_grain_number_refine_sched.size() > 0);
   int counter = 0;
   while (counter < max_iteration_count) {
      int changed = 0;

      for (int ln = 0; ln <= maxln; ln++) {
         std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

         // Fill patch boundaries
         assert(d_grain_number_refine_sched[ln]);

         bool do_physical_boundary_fill = false;
         d_grain_number_refine_sched[ln]->fillData(time,
                                                   do_physical_boundary_fill);

         for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
              p++) {

            std::shared_ptr<hier::Patch> patch(*p);
            const hier::Box& pbox(patch->getBox());

            std::shared_ptr<pdat::CellData<int> > g(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                    patch->getPatchData(d_grain_number_id)));

            std::shared_ptr<pdat::CellData<double> > w(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(weight_id)));

#if 0  // optimized loop...       
            int imin = pbox.lower(0);
            int imax = pbox.upper(0);
            int jmin = pbox.lower(1);
            int jmax = pbox.upper(1);
            int kmin = 0;
            int kmax = 0;
#if (NDIM == 3)
            kmin = pbox.lower(2);
            kmax = pbox.upper(2);
#endif

            const hier::Box& gbox(g->getGhostBox());
            int imin_g = gbox.lower(0);
            int jmin_g = gbox.lower(1);
            int kmin_g = 0;
            int inc[3]={1,gbox.numberCells(0),0};
#if (NDIM == 3)
            kmin_g = gbox.lower(2);
            inc[2] = inc[1] * gbox.numberCells(1);
#endif

            const hier::Box& w_gbox(w->getGhostBox());
            int imin_w = w_gbox.lower(0);
            int jmin_w = w_gbox.lower(1);
            int kmin_w = 0;
            int incw[3]={1,w_gbox.numberCells(0),0};
#if (NDIM == 3)
            kmin_w = w_gbox.lower(2);
            incw[2] = incw[1] * w_gbox.numberCells(1);
#endif

            int* g_temp = g->getPointer();
            double* w_temp = w->getPointer();
            
            for(int kk=kmin;kk<=kmax;++kk)
            for(int jj=jmin;jj<=jmax;++jj)
            for(int ii=imin;ii<=imax;++ii)
            {
               const int idx = (ii - imin_g) 
                             + (jj - jmin_g) * inc[1] 
                             + (kk - kmin_g) * inc[2];
               const int idw = (ii - imin_w) 
                             + (jj - jmin_w) * incw[1] 
                             + (kk - kmin_w) * incw[2];
               int n = g_temp[idx];
               double v = w_temp[idw];
               if ( n >= 0 && v > 0. ) {
                  for ( int dd = 0; dd < NDIM; ++dd ) {
                     const int idxm = idx - inc[dd];
                     int nm = g_temp[idxm];
                     if ( nm >= 0 && nm < n ) {
                        g_temp[idx] = nm;
                        changed = 1;
                     }
                     
                     const int idxp = idx + inc[dd];
                     int np = g_temp[idxp];
                     if ( np >= 0 && np < n ) {
                        g_temp[idx] = np;
                        changed = 1;
                     }
                  }
               }
            }

#else
            pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
            for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
                 ic != iend; ++ic) {
               const pdat::CellIndex icell = *ic;
               const int n = (*g)(icell);
               const double v = (*w)(icell);
               if (n >= 0 && v > 0.) {
                  for (tbox::Dimension::dir_t dd = 0; dd < NDIM; dd++) {
                     pdat::CellIndex cm(icell);
                     cm[dd]--;
                     // std::cout<<"CellIndex:"<<*ic<<", val="<<(*g)(icell)<<",
                     // (*g)(cm)="<<(*g)(cm)<<endl;
                     int nm = (*g)(cm);
                     if (nm >= 0 && nm < n) {
                        (*g)(icell) = nm;
                        changed = 1;
                     }

                     pdat::CellIndex cp(icell);
                     cp[dd]++;

                     int np = (*g)(cp);
                     if (np >= 0 && np < n) {
                        (*g)(icell) = np;
                        changed = 1;
                     }
                  }
               }
            }
#endif
         }
         if (ln > 0) {
            d_grain_number_coarsen_sched[ln - 1]->coarsenData();
         }
      }

      counter++;

      int any_changed;
      mpi.Allreduce(&changed, &any_changed, 1, MPI_INT, MPI_MAX);

      if (any_changed == 0) break;

   }  // while ( counter < max_iteration_count )

   if (counter >= max_iteration_count) {
      tbox::pout << "WARNING: grain numbering did not converge" << std::endl;
      tbox::pout << "counter = " << counter
                 << ", max_iteration_count = " << max_iteration_count
                 << std::endl;
   } else {
      tbox::pout << "Grain numbering converged in " << counter << " iterations"
                 << std::endl;
   }
   // return;

   //------------------------------------------------------------
   // Make local sorted list of grain numbers

   std::set<int> local_grain_list;

   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {

         std::shared_ptr<hier::Patch> patch = *p;
         const hier::Box& pbox = patch->getBox();

         std::shared_ptr<pdat::CellData<int> > g(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                 patch->getPatchData(d_grain_number_id)));

         std::shared_ptr<pdat::CellData<double> > w(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(weight_id)));

         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend;
              ++i) {
            pdat::CellIndex cell = *i;
            int n = (*g)(cell);
            double v = (*w)(cell);

            if (n >= 0 && v > 0.) {
               local_grain_list.insert(n);
            }
         }
      }
   }

   //------------------------------------------------------------
   // Make global sorted list of grain numbers

   int size_in = (int)local_grain_list.size();
   int* local_grain_numbers = new int[size_in];

   int ii = 0;
   for (std::set<int>::const_iterator it = local_grain_list.begin();
        it != local_grain_list.end(); it++) {
      local_grain_numbers[ii] = *it;
      ii++;
   }

   int size_out = sumReduction(size_in);
   std::vector<int> grain_numbers(size_out, -1);

   allGatherv(local_grain_numbers, size_in, &grain_numbers[0], size_out);

   std::set<int> grain_list;
   for (ii = 0; ii < size_out; ii++) {
      grain_list.insert(grain_numbers[ii]);
   }

   delete[] local_grain_numbers;
   grain_numbers.clear();

   //------------------------------------------------------------
   // Renumber

   std::map<int, int> grain_number_map;

   ii = 0;
   for (std::set<int>::const_iterator it = grain_list.begin();
        it != grain_list.end(); it++) {
      grain_number_map[*it] = ii;
      ii++;
   }
   d_number_of_grains = ii;

   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {

         std::shared_ptr<hier::Patch> patch = *p;
         const hier::Box& pbox = patch->getBox();

         std::shared_ptr<pdat::CellData<int> > g(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                 patch->getPatchData(d_grain_number_id)));

         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend;
              ++i) {
            pdat::CellIndex cell = *i;
            int n = (*g)(cell);

            if (n >= 0) {
               int gn = grain_number_map[n];
               (*g)(cell) = gn;
            }
         }
      }
   }

   for (int ln = maxln - 1; ln >= 0; ln--) {
      d_grain_number_coarsen_sched[ln]->coarsenData();
   }

   t_findAndNumberGrains_timer->stop();
}

//-----------------------------------------------------------------------

// This method relies on the weight_id to be correctly set with the
// cell volumes.

void Grains::computeGrainVolumes(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int weight_id)
{
   assert(weight_id >= 0);
   assert(d_grain_number_id >= 0);

   int maxln = hierarchy->getFinestLevelNumber();

   // Compute volume of each grain

   std::map<int, double> map_lcl_grain_v;

   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {

         std::shared_ptr<hier::Patch> patch = *p;
         const hier::Box& pbox = patch->getBox();

         std::shared_ptr<pdat::CellData<int> > grain(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                 patch->getPatchData(d_grain_number_id)));

         std::shared_ptr<pdat::CellData<double> > wgt(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(weight_id)));

         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend;
              ++i) {
            pdat::CellIndex cell = *i;
            int n = (*grain)(cell);
            double v = (*wgt)(cell);

            if (n >= 0 && v > 0.) {
               map_lcl_grain_v[n] += v;
            }
         }
      }
   }

   std::map<int, double> map_gbl_grain_v;
   listLocalToGlobal(map_lcl_grain_v, map_gbl_grain_v);

   for (std::map<int, double>::const_iterator it = map_gbl_grain_v.begin();
        it != map_gbl_grain_v.end(); it++) {
      tbox::pout << "Volume of grain " << it->first << " = " << it->second
                 << std::endl;
   }

   // If a cell is part of a grain, assign it the grain volume for vis

   if (d_visit_output) {
      for (int ln = 0; ln <= maxln; ln++) {
         std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

         for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
              p++) {

            std::shared_ptr<hier::Patch> patch = *p;
            const hier::Box& pbox = patch->getBox();

            std::shared_ptr<pdat::CellData<int> > g(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                    patch->getPatchData(d_grain_number_id)));

            std::shared_ptr<pdat::CellData<double> > v(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_grain_volume_id)));

            pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
            for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox));
                 i != iend; ++i) {
               pdat::CellIndex cell = *i;
               int n = (*g)(cell);


               if (n >= 0) {
                  (*v)(cell) = map_gbl_grain_v[n];
               } else {
                  (*v)(cell) = 0.0;
               }
            }
         }
      }
   }
}

//-----------------------------------------------------------------------

// This method relies on the weight_id to be correctly set with the
// cell volumes.

void Grains::computeGrainConcentrations(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time,
    const int conc_id, const int weight_id)
{
   assert(conc_id >= 0);
   assert(weight_id >= 0);
   assert(d_grain_number_id >= 0);

   int maxln = hierarchy->getFinestLevelNumber();

   // Compute volume average concentration of each grain

   std::map<int, double> map_lcl_grain_c;
   std::map<int, double> map_lcl_grain_w;

   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {

         std::shared_ptr<hier::Patch> patch = *p;
         const hier::Box& pbox = patch->getBox();

         std::shared_ptr<pdat::CellData<int> > grain(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                 patch->getPatchData(d_grain_number_id)));

         std::shared_ptr<pdat::CellData<double> > wgt(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(weight_id)));

         std::shared_ptr<pdat::CellData<double> > conc(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(conc_id)));

         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend;
              ++i) {
            pdat::CellIndex cell = *i;
            int n = (*grain)(cell);
            double w = (*wgt)(cell);
            double c = (*conc)(cell);

            if (n >= 0 && w > 0.) {
               map_lcl_grain_c[n] += c * w;
               map_lcl_grain_w[n] += w;
            }
         }
      }
   }

   std::map<int, double> map_gbl_grain_c;
   std::map<int, double> map_gbl_grain_w;
   listLocalToGlobal(map_lcl_grain_c, map_gbl_grain_c);
   listLocalToGlobal(map_lcl_grain_w, map_gbl_grain_w);

   double c_all = 0.0;
   double w_all = 0.0;

   std::map<int, double> map_gbl_grain_cavg;
   for (std::map<int, double>::const_iterator it = map_gbl_grain_c.begin();
        it != map_gbl_grain_c.end(); it++) {
      int n = it->first;
      double cavg = map_gbl_grain_c[n] / map_gbl_grain_w[n];
      map_gbl_grain_cavg[n] = cavg;

      c_all += map_gbl_grain_c[n];
      w_all += map_gbl_grain_w[n];

      tbox::pout << "Concentration average of grain " << n << " = " << cavg
                 << std::endl;
   }

   if (map_gbl_grain_c.size() > 0)
      tbox::pout << "Concentration average for all grains at " << time << " = "
                 << c_all / w_all << std::endl;

   // Compute concentration deviation of each grain

   std::map<int, double> map_lcl_grain_cdev;

   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {

         std::shared_ptr<hier::Patch> patch = *p;
         const hier::Box& pbox = patch->getBox();

         std::shared_ptr<pdat::CellData<int> > grain(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                 patch->getPatchData(d_grain_number_id)));

         std::shared_ptr<pdat::CellData<double> > wgt(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(weight_id)));

         std::shared_ptr<pdat::CellData<double> > conc(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(conc_id)));

         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend;
              ++i) {
            pdat::CellIndex cell = *i;
            int n = (*grain)(cell);
            double w = (*wgt)(cell);
            double c = (*conc)(cell);

            if (n >= 0 && w > 0.) {
               double cdiff = c - map_gbl_grain_cavg[n];
               double cdev = cdiff * cdiff;

               map_lcl_grain_cdev[n] += w * cdev;
            }
         }
      }
   }

   std::map<int, double> map_gbl_grain_cdev;
   listLocalToGlobal(map_lcl_grain_cdev, map_gbl_grain_cdev);

   for (std::map<int, double>::const_iterator it = map_gbl_grain_cdev.begin();
        it != map_gbl_grain_cdev.end(); it++) {
      int n = it->first;
      double cdev = map_gbl_grain_cdev[n] / map_gbl_grain_w[n];
      //      map_gbl_grain_cdev[n] = cdev;

      tbox::pout << "Concentration deviation of grain " << n << " = " << cdev
                 << std::endl;
   }
}

//=======================================================================
//
// modifies quat_id
//
void Grains::extendGrainOrientation(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time,
    const int quat_scratch_id, const int phase_id, const int quat_id)
{
   tbox::pout << "Extending grain orientation" << std::endl;
   assert(quat_scratch_id >= 0);
   assert(phase_id >= 0);
   assert(quat_id >= 0);
   assert(d_grain_extend_isActive);

   assert(d_grain_extend_id >= 0);
   assert(d_grain_extend_scr_id >= 0);
   assert(d_grain_quat_id >= 0);
   assert(d_grain_quat_scr_id >= 0);

   t_extendGrainOrientation_timer->start();

   const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
   int maxln = hierarchy->getFinestLevelNumber();

   // allocate vars for this algorithm
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      level->allocatePatchData(d_grain_extend_id, time);
      level->allocatePatchData(d_grain_extend_scr_id, time);
      level->allocatePatchData(d_grain_quat_id, time);
      level->allocatePatchData(d_grain_quat_scr_id, time);
   }

   // Copy original quat data to local variables
   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
   cellops.copyData(d_grain_quat_id, quat_scratch_id);

   //------------------------------------------------------------

   // Find maximum iteration number

   int max_iteration_count = 0;

   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      hier::BoxContainer bl(level->getPhysicalDomain(hier::BlockId::zero()));
      hier::Box bbox = bl.getBoundingBox();

      for (tbox::Dimension::dir_t dd = 0; dd < NDIM; ++dd) {
         int width = bbox.numberCells(dd);
         if (width > max_iteration_count) {
            max_iteration_count = width;
         }
      }
   }

   // Tag grain cells with 1, non-grain with 0

   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {

         std::shared_ptr<hier::Patch> patch = *p;
         const hier::Box& pbox = patch->getBox();

         std::shared_ptr<pdat::CellData<double> > phase(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));

         std::shared_ptr<pdat::CellData<int> > flag(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                 patch->getPatchData(d_grain_extend_id)));

         flag->fillAll(0);

         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend;
              ++i) {
            pdat::CellIndex cell = *i;

            if ((*phase)(cell) >= d_grain_phase_threshold) {
               (*flag)(cell) = 1;
            }
         }
      }

      if (ln > 0) {
         d_grain_extend_coarsen_sched[ln - 1]->coarsenData();
      }
   }

   // Iterate over all cells, filling with a tagged neighbor if there
   // is one.
   int nn;
   for (nn = 0; nn < max_iteration_count; nn++) {
      int one_still_unset = 1;

      for (int ln = 0; ln <= maxln; ln++) {
         std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

         // Fill grain_extend_src with grain_extend
         d_grain_extend_refine_sched[ln]->fillData(time, false);

         // fill grain_quat_scr with grain_quat
         d_grain_quat_refine_sched[ln]->fillData(time, false);

         for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
              p++) {

            std::shared_ptr<hier::Patch> patch = *p;
            const hier::Box& pbox = patch->getBox();

            std::shared_ptr<pdat::CellData<int> > flag(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                    patch->getPatchData(d_grain_extend_id)));

            std::shared_ptr<pdat::CellData<int> > flag_scr(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                    patch->getPatchData(d_grain_extend_scr_id)));

            std::shared_ptr<pdat::CellData<double> > quat(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_grain_quat_id)));

            std::shared_ptr<pdat::CellData<double> > quat_scr(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_grain_quat_scr_id)));

            one_still_unset = 0;

            pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
            for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox));
                 i != iend; ++i) {
               pdat::CellIndex cell = *i;

               bool this_cell_still_unset = false;

               if ((*flag_scr)(cell) == 0) {
                  this_cell_still_unset = true;

                  for (int dd = 0; dd < NDIM; dd++) {
                     pdat::CellIndex cp = *i;
                     cp(dd)++;

                     pdat::CellIndex cm = *i;
                     cm(dd)--;

                     if ((*flag_scr)(cp) == 1) {
                        for (int ll = 0; ll < d_qlen; ll++) {
                           (*quat)(cell, ll) = (*quat_scr)(cp, ll);
                        }
                        (*flag)(cell) = 1;
                        this_cell_still_unset = false;
                     } else if ((*flag_scr)(cm) == 1) {
                        for (int ll = 0; ll < d_qlen; ll++) {
                           (*quat)(cell, ll) = (*quat_scr)(cm, ll);
                        }
                        (*flag)(cell) = 1;
                        this_cell_still_unset = false;
                     }

                     if (!this_cell_still_unset) break;
                  }
               }

               if (this_cell_still_unset) {
                  one_still_unset = 1;
               }
            }  // loop over cells
         }

         if (ln > 0) {
            d_grain_extend_coarsen_sched[ln - 1]->coarsenData();
            d_grain_quat_coarsen_sched[ln - 1]->coarsenData();
         }

      }  // for ( int ln )

      int any_still_unset = one_still_unset;
      mpi.AllReduce(&any_still_unset, 1, MPI_MAX);
      if (any_still_unset == 0) break;
   }

   if (nn >= max_iteration_count) {
      tbox::pout << "WARNING: grain extension did not converge" << std::endl;
   } else {
      tbox::pout << "Grain extension converged in " << nn << " iterations"
                 << std::endl;
   }

   for (int ln = maxln - 1; ln >= 0; ln--) {
      d_grain_quat_coarsen_sched[ln]->coarsenData();
   }

   cellops.copyData(quat_id, d_grain_quat_id);

   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      level->deallocatePatchData(d_grain_quat_id);
      level->deallocatePatchData(d_grain_quat_scr_id);
      level->deallocatePatchData(d_grain_extend_id);
      level->deallocatePatchData(d_grain_extend_scr_id);
   }

   t_extendGrainOrientation_timer->stop();
}

//=======================================================================

void Grains::resetHierarchyConfiguration(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarsest_level, const int finest_level)
{
   tbox::pout << "Grains::resetHierarchyConfiguration()" << std::endl;

   const int nlev = hierarchy->getNumberOfLevels();

   d_grain_number_refine_sched.resize(nlev);

   if (d_grain_diag_isActive)
      for (int ln = coarsest_level; ln <= finest_level; ln++) {
         assert(ln < nlev);

         std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

         // the last argument should be NULL
         // so that physical boundary values are -1, as initialized
         d_grain_number_refine_sched[ln] =
             d_grain_number_refine_alg->createSchedule(level, ln - 1, hierarchy,
                                                       NULL);
      }

   d_grain_number_coarsen_sched.resize(hierarchy->getNumberOfLevels());

   const int ln_beg = coarsest_level - (coarsest_level > 0);
   const int ln_end = finest_level;

   for (int ln = ln_beg; ln < ln_end; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      std::shared_ptr<hier::PatchLevel> finer_level =
          hierarchy->getPatchLevel(ln + 1);

      d_grain_number_coarsen_sched[ln] =
          d_grain_number_coarsen_alg->createSchedule(level, finer_level);
   }
   if (d_grain_extend_isActive) {
      d_grain_extend_refine_sched.resize(nlev);
      d_grain_quat_refine_sched.resize(nlev);

      for (int ln = coarsest_level; ln <= finest_level; ln++) {
         std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

         d_grain_extend_refine_sched[ln] =
             d_grain_extend_refine_alg->createSchedule(level, ln - 1, hierarchy,
                                                       NULL);
         d_grain_quat_refine_sched[ln] =
             d_grain_quat_refine_alg->createSchedule(level, ln - 1, hierarchy,
                                                     NULL);
      }

      d_grain_extend_coarsen_sched.resize(hierarchy->getNumberOfLevels());
      d_grain_quat_coarsen_sched.resize(hierarchy->getNumberOfLevels());

      for (int ln = ln_beg; ln < ln_end; ln++) {
         std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
         std::shared_ptr<hier::PatchLevel> finer_level =
             hierarchy->getPatchLevel(ln + 1);

         d_grain_extend_coarsen_sched[ln] =
             d_grain_extend_coarsen_alg->createSchedule(level, finer_level);
         d_grain_quat_coarsen_sched[ln] =
             d_grain_quat_coarsen_alg->createSchedule(level, finer_level);
      }
   }
}
