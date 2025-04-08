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
#include "FieldsWriter.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

#define HAVE_NETCDF4

#ifdef HAVE_NETCDF4
#include <netcdf>
using namespace netCDF;
#endif
#ifdef HAVE_NETCDF3
#include "netcdfcpp.h"
#endif

FieldsWriter::FieldsWriter(
    QuatModelParameters& model_parameters, std::string filename,
    int initial_conditions_level,
    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry, int phase_id,
    int phase_scratch_id, int temperature_id, int temperature_scratch_id,
    int quat_id, int quat_scratch_id, int conc_id, int conc_scratch_id,
    int eta_id, int eta_scratch_id, const int ncompositions, const int qlen,
    QuatRefinePatchStrategy* all_refine_patch_strategy)
    : d_model_parameters(model_parameters),
      d_filename(filename),
      d_initial_conditions_level(initial_conditions_level),
      d_grid_geometry(grid_geometry),
      d_phase_id(phase_id),
      d_phase_scratch_id(phase_scratch_id),
      d_temperature_id(temperature_id),
      d_temperature_scratch_id(temperature_scratch_id),
      d_quat_id(quat_id),
      d_quat_scratch_id(quat_scratch_id),
      d_conc_id(conc_id),
      d_conc_scratch_id(conc_scratch_id),
      d_eta_scratch_id(eta_scratch_id),
      d_ncompositions(ncompositions),
      d_qlen(qlen),
      d_all_refine_patch_strategy(all_refine_patch_strategy)
{
}

//=======================================================================

void FieldsWriter::writeInitialConditionsFile(
    const std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
    const double time)
{
   assert(d_initial_conditions_level >= 0);
   assert(d_initial_conditions_level < 10);

   tbox::plog << "Write initial conditions file..." << std::endl;

   // get new PatchLevel with uniform mesh at level
   // "d_initial_conditions_level"
   std::shared_ptr<hier::PatchLevel> flattened_level =
       FlattenHierarchy(patch_hierarchy, d_initial_conditions_level, time);

   // get size of uniform mesh to write
   hier::BoxContainer boxes = d_grid_geometry->getPhysicalDomain();
   assert(boxes.size() == 1);

   hier::Box bf(boxes.front());
   bf.refine(flattened_level->getRatioToLevelZero());

   int nx_prob = bf.numberCells(0);
   int ny_prob = bf.numberCells(1);
   int nz_prob = 1;
#if (NDIM > 2)
   nz_prob = bf.numberCells(2);
#endif

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   const int npp = mpi.getSize();

   for (int pp = 0; pp < npp; pp++) {

      // tbox::plog << "pp=" << pp << std::endl;
      if (mpi.getRank() == pp) {

         NcFile* f;
         NcVar* nc_phase;
         if (d_model_parameters.use_FolchPlapp())
            nc_phase = new NcVar[3];
         else
            nc_phase = new NcVar[1];
#ifdef HAVE_NETCDF3
         NcVar* nc_eta = nullptr;
         NcVar** nc_conc = new NcVar*[d_ncompositions];
         NcVar** nc_qcomp = new NcVar*[d_qlen];
         NcVar* nc_temp = nullptr;
#endif
#ifdef HAVE_NETCDF4
         NcVar nc_eta;
         NcVar* nc_conc = new NcVar[d_ncompositions];
         NcVar* nc_qcomp = new NcVar[d_qlen];
         NcVar nc_temp;
#endif

         if (pp == 0) {
#ifdef HAVE_NETCDF3
            f = new NcFile(d_filename.c_str(), NcFile::Replace);
            if (!f->is_valid()) {
               TBOX_ERROR("Cannot open file " << d_filename << std::endl);
            }
#endif
#ifdef HAVE_NETCDF4
            f = new NcFile(d_filename, NcFile::replace);
            if (f->isNull()) {
               TBOX_ERROR("Cannot open file " << d_filename << std::endl);
            } else {
               std::clog << "Open/replace file " << d_filename << std::endl;
            }
#endif

#ifdef HAVE_NETCDF3
            NcDim* nc_nx = f->add_dim("x", nx_prob);
            NcDim* nc_ny = f->add_dim("y", ny_prob);
            NcDim* nc_nz = f->add_dim("z", nz_prob);
            f->add_dim("qlen", d_qlen);

            if (d_model_parameters.use_FolchPlapp()) {
               nc_phase[0] = f->add_var("phase0", ncFloat, nc_nz, nc_ny, nc_nx);
               nc_phase[1] = f->add_var("phase1", ncFloat, nc_nz, nc_ny, nc_nx);
               nc_phase[2] = f->add_var("phase2", ncFloat, nc_nz, nc_ny, nc_nx);
            } else
               nc_phase[0] = f->add_var("phase", ncFloat, nc_nz, nc_ny, nc_nx);

            if (d_model_parameters.with_third_phase()) {
               nc_eta = f->add_var("eta", ncFloat, nc_nz, nc_ny, nc_nx);
            }

            if (d_model_parameters.with_orientation()) {
               for (int ii = 0; ii < d_qlen; ii++) {
                  std::ostringstream o;
                  o << "quat" << ii + 1;
                  nc_qcomp[ii] =
                      f->add_var(o.str().c_str(), ncFloat, nc_nz, nc_ny, nc_nx);
               }
            }

            if (d_model_parameters.with_concentration()) {
               for (int ii = 0; ii < d_ncompositions; ii++) {
                  std::ostringstream o;
                  o << "concentration";
                  if (d_ncompositions > 1) o << ii;
                  nc_conc[ii] =
                      f->add_var(o.str().c_str(), ncFloat, nc_nz, nc_ny, nc_nx);
               }
            }

            nc_temp = f->add_var("temperature", ncFloat, nc_nz, nc_ny, nc_nx);

#endif
#ifdef HAVE_NETCDF4
            // std::cout<<"add variables from PE 0..."<<endl;
            NcDim nc_nx = f->addDim("x", nx_prob);
            NcDim nc_ny = f->addDim("y", ny_prob);
            NcDim nc_nz = f->addDim("z", nz_prob);
            // f->addDim( "qlen", d_qlen );

            std::vector<NcDim> dims;
            dims.push_back(nc_nz);
            dims.push_back(nc_ny);
            dims.push_back(nc_nx);
            if (d_model_parameters.use_FolchPlapp()) {
               nc_phase[0] = f->addVar("phase0", ncFloat, dims);
               nc_phase[1] = f->addVar("phase1", ncFloat, dims);
               nc_phase[2] = f->addVar("phase2", ncFloat, dims);
            } else
               nc_phase[0] = f->addVar("phase", ncFloat, dims);
            if (nc_phase[0].isNull()) {
               TBOX_ERROR("Could add variable 'phase'" << std::endl);
            }
            if (d_model_parameters.with_third_phase()) {
               nc_eta = f->addVar("eta", ncFloat, dims);
            }

            if (d_model_parameters.with_orientation()) {
               for (int ii = 0; ii < d_qlen; ii++) {
                  std::ostringstream o;
                  o << "quat" << ii + 1;
                  nc_qcomp[ii] = f->addVar(o.str(), ncFloat, dims);
               }
            }

            if (d_model_parameters.with_concentration()) {
               for (int ii = 0; ii < d_ncompositions; ii++) {
                  std::ostringstream o;
                  o << "concentration";
                  if (d_ncompositions > 1) o << ii;
                  nc_conc[ii] = f->addVar(o.str(), ncFloat, dims);
               }
            }

            nc_temp = f->addVar("temperature", ncFloat, dims);
            // std::cout<<"variables added on PE 0..."<<endl;
#endif
         } else {  // pp!=0
#ifdef HAVE_NETCDF3
            f = new NcFile(d_filename.c_str(), NcFile::Write);
            if (!f->is_valid()) {
               TBOX_ERROR("Cannot open file " << d_filename << std::endl);
            }
#endif
#ifdef HAVE_NETCDF4
            f = new NcFile(d_filename, NcFile::write);
            if (f->isNull()) {
               TBOX_ERROR("Cannot open file " << d_filename << std::endl);
               //}else{
               //   clog<<"Open/write file
               //   "<<d_filename<<endl;
            }
#endif

#ifdef HAVE_NETCDF3
            nc_phase = f->get_var("phase");

            if (d_model_parameters.with_third_phase()) {
               nc_eta = f->get_var("eta");
            }

            if (d_model_parameters.with_orientation()) {
               for (int ii = 0; ii < d_qlen; ii++) {
                  std::ostringstream o;
                  o << "quat" << ii + 1;
                  nc_qcomp[ii] = f->get_var(o.str().c_str());
               }
            }

            if (d_model_parameters.with_concentration()) {
               for (int ii = 0; ii < d_ncompositions; ii++) {
                  std::ostringstream o;
                  o << "concentration";
                  if (d_ncompositions > 1) o << ii;
                  nc_conc[ii] = f->get_var(o.str().c_str());
               }
            }

            nc_temp = f->get_var("temperature");
#endif
#ifdef HAVE_NETCDF4
            // clog<<"add variables from PE >0..."<<endl;
            if (d_model_parameters.use_FolchPlapp()) {
               nc_phase[0] = f->getVar("phase0");
               nc_phase[1] = f->getVar("phase1");
               nc_phase[2] = f->getVar("phase2");
            } else
               nc_phase[0] = f->getVar("phase");

            if (d_model_parameters.with_third_phase()) {
               nc_eta = f->getVar("eta");
            }

            if (d_model_parameters.with_orientation()) {
               for (int ii = 0; ii < d_qlen; ii++) {
                  std::ostringstream o;
                  o << "quat" << ii + 1;
                  nc_qcomp[ii] = f->getVar(o.str());
               }
            }

            if (d_model_parameters.with_concentration()) {
               for (int ii = 0; ii < d_ncompositions; ii++) {
                  std::ostringstream o;
                  o << "concentration";
                  if (d_ncompositions > 1) o << ii;
                  nc_conc[ii] = f->getVar(o.str());
               }
            }

            nc_temp = f->getVar("temperature");
#endif
         }  // pp==0 or not

#ifdef HAVE_NETCDF3
         if (nc_phase == nullptr) {
            TBOX_ERROR("Could not create variable 'phase'" << std::endl);
         }

         if (d_model_parameters.with_third_phase() && nc_eta == nullptr) {
            TBOX_ERROR("Could not create variable 'eta'" << std::endl);
         }

         if (d_model_parameters.with_orientation()) {
            for (int ii = 0; ii < d_qlen; ii++) {
               std::ostringstream o;
               o << "quat" << ii + 1;
               if (nc_qcomp[ii] == nullptr) {
                  TBOX_ERROR("Could not create variable " << o.str()
                                                          << std::endl);
               }
            }
         }

         if (d_model_parameters.with_concentration()) {
            for (int ii = 0; ii < d_ncompositions; ii++) {
               std::ostringstream o;
               o << "concentration";
               if (d_ncompositions > 1) o << ii;
               if (nc_conc[ii] == nullptr) {
                  TBOX_ERROR("Could not create variable " << o.str()
                                                          << std::endl);
               }
            }
         }

         if (nc_temp == nullptr) {
            TBOX_ERROR("Could not create variable 'temperature'" << std::endl);
         }
#endif

#ifdef HAVE_NETCDF4
         if (nc_phase[0].isNull()) {
            TBOX_ERROR("Could not create variable 'phase'" << std::endl);
         }

         if (d_model_parameters.with_third_phase() && nc_eta.isNull()) {
            TBOX_ERROR("Could not create variable 'eta'" << std::endl);
         }

         if (d_model_parameters.with_orientation()) {
            for (int ii = 0; ii < d_qlen; ii++) {
               std::ostringstream o;
               o << "quat" << ii + 1;
               if (nc_qcomp[ii].isNull()) {
                  TBOX_ERROR("Could not create variable " << o.str()
                                                          << std::endl);
               }
            }
         }

         if (d_model_parameters.with_concentration()) {
            for (int ii = 0; ii < d_ncompositions; ii++) {
               std::ostringstream o;
               o << "concentration";
               if (d_ncompositions > 1) o << ii;
               if (nc_conc[ii].isNull()) {
                  TBOX_ERROR("Could not create variable " << o.str()
                                                          << std::endl);
               }
            }
         }

         if (nc_temp.isNull()) {
            TBOX_ERROR("Could not create variable 'temperature'" << std::endl);
         }
#endif

         // std::cout<<"Write data into variable objects..."<<endl;
         for (hier::PatchLevel::Iterator p(flattened_level->begin());
              p != flattened_level->end(); ++p) {
            std::shared_ptr<hier::Patch> patch = *p;

            std::shared_ptr<pdat::CellData<double> > phase_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_phase_id)));
            assert(phase_data);

            const hier::Box& this_b = patch->getBox();
            int nx = this_b.numberCells(0);
            int ny = this_b.numberCells(1);
            int nz = 1;
#if (NDIM > 2)
            nz = this_b.numberCells(2);
#endif
            int lowz = 0;
#if (NDIM > 2)
            lowz = this_b.lower(2);
#endif

#ifdef HAVE_NETCDF3
            nc_phase->set_cur(lowz, this_b.lower(1), this_b.lower(0));
            nc_phase->put(phase_data->getPointer(), nz, ny, nx);
#endif
#ifdef HAVE_NETCDF4
            std::vector<size_t> startp(3);
            startp[0] = lowz;
            startp[1] = this_b.lower(1);
            startp[2] = this_b.lower(0);
            std::vector<size_t> countp(3);
            countp[0] = nz;
            countp[1] = ny;
            countp[2] = nx;

            // std::cout<<"Write data into variable 'phase'"<<endl;
            // std::cout<<"nx="<<countp[0]<<", ny="<<countp[1]<<",
            // nz="<<countp[2]<<endl;
            nc_phase[0].putVar(startp, countp, phase_data->getPointer(0));
            if (d_model_parameters.use_FolchPlapp()) {
               nc_phase[1].putVar(startp, countp, phase_data->getPointer(1));
               nc_phase[2].putVar(startp, countp, phase_data->getPointer(2));
            }
            // std::cout<<"Data written into variable 'phase'"<<endl;
#endif
            if (d_model_parameters.with_third_phase()) {
               std::shared_ptr<pdat::CellData<double> > eta_data(
                   SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                                          hier::PatchData>(
                       patch->getPatchData(d_eta_id)));
               assert(eta_data);

#ifdef HAVE_NETCDF3
               nc_eta->set_cur(lowz, this_b.lower(1), this_b.lower(0));
               nc_eta->put(eta_data->getPointer(), nz, ny, nx);
#endif
#ifdef HAVE_NETCDF4
               nc_eta.putVar(startp, countp, eta_data->getPointer());
#endif
            }

            std::shared_ptr<pdat::CellData<double> > temp_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_temperature_id)));
            assert(temp_data);

#ifdef HAVE_NETCDF3
            nc_temp->set_cur(lowz, this_b.lower(1), this_b.lower(0));
            nc_temp->put(temp_data->getPointer(), nz, ny, nx);
#endif
#ifdef HAVE_NETCDF4
            // std::cout<<"Write data into variable 'temperature'"<<endl;
            nc_temp.putVar(startp, countp, temp_data->getPointer());
#endif

            if (d_model_parameters.with_orientation()) {
               std::shared_ptr<pdat::CellData<double> > quat_data(
                   SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                                          hier::PatchData>(
                       patch->getPatchData(d_quat_id)));
               assert(quat_data);

               for (int dd = 0; dd < d_qlen; dd++) {
#ifdef HAVE_NETCDF3
                  nc_qcomp[dd]->set_cur(lowz, this_b.lower(1), this_b.lower(0));
                  nc_qcomp[dd]->put(quat_data->getPointer(dd), nz, ny, nx);
#endif
#ifdef HAVE_NETCDF4
                  nc_qcomp[dd].putVar(startp, countp,
                                      quat_data->getPointer(dd));
#endif
               }
            }

            if (d_model_parameters.with_concentration()) {
               std::shared_ptr<pdat::CellData<double> > conc_data(
                   SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                                          hier::PatchData>(
                       patch->getPatchData(d_conc_id)));
               assert(conc_data);

               for (int dd = 0; dd < d_ncompositions; dd++) {
#ifdef HAVE_NETCDF3
                  nc_conc[dd]->set_cur(lowz, this_b.lower(1), this_b.lower(0));
                  nc_conc[dd]->put(conc_data->getPointer(dd), nz, ny, nx);
#endif
#ifdef HAVE_NETCDF4
                  nc_conc[dd].putVar(startp, countp, conc_data->getPointer(dd));
#endif
               }
            }
         }
         // std::cout<<"Close file..."<<endl;
#ifdef HAVE_NETCDF3
         f->close();
#endif
#ifdef HAVE_NETCDF4
         delete f;
#endif
         delete[] nc_qcomp;
         delete[] nc_conc;
      }
      mpi.Barrier();
   }  // for ( int pp ...
}

//=======================================================================

std::shared_ptr<hier::PatchLevel> FieldsWriter::FlattenHierarchy(
    const std::shared_ptr<hier::PatchHierarchy> src_hierarchy,
    const int level_number, const double time)
{
   tbox::plog << "FieldsWriter::FlattenHierarchy..." << std::endl;

   assert(level_number >= 0);
   assert(level_number < 10);

   // early return for single level case
   if (level_number == 0) return src_hierarchy->getPatchLevel(level_number);

   // It will be assumed that the levels of src_hierarchy are
   // synchronized at this point.

   hier::IntVector ratio(tbox::Dimension(NDIM), 1);
   for (int l = level_number; l > 0; l--)
      ratio *= src_hierarchy->getRatioToCoarserLevel(l);

   tbox::plog << "ratio: " << ratio << std::endl;

   // Compute physical domain box array describing the index space of the
   // physical domain managed by this geometry object. If any entry of ratio
   // std::vector is negative, the index space is coarsened with respect to
   // the physical domain description. Otherwise, the index space is refined.

   // get boxes corresponding to level 0 of hierarchy, then refine them to
   // "level_number"
   hier::BoxLevel layer0(ratio, d_grid_geometry);
   hier::BoxContainer boxes;
   std::shared_ptr<hier::PatchLevel> zero_level =
       src_hierarchy->getPatchLevel(0);
   // iterate over patches
   for (hier::PatchLevel::Iterator p(zero_level->begin());
        p != zero_level->end(); ++p) {
      const hier::Box& box = (*p)->getBox();
      boxes.pushBack(box);
   }

   // hier::BoxContainer boxes ( zero_level->getBoxes() );
   boxes.refine(ratio);
   // boxes.print(cout);

   hier::BoxContainer::const_iterator boxes_itr = boxes.begin();
   for (int ib = 0; ib < boxes.size(); ib++, boxes_itr++) {
      layer0.addBox(
          hier::Box(*boxes_itr, hier::LocalId(ib), layer0.getMPI().getRank()));
   }

   std::shared_ptr<hier::PatchLevel> src_level =
       src_hierarchy->getPatchLevel(level_number);
   std::shared_ptr<hier::PatchLevel> flattened_level(new hier::PatchLevel(
       layer0, d_grid_geometry,
       hier::VariableDatabase::getDatabase()->getPatchDescriptor()));

   // allocate data on newly created level
   flattened_level->allocatePatchData(d_phase_id);
   flattened_level->allocatePatchData(d_phase_scratch_id);
   flattened_level->setTime(time, d_phase_id);
   flattened_level->setTime(time, d_phase_scratch_id);

   if (d_model_parameters.with_third_phase()) {
      flattened_level->allocatePatchData(d_eta_id);
      flattened_level->allocatePatchData(d_eta_scratch_id);
      flattened_level->setTime(time, d_eta_id);
      flattened_level->setTime(time, d_eta_scratch_id);
   }

   if (d_temperature_id >= 0) {
      flattened_level->allocatePatchData(d_temperature_id);
      flattened_level->allocatePatchData(d_temperature_scratch_id);
      flattened_level->setTime(time, d_temperature_id);
      flattened_level->setTime(time, d_temperature_scratch_id);
   }
   if (d_model_parameters.with_orientation()) {
      flattened_level->allocatePatchData(d_quat_id);
      flattened_level->allocatePatchData(d_quat_scratch_id);
      flattened_level->setTime(time, d_quat_id);
      flattened_level->setTime(time, d_quat_scratch_id);
   }
   if (d_model_parameters.with_concentration()) {
      flattened_level->allocatePatchData(d_conc_id);
      flattened_level->allocatePatchData(d_conc_scratch_id);
      flattened_level->setTime(time, d_conc_id);
      flattened_level->setTime(time, d_conc_scratch_id);
   }

   // fill new level with data
   tbox::plog << "Fill new level with data..." << std::endl;
   std::shared_ptr<hier::Variable> variable;
   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
   if (d_model_parameters.with_phase())
      vdb->mapIndexToVariable(d_phase_id, variable);
   else
      vdb->mapIndexToVariable(d_temperature_id, variable);

   std::shared_ptr<hier::RefineOperator> cell_refine_op;
   cell_refine_op =
       d_grid_geometry->lookupRefineOperator(variable, "LINEAR_REFINE");

   std::shared_ptr<xfer::RefineAlgorithm> curr_to_curr_refine_alg;
   curr_to_curr_refine_alg.reset(new xfer::RefineAlgorithm());
   if (d_model_parameters.with_phase()) {
      curr_to_curr_refine_alg->registerRefine(d_phase_id,  // destination
                                              d_phase_id,  // source
                                              d_phase_scratch_id,  // temporary
                                              cell_refine_op);
   }
   if (d_temperature_id >= 0) {
      assert(d_temperature_scratch_id >= 0);
      curr_to_curr_refine_alg->registerRefine(
          d_temperature_id,          // destination
          d_temperature_id,          // source
          d_temperature_scratch_id,  // temporary
          cell_refine_op);
   }
   if (d_model_parameters.with_third_phase()) {
      curr_to_curr_refine_alg->registerRefine(d_eta_id,          // destination
                                              d_eta_id,          // source
                                              d_eta_scratch_id,  // temporary
                                              cell_refine_op);
   }
   if (d_model_parameters.with_orientation()) {
      assert(d_quat_scratch_id >= 0);
      // we may need a different refine operator to take symmetry into
      // account here
      curr_to_curr_refine_alg->registerRefine(d_quat_id,          // destination
                                              d_quat_id,          // source
                                              d_quat_scratch_id,  // temporary
                                              cell_refine_op);
   }

   if (d_model_parameters.with_concentration()) {
      assert(d_conc_scratch_id >= 0);
      curr_to_curr_refine_alg->registerRefine(d_conc_id,          // destination
                                              d_conc_id,          // source
                                              d_conc_scratch_id,  // temporary
                                              cell_refine_op);
   }

   std::shared_ptr<xfer::RefineSchedule> schedule(
       curr_to_curr_refine_alg->createSchedule(flattened_level, src_level,
                                               level_number - 1, src_hierarchy,
                                               d_all_refine_patch_strategy));

   schedule->fillData(time);

   tbox::plog << "Flatten done..." << std::endl;

   return flattened_level;
}
