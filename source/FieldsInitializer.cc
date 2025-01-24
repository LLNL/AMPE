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
#include "FieldsInitializer.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"


#ifdef HAVE_NETCDF4
using namespace netCDF;
#endif

FieldsInitializer::FieldsInitializer(
    std::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const hier::IntVector& ratio_of_init_to_coarsest, const int verbosity)
    : d_grid_geometry(grid_geometry),
      d_ratio_of_init_to_coarsest(ratio_of_init_to_coarsest),
      d_verbosity(verbosity),
      d_phase_id(-1),
      d_eta_id(-1),
      d_temperature_id(-1),
      d_quat_id(-1),
      d_conc_id(-1),
      d_use_uniform_q_value(false),
      d_use_uniform_c_value(false),
      d_use_uniform_t_value(false)
{
}

void FieldsInitializer::registerFieldsIds(const int phase_id, const int eta_id,
                                          const int temperature_id,
                                          const int quat_id, const int qlen,
                                          const int conc_id,
                                          const int ncompositions,
                                          const int nphases)
{
   tbox::plog << "FieldsInitializer::registerFieldsIds with " << nphases
              << " phases" << std::endl;
   tbox::plog << "qlen = " << qlen << std::endl;
   d_phase_id = phase_id, d_eta_id = eta_id, d_temperature_id = temperature_id,
   d_quat_id = quat_id, d_conc_id = conc_id;

   d_qlen = qlen;
   d_ncompositions = ncompositions;
   d_nphases = nphases;
}

void FieldsInitializer::setQvalue(const std::vector<float>& qvalue)
{
   d_qvalue = qvalue;

   d_use_uniform_q_value = true;
}

void FieldsInitializer::setCvalue(const std::vector<float>& cvalue)
{
   d_cvalue = cvalue;

   d_use_uniform_c_value = true;
}

void FieldsInitializer::setTvalue(const float tvalue)
{
   d_tvalue = tvalue;

   d_use_uniform_t_value = true;
}

void FieldsInitializer::initializeLevelFromData(
    std::shared_ptr<hier::PatchLevel> level,
    const std::string& init_data_filename, const int slice_index)
{
   int ln = level->getLevelNumber();
   if (ln > 0 && d_verbosity > 0) {
      tbox::pout << "QuatModel::initializeLevelFromData(), lev# = " << ln
                 << std::endl;
   }

   std::unique_ptr<NcFile> ncf;
#ifdef HAVE_NETCDF3
   // We take care of NetCDF error checking and messages
   NcError ncerr(NcError::silent_nonfatal);
   // NcError ncerr( NcError::verbose_fatal );

   ncf.reset(new NcFile(init_data_filename.c_str()));
   if (!ncf->is_valid()) {
      TBOX_ERROR("Cannot open file " << init_data_filename << std::endl);
   }
#endif
#ifdef HAVE_NETCDF4
   if (!init_data_filename.empty()) {
      ncf.reset(new NcFile(init_data_filename, NcFile::read));
      if (ncf->isNull()) {
         TBOX_ERROR("Cannot open file " << init_data_filename << std::endl);
      }
   }
   int nvar = init_data_filename.empty() ? 0 : ncf->getVarCount();
   tbox::plog << "Number of variables in NcFile: " << nvar << std::endl;
#endif
   tbox::plog << "Opened NetCDF file " << init_data_filename << std::endl;

   size_t qlen_file = 0;
   NcVar* ncPhase = nullptr;
#ifdef HAVE_NETCDF3
   tbox::plog << "Try to read phase values..." << std::endl;
   ncPhase = ncf->get_var("phase");
   if (ncPhase == nullptr) {
      TBOX_ERROR("Could not read variable 'phase' from input data"
                 << std::endl);
   }
   // assert( ncPhase->type() == ncFloat );

   NcVar* ncEta = nullptr;
   if (d_eta_id >= 0) {
      ncEta = ncf->get_var("eta");
      if (ncEta == nullptr) {
         TBOX_ERROR("Could not read variable 'eta' from input data"
                    << std::endl);
      }
      // assert( ncEta->type() == ncFloat );
   }

   NcVar* ncTemp = nullptr;
   if (readT() && d_temperature_id >= 0) {
      ncTemp = ncf->get_var("temperature");
      if (ncTemp == nullptr) {
         TBOX_ERROR("Could not read variable 'temperature' "
                    << "from input data" << std::endl);
      }
      // assert( ncTemp->type() == ncFloat );
   }

   NcDim* ncQlen;
   if (readQ()) {
      ncQlen = ncf->get_dim("qlen");
      if (ncQlen == nullptr) {
         TBOX_ERROR("Could not read variable 'qlen' "
                    << "from input data" << std::endl);
      }
      qlen_file = ncQlen->size();
   }
#endif
#ifdef HAVE_NETCDF4
   tbox::plog << "Reading fields..." << std::endl;
   if (d_phase_id >= 0) {
      ncPhase = new NcVar[d_nphases];
      tbox::plog << "Reading phase values..." << std::endl;
      for (int i = 0; i < d_nphases; i++) {
         std::ostringstream o;
         o << "phase";
         ncPhase[i] = ncf->getVar(o.str());
         if (ncPhase[i].isNull()) {
            // add phase id to name
            o << i;
            ncPhase[i] = ncf->getVar(o.str());
            if (ncPhase[i].isNull() && (i == 0 || i < d_nphases))
               TBOX_ERROR("Could not read variable "
                          << o.str() << " from input data" << std::endl);
         }
      }
   }
   NcVar ncEta;
   if (d_eta_id >= 0) {
      ncEta = ncf->getVar("eta");
      if (ncEta.isNull())
         TBOX_ERROR("Could not read variable 'eta' from input data"
                    << std::endl);
   }

   NcVar ncTemp;
   if (readT() && d_temperature_id >= 0) {
      tbox::plog << "Reading temperature values..." << std::endl;
      ncTemp = ncf->getVar("temperature");
      if (ncTemp.isNull())
         TBOX_ERROR("Could not read variable 'temperature' "
                    << "from input data" << std::endl);
   }
   if (readQ())
      for (int ii = 0; ii < d_qlen; ii++) {
         tbox::plog << "Reading quaternion values..." << std::endl;
         std::ostringstream o;
         o << "quat" << ii + 1;
         NcVar ncv = ncf->getVar(o.str());
         if (ncv.isNull()) break;
         qlen_file++;
      }
#endif

   size_t nx_prob;
   size_t ny_prob;
   size_t nz_prob;
   getDomainSizes(nx_prob, ny_prob, nz_prob);

#ifdef HAVE_NETCDF3
   int nx_file = ncPhase->get_dim(2)->size();
   int ny_file = ncPhase->get_dim(1)->size();
   int nz_file = ncPhase->get_dim(0)->size();
#endif
#ifdef HAVE_NETCDF4
   size_t nx_file = nx_prob;
   size_t ny_file = ny_prob;
   size_t nz_file = nz_prob;
   std::vector<NcDim> dims;
   if (d_phase_id >= 0) {
      dims.push_back(ncPhase[0].getDim(2));
      dims.push_back(ncPhase[0].getDim(1));
      dims.push_back(ncPhase[0].getDim(0));
   } else if (readT()) {
      if (ncTemp.isNull())
         TBOX_ERROR("Field 'temperature' was not read" << std::endl);
      dims.push_back(ncTemp.getDim(2));
      dims.push_back(ncTemp.getDim(1));
      dims.push_back(ncTemp.getDim(0));
   }
   if (!dims.empty()) {
      nx_file = dims[0].getSize();
      ny_file = dims[1].getSize();
      nz_file = dims[2].getSize();
   }
#endif

   size_t islice = 0;
#if (NDIM == 2)
   if (slice_index < 0) {
      islice = nz_file / 2;
   }
   if (d_verbosity > 0 && nz_file > 1) {
      tbox::pout << "Using initial data slice index " << islice << std::endl;
   }
#endif

   checkInputFileDimensions(nx_file, ny_file, nz_file, qlen_file);

#ifdef HAVE_NETCDF3
   NcVar** ncQuatComponents = nullptr;
   if (readQ()) {
      ncQuatComponents = new NcVar*[d_qlen];
      for (int ii = 0; ii < d_qlen; ii++) {
         std::ostringstream o;
         o << "quat" << ii + 1;
         ncQuatComponents[ii] = ncf->get_var(o.str().c_str());
         if (ncQuatComponents[ii] == nullptr) {
            TBOX_ERROR("Could not read variable "
                       << o.str() << " from input data" << std::endl);
         }
      }
      // assert( ncQuatComponents[0]->type() == ncFloat );
   }
#endif
#ifdef HAVE_NETCDF4
   NcVar* ncQuatComponents = nullptr;
   if (readQ()) {
      ncQuatComponents = new NcVar[d_qlen];
      for (int ii = 0; ii < d_qlen; ii++) {
         std::ostringstream o;
         o << "quat" << ii + 1;
         ncQuatComponents[ii] = ncf->getVar(o.str());
         if (ncQuatComponents[ii].isNull())
            TBOX_ERROR("Could not read variable "
                       << o.str() << " from input data" << std::endl);
      }
   }
#endif

#ifdef HAVE_NETCDF3
   NcVar** ncConcComponents = nullptr;
   if (readC()) {
      tbox::pout << "With " << d_ncompositions << " composition fields"
                 << std::endl;
      ncConcComponents = new NcVar*[d_ncompositions];
      for (int ii = 0; ii < d_ncompositions; ii++) {
         std::ostringstream o;
         o << "concentration";
         if (d_ncompositions > 1) o << ii;
         ncConcComponents[ii] = ncf->get_var(o.str().c_str());
         if (ncConcComponents[ii] == nullptr && d_ncompositions == 1) {
            o << 0;
            ncConcComponents[ii] = ncf->get_var(o.str().c_str());
         }
         if (ncConcComponents[ii] == nullptr) {
            TBOX_ERROR("Could not read variable "
                       << o.str() << " from input data" << std::endl);
         }
      }
      // assert( ncConcComponents[0]->type() == ncFloat );
   }
#endif
#ifdef HAVE_NETCDF4
   NcVar* ncConcComponents = nullptr;
   if (readC()) {
      tbox::pout << "With " << d_ncompositions << " composition fields"
                 << std::endl;
      ncConcComponents = new NcVar[d_ncompositions];
      for (int ii = 0; ii < d_ncompositions; ii++) {
         std::ostringstream o;
         o << "concentration";
         if (d_ncompositions > 1) o << ii;
         ncConcComponents[ii] = ncf->getVar(o.str());
         if (ncConcComponents[ii].isNull() && d_ncompositions == 1) {
            o << 0;
            ncConcComponents[ii] = ncf->getVar(o.str());
         }
         if (ncConcComponents[ii].isNull())
            TBOX_ERROR("Could not read variable "
                       << o.str() << " from input data" << std::endl);
      }
   }
#endif

   for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

      const hier::Box& patch_box = (*p)->getBox();
#ifdef HAVE_NETCDF3
      if (ncPhase->type() == ncFloat) {
#endif
#ifdef HAVE_NETCDF4
         bool float_flag = true;
         if (d_phase_id >= 0) {
            NcType type = ncPhase[0].getType();
            if (type.getTypeClassName() == "nc_DOUBLE") float_flag = false;
         } else {
            if (readT()) {
               NcType type = ncTemp.getType();
               if (type.getTypeClassName() == "nc_DOUBLE") float_flag = false;
            }
         }

         if (float_flag) {
#endif
            float* vals = new float[patch_box.size()];

            initializePatchFromData(*p, islice, ncPhase, ncEta, ncTemp,
                                    ncQuatComponents, ncConcComponents, vals);
            delete[] vals;
         } else {
            double* vals = new double[patch_box.size()];

            initializePatchFromData(*p, islice, ncPhase, ncEta, ncTemp,
                                    ncQuatComponents, ncConcComponents, vals);

            delete[] vals;
         }

      }  // end loop over patches

#ifdef HAVE_NETCDF3
      ncf->close();  // NcVar memory deletion handled by this call
#endif

      if (readQ()) {
         delete[] ncQuatComponents;
      }
      if (readC()) {
         delete[] ncConcComponents;
      }
   }

   //=======================================================================

   template <typename T>
   void FieldsInitializer::initializePatchFromData(
       std::shared_ptr<hier::Patch> patch, size_t islice,
#ifdef HAVE_NETCDF3
       NcVar * ncPhase, NcVar * ncEta, NcVar * ncTemp,
       NcVar * *ncQuatComponents, NcVar * *ncConcComponents,
#endif
#ifdef HAVE_NETCDF4
       NcVar * ncPhase, NcVar & ncEta, NcVar & ncTemp, NcVar * ncQuatComponents,
       NcVar * ncConcComponents,
#endif
       T * vals)
   {
      const hier::Box& patch_box = patch->getBox();

      int nx = patch_box.numberCells(0);
      int ny = patch_box.numberCells(1);
      int nz = 1;
      int x_lower = patch_box.lower(0);
      int y_lower = patch_box.lower(1);
      int z_lower = static_cast<int>(islice);
      assert(x_lower >= 0);
      assert(y_lower >= 0);

#if (NDIM == 3)
      nz = patch_box.numberCells(2);
      z_lower = patch_box.lower(2);
#endif
#ifdef HAVE_NETCDF4
      std::vector<size_t> startp(3);
      startp[0] = z_lower;
      startp[1] = y_lower;
      startp[2] = x_lower;
      std::vector<size_t> countp(3);
      countp[0] = nz;
      countp[1] = ny;
      countp[2] = nx;
#endif

      // initialize phase
      if (d_phase_id >= 0) {
         std::shared_ptr<pdat::CellData<double> > phase_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_phase_id)));
         assert(phase_data);

         for (int depth = 0; depth < d_nphases; depth++) {
            tbox::plog << "Read order parameter " << depth << std::endl;
#ifdef HAVE_NETCDF3
            ncPhase->set_cur(z_lower, y_lower, x_lower);
            if (!ncPhase->get(vals, nz, ny, nx)) {
               TBOX_ERROR("Could not read 'phase' data from input data"
                          << std::endl);
            }
#endif
#ifdef HAVE_NETCDF4
            ncPhase[depth].getVar(startp, countp, vals);
#endif

            pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
            for (pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
                 i != iend; ++i) {
               const pdat::CellIndex ccell = *i;
               int ix = ccell(0) - x_lower;
               int iy = ccell(1) - y_lower;
#if (NDIM == 2)
               int idx = nx * iy + ix;
#else
            int iz = ccell(2) - z_lower;
            int idx = nx * ny * iz + nx * iy + ix;
#endif
               (*phase_data)(ccell, depth) = vals[idx];
            }
         }
         // if number of order parameters > d_nphases, set value of last one as
         // 1. - sum of others
         if (phase_data->getDepth() > d_nphases) {
            tbox::plog << "Set value of last order parameter to 1 minus sum of "
                          "the others"
                       << std::endl;
            pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
            for (pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
                 i != iend; ++i) {
               const pdat::CellIndex ccell = *i;
               double val = 0.;
               for (int depth = 0; depth < d_nphases; depth++)
                  val += (*phase_data)(ccell, depth);
               (*phase_data)(ccell, d_nphases) = 1. - val;
            }
         }
      }

      // initialize eta
      if (d_eta_id >= 0) {
         std::shared_ptr<pdat::CellData<double> > eta_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_eta_id)));
         assert(eta_data);

#ifdef HAVE_NETCDF3
         ncEta->set_cur(z_lower, y_lower, x_lower);
         if (!ncEta->get(vals, nz, ny, nx)) {
            TBOX_ERROR("Could not read 'eta' data from input data"
                       << std::endl);
         }
#endif
#ifdef HAVE_NETCDF4
         ncEta.getVar(startp, countp, vals);
#endif

         pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
              i != iend; ++i) {
            const pdat::CellIndex ccell = *i;
            int ix = ccell(0) - x_lower;
            int iy = ccell(1) - y_lower;
#if (NDIM == 2)
            int idx = nx * iy + ix;
#else
         int iz = ccell(2) - z_lower;
         int idx = nx * ny * iz + nx * iy + ix;
#endif
            (*eta_data)(ccell) = vals[idx];
         }
      }

      // initialize temperature
#ifdef HAVE_NETCDF3
      if (readT() && ncTemp != nullptr) {
#endif
#ifdef HAVE_NETCDF4
         if (readT() && !ncTemp.isNull()) {
#endif
            std::shared_ptr<pdat::CellData<double> > temp_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_temperature_id)));
            assert(temp_data);

#ifdef HAVE_NETCDF3
            ncTemp->set_cur(z_lower, y_lower, x_lower);
            if (!ncTemp->get(vals, nz, ny, nx)) {
               TBOX_ERROR("Could not read 'temperature' data from input data"
                          << std::endl);
            }
#endif
#ifdef HAVE_NETCDF4
            ncTemp.getVar(startp, countp, vals);
#endif
            pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
            for (pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
                 i != iend; ++i) {
               const pdat::CellIndex ccell = *i;
               int ix = ccell(0) - x_lower;
               int iy = ccell(1) - y_lower;
#if (NDIM == 2)
               int idx = nx * iy + ix;
#else
      int iz = ccell(2) - z_lower;
      int idx = nx * ny * iz + nx * iy + ix;
#endif
               (*temp_data)(ccell) = vals[idx];
            }
         } else if (d_temperature_id >= 0 && !readT()) {
            std::shared_ptr<pdat::CellData<double> > temp_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_temperature_id)));
            assert(temp_data);

            temp_data->fill(d_tvalue);
         }

         // initialize quaternion
         if (d_quat_id >= 0) {
            std::shared_ptr<pdat::CellData<double> > quat_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_quat_id)));
            assert(quat_data);

            for (int qq = 0; qq < d_qlen; qq++)
               if (readQ()) {
#ifdef HAVE_NETCDF3
                  NcVar* ncQuat = ncQuatComponents[qq];
                  ncQuat->set_cur(z_lower, y_lower, x_lower);

                  if (!ncQuat->get(vals, nz, ny, nx)) {
                     TBOX_ERROR("Could not read " << ncQuat->name()
                                                  << " data from input data"
                                                  << std::endl);
                  }
#endif
#ifdef HAVE_NETCDF4
                  NcVar& ncQuat = ncQuatComponents[qq];
                  ncQuat.getVar(startp, countp, vals);
#endif

                  pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
                  for (pdat::CellIterator i(
                           pdat::CellGeometry::begin(patch_box));
                       i != iend; ++i) {
                     const pdat::CellIndex ccell = *i;
                     int ix = ccell(0) - x_lower;
                     int iy = ccell(1) - y_lower;
#if (NDIM == 2)
                     int idx = nx * iy + ix;
#else
            int iz = ccell(2) - z_lower;
            int idx = nx * ny * iz + nx * iy + ix;
#endif
                     (*quat_data)(ccell, qq) = vals[idx];
                  }
               } else {
                  if (static_cast<int>(d_qvalue.size()) != d_qlen) {
                     TBOX_ERROR("need " << d_qlen << " values to specify q");
                  }
                  pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
                  for (pdat::CellIterator i(
                           pdat::CellGeometry::begin(patch_box));
                       i != iend; ++i) {
                     const pdat::CellIndex ccell = *i;
                     (*quat_data)(ccell, qq) = d_qvalue[qq];
                  }
               }
         }

         // initialize concentration
         if (d_ncompositions > 0) {
            assert(d_conc_id >= 0);
            std::shared_ptr<pdat::CellData<double> > conc_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_conc_id)));
            assert(conc_data);

            for (int cc = 0; cc < d_ncompositions; cc++)
               if (readC()) {
#ifdef HAVE_NETCDF3
                  NcVar* ncConc = ncConcComponents[cc];
                  assert(ncConc != 0);
                  ncConc->set_cur(z_lower, y_lower, x_lower);
                  if (!ncConc->get(vals, nz, ny, nx)) {
                     TBOX_ERROR("Could not read " << ncConc->name()
                                                  << " data from input data"
                                                  << std::endl);
                  }
#endif
#ifdef HAVE_NETCDF4
                  NcVar& ncConc = ncConcComponents[cc];
                  ncConc.getVar(startp, countp, vals);
#endif

                  pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
                  for (pdat::CellIterator i(
                           pdat::CellGeometry::begin(patch_box));
                       i != iend; ++i) {
                     const pdat::CellIndex ccell = *i;
                     int ix = ccell(0) - x_lower;
                     int iy = ccell(1) - y_lower;
#if (NDIM == 2)
                     int idx = nx * iy + ix;
#else
            int iz = ccell(2) - z_lower;
            int idx = nx * ny * iz + nx * iy + ix;
#endif
                     if (vals[idx] > 1.)
                        std::cerr << idx << ", vals[idx]=" << vals[idx]
                                  << std::endl;
                     assert(vals[idx] <= 1.);
                     assert(vals[idx] >= 0.);
                     (*conc_data)(ccell, cc) = vals[idx];
                  }
               } else {
                  if (static_cast<int>(d_cvalue.size()) != d_ncompositions) {
                     TBOX_ERROR("Expect " << d_ncompositions
                                          << " to specify alloy concentration");
                  }
                  pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
                  for (pdat::CellIterator i(
                           pdat::CellGeometry::begin(patch_box));
                       i != iend; ++i) {
                     const pdat::CellIndex ccell = *i;
                     (*conc_data)(ccell, cc) = d_cvalue[cc];
                  }
               }
#ifdef DEBUG_CHECK_ASSERTIONS
            math::PatchCellDataBasicOps<double> mathops;
            double maxc = mathops.max(conc_data, patch_box);
            assert(maxc == maxc);
#endif
         }

#ifdef DEBUG_CHECK_ASSERTIONS
         // math::PatchCellDataBasicOps<double> mathops;
         // tbox::pout<<"max. q="<<mathops.max(quat_data,patch_box)<<endl;
         // tbox::pout<<"min. q="<<mathops.min(quat_data,patch_box)<<endl;
         // tbox::pout<<"max. p="<<mathops.max(phase_data,patch_box)<<endl;
         // tbox::pout<<"min. p="<<mathops.min(phase_data,patch_box)<<endl;
#endif
      }

      //=======================================================================

      void FieldsInitializer::getDomainSizes(size_t & nx_prob, size_t & ny_prob,
                                             size_t & nz_prob)
      {
         hier::Box domain_box = d_grid_geometry->getPhysicalDomain().front();
         domain_box.refine(d_ratio_of_init_to_coarsest);

         nx_prob = domain_box.numberCells(0);
         ny_prob = domain_box.numberCells(1);
#if (NDIM == 3)
         nz_prob = domain_box.numberCells(2);
#else
   nz_prob = 1;
#endif
      }

      //=======================================================================

      void FieldsInitializer::checkInputFileDimensions(const size_t nx_file,
                                                       const size_t ny_file,
                                                       const size_t nz_file,
                                                       const size_t qlen_file)
      {
         size_t nx_prob;
         size_t ny_prob;
         size_t nz_prob;
         getDomainSizes(nx_prob, ny_prob, nz_prob);

         if (nx_file != nx_prob || ny_file != ny_prob || nz_file != nz_prob) {
            TBOX_ERROR("Phase input data dimensions are incorrect"
                       << ", nx_file=" << nx_file << ", ny_file=" << ny_file
                       << ", nz_file=" << nz_file << ", nx_prob=" << nx_prob
                       << ", ny_prob=" << ny_prob << ", nz_prob=" << nz_prob
                       << std::endl);
         }
         if (readQ() && (static_cast<int>(qlen_file) != d_qlen)) {
            TBOX_ERROR("Phase input data dimensions are incorrect"
                       << ", qlen_file=" << qlen_file << ", QLEN=" << d_qlen
                       << std::endl);
         }
      }
