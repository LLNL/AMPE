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
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

#define HAVE_NETCDF4
// #define HAVE_NETCDF3

#ifdef HAVE_NETCDF3
#include "netcdfcpp.h"
#endif
#ifdef HAVE_NETCDF4
#include <netcdf>
#endif

using namespace SAMRAI;

class FieldsInitializer
{
 public:
   FieldsInitializer(
       std::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
       const hier::IntVector& ratio_of_init_to_coarsest, const int verbosity);

   void registerFieldsIds(const int phase_id, const int eta_id,
                          const int temperature_id, const int quat_id,
                          const int qlen, const int conc_id,
                          const int ncompositions, const int nphases);

   void initializeLevelFromData(std::shared_ptr<hier::PatchLevel> level,
                                const std::string& filename,
                                const int slice_index);

   template <typename T>
   void initializePatchFromData(
       std::shared_ptr<hier::Patch> patch, size_t islice,
#ifdef HAVE_NETCDF3
       NcVar* ncPhase, NcVar* ncEta, NcVar* ncTemp, NcVar** ncQuatComponents,
       NcVar** ncConcComponents,
#endif
#ifdef HAVE_NETCDF4
       netCDF::NcVar* ncPhase, netCDF::NcVar& ncEta, netCDF::NcVar& ncTemp,
       netCDF::NcVar* ncQuatComponents, netCDF::NcVar* ncConcComponents,
#endif
       T* vals);

   void setQvalue(const std::vector<float>& qvalue);
   void setCvalue(const std::vector<float>& cvalue);
   void setTvalue(const float tvalue);

 private:
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   const hier::IntVector d_ratio_of_init_to_coarsest;

   const int d_verbosity;

   int d_phase_id;
   int d_eta_id;
   int d_temperature_id;
   int d_quat_id;
   int d_conc_id;

   int d_qlen;
   int d_ncompositions;
   int d_nphases;

   /*!
    * Uniform initial values for fields
    */
   std::vector<float> d_qvalue;
   std::vector<float> d_cvalue;
   float d_tvalue;

   bool d_use_uniform_q_value;
   bool d_use_uniform_c_value;
   bool d_use_uniform_t_value;

   bool readQ() const { return (d_qlen > 0 && !d_use_uniform_q_value); }
   bool readC() const
   {
      return (d_ncompositions > 0 && !d_use_uniform_c_value);
   }
   bool readT() const { return (!d_use_uniform_t_value); }

   void checkInputFileDimensions(const size_t nx_file, const size_t ny_file,
                                 const size_t nz_file, const size_t qlen_file);

   void getDomainSizes(size_t&, size_t&, size_t&);
};
