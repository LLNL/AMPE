#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <boost/make_shared.hpp>

#define HAVE_NETCDF4

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
   boost::shared_ptr<geom::CartesianGridGeometry >& grid_geometry,
   const hier::IntVector& ratio_of_init_to_coarsest,
   const int verbosity);

void registerFieldsIds(
   const int phase_id,
   const int eta_id,
   const int temperature_id,
   const int quat_id, const int qlen,
   const int conc_id, const int ncompositions);

void initializeLevelFromData(
   boost::shared_ptr< hier::PatchLevel > level,
   const std::string& filename,
   const int slice_index);

template <typename T>
void initializePatchFromData(
   boost::shared_ptr<hier::Patch > patch,
   size_t islice,
#ifdef HAVE_NETCDF3
   netCDF::NcVar* ncPhase,
   netCDF::NcVar* ncEta,
   netCDF::NcVar* ncTemp,
   netCDF::NcVar** ncQuatComponents,
   netCDF::NcVar** ncConcComponents,
#endif
#ifdef HAVE_NETCDF4
   netCDF::NcVar& ncPhase,
   netCDF::NcVar& ncEta,
   netCDF::NcVar& ncTemp,
   netCDF::NcVar* ncQuatComponents,
   netCDF::NcVar* ncConcComponents,
#endif
   T* vals);

   void setQvalue(const std::vector<float>& qvalue);
   void setCvalue(const std::vector<float>& cvalue);

private:

   boost::shared_ptr<geom::CartesianGridGeometry > d_grid_geometry;

   const hier::IntVector d_ratio_of_init_to_coarsest;

   const int d_verbosity;

   int d_phase_id;
   int d_eta_id;
   int d_temperature_id;
   int d_quat_id;
   int d_conc_id;

   int d_qlen;
   int d_ncompositions;

   /*!
    * Uniform initial values for fields
    */
   std::vector<float> d_qvalue;
   std::vector<float> d_cvalue;

   bool d_use_uniform_q_value;
   bool d_use_uniform_c_value;

   bool readQ()const{ return (d_qlen>0 && !d_use_uniform_q_value); }
   bool readC()const{ return (d_ncompositions>0 && !d_use_uniform_c_value); }

   void checkInputFileDimensions(
      const size_t nx_file, const size_t ny_file, const size_t nz_file,
      const size_t qlen_file );
};

