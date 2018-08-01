#include "FieldsInitializer.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"

using namespace std;
using namespace netCDF;

FieldsInitializer::FieldsInitializer(
   boost::shared_ptr<geom::CartesianGridGeometry >& grid_geometry,
   const hier::IntVector& ratio_of_init_to_coarsest,
   const int verbosity)
   :
   d_grid_geometry(grid_geometry),
   d_ratio_of_init_to_coarsest(ratio_of_init_to_coarsest),
   d_verbosity(verbosity),
   d_phase_id(-1),
   d_eta_id(-1),
   d_temperature_id(-1),
   d_quat_id(-1),
   d_conc_id(-1),
   d_use_uniform_q_value(false),
   d_use_uniform_c_value(false)
{
}

void FieldsInitializer::registerFieldsIds(
   const int phase_id,
   const int eta_id,
   const int temperature_id,
   const int quat_id, const int qlen,
   const int conc_id, const int ncompositions)
{
   d_phase_id       = phase_id,
   d_eta_id         = eta_id,
   d_temperature_id = temperature_id,
   d_quat_id        = quat_id,
   d_conc_id        = conc_id;

   d_qlen = qlen;
   d_ncompositions = ncompositions;
}

void FieldsInitializer::setQvalue(const vector<float>& qvalue)
{
   d_qvalue = qvalue;

   d_use_uniform_q_value = true;
}

void FieldsInitializer::setCvalue(const vector<float>& cvalue)
{
   d_cvalue = cvalue;

   d_use_uniform_c_value = true;
}

void FieldsInitializer::initializeLevelFromData(
   boost::shared_ptr< hier::PatchLevel > level,
   const string& init_data_filename,
   const int slice_index)
{
   int ln = level->getLevelNumber();
   if ( ln > 0 && d_verbosity > 0 ) {
      tbox::pout << "QuatModel::initializeLevelFromData(), lev# = "
                 << ln << endl;
   }

#ifdef HAVE_NETCDF3
   // We take care of NetCDF error checking and messages
   NcError ncerr( NcError::silent_nonfatal );
   //NcError ncerr( NcError::verbose_fatal );

   NcFile ncf( init_data_filename.c_str() );
   if ( ! ncf.is_valid() ) {
      TBOX_ERROR( "Cannot open file " << init_data_filename << endl );
   }
#endif
#ifdef HAVE_NETCDF4
   NcFile ncf( init_data_filename, NcFile::read );
   if( ncf.isNull() ){
      TBOX_ERROR( "Cannot open file " << init_data_filename << endl );
   }
   int nvar=ncf.getVarCount();
   tbox::plog << "Number of variables in NcFile: "<<nvar<<endl;
#endif
   tbox::plog<<"Opened NetCDF file "<<init_data_filename<<endl;

   size_t qlen_file = 999;
#ifdef HAVE_NETCDF3
   NcVar* ncPhase = ncf.get_var( "phase" );
   if ( ncPhase == NULL ) {
      TBOX_ERROR( "Could not read variable 'phase' from input data" << endl );
   }
   //assert( ncPhase->type() == ncFloat );

   NcVar* ncEta=NULL;
   if ( d_eta_id>=0 ) {
      ncEta = ncf.get_var( "eta" );
      if ( ncEta == NULL ) {
         TBOX_ERROR( "Could not read variable 'eta' from input data" << endl );
      }
      //assert( ncEta->type() == ncFloat );
   }

   NcVar* ncTemp = NULL;
   if ( d_model_parameters.isTemperatureConstant() ) {
      ncTemp = ncf.get_var( "temperature" );
      if ( ncTemp == NULL ) {
         TBOX_ERROR( "Could not read variable 'temperature' " <<
                     "from input data" << endl );
      }
      //assert( ncTemp->type() == ncFloat );
   }

   NcDim* ncQlen;
   if ( readQ() ) {
      ncQlen = ncf.get_dim( "qlen" );
      if ( ncQlen == NULL ) {
         TBOX_ERROR( "Could not read variable 'qlen' " <<
                     "from input data" << endl );
      }
      qlen_file = ncQlen->size();
   }
#endif
#ifdef HAVE_NETCDF4
   NcVar ncPhase;
   if ( d_phase_id>=0 ){
      ncPhase = ncf.getVar( "phase" );
      if(ncPhase.isNull())
         TBOX_ERROR( "Could not read variable 'phase' from input data" << endl );
   }
   NcVar ncEta;
   if ( d_eta_id>=0 ) {
      ncEta = ncf.getVar( "eta" );
      if(ncEta.isNull())
         TBOX_ERROR( "Could not read variable 'eta' from input data" << endl );
   }

   NcVar ncTemp;
   if ( d_temperature_id>=0 ) {
      ncTemp = ncf.getVar( "temperature" );
      if(ncTemp.isNull())
         TBOX_ERROR( "Could not read variable 'temperature' " <<
                     "from input data" << endl );
   }
   for ( int ii = 0; ii < d_qlen; ii++ ) {
      std::ostringstream o;
      o << "quat" << ii+1;
      NcVar ncv=ncf.getVar( o.str() );
      if( ncv.isNull() )break;
      qlen_file++;
   }
#endif

#ifdef HAVE_NETCDF3
   int nx_file = ncPhase->get_dim(2)->size();
   int ny_file = ncPhase->get_dim(1)->size();
   int nz_file = ncPhase->get_dim(0)->size();
#endif
#ifdef HAVE_NETCDF4
   vector<NcDim> dims;
   if ( d_phase_id>=0 ){
      dims.push_back(ncPhase.getDim(2));
      dims.push_back(ncPhase.getDim(1));
      dims.push_back(ncPhase.getDim(0));
   }else{
      dims.push_back(ncTemp.getDim(2));
      dims.push_back(ncTemp.getDim(1));
      dims.push_back(ncTemp.getDim(0));
   }
   size_t nx_file = dims[0].getSize();
   size_t ny_file = dims[1].getSize();
   size_t nz_file = dims[2].getSize();
#endif

   size_t islice = 0;
#if (NDIM == 2)
   if ( slice_index < 0 ) {
      islice = nz_file / 2;
   }
   if ( d_verbosity>0 && nz_file > 1 ) {
      tbox::pout << "Using initial data slice index " << islice << endl;
   }
#endif

   checkInputFileDimensions( nx_file, ny_file, nz_file, qlen_file );

#ifdef HAVE_NETCDF3
   NcVar** ncQuatComponents=NULL;
   if ( readQ() ) {
      ncQuatComponents = new NcVar*[d_qlen];
      for ( int ii = 0; ii < d_qlen; ii++ ) {
         std::ostringstream o;
         o << "quat" << ii+1;
         ncQuatComponents[ii] = ncf.get_var( o.str().c_str() );
         if ( ncQuatComponents[ii] == NULL ) {
            TBOX_ERROR( "Could not read variable " << o.str() <<
                        " from input data" << endl );
         }
      }
      //assert( ncQuatComponents[0]->type() == ncFloat );
   }
#endif
#ifdef HAVE_NETCDF4
   NcVar* ncQuatComponents=NULL;
   if ( readQ() ) {
      ncQuatComponents = new NcVar[d_qlen];
      for ( int ii = 0; ii < d_qlen; ii++ ) {
         std::ostringstream o;
         o << "quat" << ii+1;
         ncQuatComponents[ii] = ncf.getVar( o.str() );
         if(ncQuatComponents[ii].isNull())
            TBOX_ERROR( "Could not read variable " << o.str() <<
                        " from input data" << endl );
      }
   }
#endif

#ifdef HAVE_NETCDF3
  NcVar** ncConcComponents=NULL;
   if ( readC() ){
      tbox::pout << "With "<<d_ncompositions<<" composition fields"<<endl;
      ncConcComponents = new NcVar*[ncompositions];
      for ( int ii = 0; ii < d_ncompositions; ii++ ) {
         std::ostringstream o;
         o << "concentration";
         if( d_ncompositions>1 )o << ii;
         ncConcComponents[ii] = ncf.get_var( o.str().c_str() );
         if ( ncConcComponents[ii] == NULL && d_ncompositions==1) {
            o << 0;
            ncConcComponents[ii] = ncf.get_var( o.str().c_str() );
         }
         if ( ncConcComponents[ii] == NULL ) {
            TBOX_ERROR( "Could not read variable " << o.str() <<
                        " from input data" << endl );
         }
      }
      //assert( ncConcComponents[0]->type() == ncFloat );
   }
#endif
#ifdef HAVE_NETCDF4
   NcVar* ncConcComponents=NULL;
   if ( readC() ){
      tbox::pout << "With "<<d_ncompositions<<" composition fields"<<endl;
      ncConcComponents = new NcVar[d_ncompositions];
      for ( int ii = 0; ii < d_ncompositions; ii++ ) {
         std::ostringstream o;
         o << "concentration";
         if( d_ncompositions>1 )o << ii;
         ncConcComponents[ii] = ncf.getVar( o.str() );
         if ( ncConcComponents[ii].isNull() && d_ncompositions==1) {
            o << 0;
            ncConcComponents[ii] = ncf.getVar( o.str() );
         }
         if ( ncConcComponents[ii].isNull())
            TBOX_ERROR( "Could not read variable " << o.str() <<
                        " from input data" << endl );
      }
   }
#endif

   for ( hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p ){

      const hier::Box& patch_box = (*p)->getBox();
#ifdef HAVE_NETCDF3
      if( ncPhase->type() == ncFloat ){
#endif
#ifdef HAVE_NETCDF4
      NcType type = (d_phase_id>=0) ? ncPhase.getType() : ncTemp.getType();
      if( type.getTypeClassName() == "nc_FLOAT" ){
#endif
         float* vals = new float[patch_box.size()];

         initializePatchFromData(*p,islice, ncPhase,
            ncEta, ncTemp, ncQuatComponents, ncConcComponents, vals);
         delete[] vals;
      }else{
         double* vals = new double[patch_box.size()];

         initializePatchFromData(*p,islice, ncPhase,
             ncEta, ncTemp, ncQuatComponents, ncConcComponents, vals);

         delete[] vals;
      }

   }  // for ( hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p )

#ifdef HAVE_NETCDF3
   ncf.close();  // NcVar memory deletion handled by this call
#endif

   if ( readQ() ) {
      delete[] ncQuatComponents;
   }
   if ( readC() ){
      delete[] ncConcComponents;
   }
}

//=======================================================================

template <typename T>
void FieldsInitializer::initializePatchFromData(
   boost::shared_ptr<hier::Patch > patch,
   size_t islice,
#ifdef HAVE_NETCDF3
   NcVar* ncPhase,
   NcVar* ncEta,
   NcVar* ncTemp,
   NcVar** ncQuatComponents,
   NcVar** ncConcComponents,
#endif
#ifdef HAVE_NETCDF4
   NcVar& ncPhase,
   NcVar& ncEta,
   NcVar& ncTemp,
   NcVar* ncQuatComponents,
   NcVar* ncConcComponents,
#endif
   T* vals)
{
      const hier::Box& patch_box = patch->getBox();

      int nx = patch_box.numberCells( 0 );
      int ny = patch_box.numberCells( 1 );
      int nz = 1;
      int x_lower = patch_box.lower( 0 );
      int y_lower = patch_box.lower( 1 );
      int z_lower = static_cast<int>(islice);
      assert( x_lower>=0 );
      assert( y_lower>=0 );

#if (NDIM == 3)
      nz = patch_box.numberCells( 2 );
      z_lower = patch_box.lower( 2 );
#endif
#ifdef HAVE_NETCDF4
      vector<size_t> startp(3);
      startp[0]=z_lower;
      startp[1]=y_lower;
      startp[2]=x_lower;
      vector<size_t> countp(3);
      countp[0]=nz;
      countp[1]=ny;
      countp[2]=nx;
#endif

      // initialize phase
      if ( d_phase_id>=0 ) {
         boost::shared_ptr< pdat::CellData<double> > phase_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( d_phase_id) ) );
         assert( phase_data );

#ifdef HAVE_NETCDF3
         ncPhase->set_cur( z_lower, y_lower, x_lower );
         if ( ! ncPhase->get( vals, nz, ny, nx ) ) {
            TBOX_ERROR( "Could not read 'phase' data from input data" << endl );
         }
#endif
#ifdef HAVE_NETCDF4
         ncPhase.getVar(startp, countp, vals);
#endif

         pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
         for ( pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
                                  i!=iend; ++i ) {
            const pdat::CellIndex ccell = *i;
            int ix = ccell(0) - x_lower;
            int iy = ccell(1) - y_lower;
#if (NDIM == 2)
            int idx = nx * iy + ix;
#else
            int iz = ccell(2) - z_lower;
            int idx = nx * ny * iz + nx * iy + ix;
#endif
            (*phase_data)(ccell) = vals[idx];
         }
      }

      // initialize eta
      if ( d_eta_id>=0 ) {
         boost::shared_ptr< pdat::CellData<double> > eta_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( d_eta_id) ) );
         assert( eta_data );

#ifdef HAVE_NETCDF3
         ncEta->set_cur( z_lower, y_lower, x_lower );
         if ( ! ncEta->get( vals, nz, ny, nx ) ) {
            TBOX_ERROR( "Could not read 'eta' data from input data" << endl );
         }
#endif
#ifdef HAVE_NETCDF4
         ncEta.getVar(startp, countp, vals);
#endif

         pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
         for ( pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
                                  i!=iend; ++i ) {
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
      if ( ncTemp != NULL ) {
#endif
#ifdef HAVE_NETCDF4
      if ( !ncTemp.isNull() ){
#endif
         boost::shared_ptr< pdat::CellData<double> > temp_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_temperature_id) ) );
         assert( temp_data );

#ifdef HAVE_NETCDF3
         ncTemp->set_cur( z_lower, y_lower, x_lower );
         if ( ! ncTemp->get( vals, nz, ny, nx ) ) {
            TBOX_ERROR( "Could not read 'temperature' data from input data" << endl );
         }
#endif
#ifdef HAVE_NETCDF4
         ncTemp.getVar(startp, countp, vals);
#endif
         pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
         for ( pdat::CellIterator i(pdat::CellGeometry::begin(patch_box)); i!=iend; ++i ) {
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
      }

      // initialize quaternion
      if ( d_quat_id>=0 ){
         boost::shared_ptr< pdat::CellData<double> > quat_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( d_quat_id) ) );
         assert( quat_data );

         for ( int qq = 0; qq < d_qlen; qq++ )
         if( readQ() ){
#ifdef HAVE_NETCDF3
            NcVar* ncQuat = ncQuatComponents[qq];
            ncQuat->set_cur( z_lower, y_lower, x_lower );

            if ( ! ncQuat->get( vals, nz, ny, nx ) ) {
               TBOX_ERROR( "Could not read " << ncQuat->name() <<
                           " data from input data" << endl );
            }
#endif
#ifdef HAVE_NETCDF4
            NcVar& ncQuat = ncQuatComponents[qq];
            ncQuat.getVar(startp, countp, vals);
#endif

            pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
            for ( pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
                                     i!=iend; ++i ) {
               const pdat::CellIndex ccell = *i;
               int ix = ccell(0) - x_lower;
               int iy = ccell(1) - y_lower;
#if (NDIM == 2)
               int idx = nx * iy + ix;
#else
               int iz = ccell(2) - z_lower;
               int idx = nx * ny * iz + nx * iy + ix;
#endif
               (*quat_data)(ccell,qq) = vals[idx];
            }
         }else{
            assert( d_qvalue.size()==d_qlen );
            pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
            for ( pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
                                     i!=iend; ++i ) {
               const pdat::CellIndex ccell = *i;
               (*quat_data)(ccell,qq) = d_qvalue[qq];
            }
         }
      }

      // initialize concentration
      if ( d_ncompositions>0 ){
         assert( d_conc_id>=0 );
         boost::shared_ptr< pdat::CellData<double> > conc_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( d_conc_id) ) );
         assert( conc_data );

         for ( int cc = 0; cc < d_ncompositions; cc++ )
         if( readC() ){
#ifdef HAVE_NETCDF3
            NcVar* ncConc = ncConcComponents[cc];
            assert( ncConc!=0 );
            ncConc->set_cur( z_lower, y_lower, x_lower );
            if ( ! ncConc->get( vals, nz, ny, nx ) ) {
               TBOX_ERROR( "Could not read " << ncConc->name() <<
                           " data from input data" << endl );
            }
#endif
#ifdef HAVE_NETCDF4
            NcVar& ncConc = ncConcComponents[cc];
            ncConc.getVar(startp, countp, vals);
#endif

            pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
            for ( pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
                                     i!=iend; ++i ) {
               const pdat::CellIndex ccell = *i;
               int ix = ccell(0) - x_lower;
               int iy = ccell(1) - y_lower;
#if (NDIM == 2)
               int idx = nx * iy + ix;
#else
               int iz = ccell(2) - z_lower;
               int idx = nx * ny * iz + nx * iy + ix;
#endif
               if( vals[idx]>1. )cerr<<idx<<", vals[idx]="<<vals[idx]<<endl;
               assert( vals[idx]<=1. );
               assert( vals[idx]>=0. );
               (*conc_data)(ccell,cc) = vals[idx];
            }
         }else{
            assert( d_cvalue.size()==d_ncompositions );
            pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
            for ( pdat::CellIterator i(pdat::CellGeometry::begin(patch_box));
                                     i!=iend; ++i ) {
               const pdat::CellIndex ccell = *i;
               (*conc_data)(ccell,cc) = d_cvalue[cc];
            }
         }
#ifdef DEBUG_CHECK_ASSERTIONS
         math::PatchCellDataBasicOps<double> mathops;
         double maxc=mathops.max(conc_data, patch_box);
         assert( maxc==maxc );
#endif
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      //math::PatchCellDataBasicOps<double> mathops;
      //tbox::pout<<"max. q="<<mathops.max(quat_data,patch_box)<<endl;
      //tbox::pout<<"min. q="<<mathops.min(quat_data,patch_box)<<endl;
      //tbox::pout<<"max. p="<<mathops.max(phase_data,patch_box)<<endl;
      //tbox::pout<<"min. p="<<mathops.min(phase_data,patch_box)<<endl;
#endif
}

//=======================================================================

void FieldsInitializer::checkInputFileDimensions(
   const size_t nx_file, const size_t ny_file, const size_t nz_file,
   const size_t qlen_file )
{
   hier::Box domain_box =
      d_grid_geometry->getPhysicalDomain().front();
   domain_box.refine( d_ratio_of_init_to_coarsest );

   size_t nx_prob = domain_box.numberCells( 0 );
   size_t ny_prob = domain_box.numberCells( 1 );
   size_t nz_prob = nz_file;
#if (NDIM == 3)
   nz_prob = domain_box.numberCells( 2 );
#endif

   if ( nx_file != nx_prob ||
        ny_file != ny_prob ||
        nz_file != nz_prob ||
        static_cast<int>(qlen_file) < d_qlen ) {
      TBOX_ERROR(
         "Phase input data dimensions are incorrect"
         << ", nx_file=" << nx_file
         << ", ny_file=" << ny_file
         << ", nx_prob=" << nx_prob
         << ", ny_prob=" << ny_prob
         << ", nz_file=" << nz_file
         << ", nz_prob=" << nz_prob
         << ", qlen_file=" << qlen_file
         << ", QLEN=" << d_qlen
         << endl );
   }
}

