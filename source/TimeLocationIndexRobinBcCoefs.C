/*************************************************************************
 *
 * Adapted from SAMRAI/solv/LocationIndexRobinBcCoefs class
 *
 ************************************************************************/
#include <stdlib.h>

#include "TimeLocationIndexRobinBcCoefs.h"

#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include IOMANIP_HEADER_FILE

using namespace std;

/*
 ************************************************************************
 * Constructor using database
 ************************************************************************
 */

TimeLocationIndexRobinBcCoefs::TimeLocationIndexRobinBcCoefs(
   const tbox::Dimension& dim,
   const std::string& object_name,
   const boost::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(object_name)
{
   TBOX_ASSERT(input_db);

   getFromInput(input_db);
}

/*
 ************************************************************************
 * Destructor
 ************************************************************************
 */

TimeLocationIndexRobinBcCoefs::~TimeLocationIndexRobinBcCoefs()
{
}

/*
 ********************************************************************
 * Set state from input database
 ********************************************************************
 */

void
TimeLocationIndexRobinBcCoefs::getFromInput(
   const boost::shared_ptr<tbox::Database>& input_db)
{
   if (!input_db) {
      return;
   }

   tbox::plog<<"TimeLocationIndexRobinBcCoefs::getFromInput()"<<endl;

   //loop over faces
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      std::string name = "boundary_" + tbox::Utilities::intToString(i);
      if (input_db->isString(name)) {
         std::vector<std::string> specs = input_db->getStringVector(name);
         if (specs[0] == "file"){
            boost::shared_ptr<tbox::MemoryDatabase> bc_db( new tbox::MemoryDatabase("bc_db") );
            tbox::plog<<"Parse BC input file "<<specs[1]<<endl;
            tbox::InputManager::getManager()->parseInputFile( specs[1], bc_db);
            string type = bc_db->getString( "type" );
            //Dirichlet case
            if (type == "value") {
               int j=0;
               bool flag=true;
               do{
                  std::string timestring = "time_"+tbox::Utilities::intToString(j);
                  if( bc_db->keyExists(timestring) ){
                     double tmp[2];
                     bc_db->getDoubleArray(timestring,&tmp[0],2);
                     d_a_map[i].push_back( 1.0 );
                     d_b_map[i].push_back( 0.0 );
                     d_t_map[i].push_back( tmp[0] );
                     d_g_map[i].push_back( tmp[1] );
                     j++;
                  }else{
                     flag=false;
                  }
               } while (flag);
               TBOX_ASSERT( j>0 );
            //Neumann                           
            } else if (type == "slope") {
               int j=0;
               bool flag=true;
               do{
                  std::string timestring = "time_"+tbox::Utilities::intToString(j);
                  if( bc_db->keyExists(timestring) ){
                     double tmp[2];
                     bc_db->getDoubleArray(timestring,&tmp[0],2);
                     d_a_map[i].push_back( 0.0 );
                     d_b_map[i].push_back( 1.0 );
                     d_t_map[i].push_back( tmp[0] );
                     d_g_map[i].push_back( tmp[1] );
                     j++;
                  }else{
                     flag=false;
                  }
               }while (flag);
               tbox::plog<<"Read "<<j<<" slope values"<<endl;
               TBOX_ASSERT( j>0 );
            } else {
               TBOX_ERROR(d_object_name << ": Bad boundary specifier\n"
                                        << "'" << specs[0] << "'.  Use either 'value'\n"
                                        << "'slope'.\n");
            }
         } else {
            TBOX_ERROR(d_object_name << ": Missing file specifying boundary conditions\n");
         }
      } else {
         TBOX_ERROR(d_object_name << ": Missing boundary specifier.\n");
      }
   }
}

/*
 ************************************************************************
 * Set the bc coefficients to their mapped values
 * using a linear interpolation between the times right before and 
 * right after fill_time
 ************************************************************************
 */

void
TimeLocationIndexRobinBcCoefs::setBcCoefs(
   const boost::shared_ptr<pdat::ArrayData<double> >& acoef_data,
   const boost::shared_ptr<pdat::ArrayData<double> >& bcoef_data,
   const boost::shared_ptr<pdat::ArrayData<double> >& gcoef_data,
   const boost::shared_ptr<hier::Variable>& variable,
   const hier::Patch& patch,
   const hier::BoundaryBox& bdry_box,
   double fill_time) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(d_dim, patch, bdry_box);

   NULL_USE(variable);
   NULL_USE(patch);

   int location = bdry_box.getLocationIndex();
   TBOX_ASSERT(location >= 0 && location < 2 * d_dim.getValue());
   TBOX_ASSERT(d_t_map[location].size()>0);

   static int previous_time_slot = 0;
   static int next_time_slot = 1;

   //specify how far from previous interval we should search for current interval
   int search_range = 10;

   previous_time_slot+=search_range;
   previous_time_slot=min(previous_time_slot,(int)d_t_map[location].size()-1);
   while( d_t_map[location][previous_time_slot]>fill_time ){
      previous_time_slot--;
      assert( previous_time_slot>=0 );
   }

   next_time_slot-=search_range;
   next_time_slot=max(1,next_time_slot);
   while( d_t_map[location][next_time_slot]<fill_time )next_time_slot++;
   //tbox::plog<<"d_t_map[location][next_time_slot]="<<d_t_map[location][next_time_slot]<<endl;
   //tbox::plog<<"fill_time="<<fill_time<<endl;
   //tbox::plog<<"previous_time_slot="<<previous_time_slot<<endl;
   //tbox::plog<<"next_time_slot="<<next_time_slot<<endl;
   if( (next_time_slot-previous_time_slot)!=1 ){
      tbox::plog<<"previous_time_slot="<<previous_time_slot<<endl;
      tbox::plog<<"next_time_slot="<<next_time_slot<<endl;
      TBOX_ERROR(d_object_name << ": May need larger search range for time index"<<endl);
   }
   const double wl=fill_time-d_t_map[location][previous_time_slot];
   const double wr=fill_time-d_t_map[location][next_time_slot];

   if (acoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *acoef_data);

      const double adata=wl*d_a_map[location][previous_time_slot]
                        +wr*d_a_map[location][next_time_slot];
      acoef_data->fill(adata);
   }
   if (bcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *bcoef_data);

      const double bdata=wl*d_b_map[location][previous_time_slot]
                        +wr*d_b_map[location][next_time_slot];
      bcoef_data->fill(bdata);
   }
   if (gcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *gcoef_data);

      const double gdata=wl*d_g_map[location][previous_time_slot]
                        +wr*d_g_map[location][next_time_slot];
      gcoef_data->fill(gdata);
   }
   
}

void
TimeLocationIndexRobinBcCoefs::rescaleGcoefficients(const double factor)
{
   //loop over faces
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      vector<double>::const_iterator it=d_g_map[i].begin();
      for( vector<double>::iterator it=d_g_map[i].begin();
           it!=d_g_map[i].end();
           ++it)
      {
         (*it)*=factor;
      }
   }
}

hier::IntVector
TimeLocationIndexRobinBcCoefs::numberOfExtensionsFillable() const
{
   /*
    * Return some really big number.  We have no limits.
    */
   return hier::IntVector(d_dim, 1 << (sizeof(int) - 1));
}

/*
 ************************************************************************
 * Assignment operator
 ************************************************************************
 */

TimeLocationIndexRobinBcCoefs&
TimeLocationIndexRobinBcCoefs::operator = (
   const TimeLocationIndexRobinBcCoefs& r)
{
   d_object_name = r.d_object_name;
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      d_a_map[i] = r.d_a_map[i];
      d_b_map[i] = r.d_b_map[i];
      d_g_map[i] = r.d_g_map[i];
   }
   return *this;
}

