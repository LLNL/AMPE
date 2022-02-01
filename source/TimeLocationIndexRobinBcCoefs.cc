/*************************************************************************
 * Adapted from SAMRAI/solv/LocationIndexRobinBcCoefs.C
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 ************************************************************************/
#include <stdlib.h>

#include "TimeLocationIndexRobinBcCoefs.h"

#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include IOMANIP_HEADER_FILE


/*
 ************************************************************************
 * Constructor using database
 ************************************************************************
 */

TimeLocationIndexRobinBcCoefs::TimeLocationIndexRobinBcCoefs(
    const tbox::Dimension& dim, const std::string& object_name,
    const std::shared_ptr<tbox::Database>& input_db)
    : d_dim(dim), d_object_name(object_name)
{
   TBOX_ASSERT(input_db);

   getFromInput(input_db);
}

/*
 ************************************************************************
 * Destructor
 ************************************************************************
 */

TimeLocationIndexRobinBcCoefs::~TimeLocationIndexRobinBcCoefs() {}

/*
 ********************************************************************
 * Set state from input database
 ********************************************************************
 */

void TimeLocationIndexRobinBcCoefs::getFromInput(
    const std::shared_ptr<tbox::Database>& input_db)
{
   if (!input_db) {
      return;
   }

   tbox::plog << "TimeLocationIndexRobinBcCoefs::getFromInput()" << std::endl;

   // loop over faces
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      std::string name = "boundary_" + std::to_string(i);
      if (input_db->isString(name)) {
         std::vector<std::string> specs = input_db->getStringVector(name);
         if (specs[0] == "file") {
            std::shared_ptr<tbox::MemoryDatabase> bc_db(
                new tbox::MemoryDatabase("bc_db"));
            tbox::plog << "Parse BC input file " << specs[1] << std::endl;
            tbox::InputManager::getManager()->parseInputFile(specs[1], bc_db);
            std::string type = bc_db->getString("type");
            // Dirichlet case
            if (type == "value") {
               int j = 0;
               bool flag = true;
               do {
                  std::string timestring = "time_" + std::to_string(j);
                  if (bc_db->keyExists(timestring)) {
                     double tmp[2];
                     bc_db->getDoubleArray(timestring, &tmp[0], 2);
                     d_a_map[i].push_back(1.0);
                     d_b_map[i].push_back(0.0);
                     d_t_map[i].push_back(tmp[0]);
                     d_g_map[i].push_back(tmp[1]);
                     j++;
                  } else {
                     flag = false;
                  }
               } while (flag);
               TBOX_ASSERT(j > 0);
               // Neumann
            } else if (type == "slope") {
               int j = 0;
               bool flag = true;
               do {
                  std::string timestring = "time_" + std::to_string(j);
                  if (bc_db->keyExists(timestring)) {
                     double tmp[2];
                     bc_db->getDoubleArray(timestring, &tmp[0], 2);
                     d_a_map[i].push_back(0.0);
                     d_b_map[i].push_back(1.0);
                     d_t_map[i].push_back(tmp[0]);
                     d_g_map[i].push_back(tmp[1]);
                     j++;
                  } else {
                     flag = false;
                  }
               } while (flag);
               tbox::plog << "Read " << j << " slope values" << std::endl;
               TBOX_ASSERT(j > 0);
            } else {
               TBOX_ERROR(d_object_name << ": Bad boundary specifier\n"
                                        << "'" << specs[0]
                                        << "'.  Use either 'value'\n"
                                        << "'slope'.\n");
            }
         } else {
            TBOX_ERROR(d_object_name << ": Missing file specifying boundary "
                                        "conditions\n");
         }
      } else {
         TBOX_ERROR(d_object_name << ": Missing boundary specifier.\n");
      }
   }
}

/*
 ************************************************************************
 * Set the bc coefficients to their std::mapped values
 * using a linear interpolation between the times right before and
 * right after fill_time
 ************************************************************************
 */

void TimeLocationIndexRobinBcCoefs::setBcCoefs(
    const std::shared_ptr<pdat::ArrayData<double> >& acoef_data,
    const std::shared_ptr<pdat::ArrayData<double> >& bcoef_data,
    const std::shared_ptr<pdat::ArrayData<double> >& gcoef_data,
    const std::shared_ptr<hier::Variable>& variable, const hier::Patch& patch,
    const hier::BoundaryBox& bdry_box, double fill_time) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY2(d_dim, patch, bdry_box);

   NULL_USE(variable);
   NULL_USE(patch);

   int location = bdry_box.getLocationIndex();
   TBOX_ASSERT(location >= 0 && location < 2 * d_dim.getValue());
   TBOX_ASSERT(d_t_map[location].size() > 0);

   static int prev_time_slot[6] = {0, 0, 0, 0, 0, 0};
   static int next_time_slot[6] = {1, 1, 1, 1, 1, 1};

   // specify how far from previous interval we should search for current
   // interval. This is necessary since this function may be called for a
   // time value smaller than in previous call
   int search_range = 10;

   prev_time_slot[location] -= search_range;
   prev_time_slot[location] = std::max(prev_time_slot[location], 1);
   // loop over time slot until we reach a time larger than fill_time
   while (d_t_map[location][prev_time_slot[location]] <= fill_time) {
      prev_time_slot[location]++;
      assert(prev_time_slot[location] > 0);
   }
   prev_time_slot[location]--;
   assert(prev_time_slot[location] >= 0);

   next_time_slot[location] = prev_time_slot[location] + 1;

   // tbox::plog<<"d_t_map[location][next_time_slot[location]]="
   //          <<d_t_map[location][next_time_slot[location]]<<endl;
   // tbox::pout<<"fill_time="<<fill_time<<endl;
   // tbox::pout<<"prev_time_slot="<<prev_time_slot[location]<<endl;
   // tbox::pout<<"next_time_slot="<<next_time_slot[location]<<endl;
   if (d_t_map[location][prev_time_slot[location]] > fill_time ||
       d_t_map[location][next_time_slot[location]] < fill_time) {
      tbox::plog << "fill_time = " << fill_time << std::endl;
      tbox::plog << "previous_time = "
                 << d_t_map[location][prev_time_slot[location]] << std::endl;
      tbox::plog << "next_time = "
                 << d_t_map[location][next_time_slot[location]] << std::endl;
      TBOX_ERROR(d_object_name << ": May need larger search range for time "
                                  "index"
                               << std::endl);
   }

   const int ntime_slot = next_time_slot[location];
   const int ptime_slot = prev_time_slot[location];
   TBOX_ASSERT(fill_time <= d_t_map[location][ntime_slot]);
   TBOX_ASSERT(fill_time >= d_t_map[location][ptime_slot]);

   const double dtinv =
       1. / (d_t_map[location][ntime_slot] - d_t_map[location][ptime_slot]);
   const double wl = dtinv * (fill_time - d_t_map[location][ptime_slot]);
   const double wr = dtinv * (d_t_map[location][ntime_slot] - fill_time);

   // tbox::pout<<"wl="<<wl<<endl;
   // tbox::pout<<"wr="<<wr<<endl;
   if (acoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *acoef_data);

      const double adata = wr * d_a_map[location][ptime_slot] +
                           wl * d_a_map[location][ntime_slot];
      acoef_data->fill(adata);
      // tbox::pout<<"adata="<<adata<<endl;
   }
   if (bcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *bcoef_data);

      const double bdata = wr * d_b_map[location][ptime_slot] +
                           wl * d_b_map[location][ntime_slot];
      bcoef_data->fill(bdata);
      // tbox::pout<<"bdata="<<bdata<<endl;
   }
   if (gcoef_data) {
      TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *gcoef_data);

      const double gdata = wr * d_g_map[location][ptime_slot] +
                           wl * d_g_map[location][ntime_slot];
      gcoef_data->fill(gdata);
      // tbox::pout<<"gdata="<<gdata<<endl;
   }
}

void TimeLocationIndexRobinBcCoefs::rescaleGcoefficients(const double factor)
{
   // loop over faces
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      std::vector<double>::const_iterator it = d_g_map[i].begin();
      for (std::vector<double>::iterator it = d_g_map[i].begin();
           it != d_g_map[i].end(); ++it) {
         (*it) *= factor;
      }
   }
}

hier::IntVector TimeLocationIndexRobinBcCoefs::numberOfExtensionsFillable()
    const
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

TimeLocationIndexRobinBcCoefs& TimeLocationIndexRobinBcCoefs::operator=(
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
