/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program for testing SAMRAI data access
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#include <stdlib.h>

using namespace std;
using namespace SAMRAI;

int main(
   int argc,
   char* argv[])
{
   int error_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      /*
       * Process command line arguments.
       */
      string input_filename;

      if (argc != 2) {
         tbox::pout << "USAGE:  " << argv[0] << " <input filename> " << endl;
         exit(-1);
      } else {
         input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.
       */
      boost::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /**************************************************************************
      * Read input data and setup objects.
      **************************************************************************/

      /*
       * Retreive "Main" section of input db.
       */
      boost::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim((unsigned short)main_db->getInteger("dim"));

      const std::string log_fn =
         "iteratortest" + tbox::Utilities::intToString(dim.getValue()) + "d.log";
      tbox::PIO::logAllNodes(log_fn);

      /*
       * Regular pointer tests.
       */

      hier::Index box_lower(dim, 0);
      hier::Index box_upper(dim);

      for (int d = 0; d < dim.getValue(); d++) {
         box_upper(d) = (d + 4) * 3;
      }

      hier::Box box(box_lower, box_upper, hier::BlockId(0));

      pdat::CellData<double> cell_data(box, 1, hier::IntVector(dim, 0));
      pdat::FaceData<double> face_data(box, 1, hier::IntVector(dim, 0));
      pdat::NodeData<double> node_data(box, 1, hier::IntVector(dim, 0));
      pdat::EdgeData<double> edge_data(box, 1, hier::IntVector(dim, 0));
      pdat::SideData<double> side_data(box, 1, hier::IntVector(dim, 0));

      /*
       * The tests of the iterators first fill the patch data by directly
       * accessing the pointer to the double data.  The iterators should access
       * the data in exactly the same order that it was filled.
       */

      /*
       * The tests of the index classes loop over the box with a CellIterator.
       * The CellIterator gives each cell-centered index within the box.
       * For each cell-centered index, the tests loop over every possible
       * set of parameters for the datatype-specific index class (EdgeIndex,
       * SideIndex, etc.).  Then the test checks that accessing the data using
       * the index class retrieves data from the correct point in the data array.
       */

      /*
       * Test CellIterator
       */

      double* cell_ptr = cell_data.getPointer();

      int cell_data_size = box.size();

      for (int i = 0; i < cell_data_size; i++) {
         cell_ptr[i] = (double)i;
      }

      int j = 0;
      pdat::CellIterator ciend(box, false);
      for (pdat::CellIterator ci(box, true); ci != ciend; ++ci) {
         if (!tbox::MathUtilities<double>::equalEps(cell_data(*ci),
                cell_ptr[j])) {
            tbox::perr << "FAILED: - CellIterator test" << std::endl;
            ++error_count;
         }
         j++;
      }

      /*
       * Test FaceIterator
       */

      double* face_ptr[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

      for (int axis = 0; axis < dim.getValue(); axis++) {

         face_ptr[axis] = face_data.getPointer(axis);

         hier::Box face_box = pdat::FaceGeometry::toFaceBox(box,
               axis);

         int face_data_size = face_box.size();

         for (int i = 0; i < face_data_size; i++) {
            face_ptr[axis][i] = (double)i;
         }

         j = 0;
         pdat::FaceIterator fiend(box, axis, false);
         for (pdat::FaceIterator fi(box, axis, true); fi != fiend; ++fi) {
            if (!tbox::MathUtilities<double>::equalEps(face_data(*fi),
                   face_ptr[axis][j])) {
               tbox::perr << "FAILED: - FaceIterator test" << std::endl;
               ++error_count;
            }
            j++;
         }
      }

      /*
       * Test FaceIndex
       */

      pdat::CellIterator ifcend(box, false);
      for (pdat::CellIterator ifc(box, true); ifc != ifcend; ++ifc) {

         for (int a = 0; a < dim.getValue(); a++) {

            hier::Box face_box = pdat::FaceGeometry::toFaceBox(box,
                  a);
            hier::Index flo = face_box.lower();
            hier::Index fhi = face_box.upper();
            for (int f = 0; f < 2; f++) {

               pdat::FaceIndex findx(*ifc, a, f);

               int offset = 0;
               for (int i = dim.getValue() - 1; i > 0; i--) {
                  offset = (fhi(i - 1) - flo(i - 1) + 1)
                     * (findx(i) - flo(i) + offset);
               }
               offset += findx(0) - flo(0);

               if (!tbox::MathUtilities<double>::equalEps(face_data(findx),
                      face_ptr[a][offset])) {
                  tbox::perr << "FAILED: - FaceIndex test" << std::endl;
                  ++error_count;
               }
            }
         }
      }

      /*
       * Test NodeIterator
       */

      double* node_ptr = node_data.getPointer();

      hier::Box node_box = pdat::NodeGeometry::toNodeBox(box);

      hier::Index nlo = node_box.lower();
      hier::Index nhi = node_box.upper();

      int node_data_size = node_box.size();

      for (int i = 0; i < node_data_size; i++) {
         node_ptr[i] = (double)i;
      }

      j = 0;
      pdat::NodeIterator niend(box, false);
      for (pdat::NodeIterator ni(box, true); ni != niend; ++ni) {
         if (!tbox::MathUtilities<double>::equalEps(node_data(*ni),
                node_ptr[j])) {
            tbox::perr << "FAILED: - NodeIterator test" << std::endl;
            ++error_count;
         }
         j++;
      }

      /*
       * Test NodeIndex
       */

      pdat::CellIterator incend(box, false);
      for (pdat::CellIterator inc(box, true); inc != incend; ++inc) {

         hier::IntVector corner(dim, 0);

         bool all_corners_complete = false;

         while (!all_corners_complete) {

            pdat::NodeIndex nindx(*inc, corner);

            int offset = 0;
            for (int i = dim.getValue() - 1; i > 0; i--) {
               offset = (nhi(i - 1) - nlo(i - 1) + 1)
                  * (nindx(i) - nlo(i) + offset);
            }
            offset += nindx(0) - nlo(0);

            if (!tbox::MathUtilities<double>::equalEps(node_data(nindx),
                   node_ptr[offset])) {
               tbox::perr << "FAILED: - NodeIndex test" << std::endl;
               ++error_count;
            }

            int u;
            for (u = 0; u < dim.getValue(); u++) {
               if (corner(u) == 1) {
                  corner(u) = 0;
               } else {
                  corner(u)++;
                  break;
               }
            }
            if (u == dim.getValue()) {
               all_corners_complete = true;
            }
         }
      }

      /*
       * Test EdgeIterator
       */

      double* edge_ptr[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

      for (int axis = 0; axis < dim.getValue(); axis++) {

         edge_ptr[axis] = edge_data.getPointer(axis);

         hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box,
               axis);
         int edge_data_size = edge_box.size();

         for (int i = 0; i < edge_data_size; i++) {
            edge_ptr[axis][i] = (double)i;
         }

         j = 0;
         pdat::EdgeIterator eiend(box, axis, false);
         for (pdat::EdgeIterator ei(box, axis, true); ei != eiend; ++ei) {
            if (!tbox::MathUtilities<double>::equalEps(edge_data(*ei),
                   edge_ptr[axis][j])) {
               tbox::perr << "FAILED: - EdgeIterator test" << std::endl;
               ++error_count;
            }
            j++;
         }
      }

      /*
       * Test EdgeIndex
       */

      pdat::CellIterator iecend(box, false);
      for (pdat::CellIterator iec(box, true); iec != iecend; ++iec) {
         for (int a = 0; a < dim.getValue(); a++) {

            hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box,
                  a);
            hier::Index elo = edge_box.lower();
            hier::Index ehi = edge_box.upper();

            for (int f = 0; f < (1 << (dim.getValue() - 1)); f++) {
               pdat::EdgeIndex eindx(*iec, a, f);

               int offset = 0;
               for (int i = dim.getValue() - 1; i > 0; i--) {
                  offset = (ehi(i - 1) - elo(i - 1) + 1)
                     * (eindx(i) - elo(i) + offset);
               }
               offset += eindx(0) - elo(0);

               if (!tbox::MathUtilities<double>::equalEps(edge_data(eindx),
                      edge_ptr[a][offset])) {
                  tbox::perr << "FAILED: - EdgeIndex test" << std::endl;
                  ++error_count;
               }
            }
         }
      }

      /*
       * Test SideIterator
       */

      double* side_ptr[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

      for (int axis = 0; axis < dim.getValue(); axis++) {

         side_ptr[axis] = side_data.getPointer(axis);

         hier::Box side_box = pdat::SideGeometry::toSideBox(box,
               axis);
         int side_data_size = side_box.size();

         for (int i = 0; i < side_data_size; i++) {
            side_ptr[axis][i] = (double)i;
         }

         j = 0;
         pdat::SideIterator siend(box, axis, false);
         for (pdat::SideIterator si(box, axis, true); si != siend; ++si) {
            if (!tbox::MathUtilities<double>::equalEps(side_data(*si),
                   side_ptr[axis][j])) {
               tbox::perr << "FAILED: - SideIterator test" << std::endl;
               ++error_count;
            }
            j++;
         }
      }

      /*
       * Test SideIndex
       */

      pdat::CellIterator iscend(box, false);
      for (pdat::CellIterator isc(box, true); isc != iscend; ++isc) {
         for (int a = 0; a < dim.getValue(); a++) {

            hier::Box side_box = pdat::SideGeometry::toSideBox(box,
                  a);
            hier::Index slo = side_box.lower();
            hier::Index shi = side_box.upper();

            for (int f = 0; f < 2; f++) {
               pdat::SideIndex sindx(*isc, a, f);

               int offset = 0;
               for (int i = dim.getValue() - 1; i > 0; i--) {
                  offset = (shi(i - 1) - slo(i - 1) + 1)
                     * (sindx(i) - slo(i) + offset);
               }
               offset += sindx(0) - slo(0);

               if (!tbox::MathUtilities<double>::equalEps(side_data(sindx),
                      side_ptr[a][offset])) {
                  tbox::perr << "FAILED: - SideIndex test" << std::endl;
                  ++error_count;
               }
            }
         }
      }

      if (error_count == 0) {
         tbox::pout << "\nPASSED:  dataaccess" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return error_count;
}
