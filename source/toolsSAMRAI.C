#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include <boost/make_shared.hpp>

using namespace SAMRAI;


void copyDepthSideData(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int dst_id, const int dst_depth,
   const int src_id, const int src_depth)
{
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      boost::shared_ptr<hier::PatchLevel> level(
         hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr< pdat::SideData<double> > dst (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(
               p->getPatchData( dst_id ) ) );
         TBOX_ASSERT( dst );

         boost::shared_ptr< pdat::SideData<double> > src (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(
               p->getPatchData( src_id) ) );
         TBOX_ASSERT( src );

         dst->copyDepth(dst_depth,*src,src_depth);
      }
   }
}

void copyDepthCellData(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int dst_id, const int dst_depth,
   const int src_id, const int src_depth)
{
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      boost::shared_ptr<hier::PatchLevel> level(
         hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr< pdat::CellData<double> > dst (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               p->getPatchData( dst_id ) ) );
         TBOX_ASSERT( dst );

         boost::shared_ptr< pdat::CellData<double> > src (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               p->getPatchData( src_id) ) );
         TBOX_ASSERT( src );

         dst->copyDepth(dst_depth,*src,src_depth);

      }
   }
}


