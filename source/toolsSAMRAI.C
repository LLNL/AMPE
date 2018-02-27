#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"

#include <boost/make_shared.hpp>

using namespace SAMRAI;



void copyDepthSideData(
   const boost::shared_ptr<pdat::SideData<double> >& dst,
   const int dst_depth,
   const boost::shared_ptr<pdat::SideData<double> >& src,
   const int src_depth)
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getGhostBox().lower()==src->getGhostBox().lower());

   const hier::Box& pbox(dst->getGhostBox());

   double* ptr_src_x = src->getPointer( 0, src_depth );
   double* ptr_src_y = src->getPointer( 1, src_depth );
   double* ptr_src_z = NULL;
   if ( NDIM > 2 ) {
      ptr_src_z = src->getPointer( 2, src_depth );
   }

   double* ptr_dst_x = dst->getPointer( 0, dst_depth );
   double* ptr_dst_y = dst->getPointer( 1, dst_depth );
   double* ptr_dst_z = NULL;
   if ( NDIM > 2 ) {
      ptr_dst_z = dst->getPointer( 2, dst_depth );
   }

   hier::IntVector nc(pbox.numberCells());

   // X-side
   unsigned n=(nc[0]+1)*nc[1];
#if (NDIM == 3)
   n*=nc[2];
#endif
   memcpy(ptr_dst_x,ptr_src_x,n*sizeof(double));

   // Y-side
   n=nc[0]*(nc[1]+1);
#if (NDIM == 3)
   n*=nc[2];
#endif
   memcpy(ptr_dst_y,ptr_src_y,n*sizeof(double));

#if (NDIM == 3)
   // Z-side
   n=nc[0]*nc[1]*(nc[2]+1);
   memcpy(ptr_dst_x,ptr_src_x,n*sizeof(double));
#endif
}

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

         copyDepthSideData(dst,dst_depth,src,src_depth);
      }
   }
}

void copyDepthCellData(
   const boost::shared_ptr<pdat::CellData<double> >& dst,
   const int dst_depth,
   const boost::shared_ptr<pdat::CellData<double> >& src,
   const int src_depth)
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getGhostBox().lower()==src->getGhostBox().lower());

   const hier::Box& pbox(dst->getGhostBox());

   double* ptr_src = src->getPointer();
   double* ptr_dst = dst->getPointer();

   hier::IntVector nc(pbox.numberCells());
   unsigned n=nc[0]*nc[1];
#if (NDIM == 3)
   n*=nc[2];
#endif

   memcpy(ptr_dst,ptr_src,n*sizeof(double));
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

         copyDepthCellData(dst,dst_depth,src,src_depth);

      }
   }
}


