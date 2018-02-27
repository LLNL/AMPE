#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/CellData.h"

using namespace SAMRAI;


void copyDepthSideData(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int dst_id, const int dst_depth,
   const int src_id, const int src_depth);

void copyDepthCellData(
   const boost::shared_ptr<pdat::CellData<double> >& dst,
   const int dst_depth,
   const boost::shared_ptr<pdat::CellData<double> >& src,
   const int src_depth);

void copyDepthCellData(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int dst_id, const int dst_depth,
   const int src_id, const int src_depth);

