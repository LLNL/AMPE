#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

void computeQDiffs(
    std::shared_ptr<SAMRAI::pdat::CellData<double> > quat_data,
    std::shared_ptr<SAMRAI::pdat::SideData<double> > diff_data,
    const bool isSymmetryAware,
    std::shared_ptr<SAMRAI::pdat::SideData<int> > rotation_index);
