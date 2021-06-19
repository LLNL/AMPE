#ifndef included_DiffusionCoeffForQuat
#define included_DiffusionCoeffForQuat

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <cstring>

using namespace SAMRAI;

class DiffusionCoeffForQuat
{
 public:
   DiffusionCoeffForQuat(
       std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
       const double H_parameter, std::string orient_interp_func_type,
       std::string avg_func_type,
       std::shared_ptr<pdat::SideVariable<double> > quat_diffusion_var,
       const int quat_diffusion_id);

   void setup(const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void setDiffusion(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                     const int phase_scratch_id,
                     const int temperature_scratch_id);

 private:
   double d_H_parameter;
   std::string d_orient_interp_func_type;
   std::string d_avg_func_type;

   int d_quat_diffusion_id;

   xfer::CoarsenAlgorithm d_quat_diffusion_coarsen;

   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
       d_quat_diffusion_coarsen_schedule;

   void setDiffusionPatch(hier::Patch& patch, const int phase_scratch_id,
                          const int temperature_scratch_id);
};

#endif
