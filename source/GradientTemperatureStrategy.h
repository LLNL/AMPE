#ifndef included_GradientTemperatureStrategy
#define included_GradientTemperatureStrategy

#include "TemperatureStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <cassert>


// This strategy fills the temperature field from a single scalar value
// which may depend on time.

class GradientTemperatureStrategy:
   public TemperatureStrategy
{
public:
   GradientTemperatureStrategy(
      const int temperature_id,
      const int temperature_scratch_id,
      const double temperature0,
      boost::shared_ptr<tbox::Database> temperature_db );

   ~GradientTemperatureStrategy(){};

   virtual double getCurrentMaxTemperature(   
      boost::shared_ptr<hier::PatchHierarchy > patch_hierarchy,
      const double time );
   virtual double getCurrentMinTemperature(   
      boost::shared_ptr<hier::PatchHierarchy > patch_hierarchy,
      const double time );
   virtual double getCurrentAverageTemperature(
      boost::shared_ptr<hier::PatchHierarchy > patch_hierarchy,
      const double time );

   virtual void setCurrentTemperature(
      boost::shared_ptr<hier::PatchHierarchy > patch_hierarchy,
      const double time );

private:
   int d_temperature_id;
   int d_temperature_scratch_id;
   int d_weight_id;

   double d_temperature0;
   double d_dtemperaturedt;
   double d_target_temperature;

   double d_gradient[NDIM];
   double d_center[3];

   double getCurrentTemperature( const double time );

   void setCurrentTemperaturePrivatePatch(
      hier::Patch& patch,
      boost::shared_ptr< pdat::CellData<double> > cd_temp);

};

#endif
