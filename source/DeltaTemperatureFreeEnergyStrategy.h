#ifndef included_DeltaTemperatureFreeEnergyStrategy
#define included_DeltaTemperatureFreeEnergyStrategy 

#include "BiasDoubleWellFreeEnergyStrategy.h"

#include <cstring>

class DeltaTemperatureFreeEnergyStrategy:
   public BiasDoubleWellFreeEnergyStrategy
{
public:
   DeltaTemperatureFreeEnergyStrategy(const double Tm,
                                      const double latentHeat,
                                      const std::string phase_interp_func_type);

   virtual ~DeltaTemperatureFreeEnergyStrategy(){};
 
   virtual void addComponentRhsPhi(
      const double time,
      hier::Patch& patch,
      const int temperature_id,
      const int phase_id,
      const int eta_id,
      const int conc_id, 
      const int fl_id,
      const int fa_id,
      const int fb_id,
      const int rhs_id );

   void applydPhidTBlock(const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int phase_id,
      const int rhs_id,
      const double phase_mobility);

private:

   //melting temperature
   const double d_Tm;

   //latent heat
   const double d_L;

   const std::string d_phase_interp_func_type;
};

#endif
