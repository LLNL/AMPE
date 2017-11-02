#ifndef included_DeltaTemperatureFreeEnergyStrategy
#define included_DeltaTemperatureFreeEnergyStrategy 

#include "BiasDoubleWellFreeEnergyStrategy.h"

class DeltaTemperatureFreeEnergyStrategy:
   public BiasDoubleWellFreeEnergyStrategy
{
public:
   DeltaTemperatureFreeEnergyStrategy(const double Tm,
                                      const double latentHeat);

   virtual ~DeltaTemperatureFreeEnergyStrategy(){};
 
   virtual void addComponentRhsPhi(
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
   double d_Tm;

   //latent heat
   double d_L;
};

#endif
