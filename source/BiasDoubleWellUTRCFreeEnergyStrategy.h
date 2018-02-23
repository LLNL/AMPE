#ifndef included_BiasDoubleWellUTRCFreeEnergyStrategy
#define included_BiasDoubleWellUTRCFreeEnergyStrategy 

#include "BiasDoubleWellFreeEnergyStrategy.h"
#include "MeltingTemperatureStrategy.h"

/*!
 * @brief Class BiasDoubleWellFreeEnergyStrategy implements the free energy
 * used by Acharya et al., Acta Mat. 124 (2017)
 * with
 *    m(T)=(alpha/pi)tan^-1(gamma(Te-T))
 *
 */

class BiasDoubleWellUTRCFreeEnergyStrategy:
   public BiasDoubleWellFreeEnergyStrategy
{
public:
   BiasDoubleWellUTRCFreeEnergyStrategy(
      const double alpha,
      const double gamma,
      MeltingTemperatureStrategy* meltingTstrat );

   ~BiasDoubleWellUTRCFreeEnergyStrategy(){};
 
   void addComponentRhsPhi(
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

private:

   double d_alpha;
   double d_gamma;
   MeltingTemperatureStrategy* d_meltingTstrat;
};

#endif
