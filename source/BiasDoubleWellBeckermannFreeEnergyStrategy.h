#ifndef included_BiasDoubleWellBeckermannFreeEnergyStrategy
#define included_BiasDoubleWellBeckermannFreeEnergyStrategy 

#include "BiasDoubleWellFreeEnergyStrategy.h"
#include "MeltingTemperatureStrategy.h"

/*!
 * @brief Class BiasDoubleWellFreeEnergyStrategy implements the free energy
 * used by Beckermann et al., J. Comput. Phys. 154 (1999)
 * with
 *    m(T)=alpha*(Te-T)
 *
 */

class BiasDoubleWellBeckermannFreeEnergyStrategy:
   public BiasDoubleWellFreeEnergyStrategy
{
public:
   BiasDoubleWellBeckermannFreeEnergyStrategy(
      const double alpha,
      MeltingTemperatureStrategy* meltingTstrat );

   ~BiasDoubleWellBeckermannFreeEnergyStrategy(){};
 
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
   MeltingTemperatureStrategy* d_meltingTstrat;
};

#endif
